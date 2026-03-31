# SINGLE-SAMPLE ANALYSIS — AD  (script v5)
# v5 = v4 pipeline + Seurat v5 RNA fix: JoinLayers(Assay5) then as(., "Assay") so CellChat/NeuronChat
# see a classic v4-like Assay with @data / @counts slots (reduces "data layer" warnings).
# v4 layer: Zhao-style coarse plots (celltype group=) when meta.data$celltype exists.
# Extends v2 with: net_aggregation, netVisual_*, heatmaps, pathway loops, rankNet,
# netAnalysis_river_Neuron, pattern pheatmap/write.table exports, functional + structural
# manifold plots, lig_tar_heatmap for all interactions, CSV/txt side outputs.
# Inference uses subcelltype (group.by); celltype is only for optional coarse visualization.
# Outputs: logs, figures, tables, RDS use *_v5 so they do not overwrite v1/v2/v3/v4.

library(uwot)
library(NeuronChat)
library(CellChat)
library(Seurat)
library(ggplot2)
library(patchwork)
library(ggalluvial)
library(ComplexHeatmap)
library(circlize)
library(future)
library(Matrix)
library(pheatmap)
library(dplyr)
library(tidyr)
library(tibble)

options(stringsAsFactors = FALSE)
plan("sequential")

# ── Sample / paths ────────────────────────────────────────────────
sample_name <- "18_merged_AG_031826.28celltypes.AD"
base_dir <- "/mnt/nfs/CX000008_DS1/projects/jaeyeon/fastq_file_Dyslexia_r1/zanalysis_bogdan/samples.merged_AG_Harmony_res0.x_anno_031826"

rds_seurat <- file.path(base_dir, "18_merged_AG_031826.28celltypes.AD.rds")
out_dir <- file.path(base_dir, "18_merged_AG_031826.28celltypes.AD_neuronchat")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
fig_dir <- file.path(out_dir, paste0("results_NeuronChat_", sample_name, "_v5_figures"))
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

exports_dir <- file.path(out_dir, paste0("exports_NeuronChat_", sample_name, "_v5"))
dir.create(exports_dir, showWarnings = FALSE, recursive = TRUE)

tmpdir <- "/mnt/nfs/CX000008_DS1/projects/btanasa/tmp"
Sys.setenv(TMPDIR = tmpdir)
options(future.globals.maxSize = 200 * 1024^3)

safe_name <- function(x) gsub("[^A-Za-z0-9_\\-]+", "_", x)

# ── Logging ───────────────────────────────────────────────────────
log_file <- file.path(out_dir,
  paste0("neuronchat_AD_28celltypes_v5_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
start_time <- Sys.time()
sink(log_file, split = TRUE)

cat("╔═══════════════════════════════════════════════════════════╗\n")
cat("║  NeuronChat Single-Sample: AD v5 (subcelltype / 28ct)   ║\n")
cat("╚═══════════════════════════════════════════════════════════╝\n")
cat("Started:", format(start_time), "\n")
cat("Figures:", fig_dir, "\n")
cat("Exports:", exports_dir, "\n\n")

on.exit({
  if (sink.number() > 0) {
    cat("\nEnded  :", format(Sys.time()), "\n")
    cat("Runtime:", round(difftime(Sys.time(), start_time, units = "mins"), 2), "min\n")
    sink()
  }
}, add = TRUE)

# ── Memory check ─────────────────────────────────────────────────
check_memory <- function(threshold = 0.90) {
  meminfo <- tryCatch(readLines("/proc/meminfo", n = 20), error = function(e) NULL)
  if (is.null(meminfo)) return(invisible(NULL))
  total <- as.numeric(gsub("[^0-9]", "", meminfo[1]))
  avail <- as.numeric(gsub("[^0-9]", "", meminfo[3]))
  pct <- (total - avail) / total
  cat(sprintf("💾 Memory: %.1f%% used (%.1f / %.1f GB)\n",
    pct * 100, (total - avail) / 1024^2, total / 1024^2))
  if (pct > threshold) {
    cat("❌ Memory threshold exceeded — stopping\n")
    sink()
    stop("Memory threshold exceeded")
  }
  invisible(pct)
}

# ── Load Seurat ───────────────────────────────────────────────────
cat("=== Loading Seurat Object ===\n")
cat("File:", rds_seurat, "\n")
seurat_obj <- readRDS(rds_seurat)
cat(sprintf("Cells: %d | Genes: %d\n", ncol(seurat_obj), nrow(seurat_obj)))

if (!"subcelltype" %in% colnames(seurat_obj@meta.data)) {
  print(colnames(seurat_obj@meta.data))
  sink()
  stop("subcelltype column missing")
}

# Seurat v5: merge Assay5 layers, then classic Assay (v4-like) for NeuronChat / CellChat
DefaultAssay(seurat_obj) <- "RNA"
cat("RNA assay class (before):", paste(class(seurat_obj[["RNA"]]), collapse = ", "), "\n")
if (inherits(seurat_obj[["RNA"]], "Assay5")) {
  cat("Assay5 detected — JoinLayers(RNA) then as(., \"Assay\")\n")
  seurat_obj <- tryCatch(
    {
      seurat_obj[["RNA"]] <- JoinLayers(seurat_obj[["RNA"]])
      seurat_obj
    },
    error = function(e) {
      cat("JoinLayers(assay) failed, trying JoinLayers(object):", conditionMessage(e), "\n")
      JoinLayers(seurat_obj, assay = "RNA")
    }
  )
  seurat_obj[["RNA"]] <- as(seurat_obj[["RNA"]], "Assay")
  cat("✓ Converted Assay5 → Assay (v4-like)\n")
} else {
  cat("✓ RNA already Assay (v4-style)\n")
}
cat("RNA assay class (after):", paste(class(seurat_obj[["RNA"]]), collapse = ", "), "\n")

# ── Expression matrix (genes × cells) ─────────────────────────────
.layer_ok <- function(m) !is.null(m) && nrow(m) > 0 && ncol(m) > 0

expr_mat <- tryCatch({
  m <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")
  if (!.layer_ok(m)) stop("data layer empty")
  cat("✓ Using RNA 'data'\n")
  m
}, error = function(e) {
  cat("⚠️ RNA 'data':", conditionMessage(e), "\n")
  NULL
})

if (is.null(expr_mat)) {
  cnt <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
  if (!.layer_ok(cnt)) stop("Cannot normalize: counts missing/empty")
  seurat_obj <- NormalizeData(seurat_obj, assay = "RNA", normalization.method = "LogNormalize",
    scale.factor = 10000, verbose = FALSE)
  seurat_obj <- FindVariableFeatures(seurat_obj, assay = "RNA", selection.method = "vst",
    nfeatures = 2000, verbose = FALSE)
  seurat_obj <- ScaleData(seurat_obj, assay = "RNA", verbose = FALSE)
  expr_mat <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")
  if (!.layer_ok(expr_mat)) stop("Normalization did not produce RNA 'data'")
  cat("✓ Fallback normalization complete\n")
}

# Zhao-style group: names = subcelltype (node in x@idents), values = celltype (coarse class)
plot_group_celltype <- NULL
if ("celltype" %in% colnames(seurat_obj@meta.data)) {
  df_map <- unique(data.frame(
    subcelltype = as.character(seurat_obj$subcelltype),
    celltype    = as.character(seurat_obj$celltype),
    stringsAsFactors = FALSE
  ))
  if (any(duplicated(df_map$subcelltype))) {
    cat("⚠️  subcelltype maps to multiple celltypes — using first row per subcelltype.\n")
    df_map <- df_map[!duplicated(df_map$subcelltype), ]
  }
  plot_group_celltype <- structure(df_map$celltype, names = df_map$subcelltype)
  write.csv(df_map[order(df_map$subcelltype), ],
    file.path(exports_dir, paste0("subcelltype_to_celltype_map_", sample_name, ".csv")),
    row.names = FALSE)
}

cat("\n=== Unique subcelltype (inference / network nodes) ===\n")
u_sub <- sort(unique(as.character(seurat_obj$subcelltype)))
cat("n =", length(u_sub), "\n")
print(u_sub)
writeLines(u_sub, file.path(exports_dir, paste0("unique_subcelltypes_", sample_name, ".txt")))

cat("\n=== Unique celltype (coarse labels for plot group=) ===\n")
if ("celltype" %in% colnames(seurat_obj@meta.data)) {
  u_ct <- sort(unique(as.character(seurat_obj$celltype)))
  cat("n =", length(u_ct), "\n")
  print(u_ct)
  writeLines(u_ct, file.path(exports_dir, paste0("unique_celltypes_", sample_name, ".txt")))
} else {
  cat("(no celltype column — only subcelltype figures, no Zhao coarse group)\n")
}

Idents(seurat_obj) <- "subcelltype"
group_info <- Idents(seurat_obj)

check_memory()

# ── NeuronChat: create + run ────────────────────────────────────────
cat("\n=== createNeuronChat + run_NeuronChat ===\n")
x <- createNeuronChat(object = expr_mat, DB = "human", group.by = as.character(group_info))
x <- run_NeuronChat(x, M = 100)
cat("✓ Pathways in net:", length(x@net), "\n")

interaction_names <- names(x@net)
writeLines(interaction_names, file.path(exports_dir, paste0("all_interaction_names_", sample_name, ".txt")))

pathway_max <- sapply(x@net, function(mat) max(mat, na.rm = TRUE))
active_pathways <- names(x@net)[pathway_max > 0]
cat("Active pathways (any edge > 0):", length(active_pathways), "\n")
writeLines(active_pathways, file.path(exports_dir, paste0("active_pathways_", sample_name, ".txt")))

check_memory()

# ── Aggregated network + core visuals ─────────────────────────────
cat("\n=== net_aggregation + netVisual / heatmap_aggregated ===\n")
net_aggregated_x <- tryCatch(
  net_aggregation(x@net, method = "weight"),
  error = function(e) {
    cat("net_aggregation failed:", conditionMessage(e), "\n")
    NULL
  }
)

if (!is.null(net_aggregated_x)) {
  write.csv(net_aggregated_x, file.path(exports_dir, paste0("aggregated_network_weight_", sample_name, ".csv")),
    row.names = TRUE)
}

tryCatch({
  png(file.path(fig_dir, paste0("netVisual_circle_aggregated_", sample_name, ".png")),
    width = 12, height = 10, units = "in", res = 300)
  print(netVisual_circle_neuron(net_aggregated_x, vertex.label.cex = 1))
  dev.off()
}, error = function(e) {
  if (dev.cur() > 1) dev.off()
  cat("circle aggregated:", conditionMessage(e), "\n")
})

tryCatch({
  png(file.path(fig_dir, paste0("netVisual_chord_aggregated_", sample_name, ".png")),
    width = 12, height = 10, units = "in", res = 300)
  print(netVisual_chord_neuron(x, method = "weight", lab.cex = 1))
  dev.off()
}, error = function(e) {
  if (dev.cur() > 1) dev.off()
  cat("chord aggregated:", conditionMessage(e), "\n")
})

tryCatch({
  png(file.path(fig_dir, paste0("heatmap_aggregated_weight_", sample_name, ".png")),
    width = 12, height = 10, units = "in", res = 300)
  print(heatmap_aggregated(x, method = "weight"))
  dev.off()
}, error = function(e) {
  if (dev.cur() > 1) dev.off()
  cat("heatmap aggregated:", conditionMessage(e), "\n")
})

# Additional aggregated plots: group = subcelltype → celltype (Wei Zhao tutorial style)
if (!is.null(plot_group_celltype) && !is.null(net_aggregated_x)) {
  cat("\n=== Aggregated plots with group = celltype (per subcelltype) ===\n")
  tryCatch({
    png(file.path(fig_dir, paste0("netVisual_circle_aggregated_", sample_name, "_celltype_group.png")),
      width = 12, height = 10, units = "in", res = 300)
    print(netVisual_circle_neuron(net_aggregated_x, group = plot_group_celltype, vertex.label.cex = 1))
    dev.off()
  }, error = function(e) {
    if (dev.cur() > 1) dev.off()
    cat("circle aggregated celltype_group:", conditionMessage(e), "\n")
  })
  tryCatch({
    png(file.path(fig_dir, paste0("netVisual_chord_aggregated_", sample_name, "_celltype_group.png")),
      width = 12, height = 10, units = "in", res = 300)
    print(netVisual_chord_neuron(x, method = "weight", group = plot_group_celltype, lab.cex = 1))
    dev.off()
  }, error = function(e) {
    if (dev.cur() > 1) dev.off()
    cat("chord aggregated celltype_group:", conditionMessage(e), "\n")
  })
  tryCatch({
    png(file.path(fig_dir, paste0("heatmap_aggregated_weight_", sample_name, "_celltype_group.png")),
      width = 12, height = 10, units = "in", res = 300)
    print(heatmap_aggregated(x, method = "weight", group = plot_group_celltype))
    dev.off()
  }, error = function(e) {
    if (dev.cur() > 1) dev.off()
    cat("heatmap aggregated celltype_group:", conditionMessage(e), "\n")
  })
}

check_memory()

# ── Loop: active pathways (circle, chord, heatmaps, lig_tar) ──────
if (length(active_pathways) > 0) {
  for (pw in active_pathways) {
    pw_safe <- safe_name(pw)
    message("Pathway: ", pw)
    png(file.path(fig_dir, paste0(pw_safe, "_circle_chord_", sample_name, ".png")),
      width = 16, height = 8, units = "in", res = 300)
    par(mfrow = c(1, 2), mar = c(3, 3, 3, 3))
    try(netVisual_circle_neuron(x@net[[pw]], vertex.label.cex = 1), silent = TRUE)
    try(netVisual_chord_neuron(x, interaction_use = pw, lab.cex = 1), silent = TRUE)
    dev.off()

    png(file.path(fig_dir, paste0(pw_safe, "_heatmap_count_", sample_name, ".png")),
      width = 10, height = 9, units = "in", res = 300)
    try({
      heatmap_aggregated(x, method = "count", interaction_use = pw)
      title(main = paste0("count: ", pw))
    }, silent = TRUE)
    dev.off()

    png(file.path(fig_dir, paste0(pw_safe, "_heatmap_weight_", sample_name, ".png")),
      width = 10, height = 9, units = "in", res = 300)
    try({
      heatmap_aggregated(x, method = "weight", interaction_use = pw)
      title(main = paste0("weight: ", pw))
    }, silent = TRUE)
    dev.off()

    png(file.path(fig_dir, paste0(pw_safe, "_heatmap_single_", sample_name, ".png")),
      width = 10, height = 9, units = "in", res = 300)
    try({
      heatmap_single(x, interaction_name = pw)
      title(main = paste0("single: ", pw))
    }, silent = TRUE)
    dev.off()

    options(repr.plot.width = 30, repr.plot.height = 8, repr.plot.res = 150)
    png(file.path(fig_dir, paste0(pw_safe, "_ligand_target_heatmap_", sample_name, ".png")),
      width = 30, height = 8, units = "in", res = 300)
    try(lig_tar_heatmap(x, interaction_name = pw), silent = TRUE)
    dev.off()

    if (!is.null(plot_group_celltype)) {
      png(file.path(fig_dir, paste0(pw_safe, "_circle_chord_", sample_name, "_celltype_group.png")),
        width = 16, height = 8, units = "in", res = 300)
      par(mfrow = c(1, 2), mar = c(3, 3, 3, 3))
      try(netVisual_circle_neuron(x@net[[pw]], group = plot_group_celltype, vertex.label.cex = 1), silent = TRUE)
      try(netVisual_chord_neuron(x, interaction_use = pw, group = plot_group_celltype, lab.cex = 1), silent = TRUE)
      dev.off()

      png(file.path(fig_dir, paste0(pw_safe, "_heatmap_count_", sample_name, "_celltype_group.png")),
        width = 10, height = 9, units = "in", res = 300)
      try({
        heatmap_aggregated(x, method = "count", interaction_use = pw, group = plot_group_celltype)
        title(main = paste0("count (celltype group): ", pw))
      }, silent = TRUE)
      dev.off()

      png(file.path(fig_dir, paste0(pw_safe, "_heatmap_weight_", sample_name, "_celltype_group.png")),
        width = 10, height = 9, units = "in", res = 300)
      try({
        heatmap_aggregated(x, method = "weight", interaction_use = pw, group = plot_group_celltype)
        title(main = paste0("weight (celltype group): ", pw))
      }, silent = TRUE)
      dev.off()

      png(file.path(fig_dir, paste0(pw_safe, "_heatmap_single_", sample_name, "_celltype_group.png")),
        width = 10, height = 9, units = "in", res = 300)
      try({
        heatmap_single(x, interaction_name = pw, group = plot_group_celltype)
        title(main = paste0("single (celltype group): ", pw))
      }, silent = TRUE)
      dev.off()
    }
  }
}

check_memory()

# ── rankNet_Neuron ────────────────────────────────────────────────
cat("\n=== rankNet_Neuron ===\n")
tryCatch({
  png(file.path(fig_dir, paste0("rankNet_single_count_", sample_name, ".png")),
    width = 12, height = 8, units = "in", res = 300)
  rankNet_Neuron(x, slot.name = "net", measure = "count", mode = "single")
  dev.off()
}, error = function(e) {
  if (dev.cur() > 1) dev.off()
  cat("rankNet count:", conditionMessage(e), "\n")
})

tryCatch({
  png(file.path(fig_dir, paste0("rankNet_single_weight_", sample_name, ".png")),
    width = 12, height = 8, units = "in", res = 300)
  rankNet_Neuron(x, slot.name = "net", measure = "weight", mode = "single")
  dev.off()
}, error = function(e) {
  if (dev.cur() > 1) dev.off()
  cat("rankNet weight:", conditionMessage(e), "\n")
})

tryCatch({
  png(file.path(fig_dir, paste0("rankNet_count_vs_weight_", sample_name, ".png")),
    width = 16, height = 8, units = "in", res = 300)
  par(mfrow = c(1, 2))
  rankNet_Neuron(x, slot.name = "net", measure = "count", mode = "single")
  rankNet_Neuron(x, slot.name = "net", measure = "weight", mode = "single")
  dev.off()
}, error = function(e) {
  if (dev.cur() > 1) dev.off()
  cat("rankNet combo:", conditionMessage(e), "\n")
})

check_memory()

# ── Communication patterns (updates x@net_analysis) ───────────────
cat("\n=== identifyCommunicationPatterns_Neuron ===\n")
tryCatch({
  png(file.path(fig_dir, paste0("pattern_outgoing_", sample_name, ".png")),
    width = 7, height = 6, units = "in", res = 300)
  x <- identifyCommunicationPatterns_Neuron(x, slot.name = "net", pattern = "outgoing", k = 4, height = 18)
  dev.off()
}, error = function(e) {
  if (dev.cur() > 1) dev.off()
  cat("pattern outgoing:", conditionMessage(e), "\n")
})

tryCatch({
  png(file.path(fig_dir, paste0("pattern_incoming_", sample_name, ".png")),
    width = 7, height = 6, units = "in", res = 300)
  x <- identifyCommunicationPatterns_Neuron(x, slot.name = "net", pattern = "incoming", k = 4, height = 18)
  dev.off()
}, error = function(e) {
  if (dev.cur() > 1) dev.off()
  cat("pattern incoming:", conditionMessage(e), "\n")
})

check_memory()

# ── Pattern matrices: write.table + pheatmap ──────────────────────
get_incoming_data <- function(obj) {
  if (!is.null(obj@net_analysis$pattern$incoming$data)) return(obj@net_analysis$pattern$incoming$data)
  stop("incoming pattern$data not found")
}
get_outgoing_data <- function(obj) {
  if (!is.null(obj@net_analysis$pattern$outgoing$data)) return(obj@net_analysis$pattern$outgoing$data)
  stop("outgoing pattern$data not found")
}
get_incoming_cell_df <- function(obj) {
  if (!is.null(obj@net_analysis$pattern$incoming$pattern$cell)) return(obj@net_analysis$pattern$incoming$pattern$cell)
  if (!is.null(obj@net_analysis$cell)) return(obj@net_analysis$cell)
  stop("incoming cell df not found")
}
get_outgoing_cell_df <- function(obj) {
  if (!is.null(obj@net_analysis$pattern$outgoing$pattern$cell)) return(obj@net_analysis$pattern$outgoing$pattern$cell)
  if (!is.null(obj@net_analysis$cell)) return(obj@net_analysis$cell)
  stop("outgoing cell df not found")
}
get_incoming_signaling_df <- function(obj) {
  if (!is.null(obj@net_analysis$pattern$incoming$pattern$signaling)) return(obj@net_analysis$pattern$incoming$pattern$signaling)
  if (!is.null(obj@net_analysis$signaling)) return(obj@net_analysis$signaling)
  stop("incoming signaling df not found")
}
get_outgoing_signaling_df <- function(obj) {
  if (!is.null(obj@net_analysis$pattern$outgoing$pattern$signaling)) return(obj@net_analysis$pattern$outgoing$pattern$signaling)
  if (!is.null(obj@net_analysis$signaling)) return(obj@net_analysis$signaling)
  stop("outgoing signaling df not found")
}

cat("\n=== Pattern exports (incoming / outgoing) ===\n")
tryCatch({
  mat_in <- get_incoming_data(x)
  write.table(mat_in, file.path(exports_dir, paste0("net_analysis_incoming_data_", sample_name, ".txt")),
    sep = "\t", quote = FALSE, col.names = NA)
  png(file.path(fig_dir, paste0("net_analysis_incoming_data_heatmap_", sample_name, ".png")),
    width = 12, height = 8, units = "in", res = 300)
  pheatmap(mat_in, cluster_rows = FALSE, cluster_cols = FALSE, main = "Incoming pattern$data")
  dev.off()
}, error = function(e) cat("incoming data export:", conditionMessage(e), "\n"))

tryCatch({
  mat_out <- get_outgoing_data(x)
  write.table(mat_out, file.path(exports_dir, paste0("net_analysis_outgoing_data_", sample_name, ".txt")),
    sep = "\t", quote = FALSE, col.names = NA)
  png(file.path(fig_dir, paste0("net_analysis_outgoing_data_heatmap_", sample_name, ".png")),
    width = 12, height = 8, units = "in", res = 300)
  pheatmap(mat_out, cluster_rows = FALSE, cluster_cols = FALSE, main = "Outgoing pattern$data")
  dev.off()
}, error = function(e) cat("outgoing data export:", conditionMessage(e), "\n"))

tryCatch({
  cell_df_in <- get_incoming_cell_df(x)
  write.table(cell_df_in, file.path(exports_dir, paste0("net_analysis_incoming_cell_", sample_name, ".txt")),
    sep = "\t", quote = FALSE, row.names = FALSE)
  cell_mat_in <- cell_df_in %>%
    mutate(CellGroup = as.character(CellGroup), Pattern = as.character(Pattern)) %>%
    select(CellGroup, Pattern, Contribution) %>%
    pivot_wider(names_from = Pattern, values_from = Contribution, values_fill = 0) %>%
    column_to_rownames("CellGroup") %>% as.matrix()
  write.table(cell_mat_in, file.path(exports_dir, paste0("net_analysis_incoming_cell_matrix_", sample_name, ".txt")),
    sep = "\t", quote = FALSE, col.names = NA)
  png(file.path(fig_dir, paste0("net_analysis_incoming_cell_heatmap_", sample_name, ".png")),
    width = 7, height = 6, units = "in", res = 300)
  pheatmap(cell_mat_in, cluster_rows = FALSE, cluster_cols = FALSE, main = "Incoming cell contributions")
  dev.off()
}, error = function(e) cat("incoming cell export:", conditionMessage(e), "\n"))

tryCatch({
  cell_df_out <- get_outgoing_cell_df(x)
  write.table(cell_df_out, file.path(exports_dir, paste0("net_analysis_outgoing_cell_", sample_name, ".txt")),
    sep = "\t", quote = FALSE, row.names = FALSE)
  cell_mat_out <- cell_df_out %>%
    mutate(CellGroup = as.character(CellGroup), Pattern = as.character(Pattern)) %>%
    select(CellGroup, Pattern, Contribution) %>%
    pivot_wider(names_from = Pattern, values_from = Contribution, values_fill = 0) %>%
    column_to_rownames("CellGroup") %>% as.matrix()
  write.table(cell_mat_out, file.path(exports_dir, paste0("net_analysis_outgoing_cell_matrix_", sample_name, ".txt")),
    sep = "\t", quote = FALSE, col.names = NA)
  png(file.path(fig_dir, paste0("net_analysis_outgoing_cell_heatmap_", sample_name, ".png")),
    width = 7, height = 6, units = "in", res = 300)
  pheatmap(cell_mat_out, cluster_rows = FALSE, cluster_cols = FALSE, main = "Outgoing cell contributions")
  dev.off()
}, error = function(e) cat("outgoing cell export:", conditionMessage(e), "\n"))

tryCatch({
  sig_df_in <- get_incoming_signaling_df(x)
  write.table(sig_df_in, file.path(exports_dir, paste0("net_analysis_incoming_signaling_", sample_name, ".txt")),
    sep = "\t", quote = FALSE, row.names = FALSE)
  sig_mat_in <- sig_df_in %>%
    mutate(Signaling = as.character(Signaling), Pattern = as.character(Pattern)) %>%
    select(Signaling, Pattern, Contribution) %>%
    pivot_wider(names_from = Pattern, values_from = Contribution, values_fill = 0) %>%
    column_to_rownames("Signaling") %>% as.matrix()
  topN <- min(100L, nrow(sig_mat_in))
  ord <- order(apply(sig_mat_in, 1, max), decreasing = TRUE)
  sig_top <- sig_mat_in[ord[seq_len(topN)], , drop = FALSE]
  png(file.path(fig_dir, paste0("net_analysis_incoming_signaling_heatmap_top_", sample_name, ".png")),
    width = 10, height = 10, units = "in", res = 300)
  pheatmap(sig_top, cluster_rows = FALSE, cluster_cols = FALSE, main = "Incoming signaling (top)")
  dev.off()
}, error = function(e) cat("incoming signaling export:", conditionMessage(e), "\n"))

tryCatch({
  sig_df_out <- get_outgoing_signaling_df(x)
  write.table(sig_df_out, file.path(exports_dir, paste0("net_analysis_outgoing_signaling_", sample_name, ".txt")),
    sep = "\t", quote = FALSE, row.names = FALSE)
  sig_mat_out <- sig_df_out %>%
    mutate(Signaling = as.character(Signaling), Pattern = as.character(Pattern)) %>%
    select(Signaling, Pattern, Contribution) %>%
    pivot_wider(names_from = Pattern, values_from = Contribution, values_fill = 0) %>%
    column_to_rownames("Signaling") %>% as.matrix()
  topN <- min(100L, nrow(sig_mat_out))
  ord <- order(apply(sig_mat_out, 1, max), decreasing = TRUE)
  sig_top <- sig_mat_out[ord[seq_len(topN)], , drop = FALSE]
  png(file.path(fig_dir, paste0("net_analysis_outgoing_signaling_heatmap_top_", sample_name, ".png")),
    width = 10, height = 10, units = "in", res = 300)
  pheatmap(sig_top, cluster_rows = FALSE, cluster_cols = FALSE, main = "Outgoing signaling (top)")
  dev.off()
}, error = function(e) cat("outgoing signaling export:", conditionMessage(e), "\n"))

check_memory()

# ── River plots ───────────────────────────────────────────────────
cat("\n=== netAnalysis_river_Neuron ===\n")
tryCatch({
  png(file.path(fig_dir, paste0("netAnalysis_river_outgoing_", sample_name, ".png")),
    width = 10, height = 8, units = "in", res = 300)
  netAnalysis_river_Neuron(x, slot.name = "net", pattern = "outgoing", font.size = 2.5, cutoff.1 = 0.5, cutoff.2 = 0.5)
  dev.off()
}, error = function(e) {
  if (dev.cur() > 1) dev.off()
  cat("river outgoing:", conditionMessage(e), "\n")
})

tryCatch({
  png(file.path(fig_dir, paste0("netAnalysis_river_incoming_", sample_name, ".png")),
    width = 10, height = 8, units = "in", res = 300)
  netAnalysis_river_Neuron(x, slot.name = "net", pattern = "incoming", font.size = 2.5, cutoff.1 = 0.5, cutoff.2 = 0.5)
  dev.off()
}, error = function(e) {
  if (dev.cur() > 1) dev.off()
  cat("river incoming:", conditionMessage(e), "\n")
})

# Summary pattern matrices to exports
tryCatch({
  neta <- x@net_analysis
  if (!is.null(neta$pattern$outgoing$data)) {
    write.table(neta$pattern$outgoing$data,
      file.path(exports_dir, paste0("net_analysis_pattern_outgoing_matrix_", sample_name, ".txt")),
      sep = "\t", quote = FALSE, col.names = NA)
  }
  if (!is.null(neta$pattern$incoming$data)) {
    write.table(neta$pattern$incoming$data,
      file.path(exports_dir, paste0("net_analysis_pattern_incoming_matrix_", sample_name, ".txt")),
      sep = "\t", quote = FALSE, col.names = NA)
  }
}, error = function(e) cat("summary pattern write:", conditionMessage(e), "\n"))

check_memory()

# ── Functional similarity + embedding ────────────────────────────
cat("\n=== Functional manifold (uwot) ===\n")
tryCatch({
  x <- computeNetSimilarity_Neuron(x, slot.name = "net", type = "functional")
  x <- netEmbedding(x, slot.name = "net_analysis", type = "functional", umap.method = "uwot")
  x <- netClustering(x, slot.name = "net_analysis", type = "functional", k = 5)
  png(file.path(fig_dir, paste0("netVisual_embedding_functional_", sample_name, ".png")),
    width = 7, height = 6, units = "in", res = 300)
  netVisual_embedding_Neuron(x, type = "functional", label.size = 3.5, pathway.remove.show = FALSE)
  dev.off()
  png(file.path(fig_dir, paste0("netVisual_embeddingZoomIn_functional_", sample_name, ".png")),
    width = 12, height = 12, units = "in", res = 300)
  netVisual_embeddingZoomIn_Neuron(x, type = "functional", nCol = 2, label.size = 3)
  dev.off()
}, error = function(e) {
  if (dev.cur() > 1) dev.off()
  cat("functional manifold:", conditionMessage(e), "\n")
})

# ── Structural similarity + embedding ─────────────────────────────
cat("\n=== Structural manifold (uwot) ===\n")
tryCatch({
  x <- computeNetSimilarity_Neuron(x, slot.name = "net", type = "structural")
  x <- netEmbedding(x, slot.name = "net_analysis", type = "structural", umap.method = "uwot")
  x <- netClustering(x, slot.name = "net_analysis", type = "structural", k = 5)
  png(file.path(fig_dir, paste0("netVisual_embedding_structural_", sample_name, ".png")),
    width = 7, height = 6, units = "in", res = 300)
  netVisual_embedding_Neuron(x, type = "structural", label.size = 3.5, pathway.remove.show = FALSE)
  dev.off()
  png(file.path(fig_dir, paste0("netVisual_embeddingZoomIn_structural_", sample_name, ".png")),
    width = 12, height = 12, units = "in", res = 300)
  netVisual_embeddingZoomIn_Neuron(x, type = "structural", nCol = 2, label.size = 3)
  dev.off()
}, error = function(e) {
  if (dev.cur() > 1) dev.off()
  cat("structural manifold:", conditionMessage(e), "\n")
})

check_memory()

# ── lig_tar_heatmap: all interactions (optional long loop) ────────
cat("\n=== lig_tar_heatmap (all interactions) ===\n")
lig_tar_dir <- file.path(fig_dir, paste0("lig_tar_all_interactions_", sample_name))
dir.create(lig_tar_dir, showWarnings = FALSE, recursive = TRUE)
for (interaction_example in interaction_names) {
  if (!interaction_example %in% names(x@net)) next
  im <- x@net[[interaction_example]]
  vals <- as.vector(im)
  vals <- vals[is.finite(vals)]
  if (length(vals) < 2L || length(unique(vals)) < 2L) next
  out_file <- file.path(lig_tar_dir, paste0(safe_name(interaction_example), "_lig_tar_heatmap.png"))
  png(out_file, width = 12, height = 8, units = "in", res = 300)
  tryCatch(
    lig_tar_heatmap(x, interaction_name = interaction_example, width.vector = c(0.38, 0.35, 0.27)),
    error = function(e) cat("lig_tar", interaction_example, ":", conditionMessage(e), "\n")
  )
  dev.off()
}

check_memory()

# ── Final RDS + manifest ──────────────────────────────────────────
cat("\n=== Saving RDS ===\n")
rds_out <- file.path(out_dir, "neuronchat_18_merged_AG_031826.28celltypes.AD_at_the_end_v5.rds")
saveRDS(x, file = rds_out)
cat("✓ RDS:", rds_out, "\n")

writeLines(c(
  paste("sample:", sample_name),
  paste("rds_out:", rds_out),
  paste("fig_dir:", fig_dir),
  paste("exports_dir:", exports_dir),
  paste("ended:", format(Sys.time()))
), file.path(exports_dir, paste0("run_manifest_", sample_name, ".txt")))

cat(sprintf("Cell types in object: %d\n", length(levels(x@idents))))
print(sort(levels(x@idents)))

cat("\n╔═══════════════════════════════════════════════════════════╗\n")
cat("║  ✓ NeuronChat AD v5 COMPLETE                              ║\n")
cat("╚═══════════════════════════════════════════════════════════╝\n")

sink()
