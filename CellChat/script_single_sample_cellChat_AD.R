# ══════════════════════════════════════════════════════════════════
# CELLCHAT ANALYSIS
# Input : 18_merged_AG_031826.26celltypes.min50cells.18K.proportional.anndataR.AD_Dyslexia.rds
# Groups: AD vs Dyslexia
# ══════════════════════════════════════════════════════════════════

# ── Libraries ────────────────────────────────────────────────────
library(NMF)
library(CellChat)
library(ggalluvial)
library(patchwork)
library(Seurat)
library(future)
library(future.apply)
library(parallelly)
library(progress)
library(BiocParallel)
library(Matrix)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(gridExtra)
library(data.table)
library(cowplot)
library(foreach)
library(doParallel)
library(ggbeeswarm)
library(ComplexHeatmap)
library(circlize)
library(uwot)
options(stringsAsFactors = FALSE)

# ══════════════════════════════════════════════════════════════════
# PERFORMANCE TUNING
# ══════════════════════════════════════════════════════════════════

n_cores_total   <- parallelly::availableCores()
n_workers       <- min(n_cores_total - 1, 60)  # computeCommunProb
n_workers_light <- min(n_cores_total - 1, 8)   # lighter steps

ram_limit_gb     <- 200   # future globals max size (GB)
mem_threshold    <- 90    # stop if RAM usage exceeds this %
nboot_iterations <- 100    # 100 for publication, 50 for exploration
k_patterns       <- c(2, 3, 4, 5)  # only test k=2 and k=3
save_k_rds       <- FALSE    # TRUE = save all k, FALSE = save only best k

cat(sprintf("
══════════════════════════════════════════
  Performance settings
  Total cores    : %d
  Heavy workers  : %d  (computeCommunProb)
  Light workers  : %d  (other steps)
  RAM limit      : %d GB
  Memory stop at : %d%%
  nboot          : %d
  k patterns     : %s
══════════════════════════════════════════\n",
  n_cores_total, n_workers, n_workers_light,
  ram_limit_gb, mem_threshold, nboot_iterations,
  paste(k_patterns, collapse = ", ")))

# ── Parallelization ───────────────────────────────────────────────
if (.Platform$OS.type == "unix" &&
    !grepl("darwin", Sys.info()["sysname"], ignore.case = TRUE)) {
  BPPARAM <- MulticoreParam(workers = n_workers, RNGseed = 1337)
  cat("✓ BiocParallel: MulticoreParam (fork-based, Linux)\n")
} else {
  BPPARAM <- SnowParam(workers = n_workers, RNGseed = 1337, type = "SOCK")
  cat("✓ BiocParallel: SnowParam (socket-based)\n")
}
register(BPPARAM)
cat(sprintf("✓ BiocParallel registered: %d workers\n", n_workers))

plan(multicore, workers = n_workers)
cat(sprintf("✓ future::plan(multicore, workers = %d)\n", n_workers))

# ── Global options ────────────────────────────────────────────────
Sys.setenv(TMPDIR = "/mnt/nfs/CX000008_DS1/projects/btanasa/tmp")
options(future.globals.maxSize = ram_limit_gb * 1024^3)
options(future.seed = TRUE)
options(bitmapType = "cairo")
options(repr.plot.width = 8, repr.plot.height = 6, repr.plot.res = 160)
set.seed(1337)

if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(1)
  RhpcBLASctl::omp_set_num_threads(1)
  cat("✓ BLAS/OMP threads set to 1\n")
}

# ── ggplot theme ──────────────────────────────────────────────────
theme_set(theme_classic(base_size = 14))
update_geom_defaults("point", list(size = 1.2, alpha = 0.8))
update_geom_defaults("text",  list(size = 4))

# ══════════════════════════════════════════════════════════════════
# HELPER FUNCTIONS
# ══════════════════════════════════════════════════════════════════

is_memory_error <- function(e) {
  grepl("MEMORY THRESHOLD EXCEEDED", conditionMessage(e))
}

check_memory_usage <- function(threshold = mem_threshold, force_gc = TRUE) {
  if (force_gc) gc(verbose = FALSE, full = TRUE)
  if (file.exists("/proc/meminfo")) {
    meminfo       <- readLines("/proc/meminfo")
    mem_total     <- as.numeric(sub(".*:\\s+([0-9]+).*", "\\1",
                                    grep("^MemTotal:",     meminfo, value = TRUE)))
    mem_available <- as.numeric(sub(".*:\\s+([0-9]+).*", "\\1",
                                    grep("^MemAvailable:", meminfo, value = TRUE)))
    mem_used      <- mem_total - mem_available
    mem_pct       <- (mem_used / mem_total) * 100
    cat(sprintf("💾 Memory: %.1f%% used (%.2f GB / %.2f GB)\n",
                mem_pct, mem_used/1024^2, mem_total/1024^2))
    if (mem_pct > threshold)
      stop(sprintf(
        "❌ MEMORY THRESHOLD EXCEEDED: %.1f%% > %d%% — stopping to prevent crash",
        mem_pct, threshold))
    return(invisible(mem_pct))
  }
  invisible(NA)
}

do_gc <- function(verbose = TRUE) {
  if (verbose) cat("🧹 Running garbage collection...\n")
  gc(verbose = FALSE, full = TRUE)
  invisible(NULL)
}

safe_dev_off <- function() {
  tryCatch({
    while (grDevices::dev.cur() > 1) grDevices::dev.off()
  }, error = function(e) invisible(NULL))
}

is_gg <- function(x) inherits(x, "ggplot")

cat("\n=== Initial Memory ===\n")
check_memory_usage()

# ══════════════════════════════════════════════════════════════════
# KEY IDENTIFIERS
# ══════════════════════════════════════════════════════════════════

base_dir <- "/mnt/nfs/CX000008_DS1/projects/jaeyeon/fastq_file_Dyslexia_r1/zanalysis_bogdan/samples.merged_AG_Harmony_res0.x_anno_031826"

setwd(base_dir)

rds_filename <- "18_merged_AG_031826.28celltypes.AD.rds"
sample_name  <- sub("\\.rds$", "", rds_filename)
assay        <- "RNA"

cat(sprintf("
══════════════════════════════════════════
  Sample    : %s
  Cell type : subcelltype
  Assay     : %s
══════════════════════════════════════════\n",
  sample_name, assay))

# ── Output directory ──────────────────────────────────────────────
plot_dir <- file.path(base_dir, paste0(sample_name, "_cellchat"))
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
cat("✓ Output directory:", plot_dir, "\n")

# ── Logging ───────────────────────────────────────────────────────
log_file   <- file.path(plot_dir, paste0("cellchat_",
                                          format(Sys.time(), "%Y%m%d_%H%M%S"),
                                          ".log"))
start_time <- Sys.time()
sink(log_file, split = TRUE)

cat("╔═══════════════════════════════════════════════════════════╗\n")
cat("║  CELLCHAT ANALYSIS LOG                                    ║\n")
cat("╚═══════════════════════════════════════════════════════════╝\n")
cat("Started  :", format(start_time, "%Y-%m-%d %H:%M:%S"), "\n")
cat("Sample   :", sample_name, "\n")
cat("Log file :", log_file, "\n\n")

on.exit({
  if (sink.number() > 0) {
    cat("\nEnded  :", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
    cat("Runtime:", round(difftime(Sys.time(), start_time, units = "mins"), 2),
        "minutes\n")
    sink()
  }
}, add = TRUE)

# ══════════════════════════════════════════════════════════════════
# LOAD DATA
# ══════════════════════════════════════════════════════════════════

cat("\n=== Loading Data ===\n")
cat("File:", rds_filename, "\n")
data_human <- readRDS(rds_filename)
check_memory_usage()

cat("\nObject overview:\n")
print(data_human)

cat("\nGroups present:\n")
print(table(data_human@meta.data$group))

cat("\nCell types present (subcelltype):\n")
print(table(data_human@meta.data$subcelltype))

# ── Set assay and filter Ensembl IDs ─────────────────────────────
DefaultAssay(data_human) <- assay

genes      <- rownames(data_human[[assay]])
is_ensembl <- grepl("^ENS[A-Z]*G\\d+(\\.\\d+)?$", genes)
keep       <- genes[!is_ensembl]

cat(sprintf("\nRemoving %d Ensembl IDs, keeping %d gene symbols\n",
            sum(is_ensembl), length(keep)))

data_human <- subset(data_human, features = keep)

# ── Prepare expression matrix ─────────────────────────────────────
# Try 'data' layer first, fall back to 'normalized' (anndataR naming)
data.input <- tryCatch(
  LayerData(data_human, assay = assay, layer = "data"),
  error = function(e) NULL
)

if (is.null(data.input) || nrow(data.input) == 0) {
  cat("⚠️  'data' layer empty — using 'normalized' layer\n")
  data.input <- LayerData(data_human, assay = assay, layer = "normalized")
}

cat(sprintf("Expression matrix: %d genes x %d cells\n",
            nrow(data.input), ncol(data.input)))

# ── Metadata ──────────────────────────────────────────────────────
# FIX 1: hardcode subcelltype as labels
#         add samples column (silences CellChat v2 warning)
meta         <- data_human@meta.data
meta$labels  <- meta$subcelltype   # hardcoded — cell type identity
meta$samples <- meta$group         # AD or Dyslexia — silences CellChat warning

cat("\nCell type labels (subcelltype):\n")
print(sort(unique(meta$labels)))

cat("\nSamples (group):\n")
print(table(meta$samples))

# ── Release data_human — no longer needed ─────────────────────────
cat("\n🗑️  Releasing data_human from memory...\n")
rm(data_human, genes, keep, is_ensembl)
do_gc()
check_memory_usage()

# ══════════════════════════════════════════════════════════════════
# CREATE CELLCHAT OBJECT
# ══════════════════════════════════════════════════════════════════

cat("\n=== Creating CellChat Object ===\n")
cellchat <- createCellChat(
  object   = data.input,
  meta     = meta,
  group.by = "labels"
)

celltypes   <- sort(unique(cellchat@idents))
n_celltypes <- length(celltypes)
cat(sprintf("✓ %d cell types: %s\n\n",
            n_celltypes, paste(celltypes, collapse = ", ")))

# ── Release data.input and meta — no longer needed ────────────────
cat("🗑️  Releasing data.input and meta from memory...\n")
rm(data.input, meta)
do_gc()
check_memory_usage()

# ══════════════════════════════════════════════════════════════════
# DATABASE
# ══════════════════════════════════════════════════════════════════

cat("\n=== Setting CellChatDB ===\n")
CellChatDB  <- CellChatDB.human
cellchat@DB <- CellChatDB

cat("Interaction categories:\n")
print(unique(CellChatDB$interaction$annotation))

rm(CellChatDB)
do_gc()

cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

cat(sprintf("✓ Signaling matrix: %d genes x %d cells\n",
            nrow(cellchat@data.signaling),
            ncol(cellchat@data.signaling)))

# FIX 2: projectData → safe fallback for CellChat v2
# projectData is optional — smooths expression over PPI network
# to detect more interactions from lowly expressed ligands/receptors
cat("\n=== PPI Projection (optional) ===\n")
tryCatch({
  cellchat <- projectData(cellchat, adjMatrix = PPI.human)
  cat("✓ projectData done\n")
}, error = function(e) {
  if (is_memory_error(e)) stop(e)  # re-throw memory errors
  cat("⚠️  projectData not found — trying smoothData...\n")
  tryCatch({
    cellchat <- smoothData(cellchat, adj = PPI.human)
    cat("✓ smoothData done\n")
  }, error = function(e2) {
    if (is_memory_error(e2)) stop(e2)  # re-throw memory errors
    cat("⚠️  PPI projection skipped (not available in this CellChat version)\n")
    cat("   This is optional — analysis continues without it\n")
    cat("   Effect: fewer interactions detected for lowly expressed L-R pairs\n")
  })
})

check_memory_usage()

# ══════════════════════════════════════════════════════════════════
# COMPUTE COMMUNICATION PROBABILITIES
# ══════════════════════════════════════════════════════════════════

cat("\n=== Computing Communication Probabilities ===\n")
cat(sprintf("BiocParallel workers : %d\n", bpnworkers(BPPARAM)))
cat(sprintf("Bootstrap iterations : %d\n", nboot_iterations))
cat("Setting future::plan(sequential)\n")
cat("Reason: CellChat uses BiocParallel internally —\n")
cat("        mixing future multicore + BiocParallel causes deadlocks\n\n")

future::plan(sequential)

do_gc()
check_memory_usage()

computation_success <- FALSE
attempt             <- 1

while (!computation_success && attempt <= 3) {
  cat(sprintf("\n--- Attempt %d/3 ---\n", attempt))

  tryCatch({

    if (attempt == 1) {
      cat(sprintf("BiocParallel MulticoreParam (%d workers)...\n", n_workers))
      cellchat <- computeCommunProb(
        object  = cellchat,
        type    = "triMean",
        raw.use = TRUE,
        nboot   = nboot_iterations
      )

    } else if (attempt == 2) {
      cat("Fallback: doParallel...\n")
      cl <- makeCluster(n_workers)
      registerDoParallel(cl)
      cellchat <- computeCommunProb(
        object  = cellchat,
        type    = "triMean",
        raw.use = TRUE,
        nboot   = nboot_iterations
      )
      stopCluster(cl)
      rm(cl)

    } else {
      cat("Fallback: sequential (slow but stable)...\n")
      register(SerialParam())
      cellchat <- computeCommunProb(
        object  = cellchat,
        type    = "triMean",
        raw.use = TRUE,
        nboot   = nboot_iterations
      )
    }

    computation_success <- TRUE
    cat("✓ computeCommunProb successful!\n")

  }, error = function(e) {
    if (is_memory_error(e)) stop(e)  # re-throw memory errors
    cat("❌ Attempt", attempt, "failed:", as.character(e), "\n")
    do_gc()
    Sys.sleep(2)
  })

  attempt <- attempt + 1
}

# Restore multicore for subsequent steps
plan(multicore, workers = n_workers)
cat(sprintf("✓ Restored future::plan(multicore, workers = %d)\n", n_workers))

# ── Save after computeCommunProb ──────────────────────────────────
saveRDS(cellchat,
        file.path(plot_dir,
                  paste0("cellchat_", sample_name,
                         "_after_computeCommunProb.rds")))
cat("✓ Saved after computeCommunProb\n")

do_gc()
check_memory_usage()

# ══════════════════════════════════════════════════════════════════
# FILTER, PATHWAY PROB, AGGREGATE, CENTRALITY
# ══════════════════════════════════════════════════════════════════

cat("\n=== Filter + Aggregate ===\n")
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat, thresh = 0.05)
cellchat <- aggregateNet(cellchat, remove.isolate = TRUE, thresh = 0.05)
cellchat <- netAnalysis_computeCentrality(cellchat)

cat("\nCommunication count matrix:\n")
print(cellchat@net$count)

do_gc()
check_memory_usage()

# ── Usable pathways ───────────────────────────────────────────────
has_signal <- function(p) {
  df <- subsetCommunication(cellchat, slot.name = "netP",
                             signaling = p, thresh = 0.05)
  nrow(df) > 0
}
usable <- Filter(has_signal, cellchat@netP$pathways)
cat(sprintf("\n✓ %d / %d pathways have signal\n",
            length(usable), length(cellchat@netP$pathways)))

writeLines(cellchat@netP$pathways,
           file.path(plot_dir, "active_pathways.txt"))

# ══════════════════════════════════════════════════════════════════
# VISUALIZATIONS
# ══════════════════════════════════════════════════════════════════

# ── 1. Signaling role heatmaps ────────────────────────────────────
cat("\n=== 1. Signaling Role Heatmaps ===\n")
for (pattern in c("incoming", "outgoing")) {
  png(file.path(plot_dir,
                paste0("signaling_role_heatmap_", pattern, ".png")),
      width = 8, height = 8, units = "in", res = 300)
  netAnalysis_signalingRole_heatmap(cellchat, pattern = pattern)
  dev.off()
  cat("✓", pattern, "heatmap saved\n")
}
do_gc(verbose = FALSE)

# ── 2. Scatter plot ───────────────────────────────────────────────
cat("\n=== 2. Scatter Plot ===\n")
p_scatter <- netAnalysis_signalingRole_scatter(
  cellchat,
  x.measure = "outdeg",
  y.measure  = "indeg",
  xlabel     = "Outgoing interaction strength",
  ylabel     = "Incoming interaction strength",
  signaling  = cellchat@netP$pathways
)
ggsave(file.path(plot_dir, "signalingRole_scatter_all_pathways.png"),
       p_scatter, width = 7, height = 7, dpi = 300, bg = "white")
ggsave(file.path(plot_dir, "signalingRole_scatter_all_pathways.pdf"),
       p_scatter, width = 7, height = 7)
rm(p_scatter); do_gc(verbose = FALSE)
cat("✓ Scatter plot saved\n")

# ── 3. Signaling role networks per pathway ────────────────────────
cat("\n=== 3. Signaling Role Networks per Pathway ===\n")
pb <- progress_bar$new(
  format = "[:bar] :current/:total (:percent) ETA: :eta | :pathway",
  total  = length(usable), clear = FALSE, width = 80)

for (p in usable) {
  pb$tick(tokens = list(pathway = sprintf("%-15s", p)))
  png(file.path(plot_dir, paste0(p, "_signalingRole_network.png")),
      width = 15, height = 8, units = "in", res = 300)
  try(netAnalysis_signalingRole_network(
    cellchat, signaling = p, width = 15, height = 8, font.size = 15),
    silent = TRUE)
  dev.off()
}
do_gc(verbose = FALSE)
check_memory_usage()

# ── 4. selectK + identifyCommunicationPatterns k=2,3 ─────────────
for (pattern in c("outgoing", "incoming")) {
  cat(sprintf("\n=== 4. selectK (%s) ===\n", pattern))

  png(file.path(plot_dir, paste0("selectK_", pattern, "_patterns.png")),
      width = 14, height = 6, units = "in", res = 300)
  set.seed(12345)
  select_k <- selectK(cellchat, pattern = pattern)
  dev.off()

  best_k <- NULL
  if (!is.null(select_k) && !is.null(select_k$data) &&
      nrow(select_k$data) > 0) {
    write.csv(select_k$data,
              file.path(plot_dir,
                        paste0("selectK_", pattern, "_scores.csv")),
              row.names = FALSE)
    coph <- select_k$data[select_k$data$Measure == "Cophenetic" &
                            select_k$data$k > 1, ]
    if (nrow(coph) > 0) {
      best_k <- coph$k[which.max(coph$score)]
      cat(sprintf("★ Recommended k (%s): %d (Cophenetic = %.4f)\n",
                  pattern, best_k, max(coph$score)))
    }
    rm(coph)
  }
  rm(select_k)
  do_gc(verbose = FALSE)

  cat(sprintf("\n=== identifyCommunicationPatterns k=%s (%s) ===\n",
              paste(k_patterns, collapse = ","), pattern))

  for (k_val in k_patterns) {
    cat(sprintf("\nk = %d...\n", k_val))
    tryCatch({
      cellchat <- identifyCommunicationPatterns(
        cellchat, pattern = pattern, k = k_val)

      if (save_k_rds || (!is.null(best_k) && k_val == best_k)) {
        saveRDS(cellchat,
                file.path(plot_dir,
                          paste0("cellchat_", pattern,
                                 "_k", k_val, ".rds")))
        cat(sprintf("  ✓ RDS saved for k = %d\n", k_val))
      }

      tryCatch({
        p_river <- netAnalysis_river(cellchat, pattern = pattern)
        ggsave(file.path(plot_dir,
                         paste0("netAnalysis_river_", pattern,
                                "_k", k_val, ".png")),
               p_river, width = 10, height = 8, dpi = 300, bg = "white")
        rm(p_river)
      }, error = function(e) {
        if (is_memory_error(e)) stop(e)
        safe_dev_off()
      })

      tryCatch({
        p_dot <- netAnalysis_dot(cellchat, pattern = pattern,
                                  dot.size = c(1, 12))
        ggsave(file.path(plot_dir,
                         paste0("netAnalysis_dot_", pattern,
                                "_k", k_val, ".png")),
               p_dot, width = 14, height = 6, dpi = 300, bg = "white")
        rm(p_dot)
      }, error = function(e) {
        if (is_memory_error(e)) stop(e)
        safe_dev_off()
      })

      cat(sprintf("  ✓ k = %d done%s\n", k_val,
                  ifelse(!is.null(best_k) && k_val == best_k,
                         " ★ RECOMMENDED", "")))

    }, error = function(e) {
      if (is_memory_error(e)) stop(e)
      cat(sprintf("  ❌ k = %d failed: %s\n", k_val, as.character(e)))
    })

    do_gc(verbose = FALSE)
    check_memory_usage()
  }
}

# ── 5. Functional + Structural similarity (uwot, pure R) ──────────
cat("\n=== 5. Functional + Structural Similarity (uwot) ===\n")

for (sim_type in c("functional", "structural")) {
  cat(sprintf("\n--- %s ---\n", sim_type))

  cellchat <- computeNetSimilarity(cellchat, type = sim_type)
  check_memory_usage()

  tryCatch({
    cellchat <- netEmbedding(cellchat, type = sim_type,
                              umap.method = "uwot")
    cat("✓ Embedding done\n")

    cellchat <- netClustering(cellchat, type = sim_type,
                               do.parallel = FALSE)
    cat("✓ Clustering done\n")

    p_embed <- netVisual_embedding(
      cellchat, type = sim_type, label.size = 3.5, dot.size = c(5, 15))

    ggsave(file.path(plot_dir,
                     paste0("netVisual_embedding_", sim_type, ".png")),
           p_embed, width = 7, height = 6, dpi = 300, bg = "white")
    ggsave(file.path(plot_dir,
                     paste0("netVisual_embedding_", sim_type, ".pdf")),
           p_embed, width = 7, height = 6)
    rm(p_embed)
    cat("✓ Embedding plot saved\n")

  }, error = function(e) {
    if (is_memory_error(e)) stop(e)
    cat("⚠️  ", sim_type, "embedding failed:", as.character(e), "\n")
  })

  do_gc(verbose = FALSE)
  check_memory_usage()
}

# ── 6. Heatmaps + gene expression per pathway ─────────────────────
cat("\n=== 6. Heatmaps + Gene Expression per Pathway ===\n")
pb <- progress_bar$new(
  format = "[:bar] :current/:total (:percent) ETA: :eta | :pathway",
  total  = length(usable), clear = FALSE, width = 80)

success_heatmap <- 0
success_expr    <- 0

for (i in seq_along(usable)) {
  p <- usable[i]
  pb$tick(tokens = list(pathway = sprintf("%-15s", p)))

  if (i %% 5 == 0) check_memory_usage()

  tryCatch({
    ht <- netVisual_heatmap(cellchat, signaling = p,
                             color.heatmap = "Reds")
    png(file.path(plot_dir, paste0(p, "_netVisual_heatmap.png")),
        width = 12, height = 8, units = "in", res = 300, type = "cairo")
    print(draw(ht)); dev.off()
    pdf(file.path(plot_dir, paste0(p, "_netVisual_heatmap.pdf")),
        width = 12, height = 8)
    print(draw(ht)); dev.off()
    rm(ht)
    success_heatmap <- success_heatmap + 1
  }, error = function(e) {
    if (is_memory_error(e)) stop(e)
    safe_dev_off()
  })

  tryCatch({
    p_expr <- plotGeneExpression(cellchat, signaling = p)
    ggsave(file.path(plot_dir, paste0(p, "_plotGeneExpression.png")),
           p_expr, width = 12, height = 8, dpi = 300, bg = "white")
    rm(p_expr)
    success_expr <- success_expr + 1
  }, error = function(e) {
    if (is_memory_error(e)) stop(e)
    invisible(NULL)
  })

  do_gc(verbose = FALSE)
}

cat(sprintf("\n✓ Heatmaps       : %d/%d\n", success_heatmap, length(usable)))
cat(sprintf("✓ Gene expression : %d/%d\n", success_expr,    length(usable)))
check_memory_usage()

# ── 7. Bubble plot ────────────────────────────────────────────────
cat("\n=== 7. Bubble Plot ===\n")
tryCatch({
  p_bubble <- netVisual_bubble(cellchat)
  ggsave(file.path(plot_dir, "netVisual_bubble_all.pdf"),
         p_bubble, width = 40, height = 100,
         units = "in", limitsize = FALSE)
  ggsave(file.path(plot_dir, "netVisual_bubble_all.png"),
         p_bubble, width = 16, height = 40, dpi = 150, bg = "white")
  rm(p_bubble)
  cat("✓ Bubble plot saved (PDF + PNG)\n")
}, error = function(e) {
  if (is_memory_error(e)) stop(e)
  cat("⚠️  Bubble plot failed:", as.character(e), "\n")
})
do_gc(verbose = FALSE)
check_memory_usage()

# ── 8. Chord diagrams per pathway ─────────────────────────────────
cat("\n=== 8. Chord Diagrams per Pathway ===\n")
pb <- progress_bar$new(
  format = "[:bar] :current/:total (:percent) ETA: :eta | :pathway",
  total  = length(usable), clear = FALSE, width = 80)

for (i in seq_along(usable)) {
  p <- usable[i]
  pb$tick(tokens = list(pathway = sprintf("%-15s", p)))

  if (i %% 5 == 0) {
    do_gc(verbose = FALSE)
    check_memory_usage()
  }

  circos.clear()
  png(file.path(plot_dir, paste0(p, "_netVisual_chord_gene.png")),
      width = 20, height = 20, units = "in", res = 300)
  try({
    circos.clear()
    circos.par(gap.degree = 1, track.margin = c(0.005, 0.005))
    netVisual_chord_gene(cellchat,
                          sources.use = 1:n_celltypes,
                          targets.use = 1:n_celltypes,
                          signaling   = p)
  }, silent = TRUE)
  dev.off()
  circos.clear()
}

do_gc(verbose = FALSE)
check_memory_usage()

# ── 9. Contribution plots per pathway ─────────────────────────────
cat("\n=== 9. Contribution Plots per Pathway ===\n")
pb <- progress_bar$new(
  format = "[:bar] :current/:total (:percent) ETA: :eta | :pathway",
  total  = length(usable), clear = FALSE, width = 80)

for (i in seq_along(usable)) {
  p <- usable[i]
  pb$tick(tokens = list(pathway = sprintf("%-15s", p)))

  if (i %% 5 == 0) {
    do_gc(verbose = FALSE)
    check_memory_usage()
  }

  tryCatch({
    p_contrib <- netAnalysis_contribution(cellchat, signaling = p)
    if (is_gg(p_contrib))
      ggsave(file.path(plot_dir,
                       paste0(p, "_netAnalysis_contribution.png")),
             p_contrib, width = 6, height = 4, dpi = 300, bg = "white")
    rm(p_contrib)
  }, error = function(e) {
    if (is_memory_error(e)) stop(e)
    invisible(NULL)
  })
}

do_gc(verbose = FALSE)
check_memory_usage()

# ── 10. Combined signaling role heatmap ───────────────────────────
cat("\n=== 10. Combined Signaling Role Heatmap ===\n")
ht1 <- netAnalysis_signalingRole_heatmap(
  cellchat, pattern = "outgoing",
  font.size = 4, font.size.title = 10)
ht2 <- netAnalysis_signalingRole_heatmap(
  cellchat, pattern = "incoming",
  font.size = 4, font.size.title = 10)
ht_combined <- ht1 + ht2

png(file.path(plot_dir,
              paste0(sample_name, "_signalingRole_heatmap_combined.png")),
    width = 16, height = 16, units = "in", res = 300, type = "cairo")
print(ht_combined); dev.off()

pdf(file.path(plot_dir,
              paste0(sample_name, "_signalingRole_heatmap_combined.pdf")),
    width = 16, height = 16)
print(ht_combined); dev.off()

rm(ht1, ht2, ht_combined)
do_gc(verbose = FALSE)
cat("✓ Combined heatmap saved (PNG + PDF)\n")

# ══════════════════════════════════════════════════════════════════
# SAVE FINAL RDS
# ══════════════════════════════════════════════════════════════════

cat("\n=== Saving Final RDS ===\n")
saveRDS(cellchat,
        file.path(plot_dir,
                  paste0("cellchat_", sample_name, "_at_the_end.rds")))
cat("✓ Final RDS saved\n")

# ══════════════════════════════════════════════════════════════════
# EXPORT ALL DATA
# ══════════════════════════════════════════════════════════════════

cat("\n=== Exporting All CellChat Data ===\n")

write.csv(
  cellchat@LR$LRsig,
  file.path(plot_dir, paste0(sample_name, "_01_all_LR_pairs.csv")),
  row.names = FALSE)

all_lr_comm <- subsetCommunication(cellchat, slot.name = "net",
                                    thresh = 0.05)
write.csv(
  all_lr_comm,
  file.path(plot_dir, paste0(sample_name, "_02_all_LR_communications.csv")),
  row.names = FALSE)

all_pathway_comm <- subsetCommunication(cellchat, slot.name = "netP",
                                         thresh = 0.05)
write.csv(
  all_pathway_comm,
  file.path(plot_dir, paste0(sample_name, "_03_all_pathway_communications.csv")),
  row.names = FALSE)

write.csv(
  cellchat@net$count,
  file.path(plot_dir, paste0(sample_name,
                              "_05_communication_count_matrix.csv")))

write.csv(
  cellchat@net$weight,
  file.path(plot_dir, paste0(sample_name,
                              "_06_communication_weight_matrix.csv")))

writeLines(
  cellchat@netP$pathways,
  file.path(plot_dir, paste0(sample_name, "_07_active_pathways.txt")))

rm(all_pathway_comm)
do_gc(verbose = FALSE)

# Per-pathway CSVs
cat("\n--- Per-pathway CSVs ---\n")
pb <- progress_bar$new(
  format = "[:bar] :current/:total (:percent) | :pathway",
  total  = length(cellchat@netP$pathways), clear = FALSE, width = 80)

for (i in seq_along(cellchat@netP$pathways)) {
  pathway <- cellchat@netP$pathways[i]
  pb$tick(tokens = list(pathway = sprintf("%-15s", pathway)))

  if (i %% 20 == 0) {
    do_gc(verbose = FALSE)
    check_memory_usage()
  }

  pw_data <- subsetCommunication(cellchat, signaling = pathway,
                                  thresh = 0.05)
  if (nrow(pw_data) > 0)
    write.csv(
      pw_data,
      file.path(plot_dir,
                paste0(sample_name, "_", pathway,
                       "_pathway_interactions.csv")),
      row.names = FALSE)
  rm(pw_data)
}

cat("\n--- Export Summary ---\n")
cat(sprintf("✓ L-R pairs            : %d\n", nrow(cellchat@LR$LRsig)))
cat(sprintf("✓ Total communications : %d\n", nrow(all_lr_comm)))
cat(sprintf("✓ Active pathways      : %d\n", length(cellchat@netP$pathways)))
cat("✓ All data exported to:", plot_dir, "\n")

rm(all_lr_comm)
do_gc()
check_memory_usage()

# ══════════════════════════════════════════════════════════════════
# WRAP UP
# ══════════════════════════════════════════════════════════════════

cat("\n")
cat("╔═══════════════════════════════════════════════════════════╗\n")
cat("║                                                           ║\n")
cat("║   ✓✓✓  CELLCHAT ANALYSIS COMPLETE  ✓✓✓                  ║\n")
cat("║                                                           ║\n")
cat("╚═══════════════════════════════════════════════════════════╝\n")
cat("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Runtime  :", round(difftime(Sys.time(), start_time, units = "mins"), 2),
    "minutes\n")
cat("Log saved:", log_file, "\n\n")

sink()
