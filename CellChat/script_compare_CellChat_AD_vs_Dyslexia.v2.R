# ══════════════════════════════════════════════════════════════════
# CELLCHAT COMPARISON ANALYSIS
# Input : 18_merged_AG_031826.28celltypes.AD_at_the_end.rds
# Input2: 18_merged_AG_031826.28celltypes.Dyslexia_at_the_end.rds
# Groups: AD vs Dyslexia
# pos.dataset = "Dyslexia"
# logFC > 0 = higher in Dyslexia vs AD
# logFC < 0 = lower in Dyslexia vs AD (higher in AD)
# ══════════════════════════════════════════════════════════════════

# ══════════════════════════════════════════════════════════════════
# CONTROL FLAGS
# ══════════════════════════════════════════════════════════════════

run_AD_perspective <- FALSE   # ← set TRUE to also run AD perspective

# ── Libraries ────────────────────────────────────────────────────
library(CellChat)
library(patchwork)
library(tidyverse)
library(uwot)
library(ggrepel)
library(pheatmap)
library(ComplexHeatmap)
library(future)
library(presto)
library(Matrix)
plan("sequential")

# ══════════════════════════════════════════════════════════════════
# START LOGGING
# ══════════════════════════════════════════════════════════════════

log_file <- paste0("cellchat_comparison_AD_vs_Dyslexia_",
                   format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt")
sink(log_file, split = TRUE)

cat("=================================================\n")
cat("CellChat Comparison Analysis Log\n")
cat("AD vs Dyslexia\n")
cat("pos.dataset = Dyslexia\n")
cat("logFC > 0 = higher in Dyslexia vs AD\n")
cat("logFC < 0 = lower in Dyslexia vs AD (higher in AD)\n")
cat("Started at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("=================================================\n\n")

# ══════════════════════════════════════════════════════════════════
# MEMORY MONITORING
# ══════════════════════════════════════════════════════════════════

check_memory_usage <- function(threshold = 0.90, stop_on_exceed = TRUE) {
  meminfo <- tryCatch(readLines("/proc/meminfo", n = 20),
                      error = function(e) NULL)
  if (is.null(meminfo)) return(list(percent_used = NA, should_stop = FALSE))

  memtotal_kb     <- as.numeric(gsub("[^0-9]", "", meminfo[1]))
  memavailable_kb <- as.numeric(gsub("[^0-9]", "", meminfo[3]))
  memused_kb      <- memtotal_kb - memavailable_kb
  percent_used    <- memused_kb / memtotal_kb

  cat(sprintf("Memory: %.1f%% used (%.1f GB / %.1f GB)\n",
              percent_used * 100,
              memused_kb  / 1024^2,
              memtotal_kb / 1024^2))

  should_stop <- percent_used >= threshold
  if (should_stop && stop_on_exceed)
    cat("❌ STOPPING: Memory threshold exceeded!\n")

  return(list(percent_used = percent_used, should_stop = should_stop))
}

cat("=== INITIAL MEMORY CHECK ===\n")
mem_status <- check_memory_usage()
if (mem_status$should_stop) { sink(); stop("High memory") }

# ══════════════════════════════════════════════════════════════════
# WORKING DIRECTORY
# ══════════════════════════════════════════════════════════════════

base_dir <- "/mnt/nfs/CX000008_DS1/projects/jaeyeon/fastq_file_Dyslexia_r1/zanalysis_bogdan/samples.merged_AG_Harmony_res0.x_anno_031826"

work_dir <- file.path(base_dir, "compare_AD_vs_Dyslexia_cellchat")
dir.create(work_dir, showWarnings = FALSE, recursive = TRUE)
setwd(work_dir)
cat("Working directory:", getwd(), "\n\n")

# ══════════════════════════════════════════════════════════════════
# LOAD CELLCHAT OBJECTS
# ══════════════════════════════════════════════════════════════════

cat("=== LOADING CELLCHAT OBJECTS ===\n")

ad <- readRDS(file.path(base_dir,
  "18_merged_AG_031826.28celltypes.AD_cellchat",
  "cellchat_18_merged_AG_031826.28celltypes.AD_at_the_end.rds"))
cat("✓ Loaded: AD\n")

cat("\nAD metadata:\n")
cat("  subcelltype:\n"); print(sort(unique(ad@meta$subcelltype)))
cat("  group:\n");       print(table(ad@meta$group))

dyslexia <- readRDS(file.path(base_dir,
  "18_merged_AG_031826.28celltypes.Dyslexia_cellchat",
  "cellchat_18_merged_AG_031826.28celltypes.Dyslexia_at_the_end.rds"))
cat("✓ Loaded: Dyslexia\n")

cat("\nDyslexia metadata:\n")
cat("  subcelltype:\n"); print(sort(unique(dyslexia@meta$subcelltype)))
cat("  group:\n");       print(table(dyslexia@meta$group))

mem_status <- check_memory_usage()
if (mem_status$should_stop) { sink(); stop("High memory") }

# ══════════════════════════════════════════════════════════════════
# UPDATE AND MERGE
# ══════════════════════════════════════════════════════════════════

cat("\n=== UPDATING CELLCHAT OBJECTS ===\n")
ad  <- updateCellChat(ad)
dyslexia <- updateCellChat(dyslexia)
cat("✓ Both objects updated\n")

# ══════════════════════════════════════════════════════════════════
# CELL TYPE COMPATIBILITY
# ══════════════════════════════════════════════════════════════════

cat("\n=== CELL TYPE COMPATIBILITY CHECK ===\n")
ct_AD  <- sort(levels(ad@idents))
ct_Dyslexia <- sort(levels(dyslexia@idents))

cat(sprintf("Cell types in AD  : %d\n", length(ct_AD)))
cat(sprintf("Cell types in Dyslexia : %d\n", length(ct_Dyslexia)))

only_AD  <- setdiff(ct_AD,  ct_Dyslexia)
only_Dyslexia <- setdiff(ct_Dyslexia, ct_AD)
common_ct     <- intersect(ct_AD, ct_Dyslexia)

cat(sprintf("Only in AD  : %d : %s\n",
            length(only_AD),  paste(only_AD,  collapse = ", ")))
cat(sprintf("Only in Dyslexia : %d : %s\n",
            length(only_Dyslexia), paste(only_Dyslexia, collapse = ", ")))
cat(sprintf("Common cell types: %d\n", length(common_ct)))
cat("Common cell types:\n")
print(sort(common_ct))

cat("\n✓ Subsetting both objects to common cell types...\n")
ad_sub  <- subsetCellChat(ad,  idents.use = common_ct)
dyslexia_sub <- subsetCellChat(dyslexia, idents.use = common_ct)

cat(sprintf("✓ AD  after subset: %d cell types\n",
            length(levels(ad_sub@idents))))
cat(sprintf("✓ Dyslexia after subset: %d cell types\n",
            length(levels(dyslexia_sub@idents))))

if (length(only_AD) > 0)
  cat(sprintf("⚠️  Excluded from AD  : %s\n",
              paste(only_AD, collapse = ", ")))
if (length(only_Dyslexia) > 0)
  cat(sprintf("⚠️  Excluded from Dyslexia : %s\n",
              paste(only_Dyslexia, collapse = ", ")))

object.list <- list(AD = ad_sub, Dyslexia = dyslexia_sub)

# ══════════════════════════════════════════════════════════════════
# MERGE
# ══════════════════════════════════════════════════════════════════

cat("\n=== MERGING CELLCHAT OBJECTS ===\n")
mem_status <- check_memory_usage()
if (mem_status$should_stop) { sink(); stop("High memory") }

cellchat <- mergeCellChat(object.list, add.names = names(object.list))
cat("✓ Merged CellChat object\n")

cat("\nMerged object overview:\n")
print(cellchat)
cat("\nDatasets:\n");   print(table(cellchat@meta$datasets))
cat("\nCell types:\n"); print(sort(unique(cellchat@meta$subcelltype)))

# ── Interaction summaries ─────────────────────────────────────────
cat("\n=== INTERACTION SUMMARIES ===\n")
for (nm in names(object.list)) {
  obj <- object.list[[nm]]
  cat(sprintf("\n%s:\n", nm))
  cat(sprintf("  L-R pairs         : %d\n", nrow(obj@LR$LRsig)))
  cat(sprintf("  Total interactions: %d\n", sum(obj@net$count)))
  cat(sprintf("  (no self)         : %d\n",
              sum(obj@net$count) - sum(diag(obj@net$count))))
}

# ── Interaction count/weight heatmaps ────────────────────────────
for (nm in names(object.list)) {
  png(sprintf("%s_interaction_count_heatmap.png", nm),
      width = 8, height = 7, units = "in", res = 300)
  pheatmap::pheatmap(object.list[[nm]]@net$count,
                     main = sprintf("%s: Number of Interactions", nm))
  dev.off()

  png(sprintf("%s_interaction_weight_heatmap.png", nm),
      width = 8, height = 7, units = "in", res = 300)
  pheatmap::pheatmap(object.list[[nm]]@net$weight,
                     main = sprintf("%s: Interaction Strength", nm))
  dev.off()
}
cat("✓ Interaction heatmaps saved\n")

# ══════════════════════════════════════════════════════════════════
# COMPARE INTERACTIONS
# ══════════════════════════════════════════════════════════════════

cat("\n=== COMPARE INTERACTIONS ===\n")
gg1 <- compareInteractions(cellchat, show.legend = FALSE, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = FALSE,
                            group = c(1,2), measure = "weight")
ggsave("compareInteractions_count_weight.png",
       gg1 + gg2, width = 10, height = 5, dpi = 300)
cat("✓ compareInteractions saved\n")

# ── Differential heatmaps ─────────────────────────────────────────
tryCatch({
  gg_count <- netVisual_heatmap(cellchat)
  png("netVisual_heatmap_count.png",
      width = 6, height = 6, units = "in", res = 300)
  ComplexHeatmap::draw(gg_count); dev.off()
  cat("✓ Differential heatmap (count) saved\n")
}, error = function(e) {
  if (dev.cur() > 1) dev.off()
  cat("⚠️  netVisual_heatmap (count) failed:", as.character(e), "\n")
})

tryCatch({
  gg_weight <- netVisual_heatmap(cellchat, measure = "weight")
  png("netVisual_heatmap_weight.png",
      width = 6, height = 6, units = "in", res = 300)
  ComplexHeatmap::draw(gg_weight); dev.off()
  cat("✓ Differential heatmap (weight) saved\n")
}, error = function(e) {
  if (dev.cur() > 1) dev.off()
  cat("⚠️  netVisual_heatmap (weight) failed:", as.character(e), "\n")
})

# ── Circle plots ──────────────────────────────────────────────────
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))

png("netVisual_circle_count_comparison.png",
    width = 16, height = 8, units = "in", res = 300)
par(mfrow = c(1,2), xpd = TRUE)
for (i in seq_along(object.list)) {
  netVisual_circle(object.list[[i]]@net$count,
                   weight.scale    = TRUE,
                   label.edge      = FALSE,
                   edge.weight.max = weight.max[2],
                   edge.width.max  = 12,
                   title.name      = paste0("Interactions - ",
                                            names(object.list)[i]))
}
dev.off()
cat("✓ Circle plots saved\n")

# ── Differential interaction plots ───────────────────────────────
for (measure in c("count.merged", "weight.merged")) {
  png(sprintf("netVisual_diffInteraction_%s.png", measure),
      width = 14, height = 7, units = "in", res = 300)
  par(mfrow = c(1,2), xpd = TRUE)
  netVisual_diffInteraction(cellchat, weight.scale = TRUE,
                             measure = measure, label.edge = TRUE)
  dev.off()
}
cat("✓ Differential interaction plots saved\n")

# ── Signaling role scatter ────────────────────────────────────────
cat("\n=== SIGNALING ROLE SCATTER ===\n")

num.link <- sapply(object.list, function(x) {
  rowSums(x@net$count) + colSums(x@net$count) - diag(x@net$count)
})
weight.MinMax <- c(min(num.link), max(num.link))

gg <- list()
for (i in seq_along(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(
    object.list[[i]], title = names(object.list)[i],
    weight.MinMax = weight.MinMax)
}
ggsave("netAnalysis_signalingRole_scatter_count.png",
       patchwork::wrap_plots(plots = gg),
       width = 16, height = 6, dpi = 300)

num.link.w <- sapply(object.list, function(x) {
  rowSums(x@net$weight) + colSums(x@net$weight) - diag(x@net$weight)
})
weight.MinMax.w <- c(min(num.link.w), max(num.link.w))

gg.w <- list()
for (i in seq_along(object.list)) {
  gg.w[[i]] <- netAnalysis_signalingRole_scatter(
    object.list[[i]], title = names(object.list)[i],
    weight.MinMax = weight.MinMax.w) +
    ggplot2::scale_size_continuous(range = c(3, 10))
}
ggsave("netAnalysis_signalingRole_scatter_weight.png",
       patchwork::wrap_plots(plots = gg.w),
       width = 16, height = 6, dpi = 300)
cat("✓ Signaling role scatter plots saved\n")

# ── Signaling changes per cell type ──────────────────────────────
cat("\n=== SIGNALING CHANGES PER CELL TYPE ===\n")

celltypes_common <- intersect(
  sort(unique(ad@meta$subcelltype)),
  sort(unique(dyslexia@meta$subcelltype))
)
cat(sprintf("Cell types in common: %d\n", length(celltypes_common)))

dir.create("signalingChanges", showWarnings = FALSE)

for (ct in celltypes_common) {
  ct_safe <- gsub("[^A-Za-z0-9_]", "_", ct)
  tryCatch({
    p <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = ct)
    ggsave(file.path("signalingChanges",
                     sprintf("%s_signalingChanges.png", ct_safe)),
           p, width = 8, height = 6, dpi = 300)
  }, error = function(e) {
    cat(sprintf("⚠️  signalingChanges failed for %s: %s\n",
                ct, as.character(e)))
  })
}
cat("✓ Signaling changes plots saved to: signalingChanges/\n")

# ── Outgoing / Incoming signaling heatmaps ────────────────────────
cat("\n=== OUTGOING / INCOMING SIGNALING HEATMAPS ===\n")

pathway.union <- union(object.list[[1]]@netP$pathways,
                       object.list[[2]]@netP$pathways)
cat(sprintf("Total pathway union: %d\n", length(pathway.union)))

for (pattern in c("outgoing", "incoming", "all")) {
  col <- switch(pattern,
                outgoing = "BuGn",
                incoming = "GnBu",
                all      = "OrRd")
  tryCatch({
    ht1 <- netAnalysis_signalingRole_heatmap(
      object.list[[1]], pattern = pattern,
      signaling = pathway.union,
      title     = names(object.list)[1],
      width = 5, height = 6,
      color.heatmap = col)
    ht2 <- netAnalysis_signalingRole_heatmap(
      object.list[[2]], pattern = pattern,
      signaling = pathway.union,
      title     = names(object.list)[2],
      width = 5, height = 6,
      color.heatmap = col)
    png(sprintf("netAnalysis_signalingRole_heatmap_%s.png", pattern),
        width = 12, height = 8, units = "in", res = 300)
    ComplexHeatmap::draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
    dev.off()
    cat(sprintf("✓ Signaling heatmap (%s) saved\n", pattern))
  }, error = function(e) {
    if (dev.cur() > 1) dev.off()
    cat(sprintf("⚠️  Signaling heatmap (%s) failed: %s\n",
                pattern, as.character(e)))
  })
}

# ══════════════════════════════════════════════════════════════════
# FUNCTIONAL SIMILARITY
# ══════════════════════════════════════════════════════════════════

cat("\n=== FUNCTIONAL SIMILARITY ===\n")
output_dir_functional <- "cellchat_functional_plots"
dir.create(output_dir_functional, showWarnings = FALSE)

mem_status <- check_memory_usage()
if (mem_status$should_stop) { sink(); stop("High memory") }

cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional", umap.method = "uwot")
cellchat <- netClustering(cellchat, type = "functional", do.parallel = FALSE)

png(file.path(output_dir_functional, "cellchat_functional_embedding.png"),
    width = 12, height = 10, units = "in", res = 300)
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
dev.off()

png(file.path(output_dir_functional, "cellchat_functional_embedding_zoomin.png"),
    width = 16, height = 8, units = "in", res = 300)
netVisual_embeddingPairwiseZoomIn(cellchat, type = "functional", nCol = 2)
dev.off()
cat("✓ Functional embedding plots saved\n")

M  <- cellchat@netP$similarity$functional$matrix[["1-2"]]
DR <- cellchat@netP$similarity$functional$dr[["1-2"]]

png(file.path(output_dir_functional, "functional_similarity_heatmap.png"),
    width = 14, height = 14, units = "in", res = 300)
pheatmap(M, show_rownames = TRUE, show_colnames = TRUE,
         fontsize_row = 5, fontsize_col = 5,
         main = "Pathway Functional Similarity (AD vs Dyslexia)")
dev.off()

df <- as.data.frame(DR)
colnames(df) <- c("UMAP1", "UMAP2")
df$pathway       <- rownames(DR)
df$condition     <- sub(".*--", "", df$pathway)
df$pathway_clean <- sub("--.*", "", df$pathway)

ctl_df <- df %>% filter(condition == "AD") %>%
  select(pathway = pathway_clean, UMAP1, UMAP2)
dys_df <- df %>% filter(condition == "Dyslexia") %>%
  select(pathway = pathway_clean, UMAP1, UMAP2)

merged_func <- full_join(ctl_df, dys_df,
                          by = "pathway",
                          suffix = c("_AD", "_Dyslexia")) %>%
  mutate(shift = ifelse(
    !is.na(UMAP1_AD) & !is.na(UMAP1_Dyslexia),
    sqrt((UMAP1_AD - UMAP1_Dyslexia)^2 +
           (UMAP2_AD - UMAP2_Dyslexia)^2),
    NA_real_)) %>%
  arrange(desc(shift))

topN     <- min(500, nrow(merged_func))
to_label <- merged_func$pathway[1:topN]

p1 <- ggplot(df, aes(UMAP1, UMAP2, color = condition)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_text_repel(
    data = df %>% filter(pathway_clean %in% to_label),
    aes(label = pathway_clean),
    size = 3.5, max.overlaps = Inf, box.padding = 0.5) +
  scale_color_manual(values = c("AD"  = "#E41A1C",
                                 "Dyslexia" = "#377EB8")) +
  theme_classic(base_size = 13) +
  labs(title = sprintf("Functional Similarity - Top %d Shifted Pathways", topN),
       color = "Condition") +
  coord_fixed()

ggsave(file.path(output_dir_functional, "functional_umap_labeled.png"),
       p1, width = 12, height = 10, dpi = 300)

merged_both_func <- merged_func %>% filter(!is.na(shift))

p2 <- ggplot() +
  geom_segment(
    data = merged_both_func,
    aes(x = UMAP1_AD, y = UMAP2_AD,
        xend = UMAP1_Dyslexia, yend = UMAP2_Dyslexia),
    arrow = arrow(length = unit(0.12, "inches")),
    alpha = 0.5, color = "grey60") +
  geom_point(data = df, aes(UMAP1, UMAP2, color = condition),
             size = 2, alpha = 0.8) +
  geom_text_repel(
    data = merged_both_func[1:min(15, nrow(merged_both_func)), ],
    aes(x = UMAP1_Dyslexia, y = UMAP2_Dyslexia, label = pathway),
    size = 3.5, max.overlaps = Inf) +
  scale_color_manual(values = c("AD"  = "#E41A1C",
                                 "Dyslexia" = "#377EB8")) +
  theme_classic(base_size = 13) +
  labs(title = "Functional Shift (AD → Dyslexia)", color = "Condition") +
  coord_fixed()

ggsave(file.path(output_dir_functional, "functional_pathway_shifts.png"),
       p2, width = 12, height = 10, dpi = 300)

write.csv(merged_func,
          file.path(output_dir_functional, "functional_pathway_shifts_ranked.csv"),
          row.names = FALSE)
write.csv(M,
          file.path(output_dir_functional, "functional_similarity_matrix.csv"),
          row.names = TRUE)
write.csv(df,
          file.path(output_dir_functional, "functional_umap_coordinates.csv"),
          row.names = FALSE)
cat("✓ Functional similarity complete\n")

tryCatch({
  ranked_func <- rankSimilarity(cellchat, type = "functional", slot.name = "netP")
  png(file.path(output_dir_functional, "ranked_functional_similarity.png"),
      width = 10, height = 8, units = "in", res = 300)
  print(ranked_func); dev.off()
  cat("✓ Ranked functional similarity saved\n")
}, error = function(e) {
  if (dev.cur() > 1) dev.off()
  cat("⚠️  rankSimilarity functional failed:", as.character(e), "\n")
})

# ══════════════════════════════════════════════════════════════════
# STRUCTURAL SIMILARITY
# ══════════════════════════════════════════════════════════════════

cat("\n=== STRUCTURAL SIMILARITY ===\n")
output_dir_structural <- "cellchat_structural_plots"
dir.create(output_dir_structural, showWarnings = FALSE)

mem_status <- check_memory_usage()
if (mem_status$should_stop) { sink(); stop("High memory") }

cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural", umap.method = "uwot")
cellchat <- netClustering(cellchat, type = "structural", do.parallel = FALSE)

png(file.path(output_dir_structural, "structural_cellchat_embedding.png"),
    width = 12, height = 10, units = "in", res = 300)
netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5)
dev.off()

png(file.path(output_dir_structural, "structural_cellchat_embedding_zoomin.png"),
    width = 16, height = 8, units = "in", res = 300)
netVisual_embeddingPairwiseZoomIn(cellchat, type = "structural", nCol = 2)
dev.off()
cat("✓ Structural embedding plots saved\n")

M_struct  <- cellchat@netP$similarity$structural$matrix[["1-2"]]
DR_struct <- cellchat@netP$similarity$structural$dr[["1-2"]]

png(file.path(output_dir_structural, "similarity_heatmap.png"),
    width = 14, height = 14, units = "in", res = 300)
pheatmap(M_struct, show_rownames = TRUE, show_colnames = TRUE,
         fontsize_row = 5, fontsize_col = 5,
         main = "Structural Similarity (AD vs Dyslexia)")
dev.off()

df_s <- as.data.frame(DR_struct)
colnames(df_s) <- c("UMAP1", "UMAP2")
df_s$pathway       <- rownames(DR_struct)
df_s$condition     <- sub(".*--", "", df_s$pathway)
df_s$pathway_clean <- sub("--.*", "", df_s$pathway)

ctl_s <- df_s %>% filter(condition == "AD") %>%
  select(pathway = pathway_clean, UMAP1, UMAP2)
dys_s <- df_s %>% filter(condition == "Dyslexia") %>%
  select(pathway = pathway_clean, UMAP1, UMAP2)

merged_struct <- full_join(ctl_s, dys_s,
                            by = "pathway",
                            suffix = c("_AD", "_Dyslexia")) %>%
  mutate(shift = ifelse(
    !is.na(UMAP1_AD) & !is.na(UMAP1_Dyslexia),
    sqrt((UMAP1_AD - UMAP1_Dyslexia)^2 +
           (UMAP2_AD - UMAP2_Dyslexia)^2),
    NA_real_)) %>%
  arrange(desc(shift))

topN_s     <- min(500, nrow(merged_struct))
to_label_s <- merged_struct$pathway[1:topN_s]

p1_s <- ggplot(df_s, aes(UMAP1, UMAP2, color = condition)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_text_repel(
    data = df_s %>% filter(pathway_clean %in% to_label_s),
    aes(label = pathway_clean),
    size = 3.5, max.overlaps = Inf, box.padding = 0.5) +
  scale_color_manual(values = c("AD"  = "#E41A1C",
                                 "Dyslexia" = "#377EB8")) +
  theme_classic(base_size = 13) +
  labs(title = sprintf("Structural Similarity - Top %d Shifted Pathways", topN_s),
       color = "Condition") +
  coord_fixed()

ggsave(file.path(output_dir_structural, "umap_labeled.png"),
       p1_s, width = 12, height = 10, dpi = 300)

merged_both_s <- merged_struct %>% filter(!is.na(shift))

p2_s <- ggplot() +
  geom_segment(
    data = merged_both_s,
    aes(x = UMAP1_AD, y = UMAP2_AD,
        xend = UMAP1_Dyslexia, yend = UMAP2_Dyslexia),
    arrow = arrow(length = unit(0.12, "inches")),
    alpha = 0.5, color = "grey60") +
  geom_point(data = df_s, aes(UMAP1, UMAP2, color = condition),
             size = 2, alpha = 0.8) +
  geom_text_repel(
    data = merged_both_s[1:min(15, nrow(merged_both_s)), ],
    aes(x = UMAP1_Dyslexia, y = UMAP2_Dyslexia, label = pathway),
    size = 3.5, max.overlaps = Inf) +
  scale_color_manual(values = c("AD"  = "#E41A1C",
                                 "Dyslexia" = "#377EB8")) +
  theme_classic(base_size = 13) +
  labs(title = "Structural Shift (AD → Dyslexia)", color = "Condition") +
  coord_fixed()

ggsave(file.path(output_dir_structural, "pathway_shifts.png"),
       p2_s, width = 12, height = 10, dpi = 300)

write.csv(merged_struct,
          file.path(output_dir_structural, "structural_pathway_shifts_ranked.csv"),
          row.names = FALSE)
write.csv(M_struct,
          file.path(output_dir_structural, "structural_similarity_matrix.csv"),
          row.names = TRUE)
write.csv(df_s,
          file.path(output_dir_structural, "structural_umap_coordinates.csv"),
          row.names = FALSE)
cat("✓ Structural similarity complete\n")

tryCatch({
  ranked_struct <- rankSimilarity(cellchat, type = "structural", slot.name = "netP")
  png(file.path(output_dir_structural, "ranked_structural_similarity.png"),
      width = 10, height = 8, units = "in", res = 300)
  print(ranked_struct); dev.off()
  cat("✓ Ranked structural similarity saved\n")
}, error = function(e) {
  if (dev.cur() > 1) dev.off()
  cat("⚠️  rankSimilarity structural failed:", as.character(e), "\n")
})

# ══════════════════════════════════════════════════════════════════
# RANKNET
# ══════════════════════════════════════════════════════════════════

cat("\n=== RANKNET ===\n")

signaling.type <- c("Secreted Signaling", "ECM-Receptor",
                    "Cell-Cell Contact", "Non-protein Signaling")

ranknet_configs <- list(
  list(slot = "netP", stacked = TRUE,  file = "ranknet_weight_stacked.pdf",     h = 8),
  list(slot = "netP", stacked = FALSE, file = "ranknet_weight.pdf",             h = 8),
  list(slot = "net",  stacked = TRUE,  file = "ranknet_weight_net_stacked.pdf", h = 40),
  list(slot = "net",  stacked = FALSE, file = "ranknet_weight_net.pdf",         h = 40)
)

for (cfg in ranknet_configs) {
  tryCatch({
    gg <- rankNet(cellchat,
                  slot.name      = cfg$slot,
                  mode           = "comparison",
                  measure        = "weight",
                  stacked        = cfg$stacked,
                  do.stat        = TRUE,
                  signaling.type = signaling.type)
    ggsave(cfg$file, gg, width = 10, height = cfg$h)
    cat("✓ Saved:", cfg$file, "\n")
  }, error = function(e) {
    cat("⚠️  rankNet failed for", cfg$file, ":", as.character(e), "\n")
  })
}

# ══════════════════════════════════════════════════════════════════
# BUBBLE PLOTS
# ══════════════════════════════════════════════════════════════════

cat("\n=== BUBBLE PLOTS ===\n")

tryCatch({
  pdf("netVisual_bubble_comparison_all.pdf", width = 40, height = 40)
  netVisual_bubble(cellchat, comparison = c(1, 2), angle.x = 45)
  dev.off()
  cat("✓ Bubble plot (all) saved\n")
}, error = function(e) {
  if (dev.cur() > 1) dev.off()
  cat("⚠️  Bubble plot (all) failed:", as.character(e), "\n")
})

tryCatch({
  png("netVisual_bubble_increased_in_Dyslexia.png",
      width = 40, height = 40, units = "in", res = 200)
  netVisual_bubble(cellchat, comparison = c(1, 2),
                   max.dataset    = 2,
                   title.name     = "Increased signaling in Dyslexia vs AD",
                   angle.x        = 45,
                   remove.isolate = TRUE)
  dev.off()
  cat("✓ Bubble plot (increased in Dyslexia) saved\n")
}, error = function(e) {
  if (dev.cur() > 1) dev.off()
  cat("⚠️  Bubble plot (increased Dyslexia) failed:", as.character(e), "\n")
})

tryCatch({
  png("netVisual_bubble_decreased_in_Dyslexia.png",
      width = 40, height = 40, units = "in", res = 200)
  netVisual_bubble(cellchat, comparison = c(1, 2),
                   max.dataset    = 1,
                   title.name     = "Decreased signaling in Dyslexia (higher in AD)",
                   angle.x        = 45,
                   remove.isolate = TRUE)
  dev.off()
  cat("✓ Bubble plot (decreased in Dyslexia) saved\n")
}, error = function(e) {
  if (dev.cur() > 1) dev.off()
  cat("⚠️  Bubble plot (decreased Dyslexia) failed:", as.character(e), "\n")
})

# ══════════════════════════════════════════════════════════════════
# DIFFERENTIAL EXPRESSION — netMappingDEG
# ══════════════════════════════════════════════════════════════════

cat("\n=== DIFFERENTIAL EXPRESSION ANALYSIS ===\n")

# ── Set datasets as factor ────────────────────────────────────────
cat("--- Datasets check ---\n")
print(table(cellchat@meta$datasets, useNA = "ifany"))
cellchat@meta$datasets <- factor(as.character(cellchat@meta$datasets),
                                  levels = c("AD", "Dyslexia"))
cat("✓ datasets set as factor: AD, Dyslexia\n")

deg_success_Dyslexia <- FALSE
deg_success_AD2  <- FALSE

# ══════════════════════════════════════════════════════════════════
# RUN 1 — DYSLEXIA PERSPECTIVE (always runs)
# pos.dataset = "Dyslexia"
# logFC > 0 = higher in Dyslexia vs AD
# logFC < 0 = lower in Dyslexia vs AD (higher in AD)
# ══════════════════════════════════════════════════════════════════

cat("\n--- RUN 1: pos.dataset = Dyslexia ---\n")
cat("logFC > 0 = higher in Dyslexia vs AD\n")
cat("logFC < 0 = lower  in Dyslexia vs AD\n\n")

pos.dataset   <- "Dyslexia"
features.name <- paste0(pos.dataset, ".merged")  # "Dyslexia.merged"

tryCatch({

  # ── STEP 1: identifyOverExpressedGenes ───────────────────────
  cat("Running identifyOverExpressedGenes (Dyslexia)...\n")
  cellchat <- identifyOverExpressedGenes(
    cellchat,
    group.dataset = "datasets",
    pos.dataset   = pos.dataset,
    features.name = features.name,
    only.pos      = FALSE,
    thresh.pc     = 0.1,
    thresh.fc     = 0.05,
    thresh.p      = 0.05)

  OE_genes_Dyslexia <- cellchat@var.features[[features.name]]
  cat(sprintf("✓ OE genes in Dyslexia: %d\n", length(OE_genes_Dyslexia)))

  cat("\n===== OE Genes in Dyslexia =====\n")
  cat(sprintf("Total OE genes: %d\n", length(OE_genes_Dyslexia)))
  print(OE_genes_Dyslexia)

  write.csv(data.frame(gene = OE_genes_Dyslexia),
            "results.OE_genes_Dyslexia.csv", row.names = FALSE)
  cat("✓ Saved: results.OE_genes_Dyslexia.csv\n")

  # ── STEP 2: identifyOverExpressedInteractions ─────────────────
  # NOTE: this is a SUBSET of identifyOverExpressedGenes
  # It keeps only genes that form LR pairs in CellChatDB
  # where BOTH ligand AND receptor are OE (variable.both = TRUE)
  cat("\nRunning identifyOverExpressedInteractions (Dyslexia)...\n")
  cellchat <- identifyOverExpressedInteractions(
    cellchat,
    features.name = features.name,   # uses results from STEP 1
    variable.both = TRUE,            # BOTH ligand AND receptor must be OE
    return.object = TRUE)

  LR_OE_Dyslexia <- cellchat@var.features[[paste0(features.name, ".LR")]]
  cat(sprintf("✓ OE LR interactions in Dyslexia: %d\n", length(LR_OE_Dyslexia)))

  cat("\n===== OE LR Interactions in Dyslexia =====\n")
  cat(sprintf("Total OE LR interactions: %d\n", length(LR_OE_Dyslexia)))
  cat(sprintf("(subset of %d OE genes)\n", length(OE_genes_Dyslexia)))
  df_OE_LR_Dyslexia <- data.frame(rank        = seq_along(LR_OE_Dyslexia),
                                   interaction = LR_OE_Dyslexia)
  print(df_OE_LR_Dyslexia, row.names = FALSE)

  write.csv(df_OE_LR_Dyslexia,
            "results.overExpressed_LR_interactions_Dyslexia.csv",
            row.names = FALSE)
  write.table(LR_OE_Dyslexia,
              "results.overExpressed_LR_interactions_Dyslexia.txt",
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  cat("✓ Saved: results.overExpressed_LR_interactions_Dyslexia.csv\n")
  cat("✓ Saved: results.overExpressed_LR_interactions_Dyslexia.txt\n")

  deg_success_Dyslexia <- TRUE

}, error = function(e) {
  cat("⚠️  RUN 1 failed:", as.character(e), "\n")
})

if (deg_success_Dyslexia) {

  net_Dyslexia <- netMappingDEG(cellchat, features.name = features.name)

  cat("\nnetMappingDEG (Dyslexia) result:\n")
  cat(sprintf("  Rows    : %d\n", nrow(net_Dyslexia)))
  cat(sprintf("  Columns : %s\n", paste(colnames(net_Dyslexia), collapse = ", ")))
  print(head(net_Dyslexia, 20))

  write.csv(data.frame(net_Dyslexia),
            "results.netMappingDEG.ligand.receptor.net.Dyslexia.csv",
            row.names = FALSE)
  cat("✓ Saved: results.netMappingDEG.ligand.receptor.net.Dyslexia.csv\n")

  # ══════════════════════════════════════════════════════════════
  # FILTERING SET 1 — loose threshold logFC = 0.05
  # ══════════════════════════════════════════════════════════════
  cat("\n--- FILTERING SET 1: logFC threshold = 0.05 ---\n")

  net.up <- subsetCommunication(cellchat,
                                 net            = net_Dyslexia,
                                 datasets       = "Dyslexia",
                                 ligand.logFC   = 0.05,
                                 receptor.logFC = 0.05)

  net.down <- subsetCommunication(cellchat,
                                   net            = net_Dyslexia,
                                   datasets       = "AD",
                                   ligand.logFC   = -0.05,
                                   receptor.logFC = -0.05)

  cat(sprintf("\nnet.up   (logFC > 0.05, higher in Dyslexia): %d L-R pairs\n", nrow(net.up)))
  cat(sprintf("net.down (logFC < -0.05, lower in Dyslexia): %d L-R pairs\n", nrow(net.down)))

  cat("\n===== net.up — upregulated in Dyslexia (logFC > 0.05) =====\n")
  print(net.up)

  cat("\n===== net.down — downregulated in Dyslexia (logFC < -0.05) =====\n")
  print(net.down)

  write.csv(net.up,
            "results.netMappingDEG_upregulated_Dyslexia.datasets_Dyslexia.ligandlogFC_0.05.receptorlogFC_0.05.csv",
            row.names = FALSE)
  write.csv(net.down,
            "results.netMappingDEG_downregulated_Dyslexia.datasets_AD.ligandlogFC_-0.05.receptorlogFC_-0.05.csv",
            row.names = FALSE)
  cat("✓ Saved: results.netMappingDEG_upregulated_Dyslexia.datasets_Dyslexia.ligandlogFC_0.05.receptorlogFC_0.05.csv\n")
  cat("✓ Saved: results.netMappingDEG_downregulated_Dyslexia.datasets_AD.ligandlogFC_-0.05.receptorlogFC_-0.05.csv\n")

  # ── Extract genes — set 1 ─────────────────────────────────────
  gene.up   <- extractGeneSubsetFromPair(net.up,   cellchat)
  gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

  cat(sprintf("\nGenes upregulated   in Dyslexia (logFC > 0.05): %d\n", length(gene.up)))
  cat(sprintf("Genes downregulated in Dyslexia (logFC < -0.05): %d\n", length(gene.down)))

  cat("\n===== gene.up (logFC > 0.05) =====\n");   print(gene.up)
  cat("\n===== gene.down (logFC < -0.05) =====\n"); print(gene.down)

  write.csv(gene.up,
            "results.genes_upregulated_Dyslexia.logFC_0.05.csv",   row.names = FALSE)
  write.csv(gene.down,
            "results.genes_downregulated_Dyslexia.logFC_0.05.csv", row.names = FALSE)
  cat("✓ Saved: results.genes_upregulated_Dyslexia.logFC_0.05.csv\n")
  cat("✓ Saved: results.genes_downregulated_Dyslexia.logFC_0.05.csv\n")

  # ══════════════════════════════════════════════════════════════
  # FILTERING SET 2 — strict threshold logFC = 0.25
  # ══════════════════════════════════════════════════════════════
  cat("\n--- FILTERING SET 2: logFC threshold = 0.25 ---\n")

  net.up.strict <- subsetCommunication(cellchat,
                                        net            = net_Dyslexia,
                                        datasets       = "Dyslexia",
                                        ligand.logFC   = 0.25,
                                        receptor.logFC = 0.25)

  net.down.strict <- subsetCommunication(cellchat,
                                          net            = net_Dyslexia,
                                          datasets       = "AD",
                                          ligand.logFC   = -0.25,
                                          receptor.logFC = -0.25)

  cat(sprintf("\nnet.up.strict   (logFC > 0.25, higher in Dyslexia): %d L-R pairs\n", nrow(net.up.strict)))
  cat(sprintf("net.down.strict (logFC < -0.25, lower in Dyslexia): %d L-R pairs\n", nrow(net.down.strict)))

  cat("\n===== net.up.strict — upregulated in Dyslexia (logFC > 0.25) =====\n")
  print(net.up.strict)

  cat("\n===== net.down.strict — downregulated in Dyslexia (logFC < -0.25) =====\n")
  print(net.down.strict)

  write.csv(net.up.strict,
            "results.netMappingDEG_upregulated_Dyslexia.datasets_Dyslexia.ligandlogFC_0.25.receptorlogFC_0.25.csv",
            row.names = FALSE)
  write.csv(net.down.strict,
            "results.netMappingDEG_downregulated_Dyslexia.datasets_AD.ligandlogFC_-0.25.receptorlogFC_-0.25.csv",
            row.names = FALSE)
  cat("✓ Saved: results.netMappingDEG_upregulated_Dyslexia.datasets_Dyslexia.ligandlogFC_0.25.receptorlogFC_0.25.csv\n")
  cat("✓ Saved: results.netMappingDEG_downregulated_Dyslexia.datasets_AD.ligandlogFC_-0.25.receptorlogFC_-0.25.csv\n")

  # ── Extract genes — set 2 ─────────────────────────────────────
  gene.up.strict   <- extractGeneSubsetFromPair(net.up.strict,   cellchat)
  gene.down.strict <- extractGeneSubsetFromPair(net.down.strict, cellchat)

  cat(sprintf("\nGenes upregulated   in Dyslexia (logFC > 0.25): %d\n", length(gene.up.strict)))
  cat(sprintf("Genes downregulated in Dyslexia (logFC < -0.25): %d\n", length(gene.down.strict)))

  cat("\n===== gene.up.strict (logFC > 0.25) =====\n");   print(gene.up.strict)
  cat("\n===== gene.down.strict (logFC < -0.25) =====\n"); print(gene.down.strict)

  write.csv(gene.up.strict,
            "results.genes_upregulated_Dyslexia.logFC_0.25.csv",   row.names = FALSE)
  write.csv(gene.down.strict,
            "results.genes_downregulated_Dyslexia.logFC_0.25.csv", row.names = FALSE)
  cat("✓ Saved: results.genes_upregulated_Dyslexia.logFC_0.25.csv\n")
  cat("✓ Saved: results.genes_downregulated_Dyslexia.logFC_0.25.csv\n")

  # ══════════════════════════════════════════════════════════════
  # STRICT FILTERING — logFC=0.25 + pval<0.05 + pct>=0.10
  # ══════════════════════════════════════════════════════════════
  p_cut   <- 0.05
  lfc_cut <- 0.25
  pct_cut <- 0.10

  cat(sprintf("\nStrict cutoffs: p < %.2f | lfc > %.2f | pct >= %.2f\n",
              p_cut, lfc_cut, pct_cut))

  Dyslexia_LR_bothUp <- subset(net_Dyslexia,
    datasets         == "Dyslexia" &
    ligand.pvalues   <  p_cut      &
    receptor.pvalues <  p_cut      &
    ligand.logFC     >  lfc_cut    &
    receptor.logFC   >  lfc_cut    &
    ligand.pct.1     >= pct_cut    &
    receptor.pct.1   >= pct_cut)

  Dyslexia_LR_bothDown <- subset(net_Dyslexia,
    datasets         == "Dyslexia" &
    ligand.pvalues   <  p_cut      &
    receptor.pvalues <  p_cut      &
    ligand.logFC     < -lfc_cut    &
    receptor.logFC   < -lfc_cut    &
    ligand.pct.1     >= pct_cut    &
    receptor.pct.1   >= pct_cut)

  cat(sprintf("\nDyslexia_LR_bothUp  : %d pairs\n", nrow(Dyslexia_LR_bothUp)))
  cat(sprintf("Dyslexia_LR_bothDown: %d pairs\n", nrow(Dyslexia_LR_bothDown)))

  cat("\n===== Dyslexia_LR_bothUp =====\n")
  print(Dyslexia_LR_bothUp)
  cat("\nDyslexia_LR_bothUp — source x target:\n")
  print(table(Dyslexia_LR_bothUp$source, Dyslexia_LR_bothUp$target))

  cat("\n===== Dyslexia_LR_bothDown =====\n")
  print(Dyslexia_LR_bothDown)
  cat("\nDyslexia_LR_bothDown — source x target:\n")
  print(table(Dyslexia_LR_bothDown$source, Dyslexia_LR_bothDown$target))

  write.csv(Dyslexia_LR_bothUp,
            "results.Dyslexia_LR_bothUp.datasets_Dyslexia.ligandlogFC_0.25.receptorlogFC_0.25.pval_0.05.pct_0.10.csv",
            row.names = FALSE)
  write.csv(Dyslexia_LR_bothDown,
            "results.Dyslexia_LR_bothDown.datasets_Dyslexia.ligandlogFC_-0.25.receptorlogFC_-0.25.pval_0.05.pct_0.10.csv",
            row.names = FALSE)

  write.table(Dyslexia_LR_bothUp$interaction_name,
              "results.Dyslexia_LR_bothUp.datasets_Dyslexia.ligandlogFC_0.25.receptorlogFC_0.25.pval_0.05.pct_0.10_interactions.txt",
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(Dyslexia_LR_bothDown$interaction_name,
              "results.Dyslexia_LR_bothDown.datasets_Dyslexia.ligandlogFC_-0.25.receptorlogFC_-0.25.pval_0.05.pct_0.10_interactions.txt",
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  cat("✓ All Dyslexia strict-filtered tables saved\n")

  # ── Bubble plots — SET 1 (logFC 0.05) ────────────────────────
  cat("\n--- Bubble plots SET 1 (logFC 0.05) ---\n")
  pairLR.use.up   <- net.up[,   "interaction_name", drop = FALSE]
  pairLR.use.down <- net.down[, "interaction_name", drop = FALSE]

  tryCatch({
    gg_up <- netVisual_bubble(cellchat,
                               pairLR.use     = pairLR.use.up,
                               comparison     = c(1, 2),
                               angle.x        = 90,
                               remove.isolate = TRUE,
                               title.name     = "Up-regulated in Dyslexia vs AD (logFC > 0.05)")
    ggsave("results.pairLR_upregulated_Dyslexia.logFC_0.05_bubble.png",
           gg_up, width = 48, height = 40, dpi = 300)
    cat("✓ Saved: results.pairLR_upregulated_Dyslexia.logFC_0.05_bubble.png\n")
  }, error = function(e) {
    cat("⚠️  Up bubble plot (set1) failed:", as.character(e), "\n")
  })

  tryCatch({
    gg_down <- netVisual_bubble(cellchat,
                                 pairLR.use     = pairLR.use.down,
                                 comparison     = c(1, 2),
                                 angle.x        = 90,
                                 remove.isolate = TRUE,
                                 title.name     = "Down-regulated in Dyslexia vs AD (logFC < -0.05)")
    ggsave("results.pairLR_downregulated_Dyslexia.logFC_0.05_bubble.png",
           gg_down, width = 48, height = 40, dpi = 300)
    cat("✓ Saved: results.pairLR_downregulated_Dyslexia.logFC_0.05_bubble.png\n")
  }, error = function(e) {
    cat("⚠️  Down bubble plot (set1) failed:", as.character(e), "\n")
  })

  # ── Bubble plots — SET 2 (logFC 0.25) ────────────────────────
  cat("\n--- Bubble plots SET 2 (logFC 0.25) ---\n")
  pairLR.use.up.strict   <- net.up.strict[,   "interaction_name", drop = FALSE]
  pairLR.use.down.strict <- net.down.strict[, "interaction_name", drop = FALSE]

  tryCatch({
    gg_up_strict <- netVisual_bubble(cellchat,
                                      pairLR.use     = pairLR.use.up.strict,
                                      comparison     = c(1, 2),
                                      angle.x        = 90,
                                      remove.isolate = TRUE,
                                      title.name     = "Up-regulated in Dyslexia vs AD (logFC > 0.25)")
    ggsave("results.pairLR_upregulated_Dyslexia.logFC_0.25_bubble.png",
           gg_up_strict, width = 48, height = 40, dpi = 300)
    cat("✓ Saved: results.pairLR_upregulated_Dyslexia.logFC_0.25_bubble.png\n")
  }, error = function(e) {
    cat("⚠️  Up bubble plot (set2) failed:", as.character(e), "\n")
  })

  tryCatch({
    gg_down_strict <- netVisual_bubble(cellchat,
                                        pairLR.use     = pairLR.use.down.strict,
                                        comparison     = c(1, 2),
                                        angle.x        = 90,
                                        remove.isolate = TRUE,
                                        title.name     = "Down-regulated in Dyslexia vs AD (logFC < -0.25)")
    ggsave("results.pairLR_downregulated_Dyslexia.logFC_0.25_bubble.png",
           gg_down_strict, width = 48, height = 40, dpi = 300)
    cat("✓ Saved: results.pairLR_downregulated_Dyslexia.logFC_0.25_bubble.png\n")
  }, error = function(e) {
    cat("⚠️  Down bubble plot (set2) failed:", as.character(e), "\n")
  })

  # ── Enrichment scores — SET 1 (logFC 0.05) ───────────────────
  cat("\n--- Enrichment scores SET 1 (logFC 0.05) ---\n")
  tryCatch({
    png("results.enrichmentScore_upregulated_Dyslexia.logFC_0.05.png",
        width = 7, height = 6, units = "in", res = 300)
    computeEnrichmentScore(net.up, species = "human", variable.both = TRUE)
    dev.off()
    cat("✓ Saved: results.enrichmentScore_upregulated_Dyslexia.logFC_0.05.png\n")
  }, error = function(e) {
    if (dev.cur() > 1) dev.off()
    cat("⚠️  Enrichment score (up set1) failed:", as.character(e), "\n")
  })

  tryCatch({
    png("results.enrichmentScore_downregulated_Dyslexia.logFC_0.05.png",
        width = 7, height = 6, units = "in", res = 300)
    computeEnrichmentScore(net.down, species = "human", variable.both = TRUE)
    dev.off()
    cat("✓ Saved: results.enrichmentScore_downregulated_Dyslexia.logFC_0.05.png\n")
  }, error = function(e) {
    if (dev.cur() > 1) dev.off()
    cat("⚠️  Enrichment score (down set1) failed:", as.character(e), "\n")
  })

  # ── Enrichment scores — SET 2 (logFC 0.25) ───────────────────
  cat("\n--- Enrichment scores SET 2 (logFC 0.25) ---\n")
  tryCatch({
    png("results.enrichmentScore_upregulated_Dyslexia.logFC_0.25.png",
        width = 7, height = 6, units = "in", res = 300)
    computeEnrichmentScore(net.up.strict, species = "human", variable.both = TRUE)
    dev.off()
    cat("✓ Saved: results.enrichmentScore_upregulated_Dyslexia.logFC_0.25.png\n")
  }, error = function(e) {
    if (dev.cur() > 1) dev.off()
    cat("⚠️  Enrichment score (up set2) failed:", as.character(e), "\n")
  })

  tryCatch({
    png("results.enrichmentScore_downregulated_Dyslexia.logFC_0.25.png",
        width = 7, height = 6, units = "in", res = 300)
    computeEnrichmentScore(net.down.strict, species = "human", variable.both = TRUE)
    dev.off()
    cat("✓ Saved: results.enrichmentScore_downregulated_Dyslexia.logFC_0.25.png\n")
  }, error = function(e) {
    if (dev.cur() > 1) dev.off()
    cat("⚠️  Enrichment score (down set2) failed:", as.character(e), "\n")
  })

} # end if deg_success_Dyslexia

# ══════════════════════════════════════════════════════════════════
# RUN 2 — AD PERSPECTIVE (optional)
# ══════════════════════════════════════════════════════════════════

if (run_AD_perspective) {

  cat("\n--- RUN 2: pos.dataset = AD ---\n")
  cat("logFC > 0 = higher in AD vs Dyslexia\n")
  cat("logFC < 0 = lower  in AD vs Dyslexia\n\n")

  pos.dataset   <- "AD"
  features.name <- paste0(pos.dataset, ".merged")  # "AD.merged"

  tryCatch({

    # ── STEP 1: identifyOverExpressedGenes ─────────────────────
    cat("Running identifyOverExpressedGenes (AD)...\n")
    cellchat <- identifyOverExpressedGenes(
      cellchat,
      group.dataset = "datasets",
      pos.dataset   = pos.dataset,
      features.name = features.name,
      only.pos      = FALSE,
      thresh.pc     = 0.1,
      thresh.fc     = 0.05,
      thresh.p      = 0.05)

    OE_genes_AD2 <- cellchat@var.features[[features.name]]
    cat(sprintf("✓ OE genes in AD: %d\n", length(OE_genes_AD2)))

    cat("\n===== OE Genes in AD =====\n")
    cat(sprintf("Total OE genes: %d\n", length(OE_genes_AD2)))
    print(OE_genes_AD2)

    write.csv(data.frame(gene = OE_genes_AD2),
              "results.OE_genes_AD2.csv", row.names = FALSE)
    cat("✓ Saved: results.OE_genes_AD2.csv\n")

    # ── STEP 2: identifyOverExpressedInteractions ───────────────
    cat("\nRunning identifyOverExpressedInteractions (AD)...\n")
    cellchat <- identifyOverExpressedInteractions(
      cellchat,
      features.name = features.name,
      variable.both = TRUE,
      return.object = TRUE)

    LR_OE_AD2 <- cellchat@var.features[[paste0(features.name, ".LR")]]
    cat(sprintf("✓ OE LR interactions in AD: %d\n", length(LR_OE_AD2)))

    cat("\n===== OE LR Interactions in AD =====\n")
    cat(sprintf("Total OE LR interactions: %d\n", length(LR_OE_AD2)))
    cat(sprintf("(subset of %d OE genes)\n", length(OE_genes_AD2)))
    df_OE_LR_AD2 <- data.frame(rank        = seq_along(LR_OE_AD2),
                                    interaction = LR_OE_AD2)
    print(df_OE_LR_AD2, row.names = FALSE)

    write.csv(df_OE_LR_AD2,
              "results.overExpressed_LR_interactions_AD.csv",
              row.names = FALSE)
    write.table(LR_OE_AD2,
                "results.overExpressed_LR_interactions_AD.txt",
                quote = FALSE, row.names = FALSE, col.names = FALSE)
    cat("✓ Saved: results.overExpressed_LR_interactions_AD.csv\n")
    cat("✓ Saved: results.overExpressed_LR_interactions_AD.txt\n")

    deg_success_AD2 <- TRUE

  }, error = function(e) {
    cat("⚠️  RUN 2 failed:", as.character(e), "\n")
  })

  if (deg_success_AD2) {

    net_AD2 <- netMappingDEG(cellchat, features.name = features.name)

    cat("\nnetMappingDEG (AD) result:\n")
    cat(sprintf("  Rows    : %d\n", nrow(net_AD2)))
    cat(sprintf("  Columns : %s\n", paste(colnames(net_AD2), collapse = ", ")))
    print(head(net_AD2, 20))

    write.csv(data.frame(net_AD2),
              "results.netMappingDEG.ligand.receptor.net.AD.csv",
              row.names = FALSE)
    cat("✓ Saved: results.netMappingDEG.ligand.receptor.net.AD.csv\n")

    net.up.H <- subsetCommunication(cellchat,
                                     net            = net_AD2,
                                     datasets       = "AD",
                                     ligand.logFC   = 0.05,
                                     receptor.logFC = 0.05)

    net.down.H <- subsetCommunication(cellchat,
                                       net            = net_AD2,
                                       datasets       = "Dyslexia",
                                       ligand.logFC   = -0.05,
                                       receptor.logFC = -0.05)

    cat(sprintf("\nnet.up.H   (higher in AD)): %d L-R pairs\n", nrow(net.up.H)))
    cat(sprintf("net.down.H (lower  in AD)): %d L-R pairs\n", nrow(net.down.H)))

    cat("\n===== net.up.H — upregulated in AD =====\n")
    print(net.up.H)
    cat("\n===== net.down.H — downregulated in AD =====\n")
    print(net.down.H)

    write.csv(net.up.H,
              "results.netMappingDEG_upregulated_AD.csv",   row.names = FALSE)
    write.csv(net.down.H,
              "results.netMappingDEG_downregulated_AD.csv", row.names = FALSE)
    cat("✓ Saved: results.netMappingDEG_upregulated_AD.csv\n")
    cat("✓ Saved: results.netMappingDEG_downregulated_AD.csv\n")

    gene.up.H   <- extractGeneSubsetFromPair(net.up.H,   cellchat)
    gene.down.H <- extractGeneSubsetFromPair(net.down.H, cellchat)

    cat(sprintf("\nGenes upregulated   in AD: %d\n", length(gene.up.H)))
    cat(sprintf("Genes downregulated in AD: %d\n", length(gene.down.H)))

    cat("\n===== gene.up.H =====\n");   print(gene.up.H)
    cat("\n===== gene.down.H =====\n"); print(gene.down.H)

    write.csv(gene.up.H,
              "results.genes_upregulated_AD.csv",   row.names = FALSE)
    write.csv(gene.down.H,
              "results.genes_downregulated_AD.csv", row.names = FALSE)

    p_cut   <- 0.05
    lfc_cut <- 0.25
    pct_cut <- 0.10

    AD_LR_bothUp2 <- subset(net_AD2,
      datasets         == "AD"  &
      ligand.pvalues   <  p_cut      &
      receptor.pvalues <  p_cut      &
      ligand.logFC     >  lfc_cut    &
      receptor.logFC   >  lfc_cut    &
      ligand.pct.1     >= pct_cut    &
      receptor.pct.1   >= pct_cut)

    AD_LR_bothDown2 <- subset(net_AD2,
      datasets         == "AD"  &
      ligand.pvalues   <  p_cut      &
      receptor.pvalues <  p_cut      &
      ligand.logFC     < -lfc_cut    &
      receptor.logFC   < -lfc_cut    &
      ligand.pct.1     >= pct_cut    &
      receptor.pct.1   >= pct_cut)

    cat(sprintf("\nAD_LR_bothUp2  : %d pairs\n", nrow(AD_LR_bothUp2)))
    cat(sprintf("AD_LR_bothDown2: %d pairs\n", nrow(AD_LR_bothDown2)))

    cat("\n===== AD_LR_bothUp2 =====\n")
    print(AD_LR_bothUp2)
    cat("\nAD_LR_bothUp2 — source x target:\n")
    print(table(AD_LR_bothUp2$source, AD_LR_bothUp2$target))

    cat("\n===== AD_LR_bothDown2 =====\n")
    print(AD_LR_bothDown2)
    cat("\nAD_LR_bothDown2 — source x target:\n")
    print(table(AD_LR_bothDown2$source, AD_LR_bothDown2$target))

    write.csv(AD_LR_bothUp2,
              "results.AD_LR_bothUp2.csv",   row.names = FALSE)
    write.csv(AD_LR_bothDown2,
              "results.AD_LR_bothDown2.csv", row.names = FALSE)

    write.table(AD_LR_bothUp2$interaction_name,
                "results.AD_LR_bothUp2_interactions.txt",
                quote = FALSE, row.names = FALSE, col.names = FALSE)
    write.table(AD_LR_bothDown2$interaction_name,
                "results.AD_LR_bothDown2_interactions.txt",
                quote = FALSE, row.names = FALSE, col.names = FALSE)
    cat("✓ All AD strict-filtered tables saved\n")

    pairLR.use.up.H   <- net.up.H[,   "interaction_name", drop = FALSE]
    pairLR.use.down.H <- net.down.H[, "interaction_name", drop = FALSE]

    tryCatch({
      gg_up_H <- netVisual_bubble(cellchat,
                                   pairLR.use     = pairLR.use.up.H,
                                   comparison     = c(1, 2),
                                   angle.x        = 90,
                                   remove.isolate = TRUE,
                                   title.name     = "Up-regulated in AD vs Dyslexia")
      ggsave("results.pairLR_upregulated_AD_bubble.png",
             gg_up_H, width = 48, height = 40, dpi = 300)
      cat("✓ Saved: results.pairLR_upregulated_AD_bubble.png\n")
    }, error = function(e) {
      cat("⚠️  Up bubble plot (AD) failed:", as.character(e), "\n")
    })

    tryCatch({
      gg_down_H <- netVisual_bubble(cellchat,
                                     pairLR.use     = pairLR.use.down.H,
                                     comparison     = c(1, 2),
                                     angle.x        = 90,
                                     remove.isolate = TRUE,
                                     title.name     = "Down-regulated in AD vs Dyslexia")
      ggsave("results.pairLR_downregulated_AD_bubble.png",
             gg_down_H, width = 48, height = 40, dpi = 300)
      cat("✓ Saved: results.pairLR_downregulated_AD_bubble.png\n")
    }, error = function(e) {
      cat("⚠️  Down bubble plot (AD) failed:", as.character(e), "\n")
    })

    tryCatch({
      png("results.enrichmentScore_upregulated_AD.png",
          width = 7, height = 6, units = "in", res = 300)
      computeEnrichmentScore(net.up.H, species = "human", variable.both = TRUE)
      dev.off()
      cat("✓ Saved: results.enrichmentScore_upregulated_AD.png\n")
    }, error = function(e) {
      if (dev.cur() > 1) dev.off()
      cat("⚠️  Enrichment score (AD up) failed:", as.character(e), "\n")
    })

    tryCatch({
      png("results.enrichmentScore_downregulated_AD.png",
          width = 7, height = 6, units = "in", res = 300)
      computeEnrichmentScore(net.down.H, species = "human", variable.both = TRUE)
      dev.off()
      cat("✓ Saved: results.enrichmentScore_downregulated_AD.png\n")
    }, error = function(e) {
      if (dev.cur() > 1) dev.off()
      cat("⚠️  Enrichment score (AD down) failed:", as.character(e), "\n")
    })

  } # end if deg_success_AD2

} else {
  cat("\n--- RUN 2 (AD perspective) SKIPPED ---\n")
  cat("    set run_AD_perspective <- TRUE to enable\n")
}

# ══════════════════════════════════════════════════════════════════
# PATHWAY-SPECIFIC VISUALIZATIONS
# ══════════════════════════════════════════════════════════════════

cat("\n=== PATHWAY-SPECIFIC VISUALIZATIONS ===\n")

generate_pathway_plots <- function(pathway_list, output_label,
                                   object.list, cellchat,
                                   output_dir = ".") {
  cat(sprintf("\n--- %s (%d pathways) ---\n",
              output_label, length(pathway_list)))

  category_dir <- file.path(output_dir, output_label)
  dir.create(category_dir, showWarnings = FALSE, recursive = TRUE)

  cellchat@meta$datasets <- factor(cellchat@meta$datasets,
                                    levels = c("AD", "Dyslexia"))
  success_count <- 0

  for (pathway in pathway_list) {

    tryCatch({
      weight.max <- getMaxWeight(object.list, slot.name = "netP",
                                  attribute = pathway)
      png(file.path(category_dir,
                    sprintf("%s_%s_circle.png", pathway, output_label)),
          width = 14, height = 7, units = "in", res = 300)
      par(mfrow = c(1,2), xpd = TRUE, mar = c(2,2,4,2))
      for (i in seq_along(object.list)) {
        netVisual_aggregate(object.list[[i]], signaling = pathway,
                            layout = "circle",
                            edge.weight.max = weight.max[1],
                            edge.width.max  = 10,
                            signaling.name  = paste(pathway,
                                                    names(object.list)[i]))
      }
      dev.off()
    }, error = function(e) { if (dev.cur() > 1) dev.off() })

    tryCatch({
      ht <- list()
      for (i in seq_along(object.list)) {
        ht[[i]] <- netVisual_heatmap(object.list[[i]],
                                      signaling     = pathway,
                                      color.heatmap = "Reds",
                                      title.name    = paste(pathway,
                                                            names(object.list)[i]))
      }
      for (i in seq_along(object.list)) {
        png(file.path(category_dir,
                      sprintf("%s_%s_heatmap_%s.png",
                              pathway, output_label, names(object.list)[i])),
            width = 6, height = 6, units = "in", res = 300)
        ComplexHeatmap::draw(ht[[i]]); dev.off()
      }
    }, error = function(e) { if (dev.cur() > 1) dev.off() })

    tryCatch({
      png(file.path(category_dir,
                    sprintf("%s_%s_chord.png", pathway, output_label)),
          width = 14, height = 7, units = "in", res = 300)
      par(mfrow = c(1,2), xpd = TRUE, mar = c(2,2,4,2))
      for (i in seq_along(object.list)) {
        netVisual_aggregate(object.list[[i]], signaling = pathway,
                            layout = "chord",
                            signaling.name = paste(pathway, names(object.list)[i]))
      }
      dev.off()
    }, error = function(e) { if (dev.cur() > 1) dev.off() })

    tryCatch({
      p <- plotGeneExpression(cellchat, signaling = pathway,
                               split.by = "datasets",
                               colors.ggplot = TRUE, type = "violin")
      ggsave(file.path(category_dir,
                       sprintf("%s_%s_violin.png", pathway, output_label)),
             p, width = 7, height = 5, dpi = 300)
    }, error = function(e) invisible(NULL))

    success_count <- success_count + 1
    if (success_count %% 5 == 0) gc(verbose = FALSE)
  }

  cat(sprintf("  ✓ %d / %d pathways processed\n",
              success_count, length(pathway_list)))
}

detected_pathways <- unique(c(ad@netP$pathways, dyslexia@netP$pathways))
cat(sprintf("Total detected pathways: %d\n", length(detected_pathways)))

pathway_plots_dir <- "pathway_specific_plots"
dir.create(pathway_plots_dir, showWarnings = FALSE)

generate_pathway_plots(detected_pathways, "all_pathways",
                       object.list, cellchat, pathway_plots_dir)

if (deg_success_Dyslexia) {
  # ── SET 1 pathways (logFC 0.05) ───────────────────────────────
  if (nrow(net.up) > 0)
    generate_pathway_plots(unique(net.up$pathway_name),
                           "net_up_Dyslexia.logFC_0.05",
                           object.list, cellchat, pathway_plots_dir)
  if (nrow(net.down) > 0)
    generate_pathway_plots(unique(net.down$pathway_name),
                           "net_down_Dyslexia.logFC_0.05",
                           object.list, cellchat, pathway_plots_dir)
  # ── SET 2 pathways (logFC 0.25) ───────────────────────────────
  if (nrow(net.up.strict) > 0)
    generate_pathway_plots(unique(net.up.strict$pathway_name),
                           "net_up_Dyslexia.logFC_0.25",
                           object.list, cellchat, pathway_plots_dir)
  if (nrow(net.down.strict) > 0)
    generate_pathway_plots(unique(net.down.strict$pathway_name),
                           "net_down_Dyslexia.logFC_0.25",
                           object.list, cellchat, pathway_plots_dir)
  # ── STRICT filtered pathways ──────────────────────────────────
  if (nrow(Dyslexia_LR_bothUp) > 0)
    generate_pathway_plots(unique(Dyslexia_LR_bothUp$pathway_name),
                           "Dyslexia_LR_bothUp.logFC_0.25.pval_0.05.pct_0.10",
                           object.list, cellchat, pathway_plots_dir)
  if (nrow(Dyslexia_LR_bothDown) > 0)
    generate_pathway_plots(unique(Dyslexia_LR_bothDown$pathway_name),
                           "Dyslexia_LR_bothDown.logFC_0.25.pval_0.05.pct_0.10",
                           object.list, cellchat, pathway_plots_dir)
}

if (run_AD_perspective && deg_success_AD2) {
  if (nrow(net.up.H) > 0)
    generate_pathway_plots(unique(net.up.H$pathway_name),
                           "net_up_AD",
                           object.list, cellchat, pathway_plots_dir)
  if (nrow(net.down.H) > 0)
    generate_pathway_plots(unique(net.down.H$pathway_name),
                           "net_down_AD",
                           object.list, cellchat, pathway_plots_dir)
  if (nrow(AD_LR_bothUp2) > 0)
    generate_pathway_plots(unique(AD_LR_bothUp2$pathway_name),
                           "AD2_LR_bothUp",
                           object.list, cellchat, pathway_plots_dir)
  if (nrow(AD_LR_bothDown2) > 0)
    generate_pathway_plots(unique(AD_LR_bothDown2$pathway_name),
                           "AD2_LR_bothDown",
                           object.list, cellchat, pathway_plots_dir)
}

cat("✓ All pathway-specific plots complete\n")

# ══════════════════════════════════════════════════════════════════
# SAVE FINAL OBJECTS
# ══════════════════════════════════════════════════════════════════

cat("\n=== SAVING FINAL OBJECTS ===\n")
saveRDS(cellchat,    "cellchat_AD_vs_Dyslexia_merged.rds")
saveRDS(object.list, "cellchat_object.list_AD_Dyslexia.rds")
cat("✓ Saved: cellchat_AD_vs_Dyslexia_merged.rds\n")
cat("✓ Saved: cellchat_object.list_AD_Dyslexia.rds\n")

# ══════════════════════════════════════════════════════════════════
# SAVE DATA.SIGNALING MATRICES
# ══════════════════════════════════════════════════════════════════

cat("\n=== SAVING DATA.SIGNALING MATRICES ===\n")
mat_all <- cellchat@data.signaling
cond    <- cellchat@meta$datasets
names(cond) <- rownames(cellchat@meta)
cond_ds <- cond[colnames(mat_all)]

cells_AD  <- colnames(mat_all)[cond_ds == "AD"]
cells_Dyslexia <- colnames(mat_all)[cond_ds == "Dyslexia"]

cat(sprintf("Cells AD  : %d\n", length(cells_AD)))
cat(sprintf("Cells Dyslexia : %d\n", length(cells_Dyslexia)))

mat_AD  <- mat_all[, cells_AD,  drop = FALSE]
mat_Dyslexia <- mat_all[, cells_Dyslexia, drop = FALSE]

cat("\nTop genes by mean (AD):\n")
print(head(sort(Matrix::rowMeans(mat_AD),  decreasing = TRUE), 20))
cat("\nTop genes by mean (Dyslexia):\n")
print(head(sort(Matrix::rowMeans(mat_Dyslexia), decreasing = TRUE), 20))

saveRDS(mat_all,      "cellchat_data.signaling.rds")
saveRDS(mat_AD,  "cellchat_data.signaling_AD.rds")
saveRDS(mat_Dyslexia, "cellchat_data.signaling_Dyslexia.rds")

writeMM(mat_all, "cellchat_data.signaling.mtx")
write.table(rownames(mat_all),
            "cellchat_data.signaling_genes.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(colnames(mat_all),
            "cellchat_data.signaling_cells.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

saveRDS(list(data.signaling = mat_all,
             genes    = rownames(mat_all),
             cells    = colnames(mat_all),
             datasets = cellchat@meta$datasets,
             cond_ds  = cond_ds),
        "cellchat_data.signaling_with_metadata.rds")
cat("✓ All data.signaling matrices saved\n")

# ══════════════════════════════════════════════════════════════════
# WRAP UP
# ══════════════════════════════════════════════════════════════════

cat("\n")
cat("╔═══════════════════════════════════════════════════════════╗\n")
cat("║                                                           ║\n")
cat("║  ✓✓✓  AD vs DYSLEXIA COMPARISON COMPLETE  ✓✓✓      ║\n")
cat("║                                                           ║\n")
cat("╚═══════════════════════════════════════════════════════════╝\n")
cat("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Log saved:", log_file, "\n\n")

sink()
