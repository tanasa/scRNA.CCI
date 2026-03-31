# ══════════════════════════════════════════════════════════════════
# CELLCHAT COMPARISON ANALYSIS
# Input : 18_merged_AG_031826.28celltypes.Healthy.rds
#         18_merged_AG_031826.28celltypes.AD.rds
# Groups: Healthy vs AD
# pos.dataset = "AD"
# logFC > 0 means higher in AD vs Healthy
# logFC < 0 means lower in AD vs Healthy (higher in Healthy)
# ══════════════════════════════════════════════════════════════════

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

log_file <- paste0("cellchat_comparison_Healthy_vs_AD_",
                   format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt")
sink(log_file, split = TRUE)

cat("=================================================\n")
cat("CellChat Comparison Analysis Log\n")
cat("Healthy vs AD\n")
cat("pos.dataset = AD\n")
cat("logFC > 0 = higher in AD vs Healthy\n")
cat("logFC < 0 = lower in AD vs Healthy (higher in Healthy)\n")
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

work_dir <- file.path(base_dir, "compare_Healthy_vs_AD_cellchat")
dir.create(work_dir, showWarnings = FALSE, recursive = TRUE)
setwd(work_dir)
cat("Working directory:", getwd(), "\n\n")

# ══════════════════════════════════════════════════════════════════
# LOAD CELLCHAT OBJECTS
# ══════════════════════════════════════════════════════════════════

cat("=== LOADING CELLCHAT OBJECTS ===\n")

healthy <- readRDS(file.path(base_dir,
  "18_merged_AG_031826.28celltypes.Healthy_cellchat",
  "cellchat_18_merged_AG_031826.28celltypes.Healthy_at_the_end.rds"))
cat("✓ Loaded: Healthy\n")

cat("\nHealthy metadata:\n")
cat("  subcelltype:\n"); print(sort(unique(healthy@meta$subcelltype)))
cat("  group:\n");       print(table(healthy@meta$group))

ad <- readRDS(file.path(base_dir,
  "18_merged_AG_031826.28celltypes.AD_cellchat",
  "cellchat_18_merged_AG_031826.28celltypes.AD_at_the_end.rds"))
cat("✓ Loaded: AD\n")

cat("\nAD metadata:\n")
cat("  subcelltype:\n"); print(sort(unique(ad@meta$subcelltype)))
cat("  group:\n");       print(table(ad@meta$group))

mem_status <- check_memory_usage()
if (mem_status$should_stop) { sink(); stop("High memory") }

# ══════════════════════════════════════════════════════════════════
# UPDATE AND MERGE
# ══════════════════════════════════════════════════════════════════

cat("\n=== UPDATING CELLCHAT OBJECTS ===\n")
healthy <- updateCellChat(healthy)
ad      <- updateCellChat(ad)
cat("✓ Both objects updated\n")

# ══════════════════════════════════════════════════════════════════
# CELL TYPE COMPATIBILITY — subset to common cell types
# This avoids liftCellChat which creates empty cell types causing:
#   - non-conformable arrays in netVisual_heatmap
#   - Must have at least 2 groups in wilcoxauc/identifyOverExpressedGenes
# We keep only cell types present in BOTH conditions.
# Rare condition-specific cell types are excluded from comparison.
# ══════════════════════════════════════════════════════════════════

cat("\n=== CELL TYPE COMPATIBILITY CHECK ===\n")
ct_Healthy <- sort(levels(healthy@idents))
ct_AD      <- sort(levels(ad@idents))

cat(sprintf("Cell types in Healthy : %d\n", length(ct_Healthy)))
cat(sprintf("Cell types in AD      : %d\n", length(ct_AD)))

only_Healthy <- setdiff(ct_Healthy, ct_AD)
only_AD      <- setdiff(ct_AD,      ct_Healthy)
common_ct    <- intersect(ct_Healthy, ct_AD)

cat(sprintf("Only in Healthy  : %d : %s\n",
            length(only_Healthy), paste(only_Healthy, collapse = ", ")))
cat(sprintf("Only in AD       : %d : %s\n",
            length(only_AD),      paste(only_AD,      collapse = ", ")))
cat(sprintf("Common cell types: %d\n", length(common_ct)))
cat("Common cell types:\n")
print(sort(common_ct))

# Subset both objects to common cell types only
cat("\n✓ Subsetting both objects to common cell types...\n")
healthy_sub <- subsetCellChat(healthy, idents.use = common_ct)
ad_sub      <- subsetCellChat(ad,      idents.use = common_ct)

cat(sprintf("✓ Healthy after subset: %d cell types\n",
            length(levels(healthy_sub@idents))))
cat(sprintf("✓ AD      after subset: %d cell types\n",
            length(levels(ad_sub@idents))))

if (length(only_Healthy) > 0)
  cat(sprintf("⚠️  Excluded from Healthy (not in AD): %s\n",
              paste(only_Healthy, collapse = ", ")))
if (length(only_AD) > 0)
  cat(sprintf("⚠️  Excluded from AD (not in Healthy): %s\n",
              paste(only_AD, collapse = ", ")))

# Healthy first, AD second
# pos.dataset = "AD" → logFC > 0 means higher in AD vs Healthy
object.list <- list(Healthy = healthy_sub, AD = ad_sub)

cat("\n=== MERGING CELLCHAT OBJECTS ===\n")
mem_status <- check_memory_usage()
if (mem_status$should_stop) { sink(); stop("High memory") }

# No cell.prefix needed — same cell types in both objects
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
cat("✓ Merged CellChat object\n")

cat("\nMerged object overview:\n")
print(cellchat)
cat("\nDatasets:\n");    print(table(cellchat@meta$datasets))
cat("\nCell types:\n");  print(sort(unique(cellchat@meta$subcelltype)))

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
  cat("   Script continues...\n")
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
  cat("   Script continues...\n")
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

celltypes_Healthy <- sort(unique(healthy@meta$subcelltype))
celltypes_AD      <- sort(unique(ad@meta$subcelltype))
celltypes_common  <- intersect(celltypes_Healthy, celltypes_AD)

cat(sprintf("Cell types in Healthy : %d\n", length(celltypes_Healthy)))
cat(sprintf("Cell types in AD      : %d\n", length(celltypes_AD)))
cat(sprintf("Cell types in common  : %d\n", length(celltypes_common)))

dir.create("signalingChanges", showWarnings = FALSE)

for (ct in celltypes_common) {
  ct_safe <- gsub("[^A-Za-z0-9_]", "_", ct)
  tryCatch({
    p <- netAnalysis_signalingChanges_scatter(cellchat,
                                               idents.use = ct)
    ggsave(file.path("signalingChanges",
                     sprintf("%s_signalingChanges.png", ct_safe)),
           p, width = 8, height = 6, dpi = 300)
  }, error = function(e) {
    cat(sprintf("⚠️  signalingChanges failed for %s: %s\n",
                ct, as.character(e)))
  })
}
cat("✓ Signaling changes plots saved to: signalingChanges/\n")

# ══════════════════════════════════════════════════════════════════
# FUNCTIONAL SIMILARITY
# ══════════════════════════════════════════════════════════════════

cat("\n=== FUNCTIONAL SIMILARITY ===\n")
output_dir_functional <- "cellchat_functional_plots"
dir.create(output_dir_functional, showWarnings = FALSE)

mem_status <- check_memory_usage()
if (mem_status$should_stop) { sink(); stop("High memory") }

cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional",
                          umap.method = "uwot")
cellchat <- netClustering(cellchat, type = "functional",
                           do.parallel = FALSE)

png(file.path(output_dir_functional,
              "cellchat_functional_embedding.png"),
    width = 12, height = 10, units = "in", res = 300)
netVisual_embeddingPairwise(cellchat, type = "functional",
                             label.size = 3.5)
dev.off()

png(file.path(output_dir_functional,
              "cellchat_functional_embedding_zoomin.png"),
    width = 16, height = 8, units = "in", res = 300)
netVisual_embeddingPairwiseZoomIn(cellchat, type = "functional",
                                   nCol = 2)
dev.off()
cat("✓ Functional embedding plots saved\n")

M  <- cellchat@netP$similarity$functional$matrix[["1-2"]]
DR <- cellchat@netP$similarity$functional$dr[["1-2"]]

png(file.path(output_dir_functional,
              "functional_similarity_heatmap.png"),
    width = 14, height = 14, units = "in", res = 300)
pheatmap(M, show_rownames = TRUE, show_colnames = TRUE,
         fontsize_row = 5, fontsize_col = 5,
         main = "Pathway Functional Similarity (Healthy vs AD)")
dev.off()

df <- as.data.frame(DR)
colnames(df) <- c("UMAP1", "UMAP2")
df$pathway       <- rownames(DR)
df$condition     <- sub(".*--", "", df$pathway)
df$pathway_clean <- sub("--.*", "", df$pathway)

ctl_df <- df %>% filter(condition == "Healthy") %>%
  select(pathway = pathway_clean, UMAP1, UMAP2)
ad_df  <- df %>% filter(condition == "AD") %>%
  select(pathway = pathway_clean, UMAP1, UMAP2)

merged_func <- full_join(ctl_df, ad_df,
                          by = "pathway",
                          suffix = c("_Healthy", "_AD")) %>%
  mutate(shift = ifelse(
    !is.na(UMAP1_Healthy) & !is.na(UMAP1_AD),
    sqrt((UMAP1_Healthy - UMAP1_AD)^2 +
           (UMAP2_Healthy - UMAP2_AD)^2),
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
  scale_color_manual(values = c("Healthy" = "#4DAF4A",
                                 "AD"      = "#E41A1C")) +
  theme_classic(base_size = 13) +
  labs(title = sprintf("Functional Similarity - Top %d Shifted Pathways",
                       topN), color = "Condition") +
  coord_fixed()

ggsave(file.path(output_dir_functional,
                 "functional_umap_labeled.png"),
       p1, width = 12, height = 10, dpi = 300)

merged_both_func <- merged_func %>% filter(!is.na(shift))

p2 <- ggplot() +
  geom_segment(
    data = merged_both_func,
    aes(x = UMAP1_Healthy, y = UMAP2_Healthy,
        xend = UMAP1_AD, yend = UMAP2_AD),
    arrow = arrow(length = unit(0.12, "inches")),
    alpha = 0.5, color = "grey60") +
  geom_point(data = df, aes(UMAP1, UMAP2, color = condition),
             size = 2, alpha = 0.8) +
  geom_text_repel(
    data = merged_both_func[1:min(15, nrow(merged_both_func)), ],
    aes(x = UMAP1_AD, y = UMAP2_AD, label = pathway),
    size = 3.5, max.overlaps = Inf) +
  scale_color_manual(values = c("Healthy" = "#4DAF4A",
                                 "AD"      = "#E41A1C")) +
  theme_classic(base_size = 13) +
  labs(title = "Functional Shift (Healthy → AD)",
       color = "Condition") +
  coord_fixed()

ggsave(file.path(output_dir_functional,
                 "functional_pathway_shifts.png"),
       p2, width = 12, height = 10, dpi = 300)

write.csv(merged_func,
          file.path(output_dir_functional,
                    "functional_pathway_shifts_ranked.csv"),
          row.names = FALSE)
write.csv(M,
          file.path(output_dir_functional,
                    "functional_similarity_matrix.csv"),
          row.names = TRUE)
write.csv(df,
          file.path(output_dir_functional,
                    "functional_umap_coordinates.csv"),
          row.names = FALSE)
cat("✓ Functional similarity complete\n")

tryCatch({
  ranked_func <- rankSimilarity(cellchat, type = "functional",
                                 slot.name = "netP")
  png(file.path(output_dir_functional,
                "ranked_functional_similarity.png"),
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
cellchat <- netEmbedding(cellchat, type = "structural",
                          umap.method = "uwot")
cellchat <- netClustering(cellchat, type = "structural",
                           do.parallel = FALSE)

png(file.path(output_dir_structural,
              "structural_cellchat_embedding.png"),
    width = 12, height = 10, units = "in", res = 300)
netVisual_embeddingPairwise(cellchat, type = "structural",
                             label.size = 3.5)
dev.off()

png(file.path(output_dir_structural,
              "structural_cellchat_embedding_zoomin.png"),
    width = 16, height = 8, units = "in", res = 300)
netVisual_embeddingPairwiseZoomIn(cellchat, type = "structural",
                                   nCol = 2)
dev.off()
cat("✓ Structural embedding plots saved\n")

M_struct  <- cellchat@netP$similarity$structural$matrix[["1-2"]]
DR_struct <- cellchat@netP$similarity$structural$dr[["1-2"]]

png(file.path(output_dir_structural, "similarity_heatmap.png"),
    width = 14, height = 14, units = "in", res = 300)
pheatmap(M_struct, show_rownames = TRUE, show_colnames = TRUE,
         fontsize_row = 5, fontsize_col = 5,
         main = "Structural Similarity (Healthy vs AD)")
dev.off()

df_s <- as.data.frame(DR_struct)
colnames(df_s) <- c("UMAP1", "UMAP2")
df_s$pathway       <- rownames(DR_struct)
df_s$condition     <- sub(".*--", "", df_s$pathway)
df_s$pathway_clean <- sub("--.*", "", df_s$pathway)

ctl_s <- df_s %>% filter(condition == "Healthy") %>%
  select(pathway = pathway_clean, UMAP1, UMAP2)
ad_s  <- df_s %>% filter(condition == "AD") %>%
  select(pathway = pathway_clean, UMAP1, UMAP2)

merged_struct <- full_join(ctl_s, ad_s,
                            by = "pathway",
                            suffix = c("_Healthy", "_AD")) %>%
  mutate(shift = ifelse(
    !is.na(UMAP1_Healthy) & !is.na(UMAP1_AD),
    sqrt((UMAP1_Healthy - UMAP1_AD)^2 +
           (UMAP2_Healthy - UMAP2_AD)^2),
    NA_real_)) %>%
  arrange(desc(shift))

topN_s    <- min(500, nrow(merged_struct))
to_label_s <- merged_struct$pathway[1:topN_s]

p1_s <- ggplot(df_s, aes(UMAP1, UMAP2, color = condition)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_text_repel(
    data = df_s %>% filter(pathway_clean %in% to_label_s),
    aes(label = pathway_clean),
    size = 3.5, max.overlaps = Inf, box.padding = 0.5) +
  scale_color_manual(values = c("Healthy" = "#4DAF4A",
                                 "AD"      = "#E41A1C")) +
  theme_classic(base_size = 13) +
  labs(title = sprintf("Structural Similarity - Top %d Shifted Pathways",
                       topN_s), color = "Condition") +
  coord_fixed()

ggsave(file.path(output_dir_structural, "umap_labeled.png"),
       p1_s, width = 12, height = 10, dpi = 300)

merged_both_s <- merged_struct %>% filter(!is.na(shift))

p2_s <- ggplot() +
  geom_segment(
    data = merged_both_s,
    aes(x = UMAP1_Healthy, y = UMAP2_Healthy,
        xend = UMAP1_AD, yend = UMAP2_AD),
    arrow = arrow(length = unit(0.12, "inches")),
    alpha = 0.5, color = "grey60") +
  geom_point(data = df_s, aes(UMAP1, UMAP2, color = condition),
             size = 2, alpha = 0.8) +
  geom_text_repel(
    data = merged_both_s[1:min(15, nrow(merged_both_s)), ],
    aes(x = UMAP1_AD, y = UMAP2_AD, label = pathway),
    size = 3.5, max.overlaps = Inf) +
  scale_color_manual(values = c("Healthy" = "#4DAF4A",
                                 "AD"      = "#E41A1C")) +
  theme_classic(base_size = 13) +
  labs(title = "Structural Shift (Healthy → AD)",
       color = "Condition") +
  coord_fixed()

ggsave(file.path(output_dir_structural, "pathway_shifts.png"),
       p2_s, width = 12, height = 10, dpi = 300)

write.csv(merged_struct,
          file.path(output_dir_structural,
                    "structural_pathway_shifts_ranked.csv"),
          row.names = FALSE)
write.csv(M_struct,
          file.path(output_dir_structural,
                    "structural_similarity_matrix.csv"),
          row.names = TRUE)
write.csv(df_s,
          file.path(output_dir_structural,
                    "structural_umap_coordinates.csv"),
          row.names = FALSE)
cat("✓ Structural similarity complete\n")

tryCatch({
  ranked_struct <- rankSimilarity(cellchat, type = "structural",
                                   slot.name = "netP")
  png(file.path(output_dir_structural,
                "ranked_structural_similarity.png"),
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
  list(slot = "netP", stacked = TRUE,
       file = "ranknet_weight_stacked.pdf",       h = 8),
  list(slot = "netP", stacked = FALSE,
       file = "ranknet_weight.pdf",               h = 8),
  list(slot = "net",  stacked = TRUE,
       file = "ranknet_weight_net_stacked.pdf",   h = 40),
  list(slot = "net",  stacked = FALSE,
       file = "ranknet_weight_net.pdf",           h = 40)
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

# All interactions
tryCatch({
  pdf("netVisual_bubble_comparison_all.pdf", width = 40, height = 40)
  netVisual_bubble(cellchat, comparison = c(1, 2), angle.x = 45)
  dev.off()
  cat("✓ Bubble plot (all) saved\n")
}, error = function(e) {
  if (dev.cur() > 1) dev.off()
  cat("⚠️  Bubble plot (all) failed:", as.character(e), "\n")
})

# Increased in AD
tryCatch({
  png("netVisual_bubble_increased_in_AD.png",
      width = 40, height = 40, units = "in", res = 200)
  netVisual_bubble(cellchat, comparison = c(1, 2),
                   max.dataset    = 2,
                   title.name     = "Increased signaling in AD",
                   angle.x        = 45,
                   remove.isolate = TRUE)
  dev.off()
  cat("✓ Bubble plot (increased in AD) saved\n")
}, error = function(e) {
  if (dev.cur() > 1) dev.off()
  cat("⚠️  Bubble plot (increased AD) failed:", as.character(e), "\n")
})

# Decreased in AD (increased in Healthy)
tryCatch({
  png("netVisual_bubble_decreased_in_AD.png",
      width = 40, height = 40, units = "in", res = 200)
  netVisual_bubble(cellchat, comparison = c(1, 2),
                   max.dataset    = 1,
                   title.name     = "Decreased signaling in AD (higher in Healthy)",
                   angle.x        = 45,
                   remove.isolate = TRUE)
  dev.off()
  cat("✓ Bubble plot (decreased in AD) saved\n")
}, error = function(e) {
  if (dev.cur() > 1) dev.off()
  cat("⚠️  Bubble plot (decreased AD) failed:", as.character(e), "\n")
})

# Cell-type-specific bubble plots — COMMENTED OUT (fill in subcelltype names)
# Run: unique(healthy@meta$subcelltype) to get the names first
#
# tryCatch({
#   exc_types <- c("...")  # your excitatory subcelltype names
#   png("netVisual_bubble_excitatory.png",
#       width = 40, height = 40, units = "in", res = 200)
#   netVisual_bubble(cellchat,
#                    sources.use = exc_types,
#                    targets.use = exc_types,
#                    comparison  = c(1, 2), angle.x = 45)
#   dev.off()
# }, error = function(e) {
#   if (dev.cur() > 1) dev.off()
# })
#
# tryCatch({
#   inh_types <- c("...")  # your inhibitory subcelltype names
#   png("netVisual_bubble_inhibitory.png",
#       width = 40, height = 40, units = "in", res = 200)
#   netVisual_bubble(cellchat,
#                    sources.use = inh_types,
#                    targets.use = inh_types,
#                    comparison  = c(1, 2), angle.x = 45)
#   dev.off()
# }, error = function(e) {
#   if (dev.cur() > 1) dev.off()
# })

# ══════════════════════════════════════════════════════════════════
# DIFFERENTIAL EXPRESSION — netMappingDEG
# ══════════════════════════════════════════════════════════════════

cat("\n=== DIFFERENTIAL EXPRESSION ANALYSIS ===\n")

# ── Set datasets as factor ────────────────────────────────────────
# No cell.prefix was used so datasets column is clean
cat("--- Datasets check ---\n")
print(table(cellchat@meta$datasets, useNA = "ifany"))
cellchat@meta$datasets <- factor(as.character(cellchat@meta$datasets),
                                  levels = c("Healthy", "AD"))
cat("✓ datasets set as factor: Healthy, AD\n")

pos.dataset   <- "AD"
features.name <- "AD.merged"

cat(sprintf("\nPositive dataset : %s\n", pos.dataset))
cat(sprintf("Features name    : %s\n", features.name))
cat("logFC > 0 = higher in AD vs Healthy\n")
cat("logFC < 0 = lower in AD (higher in Healthy)\n\n")

deg_success <- FALSE

tryCatch({
  cat("Running identifyOverExpressedGenes...\n")
  cellchat <- identifyOverExpressedGenes(
    cellchat,
    group.dataset = "datasets",
    pos.dataset   = pos.dataset,
    features.name = features.name,
    only.pos      = FALSE,
    thresh.pc     = 0.1,
    thresh.fc     = 0.05,
    thresh.p      = 0.05)
  deg_success <- TRUE
  cat(sprintf("✓ Overexpressed genes in AD: %d\n",
              length(cellchat@var.features[[features.name]])))
}, error = function(e) {
  cat("⚠️  identifyOverExpressedGenes failed:", as.character(e), "\n")
})

if (!deg_success) {
  cat("❌ DEG failed — skipping netMappingDEG section\n")
} else {

net <- netMappingDEG(cellchat, features.name = features.name)

cat("\nnetMappingDEG result:\n")
cat(sprintf("  Rows    : %d\n", nrow(net)))
cat(sprintf("  Columns : %s\n", paste(colnames(net), collapse = ", ")))

write.csv(data.frame(net),
          "results.netMappingDEG.ligand.receptor.net.csv",
          row.names = FALSE)
cat("✓ Full netMappingDEG table saved\n")

# ── Subsets ───────────────────────────────────────────────────────
# Higher in AD vs Healthy
net.up <- subsetCommunication(cellchat, net = net,
                               datasets       = "AD",
                               ligand.logFC   = 0.05,
                               receptor.logFC = 0.05)

# Higher in Healthy vs AD
net.down <- subsetCommunication(cellchat, net = net,
                                 datasets       = "Healthy",
                                 ligand.logFC   = -0.05,
                                 receptor.logFC = -0.05)

cat(sprintf("\nnet.up   (higher in AD)     : %d L-R pairs\n", nrow(net.up)))
cat(sprintf("net.down (higher in Healthy): %d L-R pairs\n", nrow(net.down)))

write.csv(net.up,   "netMappingDEG_upregulated_AD.csv",      row.names = FALSE)
write.csv(net.down, "netMappingDEG_upregulated_Healthy.csv", row.names = FALSE)

gene.up   <- extractGeneSubsetFromPair(net.up,   cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

cat(sprintf("\nGenes upregulated in AD     : %d\n", length(gene.up)))
cat(sprintf("Genes upregulated in Healthy: %d\n", length(gene.down)))

write.csv(gene.up,   "genes_upregulated_AD.csv",      row.names = FALSE)
write.csv(gene.down, "genes_upregulated_Healthy.csv", row.names = FALSE)

# ── Manual filtering ──────────────────────────────────────────────
p_cut   <- 0.05
lfc_cut <- 0.25   # standard publication threshold
pct_cut <- 0.10

# AD: both ligand AND receptor upregulated vs Healthy
AD_LR_bothUp <- subset(net,
  datasets         == "AD" &
  ligand.pvalues   <  p_cut &
  receptor.pvalues <  p_cut &
  ligand.logFC     >  lfc_cut &
  receptor.logFC   >  lfc_cut &
  ligand.pct.1     >= pct_cut &
  receptor.pct.1   >= pct_cut)

# AD: both ligand AND receptor downregulated vs Healthy
AD_LR_bothDown <- subset(net,
  datasets         == "AD" &
  ligand.pvalues   <  p_cut &
  receptor.pvalues <  p_cut &
  ligand.logFC     < -lfc_cut &
  receptor.logFC   < -lfc_cut &
  ligand.pct.1     >= pct_cut &
  receptor.pct.1   >= pct_cut)

# Healthy: both ligand AND receptor upregulated vs AD
Healthy_LR_bothUp <- subset(net,
  datasets         == "Healthy" &
  ligand.pvalues   <  p_cut &
  receptor.pvalues <  p_cut &
  ligand.logFC     >  lfc_cut &
  receptor.logFC   >  lfc_cut &
  ligand.pct.1     >= pct_cut &
  receptor.pct.1   >= pct_cut)

# Healthy: both ligand AND receptor downregulated vs AD
Healthy_LR_bothDown <- subset(net,
  datasets         == "Healthy" &
  ligand.pvalues   <  p_cut &
  receptor.pvalues <  p_cut &
  ligand.logFC     < -lfc_cut &
  receptor.logFC   < -lfc_cut &
  ligand.pct.1     >= pct_cut &
  receptor.pct.1   >= pct_cut)

cat(sprintf("\nAD_LR_bothUp       : %d pairs\n", nrow(AD_LR_bothUp)))
cat(sprintf("AD_LR_bothDown     : %d pairs\n", nrow(AD_LR_bothDown)))
cat(sprintf("Healthy_LR_bothUp  : %d pairs\n", nrow(Healthy_LR_bothUp)))
cat(sprintf("Healthy_LR_bothDown: %d pairs\n", nrow(Healthy_LR_bothDown)))

cat("\nAD_LR_bothUp — source x target:\n")
print(table(AD_LR_bothUp$source, AD_LR_bothUp$target))

cat("\nHealthy_LR_bothUp — source x target:\n")
print(table(Healthy_LR_bothUp$source, Healthy_LR_bothUp$target))

write.csv(AD_LR_bothUp,       "AD_LR_bothUp.csv",       row.names = FALSE)
write.csv(AD_LR_bothDown,     "AD_LR_bothDown.csv",     row.names = FALSE)
write.csv(Healthy_LR_bothUp,  "Healthy_LR_bothUp.csv",  row.names = FALSE)
write.csv(Healthy_LR_bothDown,"Healthy_LR_bothDown.csv",row.names = FALSE)

write.table(AD_LR_bothUp$interaction_name,
            "AD_LR_bothUp_interactions.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(AD_LR_bothDown$interaction_name,
            "AD_LR_bothDown_interactions.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(Healthy_LR_bothUp$interaction_name,
            "Healthy_LR_bothUp_interactions.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(Healthy_LR_bothDown$interaction_name,
            "Healthy_LR_bothDown_interactions.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)
cat("✓ All differential L-R tables saved\n")

# ── Bubble plots for up/down ──────────────────────────────────────
pairLR.use.up   <- net.up[,   "interaction_name", drop = FALSE]
pairLR.use.down <- net.down[, "interaction_name", drop = FALSE]

tryCatch({
  gg_up <- netVisual_bubble(cellchat,
                             pairLR.use     = pairLR.use.up,
                             comparison     = c(1, 2),
                             angle.x        = 90,
                             remove.isolate = TRUE,
                             title.name     = "Up-regulated in AD vs Healthy")
  ggsave("pairLR_upregulated_AD_bubble.png",
         gg_up, width = 12, height = 10, dpi = 300)
  cat("✓ Up-regulated bubble plot saved\n")
}, error = function(e) {
  cat("⚠️  Up bubble plot failed:", as.character(e), "\n")
})

tryCatch({
  gg_down <- netVisual_bubble(cellchat,
                               pairLR.use     = pairLR.use.down,
                               comparison     = c(1, 2),
                               angle.x        = 90,
                               remove.isolate = TRUE,
                               title.name     = "Up-regulated in Healthy vs AD")
  ggsave("pairLR_upregulated_Healthy_bubble.png",
         gg_down, width = 12, height = 10, dpi = 300)
  cat("✓ Healthy up-regulated bubble plot saved\n")
}, error = function(e) {
  cat("⚠️  Down bubble plot failed:", as.character(e), "\n")
})

# ── Enrichment scores ─────────────────────────────────────────────
tryCatch({
  png("enrichmentScore_upregulated_AD.png",
      width = 7, height = 6, units = "in", res = 300)
  computeEnrichmentScore(net.up, species = "human", variable.both = TRUE)
  dev.off()
  cat("✓ Enrichment score (AD up) saved\n")
}, error = function(e) {
  if (dev.cur() > 1) dev.off()
  cat("⚠️  Enrichment score failed:", as.character(e), "\n")
})

tryCatch({
  png("enrichmentScore_upregulated_Healthy.png",
      width = 7, height = 6, units = "in", res = 300)
  computeEnrichmentScore(net.down, species = "human", variable.both = TRUE)
  dev.off()
  cat("✓ Enrichment score (Healthy up) saved\n")
}, error = function(e) {
  if (dev.cur() > 1) dev.off()
  cat("⚠️  Enrichment score failed:", as.character(e), "\n")
})

} # end if (!deg_success)

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
                                    levels = c("Healthy", "AD"))
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
                              pathway, output_label,
                              names(object.list)[i])),
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
                            signaling.name = paste(pathway,
                                                   names(object.list)[i]))
      }
      dev.off()
    }, error = function(e) { if (dev.cur() > 1) dev.off() })

    tryCatch({
      p <- plotGeneExpression(cellchat, signaling = pathway,
                               split.by = "datasets",
                               colors.ggplot = TRUE, type = "violin")
      ggsave(file.path(category_dir,
                       sprintf("%s_%s_violin.png",
                               pathway, output_label)),
             p, width = 7, height = 5, dpi = 300)
    }, error = function(e) invisible(NULL))

    success_count <- success_count + 1
    if (success_count %% 5 == 0) gc(verbose = FALSE)
  }

  cat(sprintf("  ✓ %d / %d pathways processed\n",
              success_count, length(pathway_list)))
}

detected_pathways <- unique(c(healthy@netP$pathways, ad@netP$pathways))
cat(sprintf("Total detected pathways: %d\n", length(detected_pathways)))

pathway_plots_dir <- "pathway_specific_plots"
dir.create(pathway_plots_dir, showWarnings = FALSE)

generate_pathway_plots(detected_pathways, "all_pathways",
                       object.list, cellchat, pathway_plots_dir)

if (nrow(net.up) > 0)
  generate_pathway_plots(unique(net.up$pathway_name),
                         "net_up_AD", object.list, cellchat,
                         pathway_plots_dir)

if (nrow(net.down) > 0)
  generate_pathway_plots(unique(net.down$pathway_name),
                         "net_up_Healthy", object.list, cellchat,
                         pathway_plots_dir)

if (nrow(AD_LR_bothUp) > 0)
  generate_pathway_plots(unique(AD_LR_bothUp$pathway_name),
                         "AD_LR_bothUp", object.list, cellchat,
                         pathway_plots_dir)

if (nrow(AD_LR_bothDown) > 0)
  generate_pathway_plots(unique(AD_LR_bothDown$pathway_name),
                         "AD_LR_bothDown", object.list, cellchat,
                         pathway_plots_dir)

if (nrow(Healthy_LR_bothUp) > 0)
  generate_pathway_plots(unique(Healthy_LR_bothUp$pathway_name),
                         "Healthy_LR_bothUp", object.list, cellchat,
                         pathway_plots_dir)

if (nrow(Healthy_LR_bothDown) > 0)
  generate_pathway_plots(unique(Healthy_LR_bothDown$pathway_name),
                         "Healthy_LR_bothDown", object.list, cellchat,
                         pathway_plots_dir)

cat("✓ All pathway-specific plots complete\n")

# ══════════════════════════════════════════════════════════════════
# SAVE FINAL OBJECTS
# ══════════════════════════════════════════════════════════════════

cat("\n=== SAVING FINAL OBJECTS ===\n")
saveRDS(cellchat,    "cellchat_Healthy_vs_AD_merged.rds")
saveRDS(object.list, "cellchat_object.list_Healthy_AD.rds")
cat("✓ Saved: cellchat_Healthy_vs_AD_merged.rds\n")
cat("✓ Saved: cellchat_object.list_Healthy_AD.rds\n")

# ══════════════════════════════════════════════════════════════════
# SAVE DATA.SIGNALING MATRICES
# ══════════════════════════════════════════════════════════════════

cat("\n=== SAVING DATA.SIGNALING MATRICES ===\n")
mat_all  <- cellchat@data.signaling
cond     <- cellchat@meta$datasets
names(cond) <- rownames(cellchat@meta)
cond_ds  <- cond[colnames(mat_all)]

cells_Healthy <- colnames(mat_all)[cond_ds == "Healthy"]
cells_AD      <- colnames(mat_all)[cond_ds == "AD"]

cat(sprintf("Cells Healthy: %d\n", length(cells_Healthy)))
cat(sprintf("Cells AD     : %d\n", length(cells_AD)))

mat_Healthy <- mat_all[, cells_Healthy, drop = FALSE]
mat_AD      <- mat_all[, cells_AD,      drop = FALSE]

cat("\nTop genes by mean (Healthy):\n")
print(head(sort(Matrix::rowMeans(mat_Healthy), decreasing = TRUE), 20))
cat("\nTop genes by mean (AD):\n")
print(head(sort(Matrix::rowMeans(mat_AD),      decreasing = TRUE), 20))

saveRDS(mat_all,     "cellchat_data.signaling.rds")
saveRDS(mat_Healthy, "cellchat_data.signaling_Healthy.rds")
saveRDS(mat_AD,      "cellchat_data.signaling_AD.rds")

writeMM(mat_all, "cellchat_data.signaling.mtx")
write.table(rownames(mat_all), "cellchat_data.signaling_genes.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(colnames(mat_all), "cellchat_data.signaling_cells.txt",
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
cat("║   ✓✓✓  HEALTHY vs AD COMPARISON COMPLETE  ✓✓✓           ║\n")
cat("║                                                           ║\n")
cat("╚═══════════════════════════════════════════════════════════╝\n")
cat("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Log saved:", log_file, "\n\n")

sink()
