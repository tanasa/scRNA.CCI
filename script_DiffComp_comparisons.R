# ============================================================
# scDiffCom Analysis : AD vs Dyslexia 
# Author  : B. Tanasa
# Dataset : 18_merged_AG_031826.26celltypes.min50cells.18K.proportional
# Compare : AD (cond1) vs Dyslexia (cond2)
# ============================================================

# ── Working directory ────────────────────────────────────────
setwd("/mnt/nfs/CX000008_DS1/projects/jaeyeon/fastq_file_Dyslexia_r1/zanalysis_bogdan/samples.merged_AG_Harmony_res0.x_anno_031826")

# ── Libraries ────────────────────────────────────────────────
library(Seurat)
library(Matrix)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(gridExtra)
library(data.table)
library(cowplot)
library(patchwork)
library(future)
library(parallelly)
library(foreach)
library(doParallel)
library(ggbeeswarm)
library(BiocParallel)
library(scDiffCom)
library(readr)
library(forcats)
library(RColorBrewer)


# ── Parallelization (multicore / fork-based, Linux only) ────
plan(multicore, workers = min(parallelly::availableCores() - 1, 60))


# ── Global options ───────────────────────────────────────────
Sys.setenv(TMPDIR = "/mnt/nfs/CX000008_DS1/projects/btanasa/tmp")
options(future.globals.maxSize = 100 * 1024^3)   # 100 GB
options(bitmapType = "cairo")                      # crisp plots on server


# ── ggplot theme ─────────────────────────────────────────────
theme_set(theme_classic(base_size = 14))
update_geom_defaults("point", list(size = 1.2, alpha = 0.8))
update_geom_defaults("text",  list(size = 4))


# ============================================================
# Key analysis identifiers  (edit here — used throughout)
# ============================================================

rds_file        <- "18_merged_AG_031826.26celltypes.min50cells.18K.proportional.anndataR.AD_Dyslexia.rds"
analysis_name   <- "AD_vs_Dyslexia_v2"                    # used in filenames
condition_col   <- "group"                                # metadata column
cond1           <- "AD"
cond2           <- "Dyslexia"
celltype_col    <- "subcelltype"
scdiffcom_name  <- "scDiffCom_AD_vs_Dyslexia_v2"          # name stored inside the object

# ============================================================
# Output directory structure
# ============================================================

base_dir   <- "/mnt/nfs/CX000008_DS1/projects/jaeyeon/fastq_file_Dyslexia_r1/zanalysis_bogdan/samples.merged_AG_Harmony_res0.x_anno_031826"
output_dir <- file.path(base_dir, paste0("scDiffCom_results_", analysis_name))

for (subdir in c("plots", "tables", "rds_objects", "logs")) {
  dir.create(file.path(output_dir, subdir), showWarnings = FALSE, recursive = TRUE)
}

# Convenience path helpers
outfn_plot  <- function(fn) file.path(output_dir, "plots",       fn)
outfn_table <- function(fn) file.path(output_dir, "tables",      fn)
outfn_rds   <- function(fn) file.path(output_dir, "rds_objects", fn)
outfn_log   <- function(fn) file.path(output_dir, "logs",        fn)

cat("Output directory:\n", output_dir, "\n")
cat("  ├── plots/\n  ├── tables/\n  ├── rds_objects/\n  └── logs/\n\n")


# ── Start log ────────────────────────────────────────────────
sink(outfn_log(paste0(analysis_name, "_run.log")), split = TRUE)

cat("Run started  :", format(Sys.time()), "\n")
cat("R version    :", R.version$version.string, "\n")
cat("Seurat       :", as.character(packageVersion("Seurat")), "\n")
cat("scDiffCom    :", as.character(packageVersion("scDiffCom")), "\n\n")


# ============================================================
# Load Seurat object
# ============================================================

if (!file.exists(rds_file)) {
  stop("RDS file not found: ", rds_file, "\nCurrent dir: ", getwd())
}

obj <- readRDS(rds_file)

cat("Object class :", class(obj), "\n")
cat("Cells        :", ncol(obj), "| Genes:", nrow(obj), "\n")
cat("Metadata columns:\n")
print(colnames(obj@meta.data))
cat("\nCondition counts:\n")
print(table(obj[[condition_col]]))
cat("\nSubcelltype counts:\n")
print(table(obj[[celltype_col]]))


# ── UMAP overview plot ───────────────────────────────────────
# p_umap <- DimPlot(obj, reduction = "umap.harmony",
#                  group.by = celltype_col, label = TRUE, repel = TRUE) +
#  ggtitle(paste("UMAP —", analysis_name)) +
#  theme_minimal()

# ggsave(outfn_plot(paste0("UMAP_subcelltype_", analysis_name, ".pdf")),
#        plot = p_umap, width = 12, height = 9)


# ============================================================
# Seurat v5 compatibility fix
# ============================================================

cat("\nSeurat version:", as.character(packageVersion("Seurat")), "\n")

if (packageVersion("Seurat") >= "5.0.0") {
  cat("Joining layers for Seurat v5 compatibility...\n")
  obj <- JoinLayers(obj)
}

cat("Assay layers after joining:\n")
print(Layers(obj))


# ============================================================
# Run scDiffCom
# ============================================================

cat("\n=== scDiffCom : starting", format(Sys.time()), "===\n")

scdiffcom_result <- run_interaction_analysis(
  seurat_object                 = obj,
  LRI_species                   = "human",
  seurat_celltype_id            = celltype_col,
  seurat_condition_id           = list(
    column_name = condition_col,
    cond1_name  = cond1,
    cond2_name  = cond2
  ),
  iterations                    = 1000,
  scdiffcom_object_name         = scdiffcom_name,
  seurat_assay                  = "RNA",
  seurat_layer                  = "data",
  log_scale                     = TRUE,
  score_type                    = "geometric_mean",
  threshold_min_cells           = 5,
  threshold_pct                 = 0.1,
  threshold_quantile_score      = 0.2,
  threshold_p_value_specificity = 0.05,
  threshold_p_value_de          = 0.05,
  threshold_logfc               = log(1.5),
  return_distributions          = FALSE,
  seed                          = 42,
  verbose                       = TRUE
)

cat("=== scDiffCom : finished", format(Sys.time()), "===\n\n")


# ============================================================
# Save results
# ============================================================

# Full scDiffCom object
saveRDS(scdiffcom_result,
        outfn_rds(paste0(scdiffcom_name, ".rds")))
cat("Saved RDS:", outfn_rds(paste0(scdiffcom_name, ".rds")), "\n")

# Result table as CSV
result_table <- GetTableCCI(scdiffcom_result, type = "detected", simplified = FALSE)
write.csv(result_table,
          outfn_table(paste0(scdiffcom_name, "_CCI_table.csv")),
          row.names = FALSE)
cat("Saved CCI table:", outfn_table(paste0(scdiffcom_name, "_CCI_table.csv")), "\n")


# ── Close log ────────────────────────────────────────────────
cat("\nRun finished:", format(Sys.time()), "\n")
sink()

# ============================================================
# scDiffCom Analysis : AD vs Dyslexia (v1)
# PART 2 : Summaries, ORA, Plots, Tables
# Requires: scdiffcom_result, outfn_*, analysis_name, cci_detected
# ============================================================

# ============================================================
# Analysis Summary
# ============================================================

cat("\n=== Analysis Summary ===\n")

cci_detected <- GetTableCCI(scdiffcom_result,
                            type       = "detected",
                            simplified = TRUE)

cat("Total detected CCIs    :", nrow(cci_detected), "\n")
cat("UP-regulated           :", sum(cci_detected$REGULATION == "UP",   na.rm = TRUE), "\n")
cat("DOWN-regulated         :", sum(cci_detected$REGULATION == "DOWN", na.rm = TRUE), "\n")
cat("FLAT (no change)       :", sum(cci_detected$REGULATION == "FLAT", na.rm = TRUE), "\n")

head(cci_detected)
tail(cci_detected)


# ============================================================
# ORA results — inspect
# ============================================================

ORA_results <- GetTableORA(scdiffcom_result, categories = "all", simplified = TRUE)

cat("\nAvailable ORA categories:\n")
print(names(ORA_results))

head(ORA_results$LRI, 2)
head(ORA_results$LIGAND_COMPLEX, 2)


# ============================================================
# Save all ORA results by category
# ============================================================

cat("\n=== Saving ORA Results ===\n")

ora_results          <- scdiffcom_result@ora_table
available_categories <- names(ora_results)

for (category in available_categories) {
  ora_data <- ora_results[[category]]
  if (!is.null(ora_data) && nrow(ora_data) > 0) {
    write.csv(ora_data,
              outfn_table(paste0("ORA_", category, "_", analysis_name, ".csv")),
              row.names = FALSE)
    cat(sprintf("ORA %s saved (%d rows)\n", category, nrow(ora_data)))
  } else {
    cat(sprintf("No data for ORA %s\n", category))
  }
}

# Save complete ORA as RDS
saveRDS(ora_results, outfn_rds(paste0("ORA_all_categories_", analysis_name, ".rds")))
cat("All ORA results saved as single RDS\n")

# ORA summary table
ora_summary <- data.frame(
  Category  = available_categories,
  N_Enriched = sapply(available_categories, function(cat) {
    if (!is.null(ora_results[[cat]])) nrow(ora_results[[cat]]) else 0
  })
)
write.csv(ora_summary,
          outfn_table(paste0("ORA_summary_all_categories_", analysis_name, ".csv")),
          row.names = FALSE)
cat("\n=== ORA Summary ===\n")
print(ora_summary)


# ============================================================
# Volcano plot
# ============================================================

cat("\n=== Creating Volcano Plot ===\n")

p_volcano <- ggplot(
  cci_detected,
  aes(x = LOGFC, y = -log10(BH_P_VALUE_DE + 1E-2), colour = REGULATION)
) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_colour_manual(values = c("UP" = "red", "DOWN" = "blue",
                                 "FLAT" = "green", "NSC" = "grey")) +
  xlab("log(FC)") +
  ylab("-log10(Adj. p-value)") +
  labs(title    = paste("Volcano Plot — Differential CCIs"),
       subtitle = paste(cond1, "vs", cond2),
       colour   = "Regulation") +
  theme_classic(base_size = 14) +
  theme(plot.title    = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "right") +
  geom_vline(xintercept = 0,        linetype = "dashed", color = "black", alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.5)

ggsave(outfn_plot(paste0("volcano_CCIs_", analysis_name, ".png")),
       p_volcano, width = 10, height = 8, dpi = 300)
ggsave(outfn_plot(paste0("volcano_CCIs_", analysis_name, ".pdf")),
       p_volcano, width = 10, height = 8)
cat("Volcano plot saved\n")


# ── Labeled volcano (top 10 UP + top 10 DOWN by significance + effect) ──
# Selection criterion: significant (BH p < 0.05) AND ranked by a combined
# score = |LOGFC| * -log10(BH_P_VALUE_DE + 1e-2)
# This ensures labeled points are both large in effect size AND
# statistically supported — not just extreme logFC outliers.

top_interactions <- bind_rows(
  cci_detected %>%
    filter(REGULATION == "UP",   BH_P_VALUE_DE < 0.05) %>%
    mutate(rank_score = LOGFC * -log10(BH_P_VALUE_DE + 1e-2)) %>%
    arrange(desc(rank_score)) %>%
    head(10),
  cci_detected %>%
    filter(REGULATION == "DOWN", BH_P_VALUE_DE < 0.05) %>%
    mutate(rank_score = abs(LOGFC) * -log10(BH_P_VALUE_DE + 1e-2)) %>%
    arrange(desc(rank_score)) %>%
    head(10)
)

cat("Labeled interactions (UP significant)  :",
    sum(top_interactions$REGULATION == "UP"), "\n")
cat("Labeled interactions (DOWN significant):",
    sum(top_interactions$REGULATION == "DOWN"), "\n")

p_volcano_labeled <- p_volcano +
  geom_text_repel(
    data          = top_interactions,
    aes(label     = LRI),
    size          = 3,
    max.overlaps  = 20,
    box.padding   = 0.5,
    point.padding = 0.3,
    segment.color = "grey50",
    segment.size  = 0.3
  ) +
  labs(subtitle = paste("Top significant interactions labeled (BH p < 0.05) —",
                        cond1, "vs", cond2))

ggsave(outfn_plot(paste0("volcano_CCIs_labeled_", analysis_name, ".png")),
       p_volcano_labeled, width = 12, height = 10, dpi = 300)
ggsave(outfn_plot(paste0("volcano_CCIs_labeled_", analysis_name, ".pdf")),
       p_volcano_labeled, width = 12, height = 10)
cat("Labeled volcano plot saved\n")


# ── Faceted volcano (UP / DOWN only) ────────────────────────

p_volcano_facet <- ggplot(
  cci_detected %>% filter(REGULATION %in% c("UP", "DOWN")),
  aes(x = LOGFC, y = -log10(BH_P_VALUE_DE + 1E-2), colour = REGULATION)
) +
  geom_point(alpha = 0.7, size = 2) +
  scale_colour_manual(values = c("UP" = "red", "DOWN" = "blue")) +
  facet_wrap(~ REGULATION, ncol = 2, scales = "free_x") +
  xlab("log(FC)") + ylab("-log10(Adj. p-value)") +
  labs(title    = "Volcano Plots: UP vs DOWN Regulated CCIs",
       subtitle = paste(cond1, "vs", cond2)) +
  theme_classic(base_size = 14) +
  theme(plot.title       = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle    = element_text(hjust = 0.5),
        legend.position  = "none",
        strip.background = element_rect(fill = "lightgray", color = "black"),
        strip.text       = element_text(face = "bold")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.5)

ggsave(outfn_plot(paste0("volcano_CCIs_faceted_", analysis_name, ".png")),
       p_volcano_facet, width = 12, height = 6, dpi = 300)
cat("Faceted volcano plot saved\n")

# Volcano statistics table
volcano_stats <- cci_detected %>%
  group_by(REGULATION) %>%
  summarise(count        = n(),
            mean_logfc   = mean(LOGFC, na.rm = TRUE),
            median_logfc = median(LOGFC, na.rm = TRUE),
            mean_pval    = mean(-log10(BH_P_VALUE_DE + 1E-2), na.rm = TRUE),
            .groups = "drop")
write.csv(volcano_stats,
          outfn_table(paste0("volcano_statistics_", analysis_name, ".csv")),
          row.names = FALSE)


# ============================================================
# ORA plots — helper to avoid repetition
# ============================================================

plot_ora <- function(category, regulation, go_aspect = NULL,
                     max_terms = 50, width = 12, height = 10) {

  args <- list(object        = scdiffcom_result,
               category      = category,
               regulation    = regulation,
               max_terms_show = max_terms)
  if (!is.null(go_aspect)) args$GO_aspect <- go_aspect

  p <- tryCatch(do.call(PlotORA, args), error = function(e) {
    message("PlotORA skipped (", category, " / ", regulation, "): ", e$message)
    return(NULL)
  })
  if (is.null(p)) return(invisible(NULL))

  aspect_tag <- if (!is.null(go_aspect)) paste0("_", go_aspect) else ""
  fname <- paste0("ORA_", category, aspect_tag, "_", regulation, "_", analysis_name)

  p <- p +
    theme(legend.position  = c(0.85, 0.3),
          legend.key.size  = unit(0.5, "cm"),
          plot.title       = element_text(hjust = 0.5, face = "bold", size = 14),
          axis.text.y      = element_text(size = 8)) +
    labs(title = paste("ORA:", category, go_aspect, "—", regulation, "regulated"))

  print(p)
  ggsave(outfn_plot(paste0(fname, ".png")), p, width = width, height = height, dpi = 300)
  ggsave(outfn_plot(paste0(fname, ".pdf")), p, width = width, height = height)
  cat("ORA plot saved:", fname, "\n")
}

# ── LRI ──────────────────────────────────────────────────────
cat("\n=== ORA Plots : LRI ===\n")
plot_ora("LRI", "UP",   max_terms = 100)
plot_ora("LRI", "DOWN", max_terms = 100)

# ── GO Terms ─────────────────────────────────────────────────
cat("\n=== ORA Plots : GO Terms ===\n")
for (reg in c("UP", "DOWN")) {
  for (aspect in c("biological_process", "molecular_function", "cellular_component")) {
    plot_ora("GO_TERMS", reg, go_aspect = aspect, width = 14, height = 10)
  }
}

# ── KEGG ─────────────────────────────────────────────────────
cat("\n=== ORA Plots : KEGG ===\n")
plot_ora("KEGG_PWS", "UP",   width = 14)
plot_ora("KEGG_PWS", "DOWN", width = 14)

# ── Cell type pairs ───────────────────────────────────────────
cat("\n=== ORA Plots : Cell Types ===\n")
for (reg in c("UP", "DOWN")) {
  plot_ora("ER_CELLTYPES",      reg, width = 12)
  plot_ora("EMITTER_CELLTYPE",  reg, width = 10, height = 8)
  plot_ora("RECEIVER_CELLTYPE", reg, width = 10, height = 8)
}

cat("All ORA plots saved\n")


# ============================================================
# Export main CCI tables
# ============================================================

cat("\n=== Exporting CCI Tables ===\n")

# Detected (filtered)
write.csv(cci_detected,
          outfn_table(paste0("CCI_detected_", analysis_name, ".csv")),
          row.names = FALSE)
cat("Detected CCIs saved\n")

# Raw (all interactions before thresholds)
cci_raw <- GetTableCCI(scdiffcom_result, type = "raw", simplified = TRUE)
write.csv(cci_raw,
          outfn_table(paste0("CCI_raw_", analysis_name, ".csv")),
          row.names = FALSE)
cat("Raw CCIs saved\n")

# Top 100 UP / DOWN
write.csv(cci_detected %>% filter(REGULATION == "UP")   %>% arrange(desc(LOGFC)) %>% head(100),
          outfn_table(paste0("CCI_top100_UP_",   analysis_name, ".csv")), row.names = FALSE)
write.csv(cci_detected %>% filter(REGULATION == "DOWN") %>% arrange(LOGFC)       %>% head(100),
          outfn_table(paste0("CCI_top100_DOWN_", analysis_name, ".csv")), row.names = FALSE)

# Significant UP / DOWN (BH p < 0.05)
cci_UP_sig   <- cci_detected %>% filter(REGULATION == "UP",   BH_P_VALUE_DE < 0.05)
cci_DOWN_sig <- cci_detected %>% filter(REGULATION == "DOWN", BH_P_VALUE_DE < 0.05)
cci_sig_both <- bind_rows(cci_UP_sig, cci_DOWN_sig)

write.csv(cci_UP_sig,   outfn_table(paste0("CCI_UP_significant_",           analysis_name, ".csv")), row.names = FALSE)
write.csv(cci_DOWN_sig, outfn_table(paste0("CCI_DOWN_significant_",         analysis_name, ".csv")), row.names = FALSE)
write.csv(cci_sig_both, outfn_table(paste0("CCI_significant_UP_and_DOWN_",  analysis_name, ".csv")), row.names = FALSE)

cat("Significant UP  :", nrow(cci_UP_sig),
    sprintf("(%.1f%%)\n", 100 * nrow(cci_UP_sig)   / sum(cci_detected$REGULATION == "UP")))
cat("Significant DOWN:", nrow(cci_DOWN_sig),
    sprintf("(%.1f%%)\n", 100 * nrow(cci_DOWN_sig) / sum(cci_detected$REGULATION == "DOWN")))


# ============================================================
# Cell type summaries (Emitter / Receiver)
# ============================================================

cat("\n=== Cell Type Summaries ===\n")

# Split ER_CELLTYPES into EMITTER and RECEIVER
# NOTE: scDiffCom uses "-" as separator (underscores in cell type
#       names are converted to dashes internally)
cci_split <- cci_detected %>%
  separate(col    = ER_CELLTYPES,
           into   = c("EMITTER", "RECEIVER"),
           sep    = "-",
           extra  = "merge",
           remove = FALSE)

# Overall emitter / receiver counts
emitter_counts  <- cci_split %>% count(EMITTER,  sort = TRUE)
receiver_counts <- cci_split %>% count(RECEIVER, sort = TRUE)

write.csv(emitter_counts,  outfn_table(paste0("emitter_counts_",  analysis_name, ".csv")), row.names = FALSE)
write.csv(receiver_counts, outfn_table(paste0("receiver_counts_", analysis_name, ".csv")), row.names = FALSE)

# UP / DOWN emitter counts
write.csv(cci_split %>% filter(REGULATION == "UP")   %>% count(EMITTER, sort = TRUE),
          outfn_table(paste0("emitter_counts_UP_",   analysis_name, ".csv")), row.names = FALSE)
write.csv(cci_split %>% filter(REGULATION == "DOWN") %>% count(RECEIVER, sort = TRUE),
          outfn_table(paste0("receiver_counts_DOWN_", analysis_name, ".csv")), row.names = FALSE)

# Emitter bar plot
p_emitter <- ggplot(emitter_counts,
                    aes(x = reorder(EMITTER, -n), y = n)) +
  geom_col(fill = "#1976D2") +
  coord_flip() +
  labs(title = "Interactions per Emitter Cell Type",
       x = "Emitter", y = "Count") +
  theme_classic(base_size = 12)

ggsave(outfn_plot(paste0("barplot_emitter_", analysis_name, ".png")),
       p_emitter, width = 8, height = 6, dpi = 300)

# Receiver bar plot
p_receiver <- ggplot(receiver_counts,
                     aes(x = reorder(RECEIVER, -n), y = n)) +
  geom_col(fill = "#C62828") +
  coord_flip() +
  labs(title = "Interactions per Receiver Cell Type",
       x = "Receiver", y = "Count") +
  theme_classic(base_size = 12)

ggsave(outfn_plot(paste0("barplot_receiver_", analysis_name, ".png")),
       p_receiver, width = 8, height = 6, dpi = 300)

print(p_emitter)
print(p_receiver)


# ── ER_CELLTYPE summary by regulation ────────────────────────
emitter_summary <- cci_detected %>%
  group_by(ER_CELLTYPES, REGULATION) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = REGULATION, values_from = n, values_fill = 0) %>%
  mutate(TOTAL = rowSums(across(where(is.numeric)))) %>%
  arrange(desc(if ("UP" %in% names(.)) UP else TOTAL))

write.csv(emitter_summary,
          outfn_table(paste0("celltype_pair_summary_", analysis_name, ".csv")),
          row.names = FALSE)
cat("\nTop cell type pairs:\n")
print(head(emitter_summary, 10))

# ── LRI summary by regulation ─────────────────────────────────
lri_summary <- cci_detected %>%
  group_by(LRI, REGULATION) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = REGULATION, values_from = n, values_fill = 0) %>%
  mutate(TOTAL = rowSums(across(where(is.numeric)))) %>%
  arrange(desc(if ("UP" %in% names(.)) UP else TOTAL))

write.csv(lri_summary,
          outfn_table(paste0("LRI_summary_", analysis_name, ".csv")),
          row.names = FALSE)
cat("\nTop LRIs:\n")
print(head(lri_summary, 10))


# ============================================================
# CCI Score boxplots (cond1 vs cond2 by REGULATION)
# ============================================================

cat("\n=== CCI Score Boxplots ===\n")

# Column names for scores depend on condition labels
score_col1 <- paste0("CCI_SCORE_", cond1)
score_col2 <- paste0("CCI_SCORE_", cond2)

if (all(c(score_col1, score_col2) %in% colnames(cci_detected))) {

  df_long <- cci_detected %>%
    select(REGULATION, all_of(c(score_col1, score_col2))) %>%
    pivot_longer(cols = all_of(c(score_col1, score_col2)),
                 names_to  = "Condition",
                 values_to = "CCI_SCORE") %>%
    mutate(Condition  = recode(Condition,
                               !!score_col1 := cond1,
                               !!score_col2 := cond2),
           REGULATION = factor(REGULATION, levels = c("UP", "DOWN", "FLAT"))) %>%
    filter(is.finite(CCI_SCORE))

  p_box <- ggplot(df_long, aes(x = Condition, y = CCI_SCORE, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, width = 0.7) +
    geom_jitter(width = 0.15, alpha = 0.15, size = 0.6) +
    facet_wrap(~ REGULATION, nrow = 1, scales = "fixed") +
    scale_fill_manual(values = setNames(c("#1976D2", "#C62828"), c(cond1, cond2))) +
    labs(title = paste("CCI Scores by Condition —", analysis_name),
         x = NULL, y = "CCI Score") +
    coord_cartesian(ylim = c(0, 10)) +
    theme_classic(base_size = 12) +
    theme(legend.position = "none")

  ggsave(outfn_plot(paste0("boxplot_CCI_score_", analysis_name, ".png")),
         p_box, width = 12, height = 5, dpi = 300)
  print(p_box)

  # Summary table
  score_summary <- df_long %>%
    group_by(REGULATION, Condition) %>%
    summarise(n = n(), mean = mean(CCI_SCORE), median = median(CCI_SCORE),
              sd = sd(CCI_SCORE), IQR = IQR(CCI_SCORE), .groups = "drop")
  write.csv(score_summary,
            outfn_table(paste0("CCI_score_summary_", analysis_name, ".csv")),
            row.names = FALSE)
  print(score_summary)

} else {
  cat("Score columns not found — skipping boxplot\n",
      "Expected:", score_col1, "and", score_col2, "\n",
      "Available:", paste(colnames(cci_detected), collapse = ", "), "\n")
}


# ============================================================
# LOGFC distribution plots
# ============================================================

cat("\n=== LOGFC Distribution Plots ===\n")

df_logfc <- cci_detected %>%
  mutate(REGULATION = factor(REGULATION, levels = c("UP", "DOWN", "FLAT"))) %>%
  filter(is.finite(LOGFC))

# Faceted by REGULATION
p_logfc <- ggplot(df_logfc, aes(x = LOGFC, fill = REGULATION, color = REGULATION)) +
  geom_histogram(aes(y = after_stat(density)), bins = 60, alpha = 0.2, position = "identity") +
  geom_density(alpha = 0.5, linewidth = 0.7) +
  facet_wrap(~ REGULATION, nrow = 1, scales = "free_y") +
  labs(title = paste("LOGFC Distribution by Regulation —", analysis_name),
       x = "LOGFC", y = "Density") +
  coord_cartesian(ylim = c(0, 10)) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none")

ggsave(outfn_plot(paste0("LOGFC_distribution_by_regulation_", analysis_name, ".png")),
       p_logfc, width = 10, height = 4, dpi = 300)

# Overall
p_logfc_all <- ggplot(df_logfc, aes(x = LOGFC)) +
  geom_histogram(aes(y = after_stat(density)), bins = 80, alpha = 0.25) +
  geom_density(linewidth = 0.9) +
  labs(title = paste("LOGFC Distribution (all CCIs) —", analysis_name),
       x = "LOGFC", y = "Density") +
  coord_cartesian(ylim = c(0, 10)) +
  theme_classic(base_size = 12)

ggsave(outfn_plot(paste0("LOGFC_distribution_overall_", analysis_name, ".png")),
       p_logfc_all, width = 6, height = 4, dpi = 300)

# LOGFC summary
logfc_summary <- df_logfc %>%
  group_by(REGULATION) %>%
  summarise(n = n(), mean = mean(LOGFC), median = median(LOGFC),
            sd = sd(LOGFC), IQR = IQR(LOGFC),
            min = min(LOGFC), max = max(LOGFC), .groups = "drop")
write.csv(logfc_summary,
          outfn_table(paste0("LOGFC_summary_", analysis_name, ".csv")),
          row.names = FALSE)
print(logfc_summary)


# ============================================================
# Detection flags
# ============================================================

cat("\n=== Detection Flags ===\n")

flag_cols <- c("IS_CCI_DETECTED_COND1", "IS_CCI_DETECTED_COND2")
# Try the actual condition-named columns if the generic ones don't exist
if (!all(flag_cols %in% colnames(cci_detected))) {
  flag_cols <- grep("IS_CCI_DETECTED", colnames(cci_detected), value = TRUE)
}
cat("Detection flag columns found:", paste(flag_cols, collapse = ", "), "\n")

if (length(flag_cols) == 2) {
  det_both <- cci_detected %>%
    rename(FLAG1 = !!flag_cols[1], FLAG2 = !!flag_cols[2]) %>%
    count(FLAG1, FLAG2, name = "count") %>%
    mutate(label = case_when(
       FLAG1 &  FLAG2 ~ "Detected in BOTH",
       FLAG1 & !FLAG2 ~ paste0("Only ", cond1),
      !FLAG1 &  FLAG2 ~ paste0("Only ", cond2),
      TRUE            ~ "Detected in NONE"
    ))

  write.csv(det_both,
            outfn_table(paste0("detection_flags_contingency_", analysis_name, ".csv")),
            row.names = FALSE)

  p_det <- ggplot(det_both, aes(x = label, y = count, fill = label)) +
    geom_col(show.legend = FALSE) +
    geom_text(aes(label = count), vjust = -0.4) +
    labs(title = paste("CCI Detection Flags —", analysis_name), x = NULL, y = "Count") +
    theme_classic(base_size = 12)

  ggsave(outfn_plot(paste0("detection_flags_barplot_", analysis_name, ".png")),
         p_det, width = 6, height = 4, dpi = 300)
  print(p_det)
}


# ============================================================
# Close log
# ============================================================

cat("\n=== Run Complete:", format(Sys.time()), "===\n")
sink()
