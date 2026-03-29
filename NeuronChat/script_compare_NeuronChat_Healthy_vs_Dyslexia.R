# ══════════════════════════════════════════════════════════════════
# NEURONCHAT COMPARISON ANALYSIS
# Comparison : Healthy vs Dyslexia
# Cond 1     : Healthy  (reference)
# Cond 2     : Dyslexia (logFC > 0 = higher in Dyslexia vs Healthy)
# ══════════════════════════════════════════════════════════════════

# ── Libraries ────────────────────────────────────────────────────
library(NeuronChat)
library(CellChat)
library(uwot)
library(future)
library(progress)
library(Matrix)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(scales)
library(pheatmap)
library(ComplexHeatmap)
library(RColorBrewer)
library(gplots)
library(grid)
library(circlize)
options(stringsAsFactors = FALSE)

plan("sequential")

# ══════════════════════════════════════════════════════════════════
# KEY IDENTIFIERS — change these for each comparison
# ══════════════════════════════════════════════════════════════════

base_dir  <- "/mnt/nfs/CX000008_DS1/projects/jaeyeon/fastq_file_Dyslexia_r1/zanalysis_bogdan/samples.merged_AG_Harmony_res0.x_anno_031826"

name1     <- "Healthy"
name2     <- "Dyslexia"
col1      <- "#4DAF4A"   # green = Healthy
col2      <- "#377EB8"   # blue  = Dyslexia

rds1 <- file.path(base_dir,
  "18_merged_AG_031826.28celltypes.Healthy_neuronchat",
  "neuronchat_18_merged_AG_031826.28celltypes.Healthy_at_the_end.rds")

rds2 <- file.path(base_dir,
  "18_merged_AG_031826.28celltypes.Dyslexia_neuronchat",
  "neuronchat_18_merged_AG_031826.28celltypes.Dyslexia_at_the_end.rds")

work_dir   <- file.path(base_dir,
  paste0("compare_NeuronChat_", name1, "_vs_", name2))
output_dir <- file.path(work_dir, "neuronchat_diff_comp")

dir.create(work_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ══════════════════════════════════════════════════════════════════
# LOGGING
# ══════════════════════════════════════════════════════════════════

log_file   <- file.path(work_dir,
  paste0("neuronchat_compare_", name1, "_vs_", name2, "_",
         format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
start_time <- Sys.time()
sink(log_file, split = TRUE)

cat("╔═══════════════════════════════════════════════════════════╗\n")
cat(sprintf("║  NeuronChat Comparison: %s vs %s\n", name1, name2))
cat("╚═══════════════════════════════════════════════════════════╝\n")
cat("Started :", format(start_time, "%Y-%m-%d %H:%M:%S"), "\n\n")

on.exit({
  if (sink.number() > 0) {
    cat("\nEnded  :", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
    cat("Runtime:", round(difftime(Sys.time(), start_time,
                                   units = "mins"), 2), "minutes\n")
    sink()
  }
}, add = TRUE)

# ══════════════════════════════════════════════════════════════════
# HELPER FUNCTIONS
# ══════════════════════════════════════════════════════════════════

do_gc <- function() { gc(verbose = FALSE, full = TRUE); invisible(NULL) }

safe_dev_off <- function() {
  tryCatch({
    while (grDevices::dev.cur() > 1) grDevices::dev.off()
  }, error = function(e) invisible(NULL))
}

check_memory_usage <- function(threshold = 0.90) {
  meminfo <- tryCatch(readLines("/proc/meminfo", n = 20),
                      error = function(e) NULL)
  if (is.null(meminfo)) return(invisible(NA))
  memtotal_kb     <- as.numeric(gsub("[^0-9]", "", meminfo[1]))
  memavailable_kb <- as.numeric(gsub("[^0-9]", "", meminfo[3]))
  pct <- (memtotal_kb - memavailable_kb) / memtotal_kb
  cat(sprintf("💾 Memory: %.1f%% used (%.1f / %.1f GB)\n",
              pct * 100,
              (memtotal_kb - memavailable_kb) / 1024^2,
              memtotal_kb / 1024^2))
  invisible(pct)
}

# ══════════════════════════════════════════════════════════════════
# LOAD OBJECTS
# ══════════════════════════════════════════════════════════════════

cat("\n=== Loading NeuronChat Objects ===\n")
cat("RDS 1:", rds1, "\n")
cat("RDS 2:", rds2, "\n\n")

obj1 <- readRDS(rds1); cat(sprintf("✓ Loaded: %s\n", name1))
obj2 <- readRDS(rds2); cat(sprintf("✓ Loaded: %s\n", name2))

cat(sprintf("\n%s pathways : %d\n", name1, length(obj1@net)))
cat(sprintf("%s pathways : %d\n", name2, length(obj2@net)))
cat(sprintf("%s cell types: %d\n", name1, length(levels(obj1@idents))))
cat(sprintf("%s cell types: %d\n", name2, length(levels(obj2@idents))))

check_memory_usage()

# ── Fix meta slot if empty ────────────────────────────────────────
for (j in 1:2) {
  obj  <- if (j == 1) obj1 else obj2
  nm   <- if (j == 1) name1 else name2
  if (is.null(obj@meta) || nrow(obj@meta) == 0) {
    cb <- colnames(obj@data)
    obj@meta <- data.frame(cell_barcode = cb, dataset = nm,
                            sample = nm, row.names = cb,
                            stringsAsFactors = FALSE)
    cat(sprintf("✓ Meta created for %s (%d cells)\n", nm, length(cb)))
  }
  if (j == 1) obj1 <- obj else obj2 <- obj
}

# ── Cell type compatibility — subset to common cell types ─────────
# Avoids non-conformable arrays in computeNetSimilarityPairwise_Neuron
# and NULL errors in the cell-to-cell comparison table.
cat("\n=== CELL TYPE COMPATIBILITY CHECK ===\n")
ct1 <- sort(levels(obj1@idents))
ct2 <- sort(levels(obj2@idents))

cat(sprintf("Cell types in %s : %d\n", name1, length(ct1)))
cat(sprintf("Cell types in %s : %d\n", name2, length(ct2)))

only1     <- setdiff(ct1, ct2)
only2     <- setdiff(ct2, ct1)
common_ct <- intersect(ct1, ct2)

cat(sprintf("Only in %s   : %d : %s\n",
            name1, length(only1), paste(only1, collapse = ", ")))
cat(sprintf("Only in %s : %d : %s\n",
            name2, length(only2), paste(only2, collapse = ", ")))
cat(sprintf("Common cell types : %d\n", length(common_ct)))

if (length(only1) > 0 || length(only2) > 0) {
  cat("⚠️  Different cell type compositions — subsetting to common types\n")

  subset_neuronchat <- function(obj, keep_ct) {
    obj@net <- lapply(obj@net, function(mat) {
      idx <- rownames(mat) %in% keep_ct
      mat[idx, idx, drop = FALSE]
    })
    obj@idents <- factor(as.character(obj@idents)[
      as.character(obj@idents) %in% keep_ct], levels = keep_ct)
    keep_cells <- colnames(obj@data)[
      as.character(obj@idents) %in% keep_ct]
    if (length(keep_cells) > 0) {
      obj@data <- obj@data[, keep_cells, drop = FALSE]
      if (!is.null(obj@data.signaling) && ncol(obj@data.signaling) > 0)
        obj@data.signaling <- obj@data.signaling[
          , keep_cells[keep_cells %in% colnames(obj@data.signaling)],
          drop = FALSE]
    }
    obj
  }

  obj1 <- subset_neuronchat(obj1, common_ct)
  obj2 <- subset_neuronchat(obj2, common_ct)
  cat(sprintf("✓ %s after subset: %d cell types\n",
              name1, length(levels(obj1@idents))))
  cat(sprintf("✓ %s after subset: %d cell types\n",
              name2, length(levels(obj2@idents))))
} else {
  cat("✓ Same cell type compositions — no subset needed\n")
}

# cell type subset inserted below
cortex_list <- list(obj1, obj2)
names(cortex_list) <- c(name1, name2)

# ── All pathways (union) ──────────────────────────────────────────
all_pathways <- union(names(obj1@net), names(obj2@net))
cat(sprintf("\nTotal pathways (union): %d\n", length(all_pathways)))

# ══════════════════════════════════════════════════════════════════
# PART I: INDIVIDUAL PATHWAY CIRCLE PLOTS
# ══════════════════════════════════════════════════════════════════

cat(sprintf("\n=== Part I: Per-Pathway Circle Plots (%d pathways) ===\n",
            length(all_pathways)))

pb <- progress_bar$new(
  format = "[:bar] :current/:total (:percent) | :pathway",
  total = length(all_pathways), clear = FALSE, width = 80)

for (pathway in all_pathways) {
  pb$tick(tokens = list(pathway = sprintf("%-25s", pathway)))

  tryCatch({
    png(file.path(output_dir,
                  paste0(pathway, "_circle_comparison.png")),
        width = 1400, height = 700, res = 100)
    par(mfrow = c(1, 2))
    for (j in 1:2) {
      nm <- names(cortex_list)[j]
      if (pathway %in% names(cortex_list[[j]]@net)) {
        netVisual_circle_neuron(cortex_list[[j]]@net[[pathway]],
                                title.name = paste(pathway, "--", nm),
                                arrow.size = 0.5, margin = 0.3,
                                edge.width.max = 8)
      } else {
        plot.new()
        title(main = paste(pathway, "--", nm, "\n(Not detected)"))
      }
    }
    dev.off()
  }, error = function(e) { safe_dev_off() })
}
do_gc()
cat("✓ Per-pathway circle plots saved\n")

# ══════════════════════════════════════════════════════════════════
# PART II: AGGREGATED NETWORK PLOTS
# ══════════════════════════════════════════════════════════════════

cat("\n=== Part II: Aggregated Network Plots ===\n")

for (method in c("count", "weight")) {
  tryCatch({
    png(file.path(output_dir,
                  paste0("fig_aggregated_network_", method,
                         "_comparison.png")),
        width = 1400, height = 700, res = 100)
    par(mfrow = c(1, 2))
    for (j in 1:2) {
      net_agg <- net_aggregation(cortex_list[[j]]@net, method = method)
      netVisual_circle_neuron(net_agg,
                              title.name = names(cortex_list)[j],
                              arrow.size = 0.5, margin = 0.3,
                              edge.width.max = 8)
    }
    dev.off()
    cat(sprintf("✓ Aggregated network (%s) saved\n", method))
  }, error = function(e) {
    safe_dev_off()
    cat(sprintf("⚠️  Aggregated network (%s) failed: %s\n",
                method, as.character(e)))
  })
}

# ══════════════════════════════════════════════════════════════════
# PART III: MERGE + COMPARE INTERACTIONS
# ══════════════════════════════════════════════════════════════════

cat("\n=== Part III: Merge + Compare Interactions ===\n")

neuronchat_list <- tryCatch({
  mergeNeuronChat(cortex_list, add.names = names(cortex_list))
}, error = function(e) {
  cat("⚠️  mergeNeuronChat failed:", as.character(e), "\n")
  NULL
})

if (!is.null(neuronchat_list)) {
  cat("✓ Merged NeuronChat object\n")

  # compareInteractions barplot
  tryCatch({
    p1 <- compareInteractions_Neuron(neuronchat_list,
                                      measure = "count",
                                      comparison = c(1, 2),
                                      group = c(1, 2),
                                      show.legend = FALSE)
    p2 <- compareInteractions_Neuron(neuronchat_list,
                                      measure = "weight",
                                      comparison = c(1, 2),
                                      group = c(1, 2),
                                      show.legend = FALSE)
    ggsave(file.path(output_dir,
                     "fig_compareInteractions_count_weight.png"),
           p1 + p2, width = 10, height = 5, dpi = 300, bg = "white")
    cat("✓ compareInteractions saved\n")
  }, error = function(e) {
    cat("⚠️  compareInteractions failed:", as.character(e), "\n")
  })

  # rankNet_Neuron
  tryCatch({
    g1 <- rankNet_Neuron(neuronchat_list, mode = "comparison",
                          measure = "count", comparison = 1:2,
                          do.stat = FALSE, tol = 0.1,
                          stacked = FALSE, font.size = 11)
    g2 <- rankNet_Neuron(neuronchat_list, mode = "comparison",
                          measure = "weight", comparison = 1:2,
                          do.stat = FALSE, tol = 0.1,
                          stacked = FALSE, font.size = 11)
    ggsave(file.path(output_dir,
                     "fig_rankNet_count_weight.png"),
           g1 + g2, width = 14, height = 16, dpi = 300, bg = "white")
    cat("✓ rankNet plots saved\n")
  }, error = function(e) {
    cat("⚠️  rankNet failed:", as.character(e), "\n")
  })

  do_gc()
  check_memory_usage()

  # ══════════════════════════════════════════════════════════════
  # PART IV: FUNCTIONAL SIMILARITY + UMAP
  # ══════════════════════════════════════════════════════════════

  cat("\n=== Part IV: Functional Similarity + UMAP ===\n")

  tryCatch({
    neuronchat_list <- computeNetSimilarityPairwise_Neuron(
      neuronchat_list, slot.name = "net",
      type = "functional", comparison = c(1, 2))
    cat("✓ Similarity computed\n")

    neuronchat_list <- netEmbedding(neuronchat_list,
                                     slot.name = "net_analysis",
                                     type = "functional",
                                     umap.method = "uwot",
                                     comparison = c(1, 2))
    cat("✓ UMAP embedding done\n")

    neuronchat_list <- netClustering(neuronchat_list,
                                      slot.name = "net_analysis",
                                      type = "functional",
                                      comparison = c(1, 2), k = 5)
    cat("✓ Clustering done\n")

    # Embedding plot
    tryCatch({
      png(file.path(output_dir,
                    "fig_functional_embedding_pairwise.png"),
          width = 1200, height = 1000, res = 150)
      netVisual_embeddingPairwise_Neuron(neuronchat_list,
                                         slot.name = "net_analysis",
                                         type = "functional",
                                         label.size = 3.5,
                                         comparison = c(1, 2),
                                         pathway.remove.show = FALSE,
                                         pathway.labeled = FALSE)
      dev.off()
      cat("✓ Embedding plot saved\n")
    }, error = function(e) { safe_dev_off() })

    # ZoomIn plot
    tryCatch({
      png(file.path(output_dir,
                    "fig_functional_embedding_zoomin.png"),
          width = 1600, height = 800, res = 150)
      netVisual_embeddingPairwiseZoomIn_Neuron(neuronchat_list,
                                                slot.name = "net_analysis",
                                                type = "functional",
                                                label.size = 5,
                                                comparison = c(1, 2),
                                                nCol = 3)
      dev.off()
      cat("✓ ZoomIn plot saved\n")
    }, error = function(e) { safe_dev_off() })

    # Similarity matrix heatmap
    matrix_data <- neuronchat_list@net_analysis$similarity$functional$matrix[["1-2"]]
    if (!is.null(matrix_data)) {
      cat(sprintf("Similarity matrix: %d x %d\n",
                  nrow(matrix_data), ncol(matrix_data)))

      tryCatch({
        png(file.path(output_dir,
                      "fig_similarity_heatmap_full.png"),
            width = 6000, height = 6000, res = 100)
        pheatmap(matrix_data, cluster_rows = TRUE, cluster_cols = TRUE,
                 color = colorRampPalette(c("white", "red"))(100),
                 border_color = NA, fontsize_row = 10, fontsize_col = 10,
                 main = sprintf("Functional Similarity: %s vs %s",
                                name1, name2))
        dev.off()
        cat("✓ Similarity heatmap (full) saved\n")
      }, error = function(e) { safe_dev_off() })

      # Filtered heatmap
      tryCatch({
        rows_keep <- rowSums(matrix_data != 0) > 1
        cols_keep <- colSums(matrix_data != 0) > 1
        mat_filt  <- matrix_data[rows_keep, cols_keep]
        cat(sprintf("Filtered matrix: %d x %d\n",
                    nrow(mat_filt), ncol(mat_filt)))

        png(file.path(output_dir,
                      "fig_similarity_heatmap_filtered.png"),
            width = 3000, height = 3000, res = 100)
        pheatmap(mat_filt, cluster_rows = TRUE, cluster_cols = TRUE,
                 color = colorRampPalette(c("white", "red"))(100),
                 border_color = NA, fontsize_row = 10, fontsize_col = 10,
                 main = sprintf("Functional Similarity (filtered): %s vs %s",
                                name1, name2))
        dev.off()
        cat("✓ Similarity heatmap (filtered) saved\n")
      }, error = function(e) { safe_dev_off() })

      # Save similarity matrix CSV
      write.csv(matrix_data,
                file.path(output_dir,
                          "res_functional_similarity_matrix.csv"),
                row.names = TRUE)

      # ── UMAP distance analysis ──────────────────────────────────
      cat("\n--- UMAP Distance Analysis ---\n")

      embedding_data <- neuronchat_list@net_analysis$similarity$functional$dr[["1-2"]]
      cluster_info   <- neuronchat_list@net_analysis$similarity$functional$group[["1-2"]]

      if (!is.null(embedding_data) && !is.null(cluster_info)) {
        pathway_names   <- rownames(embedding_data)
        names1_paths    <- grep(paste0("--", name1, "$"),
                                pathway_names, value = TRUE)
        names2_paths    <- grep(paste0("--", name2, "$"),
                                pathway_names, value = TRUE)
        base1 <- gsub(paste0("--", name1, "$"), "", names1_paths)
        base2 <- gsub(paste0("--", name2, "$"), "", names2_paths)
        shared_paths <- intersect(base1, base2)

        distances <- data.frame(
          pathway       = character(),
          distance      = numeric(),
          cluster_cond1 = integer(),
          cluster_cond2 = integer(),
          stringsAsFactors = FALSE)

        for (pw in union(base1, base2)) {
          n1 <- paste0(pw, "--", name1)
          n2 <- paste0(pw, "--", name2)
          c1 <- if (n1 %in% rownames(embedding_data)) embedding_data[n1, ] else NA
          c2 <- if (n2 %in% rownames(embedding_data)) embedding_data[n2, ] else NA
          d  <- if (!any(is.na(c1)) && !any(is.na(c2))) {
            sqrt(sum((c1 - c2)^2))
          } else NA
          k1 <- if (n1 %in% names(cluster_info)) cluster_info[n1] else NA
          k2 <- if (n2 %in% names(cluster_info)) cluster_info[n2] else NA
          distances <- rbind(distances, data.frame(
            pathway = pw, distance = d,
            cluster_cond1 = k1, cluster_cond2 = k2))
        }

        distances <- distances[order(distances$distance,
                                     decreasing = TRUE, na.last = TRUE), ]

        write.csv(distances,
                  file.path(output_dir,
                            "res_UMAP_pathway_distances.csv"),
                  row.names = FALSE)
        cat(sprintf("✓ UMAP distances saved (%d pathways)\n",
                    nrow(distances)))

        # Bar plot — all pathways
        tryCatch({
          p_dist <- ggplot(distances[!is.na(distances$distance), ],
                           aes(x = reorder(pathway, distance),
                               y = distance)) +
            geom_bar(stat = "identity", fill = "steelblue") +
            coord_flip() +
            labs(title = sprintf("UMAP Distance: %s vs %s", name1, name2),
                 x = "Pathway", y = "Euclidean Distance") +
            theme_minimal(base_size = 12)

          ggsave(file.path(output_dir,
                           "fig_UMAP_pathway_distances.png"),
                 p_dist,
                 width = 7,
                 height = max(6, sum(!is.na(distances$distance)) * 0.12),
                 dpi = 300, bg = "white", limitsize = FALSE)
          cat("✓ UMAP distance bar plot saved\n")
        }, error = function(e) invisible(NULL))

        # Arrow UMAP plot
        tryCatch({
          emb_df <- as.data.frame(embedding_data)
          colnames(emb_df) <- c("UMAP1", "UMAP2")
          emb_df$pathway_full  <- rownames(embedding_data)
          emb_df$condition     <- ifelse(
            grepl(paste0("--", name1, "$"), emb_df$pathway_full),
            name1, name2)
          emb_df$pathway_clean <- gsub(
            paste0("--", name1, "$|--", name2, "$"), "",
            emb_df$pathway_full)
          emb_df$cluster <- sapply(emb_df$pathway_full,
                                   function(x) cluster_info[x])

          high_dist <- distances[!is.na(distances$distance) &
                                   distances$distance > 1, ]

          arrow_df <- lapply(seq_len(nrow(high_dist)), function(i) {
            pw <- high_dist$pathway[i]
            n1 <- paste0(pw, "--", name1)
            n2 <- paste0(pw, "--", name2)
            if (!n1 %in% rownames(emb_df) ||
                !n2 %in% rownames(emb_df)) return(NULL)
            data.frame(
              pathway  = pw,
              x_start  = emb_df[n1, "UMAP1"],
              y_start  = emb_df[n1, "UMAP2"],
              x_end    = emb_df[n2, "UMAP1"],
              y_end    = emb_df[n2, "UMAP2"],
              x_mid    = mean(c(emb_df[n1, "UMAP1"],
                                emb_df[n2, "UMAP1"])),
              y_mid    = mean(c(emb_df[n1, "UMAP2"],
                                emb_df[n2, "UMAP2"])),
              distance = high_dist$distance[i])
          })
          arrow_df <- do.call(rbind, Filter(Negate(is.null), arrow_df))

          if (!is.null(arrow_df) && nrow(arrow_df) > 0) {
            p_arrow <- ggplot() +
              geom_point(data = emb_df,
                         aes(UMAP1, UMAP2, color = condition,
                             shape = condition),
                         size = 3, alpha = 0.5) +
              geom_segment(data = arrow_df,
                           aes(x = x_start, y = y_start,
                               xend = x_end, yend = y_end),
                           arrow = arrow(length = unit(0.2, "cm")),
                           color = "black", alpha = 0.6) +
              geom_text_repel(data = arrow_df,
                              aes(x = x_mid, y = y_mid,
                                  label = pathway),
                              size = 3, max.overlaps = Inf) +
              scale_color_manual(values = c(col1, col2),
                                 labels = c(name1, name2)) +
              scale_shape_manual(values = c(16, 17),
                                 labels = c(name1, name2)) +
              labs(title = sprintf("Functional shift: %s → %s",
                                   name1, name2),
                   subtitle = "Arrows: pathways with UMAP distance > 1") +
              theme_bw(base_size = 13) +
              coord_fixed()

            ggsave(file.path(output_dir,
                             "fig_UMAP_arrows_high_distance.png"),
                   p_arrow, width = 12, height = 10,
                   dpi = 300, bg = "white")
            cat("✓ UMAP arrow plot saved\n")
          }
        }, error = function(e) invisible(NULL))

        # Save UMAP coordinates
        write.csv(as.data.frame(embedding_data),
                  file.path(output_dir,
                            "res_UMAP_coordinates.csv"),
                  row.names = TRUE)
      }
    }

    # Pattern cluster heatmaps
    cat("\n--- Pattern Cluster Heatmaps ---\n")
    tryCatch({
      net12    <- neuronchat_list@net[c(1, 2)]
      net_c1   <- net12[[1]]
      net_c2   <- net12[[2]]
      names(net_c1) <- paste0(names(net_c1), "--", name1)
      names(net_c2) <- paste0(names(net_c2), "--", name2)
      net12_list_p  <- append(net_c1, net_c2)

      int_grp  <- neuronchat_list@net_analysis$similarity$functional$group[["1-2"]]
      n_clust  <- length(sort(unique(int_grp)))

      # Stacked heatmap
      hmap_list <- list()
      for (j in seq_len(n_clust)) {
        mat_j <- net_aggregation(
          net12_list_p[names(int_grp[int_grp == j])],
          method = "weight")
        hmap_list[[j]] <- Heatmap(
          mat_j, name = paste("Weight", j),
          col = brewer.pal(8, "YlOrBr"),
          cluster_rows = FALSE, cluster_columns = FALSE,
          column_title = paste("Pattern Cluster", j),
          column_title_gp = gpar(fontsize = 14, fontface = "bold"),
          row_title = "Sender", column_names_rot = 45, border = TRUE)
      }

      ht_stack <- Reduce(`%v%`, hmap_list)

      png(file.path(output_dir,
                    "fig_pattern_clusters_heatmaps.png"),
          width = 2400, height = 3200, res = 150, bg = "white")
      draw(ht_stack)
      dev.off()

      pdf(file.path(output_dir,
                    "fig_pattern_clusters_heatmaps.pdf"),
          width = 12, height = 16)
      draw(ht_stack)
      dev.off()
      cat("✓ Pattern cluster heatmaps saved\n")
    }, error = function(e) {
      safe_dev_off()
      cat("⚠️  Pattern cluster heatmaps failed:", as.character(e), "\n")
    })

  }, error = function(e) {
    safe_dev_off()
    cat("⚠️  Functional similarity section failed:", as.character(e), "\n")
  })

  do_gc()
  check_memory_usage()

  # ══════════════════════════════════════════════════════════════
  # PART V: DIFFERENTIAL PATHWAY ANALYSIS
  # ══════════════════════════════════════════════════════════════

  cat("\n=== Part V: Differential Pathway Analysis ===\n")

  # ── Get raw nets from cortex_list (not suffixed) ──────────────
  net_c1 <- cortex_list[[1]]@net
  net_c2 <- cortex_list[[2]]@net

  # Determine cell types from any available pathway
  any_pw <- all_pathways[1]
  ref_net <- if (any_pw %in% names(net_c1)) net_c1 else net_c2
  celltypes <- rownames(ref_net[[any_pw]])
  cat(sprintf("Cell types (%d): %s\n", length(celltypes),
              paste(celltypes, collapse = ", ")))

  # ── Build long comparison table ───────────────────────────────
  cat("Building cell-to-cell comparison table...\n")

  cmp_df_list <- lapply(all_pathways, function(p) {
    tryCatch({
      A <- if (p %in% names(net_c1)) net_c1[[p]] else {
        m <- matrix(0, length(celltypes), length(celltypes))
        rownames(m) <- celltypes; colnames(m) <- celltypes; m }
      B <- if (p %in% names(net_c2)) net_c2[[p]] else {
        m <- matrix(0, length(celltypes), length(celltypes))
        rownames(m) <- celltypes; colnames(m) <- celltypes; m }
      ct <- intersect(rownames(A), rownames(B))
      if (length(ct) == 0) return(NULL)
      A  <- A[ct, ct, drop = FALSE]
      B  <- B[ct, ct, drop = FALSE]
      df <- expand.grid(sender = rownames(A), receiver = colnames(A),
                        stringsAsFactors = FALSE)
      df$pathway <- p
      df$cond1   <- as.vector(A)
      df$cond2   <- as.vector(B)
      df$delta   <- df$cond2 - df$cond1
      df$mean_w  <- (df$cond1 + df$cond2) / 2
      df$log2FC  <- log2((df$cond2 + 1e-6) / (df$cond1 + 1e-6))
      df
    }, error = function(e) NULL)
  })

  cmp_df_list <- Filter(Negate(is.null), cmp_df_list)
  cmp_df      <- do.call(rbind, cmp_df_list)
  rm(cmp_df_list); do_gc()

  names(cmp_df)[names(cmp_df) == "cond1"] <- name1
  names(cmp_df)[names(cmp_df) == "cond2"] <- name2

  cat(sprintf("✓ Comparison table: %d rows (%d pathways × %d cell pairs)\n",
              nrow(cmp_df), length(all_pathways), length(celltypes)^2))

  write.csv(cmp_df,
            file.path(output_dir,
                      sprintf("res_cell_to_cell_%s_vs_%s.csv",
                              name1, name2)),
            row.names = FALSE)
  cat("✓ Cell-to-cell table saved\n")

  # ── Pathway-level summary ──────────────────────────────────────
  cat("\n--- Pathway-Level Summary ---\n")

  pathway_rank <- cmp_df %>%
    filter(mean_w > 0.01) %>%
    group_by(pathway) %>%
    summarise(
      n_edges       = n(),
      log2FC_wmean  = weighted.mean(log2FC, w = mean_w, na.rm = TRUE),
      delta_sum     = sum(delta, na.rm = TRUE),
      .groups       = "drop") %>%
    arrange(desc(abs(log2FC_wmean)))

  top_up   <- pathway_rank %>%
    arrange(desc(log2FC_wmean)) %>% head(30) %>%
    mutate(direction = sprintf("UP (%s > %s)", name2, name1))
  top_down <- pathway_rank %>%
    arrange(log2FC_wmean) %>% head(30) %>%
    mutate(direction = sprintf("DOWN (%s < %s)", name2, name1))

  display_tbl <- bind_rows(top_up, top_down)

  write.csv(display_tbl,
            file.path(output_dir,
                      sprintf("res_top_pathways_%s_vs_%s.csv",
                              name1, name2)),
            row.names = FALSE)

  cat(sprintf("✓ Top pathways saved\n"))
  cat(sprintf("  UP   in %s: %d\n", name2, nrow(top_up)))
  cat(sprintf("  DOWN in %s: %d\n", name2, nrow(top_down)))

  # Bar plot
  tryCatch({
    max_abs <- max(abs(display_tbl$log2FC_wmean), na.rm = TRUE)

    p_bar <- ggplot(display_tbl,
                    aes(x = reorder(pathway, log2FC_wmean),
                        y = log2FC_wmean, fill = direction)) +
      geom_col(width = 0.7, color = "grey30", linewidth = 0.2) +
      coord_flip() +
      geom_hline(yintercept = 0, linewidth = 0.4) +
      scale_y_continuous(limits = c(-max_abs * 1.05,
                                    max_abs * 1.05)) +
      scale_fill_manual(
        values = setNames(c("#b2182b", "#2166ac"),
                          c(sprintf("UP (%s > %s)", name2, name1),
                            sprintf("DOWN (%s < %s)", name2, name1)))) +
      labs(title = sprintf("Differential NeuronChat Pathways\n%s vs %s",
                           name1, name2),
           subtitle = "Weighted mean log2FC (pathway level)",
           x = NULL,
           y = sprintf("Weighted mean log2FC (%s / %s)", name2, name1),
           fill = NULL) +
      theme_bw(base_size = 13) +
      theme(panel.grid.major.y = element_blank(),
            legend.position = "top")

    ggsave(file.path(output_dir,
                     sprintf("fig_pathway_log2FC_%s_vs_%s.png",
                             name1, name2)),
           p_bar, width = 8, height = 10,
           dpi = 300, bg = "white")
    cat("✓ Pathway log2FC bar plot saved\n")
  }, error = function(e) invisible(NULL))

  # ── Scatter plot: cond1 vs cond2 strength ─────────────────────
  cat("\n--- Pathway Strength Scatter ---\n")
  tryCatch({
    pw_strength <- cmp_df %>%
      group_by(pathway) %>%
      summarise(
        s1     = sum(.data[[name1]], na.rm = TRUE),
        s2     = sum(.data[[name2]], na.rm = TRUE),
        log2FC = log2((s2 + 1e-6) / (s1 + 1e-6)),
        .groups = "drop") %>%
      mutate(
        color_group = case_when(
          log2FC > 0 ~ sprintf("UP in %s", name2),
          log2FC < 0 ~ sprintf("DOWN in %s", name2),
          TRUE ~ "No change"),
        label_text = ifelse(s1 > 1 & s2 > 1, pathway, ""))

    colnames(pw_strength)[colnames(pw_strength) == "s1"] <- name1
    colnames(pw_strength)[colnames(pw_strength) == "s2"] <- name2

    p_scatter <- ggplot(pw_strength,
                        aes(x = .data[[name1]],
                            y = .data[[name2]])) +
      geom_point(aes(color = color_group), size = 4, alpha = 0.85) +
      geom_text_repel(aes(label = label_text),
                      size = 3.5, max.overlaps = Inf,
                      box.padding = 0.4) +
      geom_abline(slope = 1, intercept = 0,
                  linetype = "dashed", color = "grey50") +
      scale_color_manual(
        values = setNames(c("#b2182b", "#2166ac", "grey70"),
                          c(sprintf("UP in %s", name2),
                            sprintf("DOWN in %s", name2),
                            "No change"))) +
      labs(title = sprintf("Pathway Strength: %s vs %s", name1, name2),
           x = sprintf("Total strength (%s)", name1),
           y = sprintf("Total strength (%s)", name2),
           color = NULL) +
      theme_minimal(base_size = 13) +
      theme(plot.background = element_rect(fill = "white"))

    ggsave(file.path(output_dir,
                     sprintf("fig_pathway_scatter_%s_vs_%s.png",
                             name1, name2)),
           p_scatter, width = 10, height = 8,
           dpi = 300, bg = "white")
    cat("✓ Scatter plot saved\n")
  }, error = function(e) invisible(NULL))

  # ── Per-pathway bubble plots ───────────────────────────────────
  cat("\n--- Per-Pathway Bubble Plots ---\n")
  bubble_dir <- file.path(output_dir, "pathway_bubble_plots")
  dir.create(bubble_dir, showWarnings = FALSE)

  max_abs_delta <- quantile(abs(cmp_df$delta), 0.95, na.rm = TRUE)

  pb2 <- progress_bar$new(
    format = "[:bar] :current/:total (:percent) | :pathway",
    total = length(all_pathways), clear = FALSE, width = 80)

  for (pw in all_pathways) {
    pb2$tick(tokens = list(pathway = sprintf("%-25s", pw)))
    df_pw <- cmp_df %>% filter(pathway == pw)
    if (all(df_pw$mean_w == 0)) next
    tryCatch({
      p_pw <- ggplot(df_pw, aes(x = receiver, y = sender)) +
        geom_point(aes(size = mean_w, color = delta), alpha = 0.9) +
        scale_color_gradient2(low = col1, mid = "white", high = col2,
                              midpoint = 0,
                              limits = c(-max_abs_delta, max_abs_delta),
                              oob = squish,
                              name = sprintf("%s − %s", name2, name1)) +
        scale_size_continuous(name = "Mean weight", range = c(1, 8)) +
        theme_bw(base_size = 12) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = sprintf("NeuronChat: %s\n(%s − %s)",
                             pw, name2, name1),
             x = "Receiver", y = "Sender")
      ggsave(file.path(bubble_dir,
                       paste0(pw, "_bubble.png")),
             p_pw, width = 8, height = 7,
             dpi = 250, bg = "white")
    }, error = function(e) invisible(NULL))
  }
  cat("\n✓ Per-pathway bubble plots saved\n")

  # ── Difference heatmaps ───────────────────────────────────────
  cat("\n--- Difference Heatmaps per Pathway ---\n")
  heatmap_dir <- file.path(output_dir, "pathway_difference_heatmaps")
  dir.create(heatmap_dir, showWarnings = FALSE)

  pw_scores <- cmp_df %>%
    group_by(pathway) %>%
    summarise(
      s1        = sum(.data[[name1]], na.rm = TRUE),
      s2        = sum(.data[[name2]], na.rm = TRUE),
      delta_sum = sum(delta, na.rm = TRUE),
      .groups   = "drop") %>%
    filter(abs(delta_sum) > 0) %>%
    arrange(desc(abs(delta_sum))) %>%
    mutate(rank = row_number())

  write.csv(pw_scores,
            file.path(output_dir,
                      sprintf("res_pathway_delta_ranking_%s_vs_%s.csv",
                              name1, name2)),
            row.names = FALSE)

  pb3 <- progress_bar$new(
    format = "[:bar] :current/:total (:percent) | :pathway",
    total = nrow(pw_scores), clear = FALSE, width = 80)

  for (i in seq_len(nrow(pw_scores))) {
    pw       <- pw_scores$pathway[i]
    rk       <- pw_scores$rank[i]
    d_sum    <- pw_scores$delta_sum[i]
    pb3$tick(tokens = list(pathway = sprintf("%-25s", pw)))

    df_pw    <- cmp_df %>% filter(pathway == pw)
    diff_mat <- matrix(df_pw$delta,
                       nrow = length(unique(df_pw$sender)),
                       ncol = length(unique(df_pw$receiver)))
    rownames(diff_mat) <- unique(df_pw$sender)
    colnames(diff_mat) <- unique(df_pw$receiver)

    max_d <- max(abs(diff_mat), na.rm = TRUE)
    if (is.na(max_d) || max_d == 0) next

    rank_pad <- sprintf("%03d", rk)
    pw_clean <- gsub("[^A-Za-z0-9_-]", "_", pw)
    fname    <- file.path(heatmap_dir,
                          paste0("rank_", rank_pad, "_",
                                 pw_clean, ".png"))

    tryCatch({
      png(fname, width = 1600, height = 1200, res = 150, bg = "white")
      pheatmap::pheatmap(diff_mat,
               main = sprintf("Rank %d: %s\n(%s − %s | Δ = %.3f)",
                              rk, pw, name2, name1, d_sum),
               color = colorRampPalette(c(col1, "white", col2))(100),
               breaks = seq(-max_d, max_d, length.out = 101),
               cluster_rows = FALSE, cluster_cols = FALSE,
               fontsize = 11, angle_col = 45, border_color = "grey80")
      dev.off()
    }, error = function(e) { safe_dev_off() })
  }
  cat("\n✓ Difference heatmaps saved\n")

  do_gc()
  check_memory_usage()

  # ── Save pathway matrices ─────────────────────────────────────
  cat("\n--- Saving Pathway Matrices as CSV ---\n")
  matrix_dir <- file.path(output_dir, "pathway_matrices")
  dir.create(matrix_dir, showWarnings = FALSE)

  for (pw in all_pathways) {
    pw_clean <- gsub("[^A-Za-z0-9_-]", "_", pw)
    if (pw %in% names(net_c1))
      write.csv(net_c1[[pw]],
                file.path(matrix_dir,
                          paste0(pw_clean, "_", name1, ".csv")),
                row.names = TRUE)
    if (pw %in% names(net_c2))
      write.csv(net_c2[[pw]],
                file.path(matrix_dir,
                          paste0(pw_clean, "_", name2, ".csv")),
                row.names = TRUE)
  }
  cat(sprintf("✓ Pathway matrices saved (%d files)\n",
              length(all_pathways) * 2))

  # ── Save merged object ────────────────────────────────────────
  saveRDS(neuronchat_list,
          file.path(work_dir,
                    sprintf("neuronchat_merged_%s_vs_%s.rds",
                            name1, name2)))
  cat("✓ Merged NeuronChat object saved\n")
}

# ══════════════════════════════════════════════════════════════════
# WRAP UP
# ══════════════════════════════════════════════════════════════════

cat("\n")
cat("╔═══════════════════════════════════════════════════════════╗\n")
cat(sprintf("║  ✓✓✓  NeuronChat %s vs %s COMPLETE  ✓✓✓\n",
            name1, name2))
cat("╚═══════════════════════════════════════════════════════════╝\n")
cat("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Output   :", output_dir, "\n")

sink()
