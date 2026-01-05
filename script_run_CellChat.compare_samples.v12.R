# Set Python path BEFORE loading any libraries that use reticulate
# python_path <- "/mnt/nfs/CX000008_DS1/projects/btanasa/virtual_env/cellbender/bin/python"

# Verify the path exists
# if (file.exists(python_path)) {
#  Sys.setenv(RETICULATE_PYTHON = python_path)
#  cat("✓ Python environment set to:", python_path, "\n")
# } else {
#  stop("Python path not found:", python_path)
# }

# a note about the files output :

#    554 adctl_LR_bothDown.csv
#    331 adctl_LR_bothUp.csv
#    480 addys_LR_bothDown.csv
#    519 addys_LR_bothUp.csv
#  12633 ligand.receptor.net.csv
#    558 netMappingDEG.net_downregulated_DYS.csv
#    519 netMappingDEG.net_upregulated_DYS.csv
#  12633 results.netMappingDEG.ligand.receptor.net.csv

#    554 adctl_LR_bothDown.csv
#    558 netMappingDEG.net_downregulated_DYS.csv

#    519 addys_LR_bothUp.csv
#    519 netMappingDEG.net_upregulated_DYS.csv

# from another script running cellchat -v20 :

# conda activate cellbender

# Set Python path BEFORE loading any libraries that use reticulate
# Must be set before reticulate initializes (which happens on first import)
# python_path <- "~/.virtualenvs/r-reticulate/bin/python"
# python_path_expanded <- path.expand(python_path)

# if (file.exists(python_path_expanded)) {
# Sys.setenv(RETICULATE_PYTHON = python_path_expanded)
#  cat("✓ Python environment set to:", python_path_expanded, "\n")
# } else {
#  cat("⚠️  Custom Python path not found:", python_path_expanded, "\n")
#  cat("   Will use system default Python\n")
# }

# Test import
# umap_learn <- import("umap")
# print(paste("UMAP version:", umap_learn$`__version__`))

# ═══════════════════════════════════════════════════════════
# PYTHON SETUP FOR RETICULATE
# ═══════════════════════════════════════════════════════════
py <- "/mnt/nfs/CX000008_DS1/projects/btanasa/conda_envs/cellbender/bin/python"
stopifnot(file.exists(py))  # Fail early if Python not found

library(reticulate)
use_python(py, required = TRUE)
py_config()

# ═══════════════════════════════════════════════════════════
# LOAD REQUIRED LIBRARIES
# ═══════════════════════════════════════════════════════════
library(CellChat)
library(patchwork)
library(tidyverse)      # Includes dplyr, ggplot2, tidyr, etc.
library("uwot")
library(ggrepel)
library(pheatmap)
library(ComplexHeatmap)

library(future)
plan("sequential")
# Fix the multiprocess issue
# plan("multisession", workers = 4)  # or however many cores you want to use

# ═══════════════════════════════════════════════════════════
# START LOGGING
# ═══════════════════════════════════════════════════════════
log_file <- paste0("cellchat_analysis_log_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt")
sink(log_file, split = TRUE)  # split=TRUE shows output in console AND saves to file

cat("=================================================\n")
cat("CellChat Analysis Log\n")
cat("Started at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("=================================================\n\n")

# ═══════════════════════════════════════════════════════════
# MEMORY MONITORING FUNCTION
# ═══════════════════════════════════════════════════════════

check_memory_usage <- function(threshold = 0.90, stop_on_exceed = TRUE) {
  # Read memory info from /proc/meminfo (Linux)
  meminfo <- tryCatch({
    readLines("/proc/meminfo", n = 20)  # Read more lines to get Buffers and Cached
  }, error = function(e) {
    # Fallback: try using 'free' command
    tryCatch({
      free_output <- system("free -k", intern = TRUE)
      return(free_output)
    }, error = function(e2) {
      cat("⚠️  Could not read memory information\n")
      return(NULL)
    })
  })
  
  if (is.null(meminfo)) {
    return(list(percent_used = NA, should_stop = FALSE))
  }
  
  # Parse /proc/meminfo format
  if (length(meminfo) >= 3 && grepl("MemTotal", meminfo[1])) {
    memtotal_kb <- as.numeric(gsub("[^0-9]", "", meminfo[1]))
    memavailable_kb <- as.numeric(gsub("[^0-9]", "", meminfo[3]))
    
    if (is.na(memavailable_kb) || memavailable_kb == 0) {
      # If MemAvailable not available, calculate from MemFree + Buffers + Cached
      memfree_kb <- as.numeric(gsub("[^0-9]", "", meminfo[2]))
      buffers_kb <- 0
      cached_kb <- 0
      
      # Try to get buffers and cached
      if (length(meminfo) > 3) {
        for (line in meminfo[4:min(10, length(meminfo))]) {
          if (grepl("^Buffers:", line)) {
            buffers_kb <- as.numeric(gsub("[^0-9]", "", line))
          }
          if (grepl("^Cached:", line)) {
            cached_kb <- as.numeric(gsub("[^0-9]", "", line))
          }
        }
      }
      memavailable_kb <- memfree_kb + buffers_kb + cached_kb
    }
    
    memused_kb <- memtotal_kb - memavailable_kb
    percent_used <- memused_kb / memtotal_kb
    
  } else {
    # Parse 'free' command output
    if (length(meminfo) >= 2) {
      mem_line <- strsplit(meminfo[2], "\\s+")[[1]]
      if (length(mem_line) >= 3) {
        memtotal_kb <- as.numeric(mem_line[2])
        memused_kb <- as.numeric(mem_line[3])
        if (length(mem_line) >= 4) {
          memavailable_kb <- as.numeric(mem_line[4])
        } else {
          memavailable_kb <- memtotal_kb - memused_kb
        }
        percent_used <- memused_kb / memtotal_kb
      } else {
        return(list(percent_used = NA, should_stop = FALSE))
      }
    } else {
      return(list(percent_used = NA, should_stop = FALSE))
    }
  }
  
  # Format output
  memtotal_gb <- memtotal_kb / 1024 / 1024
  memused_gb <- memused_kb / 1024 / 1024
  memavailable_gb <- memavailable_kb / 1024 / 1024
  
  cat(sprintf("Memory Status: %.1f%% used (%.1f GB / %.1f GB total, %.1f GB available)\n",
              percent_used * 100, memused_gb, memtotal_gb, memavailable_gb))
  
  should_stop <- FALSE
  if (percent_used >= threshold) {
    cat(sprintf("⚠️  WARNING: Memory usage (%.1f%%) exceeds threshold (%.1f%%)\n",
                percent_used * 100, threshold * 100))
    if (stop_on_exceed) {
      cat("❌ STOPPING SCRIPT: Memory usage too high!\n")
      should_stop <- TRUE
    }
  }
  
  return(list(
    percent_used = percent_used,
    memtotal_gb = memtotal_gb,
    memused_gb = memused_gb,
    memavailable_gb = memavailable_gb,
    should_stop = should_stop
  ))
}

cat("=== R SESSION INFO ===\n")
print(sessionInfo())
cat("\n")

# Check initial memory usage
cat("=== INITIAL MEMORY CHECK ===\n")
mem_status <- check_memory_usage(threshold = 0.90, stop_on_exceed = TRUE)
if (mem_status$should_stop) {
  sink()  # Close log file
  stop("Script stopped due to high memory usage (>90%)")
}
cat("\n")

cat("=== WORKING DIRECTORY ===\n")

# Define working directory where the files are located
work_dir <- "/mnt/nfs/CX000008_DS1/projects/jaeyeon/fastq_file_Dyslexia_r1/zanalysis_bogdan/samples.merged_AG_Harmony_res0.1_anno_091625/merged_AG_Harmony_res0.1_anno_091625.to.compare.results_cellchat"
setwd(work_dir)
cat("Set to:", getwd(), "\n\n")

cat("=== FILES IN WORKING DIRECTORY ===\n")
print(list.files())
cat("\n")

cat("=== LOADING CELLCHAT OBJECTS ===\n")
ad_ctl = readRDS("cellchat_merged_AG_Harmony_res0.1_anno_091625.AD_CTL_at_the_end.rds")
cat("✓ Loaded: AD_CTL\n")

# glimpse(ad_ctl)
unique(ad_ctl@meta$orig.ident)
unique(ad_ctl@meta$group)
unique(ad_ctl@meta$celltype)

ad_dys = readRDS("cellchat_merged_AG_Harmony_res0.1_anno_091625.AD_Dyslexia_at_the_end.rds")
cat("✓ Loaded: AD_Dyslexia\n")

# Check memory after loading large objects
cat("\n=== MEMORY CHECK AFTER LOADING DATA ===\n")
mem_status <- check_memory_usage(threshold = 0.90, stop_on_exceed = TRUE)
if (mem_status$should_stop) {
  sink()  # Close log file
  stop("Script stopped due to high memory usage (>90%)")
}
cat("\n")

# glimpse(ad_dys)
unique(ad_dys@meta$orig.ident)
unique(ad_dys@meta$group)
unique(ad_dys@meta$celltype)

# ctl_ctl = readRDS("cellchat_merged_AG_Harmony_res0.1_anno_091625.Ctr_Ctr_at_the_end.rds")

# glimpse(ctl_ctl)
# unique(ctl_ctl@meta$orig.ident)
# unique(ctl_ctl@meta$group)
# unique(ctl_ctl@meta$celltype)

# ctl_ctl@data.raw
# ctl_ctl@data 
# ctl_ctl@data.signaling 
# ctl_ctl@data.scale 
# ctl_ctl@net 
# ctl_ctl@netP 
# ctl_ctl@meta
# ctl_ctl@idents
# ctl_ctl@DB 
# ctl_ctl@LR 
# ctl_ctl@var.features 

# ctl_ctl@data.raw       : empty
# ctl_ctl@data
# ctl_ctl@data.signaling
# ctl_ctl@data.scale     : empty

# unique(ctl_ctl@meta$orig.ident)
# unique(ctl_ctl@meta$group)
# unique(ctl_ctl@meta$celltype)

# ad_ctl
# ad_dys
# ctl_ctl

# ad_ctl = updateCellChat(ad_ctl)
# ad_dys = updateCellChat(ad_dys)
# ctl_ctl = updateCellChat(ctl_ctl)

# working with :
ad_ctl = updateCellChat(ad_ctl)
ad_dys = updateCellChat(ad_dys)

# Create list of objects
object.list <- list(adctl = ad_ctl, addys = ad_dys)
print(object.list)

cat("\n=== MERGING CELLCHAT OBJECTS ===\n")
# Check memory before merging (memory-intensive operation)
cat("=== MEMORY CHECK BEFORE MERGING ===\n")
mem_status <- check_memory_usage(threshold = 0.90, stop_on_exceed = TRUE)
if (mem_status$should_stop) {
  sink()  # Close log file
  stop("Script stopped due to high memory usage (>90%)")
}
cat("\n")

# Merge CellChat objects
# cellchat <- mergeCellChat(object.list, add.names = names(object.list))
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
cat("✓ Merged CellChat objects\n\n")

cellchat 

glimpse(cellchat)

glimpse(cellchat@idents)

unique(cellchat@meta$datasets)

# cellchat@net
# cellchat@netP
# slotNames(cellchat)
str(cellchat, max.level = 2)
head(cellchat@meta, 2)
unique(cellchat@meta$datasets)
table(cellchat@meta$datasets)

names(cellchat@net)  # Should show: 'adctl''addys'
names(cellchat@net$adctl)   # Shows: "count", "weight", "sum", etc.
names(cellchat@net$addys)   # Same structure

pheatmap::pheatmap(cellchat@net$adctl$count, 
                   main = "adctl: Number of Interactions")


png("adctl_interaction_count_heatmap.png", width = 8, height = 7, units = "in", res = 300)
pheatmap::pheatmap(cellchat@net$adctl$count, 
                   main = "adctl: Number of Interactions")
dev.off()

pheatmap::pheatmap(cellchat@net$addys$count, 
                   main = "addys: Number of Interactions")

png("addys_interaction_count_heatmap.png", width = 8, height = 7, units = "in", res = 300)
pheatmap::pheatmap(cellchat@net$addys$count, 
                   main = "addys: Number of Interactions")
dev.off()

pheatmap::pheatmap(cellchat@net$adctl$weight, 
                   main = "adctl: Interaction Strength")

png("adctl_interaction_weight_heatmap.png", width = 8, height = 7, units = "in", res = 300)
pheatmap::pheatmap(cellchat@net$adctl$weight, 
                   main = "adctl: Interaction Strength")
dev.off()

pheatmap::pheatmap(cellchat@net$addys$weight, 
                   main = "addys: Interaction Strength")

png("addys_interaction_weight_heatmap.png", width = 8, height = 7, units = "in", res = 300)
pheatmap::pheatmap(cellchat@net$addys$weight, 
                   main = "addys: Interaction Strength")
dev.off()

head(cellchat@meta, 2)
unique(cellchat@meta$orig.ident)
unique(cellchat@meta$datasets)

# cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("adctl", "addys"))

# Number of significant ligand-receptor pairs detected
num_LR_pairs <- nrow(ad_ctl@LR$LRsig)
print(paste("Number of significant L-R pairs:", num_LR_pairs))

# Total number of interactions (count matrix)
total_interactions <- sum(ad_ctl@net$count)
print(paste("Total number of interactions:", total_interactions))

# Or if you want unique interactions (excluding self-interactions)
total_interactions_no_diag <- sum(ad_ctl@net$count) - sum(diag(ad_ctl@net$count))
print(paste("Total interactions (no self):", total_interactions_no_diag))

# Number of significant ligand-receptor pairs detected
num_LR_pairs <- nrow(ad_dys@LR$LRsig)
print(paste("Number of significant L-R pairs:", num_LR_pairs))

# Total number of interactions (count matrix)
total_interactions <- sum(ad_dys@net$count)
print(paste("Total number of interactions:", total_interactions))

# Or if you want unique interactions (excluding self-interactions)
total_interactions_no_diag <- sum(ad_dys@net$count) - sum(diag(ad_dys@net$count))
print(paste("Total interactions (no self):", total_interactions_no_diag))

ptm = Sys.time()
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

p <- gg1 + gg2
ggsave("compareInteractions_count_weight.png", plot = p, width = 10, height = 5, units = "in", dpi = 300)

table(cellchat@meta$datasets)

# par(mfrow = c(1,2), xpd=TRUE)
# netVisual_diffInteraction(cellchat, weight.scale = T)
# netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

# cellchat
# glimpse(cellchat@idents)

# cell_counts <- table(cellchat@meta$datasets)  # Assuming 'datasets' marks your conditions
# print(cell_counts)

# netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count.merged")
# netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight.merged")

library(ComplexHeatmap)

options(repr.plot.width = 6, repr.plot.height = 6)
gg1 <- netVisual_heatmap(cellchat)

# Do heatmap based on a merged object
draw(gg1)  # Display first heatmap

gg1

png("netVisual_heatmap_count.png", width = 6, height = 6, units = "in", res = 300)
ComplexHeatmap::draw(gg1)
dev.off()

gg2 <- netVisual_heatmap(cellchat, measure = "weight")
# Do heatmap based on a merged object
# Use draw() instead of direct printing
# Display separately instead of combining
draw(gg2)  # Display second heatmap

png("netVisual_heatmap_weight.png", width = 6, height = 6, units = "in", res = 300)
ComplexHeatmap::draw(gg2)
dev.off()

object.list

# Set figure size (width, height in inches)
options(repr.plot.width = 16, repr.plot.height = 8)
options(repr.plot.res = 150)  # Resolution (DPI)

weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))

par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, 
                   weight.scale = T, 
                   label.edge= F, 
                   edge.weight.max = weight.max[2], 
                   edge.width.max = 12, 
                   title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

png("netVisual_circle_count_comparison.png", width = 16, height = 8, units = "in", res = 300)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, 
                   weight.scale = T, 
                   label.edge= F, 
                   edge.weight.max = weight.max[2], 
                   edge.width.max = 12, 
                   title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
dev.off()

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count.merged", label.edge = T)

png("netVisual_diffInteraction_count_merged.png", width = 14, height = 7, units = "in", res = 300)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count.merged", label.edge = T)
dev.off()

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight.merged", label.edge = T)

png("netVisual_diffInteraction_weight_merged.png", width = 14, height = 7, units = "in", res = 300)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight.merged", label.edge = T)
dev.off()

# analysis based on COUNT
options(repr.plot.width = 16, repr.plot.height = 6)

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], 
                                               title = names(object.list)[i], 
                                               weight.MinMax = weight.MinMax)
}

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)

p <- patchwork::wrap_plots(plots = gg)
ggsave("netAnalysis_signalingRole_scatter_count.png", plot = p, width = 16, height = 6, units = "in", dpi = 300)

# analysis based on WEIGHT (STRENGTH)
options(repr.plot.width = 16, repr.plot.height = 6)

# For WEIGHT-based analysis (communication strength):
# Modify the plot after creation to increase dot size

num.link <- sapply(object.list, function(x) {
  rowSums(x@net$weight) + colSums(x@net$weight) - diag(x@net$weight)
})
num.link

weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets

gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], 
                                               title = names(object.list)[i], 
                                               weight.MinMax = weight.MinMax)
  
  # Increase point size manually
  gg[[i]] <- gg[[i]] + 
  ggplot2::scale_size_continuous(range = c(3, 10))  # Adjust these values
}

patchwork::wrap_plots(plots = gg)

p <- patchwork::wrap_plots(plots = gg)
ggsave("netAnalysis_signalingRole_scatter_weight.png", plot = p, width = 16, height = 6, units = "in", dpi = 300)



options(repr.plot.width = 16, repr.plot.height = 6)

gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "L6_IT_Exc")
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "L6_Exc")
patchwork::wrap_plots(plots = list(gg1,gg2))

p <- patchwork::wrap_plots(plots = list(gg1,gg2))
ggsave("netAnalysis_signalingChanges_scatter_L6.png", plot = p, width = 16, height = 6, units = "in", dpi = 300)

options(repr.plot.width = 16, repr.plot.height = 6)

gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "L56_Exc")
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "L4_Exc")
patchwork::wrap_plots(plots = list(gg1,gg2))

p <- patchwork::wrap_plots(plots = list(gg1,gg2))
ggsave("netAnalysis_signalingChanges_scatter_L56_L4.png", plot = p, width = 16, height = 6, units = "in", dpi = 300)


# conda activate cellbender

# Set Python path BEFORE loading any libraries that use reticulate
# Must be set before reticulate initializes (which happens on first import)
# python_path <- "~/.virtualenvs/r-reticulate/bin/python"
# python_path_expanded <- path.expand(python_path)

# if (file.exists(python_path_expanded)) {
#  Sys.setenv(RETICULATE_PYTHON = python_path_expanded)
#  cat("✓ Python environment set to:", python_path_expanded, "\n")
# } else {
#  cat("⚠️  Custom Python path not found:", python_path_expanded, "\n")
#  cat("   Will use system default Python\n")
# }

# Set output directory for FUNCTIONAL similarity
output_dir_functional <- "cellchat_functional_plots"
dir.create(output_dir_functional, showWarnings = FALSE)
cat("✓ Created output directory:", file.path(getwd(), output_dir_functional), "\n")

library(future)

# Set global plan at the beginning
# plan("sequential")

# Force multisession no matter what downstream code tries to do
# options(future.plan = "multisession")
# plan(multisession, workers = 4)

# library(future)
# library(BiocParallel)
# plan(sequential)

cat("Identify signaling groups based on their FUNCTIONAL similarity")

# Check memory before similarity computation (very memory-intensive)
cat("\n=== MEMORY CHECK BEFORE FUNCTIONAL SIMILARITY ===\n")
mem_status <- check_memory_usage(threshold = 0.90, stop_on_exceed = TRUE)
if (mem_status$should_stop) {
  sink()  # Close log file
  stop("Script stopped due to high memory usage (>90%)")
}
cat("\n")

# Then run all steps
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")

cellchat <- netEmbedding(cellchat, 
                         type = "functional", 
                         umap.method = "uwot") 

cellchat <- netClustering(cellchat, type = "functional", do.parallel = FALSE)
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
netVisual_embeddingPairwiseZoomIn(cellchat, type = "functional", nCol = 2)

# Standard embedding plot
functional_embedding_file <- file.path(output_dir_functional, "cellchat_functional_embedding.png")
png(functional_embedding_file,
    width = 12, height = 10, units = "in", res = 300)
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
dev.off()
full_path_embedding <- file.path(getwd(), functional_embedding_file)
cat("✓ Saved:", functional_embedding_file, "\n")
cat("   Full path:", full_path_embedding, "\n")
if (file.exists(full_path_embedding)) {
  cat("   ✓ File verified to exist\n")
} else {
  cat("   ⚠️  WARNING: File not found at expected location!\n")
}

# Zoom-in version (faceted)
functional_zoomin_file <- file.path(output_dir_functional, "cellchat_functional_embedding_zoomin.png")
png(functional_zoomin_file,
    width = 16, height = 8, units = "in", res = 300)
netVisual_embeddingPairwiseZoomIn(cellchat, type = "functional", nCol = 2)
dev.off()
full_path_zoomin <- file.path(getwd(), functional_zoomin_file)
cat("✓ Saved:", functional_zoomin_file, "\n")
cat("   Full path:", full_path_zoomin, "\n")
if (file.exists(full_path_zoomin)) {
  cat("   ✓ File verified to exist\n")
} else {
  cat("   ⚠️  WARNING: File not found at expected location!\n")
}

# glimpse(cellchat@netP)
# cellchat@netP$similarity           # Pathway similarity matrices
cellchat@netP$similarity$functional  # Functional similarity specifically

# You can use that object exactly like a 2-condition signaling similarity / distance result:

# $matrix$"1-2" is a 115×115 pairwise similarity (or distance) matrix between pathways 
# from condition 1 (adctl) and condition 2 (addys).
# $dr$"1-2" is a 2D embedding (UMAP1/UMAP2) placing those pathways in a plane so you can visualize closeness.

# The row/col names like NRXN--adctl and CSPG4--addys are the “nodes” in the embedding.

# Functional similarity matrix (condition 1 vs 2)
M  <- cellchat@netP$similarity$functional$matrix[["1-2"]]   # 115 x 115

# 2D embedding coordinates for that same comparison
DR <- cellchat@netP$similarity$functional$dr[["1-2"]]       # 115 x 2

dim(M); 
dim(DR)
head(rownames(M)); 
head(rownames(DR))

pheatmap(M, 
         show_rownames = TRUE, 
         show_colnames = TRUE, 
         fontsize_row = 5,
         fontsize_col = 5,
         main = "Pathway Functional Similarity")

png(file.path(output_dir_functional, "functional_similarity_heatmap.png"),
    width = 14, height = 14, units = "in", res = 300)
pheatmap(M, 
         show_rownames = TRUE, 
         show_colnames = TRUE, 
         fontsize_row = 5,
         fontsize_col = 5,
         main = "Pathway Functional Similarity")
dev.off()




df <- as.data.frame(DR)
colnames(df) <- c("UMAP1", "UMAP2")
df$pathway <- rownames(DR)
df$condition <- sub(".*--", "", df$pathway)
df$pathway_clean <- sub("--.*", "", df$pathway)

# Calculate shifts
ctl <- df %>% filter(condition == "adctl") %>% 
  select(pathway = pathway_clean, UMAP1, UMAP2)
dys <- df %>% filter(condition == "addys") %>% 
  select(pathway = pathway_clean, UMAP1, UMAP2)

# Use full_join to keep condition-specific pathways
merged <- full_join(ctl, dys, by = "pathway", suffix = c("_ctl", "_dys")) %>%
  mutate(
    # Calculate shift only for pathways present in both conditions
    shift = ifelse(
      !is.na(UMAP1_ctl) & !is.na(UMAP1_dys),
      sqrt((UMAP1_ctl - UMAP1_dys)^2 + (UMAP2_ctl - UMAP2_dys)^2),
      NA_real_
    )
  ) %>%
  arrange(desc(shift))


# Now merged has: pathway, UMAP1_ctl, UMAP2_ctl, UMAP1_dys, UMAP2_dys, shift
# Pathways specific to one condition will have NA coordinates for the other

topN <- min(500, nrow(merged))  # Safety check: don't exceed available pathways
to_label <- merged$pathway[1:topN]

p1 <- ggplot(df, aes(UMAP1, UMAP2, color = condition)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_text_repel(
    data = df %>% filter(pathway_clean %in% to_label),
    aes(label = pathway_clean),
    size = 3.5, max.overlaps = Inf, box.padding = 0.5
  ) +
  scale_color_manual(values = c("adctl" = "#E41A1C", "addys" = "#377EB8")) +
  theme_classic(base_size = 13) +
  labs(title = sprintf("Functional Similarity - Top %d Shifted Pathways", topN),
       color = "Condition") +
  coord_fixed()

ggsave(file.path(output_dir_functional, "functional_umap_labeled.png"), 
       p1, width = 12, height = 10, dpi = 300)

# Pathway shift arrows - FUNCTIONAL SIMILARITY
# Filter to only pathways present in both conditions (where shift is not NA)
merged_both <- merged %>% filter(!is.na(shift))

p2 <- ggplot() +
  geom_segment(
    data = merged_both,
    aes(x = UMAP1_ctl, y = UMAP2_ctl, xend = UMAP1_dys, yend = UMAP2_dys),
    arrow = arrow(length = unit(0.12, "inches")),
    alpha = 0.5, color = "grey60"
  ) +
  geom_point(data = df, aes(UMAP1, UMAP2, color = condition), 
             size = 2, alpha = 0.8) +
  geom_text_repel(
    data = merged_both[1:min(15, nrow(merged_both)), ],
    aes(x = UMAP1_dys, y = UMAP2_dys, label = pathway),
    size = 3.5, max.overlaps = Inf
  ) +
  scale_color_manual(values = c("adctl" = "#E41A1C", "addys" = "#377EB8")) +
  theme_classic(base_size = 13) +
  labs(title = "Functional Shift (adctl → addys)", color = "Condition") +
  coord_fixed()

ggsave(file.path(output_dir_functional, "functional_pathway_shifts.png"), 
       p2, width = 12, height = 10, dpi = 300)

# Shift rankings
write.csv(merged, 
          file.path(output_dir_functional, "functional_pathway_shifts_ranked.csv"), 
          row.names = FALSE)

# Similarity matrix
write.csv(M, 
          file.path(output_dir_functional, "functional_similarity_matrix.csv"), 
          row.names = TRUE)

# UMAP coordinates
write.csv(df, 
          file.path(output_dir_functional, "functional_umap_coordinates.csv"), 
          row.names = FALSE)

p1

p2



cat("Identify signaling groups based on their STRUCTURAL similarity")

# Set output directory for STRUCTURAL similarity
output_dir_structural <- "cellchat_structural_plots"
dir.create(output_dir_structural, showWarnings = FALSE)
cat("✓ Created output directory:", file.path(getwd(), output_dir_structural), "\n")

# Check memory before structural similarity computation (very memory-intensive)
cat("\n=== MEMORY CHECK BEFORE STRUCTURAL SIMILARITY ===\n")
mem_status <- check_memory_usage(threshold = 0.90, stop_on_exceed = TRUE)
if (mem_status$should_stop) {
  sink()  # Close log file
  stop("Script stopped due to high memory usage (>90%)")
}
cat("\n")

cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")

cellchat <- netEmbedding(cellchat, 
                            type = "structural", 
                            umap.method = "uwot")

cellchat <- netClustering(cellchat, type = "structural", do.parallel = FALSE)

# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5)
netVisual_embeddingPairwiseZoomIn(cellchat, type = "structural", nCol = 2)

# Standard embedding plot
structural_embedding_file <- file.path(output_dir_structural, "structural_cellchat_embedding.png")
png(structural_embedding_file,
    width = 12, height = 10, units = "in", res = 300)
netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5)
dev.off()
full_path_structural <- file.path(getwd(), structural_embedding_file)
cat("✓ Saved:", structural_embedding_file, "\n")
cat("   Full path:", full_path_structural, "\n")
if (file.exists(full_path_structural)) {
  cat("   ✓ File verified to exist\n")
} else {
  cat("   ⚠️  WARNING: File not found at expected location!\n")
}

netVisual_embeddingPairwiseZoomIn(cellchat, type = "structural", nCol = 2)

# Zoom-in version (faceted)
structural_zoomin_file <- file.path(output_dir_structural, "structural_cellchat_embedding_zoomin.png")
png(structural_zoomin_file,
    width = 16, height = 8, units = "in", res = 300)
netVisual_embeddingPairwiseZoomIn(cellchat, type = "structural", nCol = 2)
dev.off()
full_path_structural_zoomin <- file.path(getwd(), structural_zoomin_file)
cat("✓ Saved:", structural_zoomin_file, "\n")
cat("   Full path:", full_path_structural_zoomin, "\n")
if (file.exists(full_path_structural_zoomin)) {
  cat("   ✓ File verified to exist\n")
} else {
  cat("   ⚠️  WARNING: File not found at expected location!\n")
}

# glimpse(cellchat@netP)
# cellchat@netP$similarity        # Pathway similarity matrices
cellchat@netP$similarity$structural  # Functional similarity specifically

M <- cellchat@netP$similarity$structural$matrix[["1-2"]]

pheatmap(M, 
         show_rownames = TRUE, 
         show_colnames = TRUE, 
         fontsize_row = 5,
         fontsize_col = 5,
         main = "Structural Similarity")

png(file.path(output_dir_structural, "similarity_heatmap.png"),
    width = 14, height = 14, units = "in", res = 300)
pheatmap(M, 
         show_rownames = TRUE, 
         show_colnames = TRUE, 
         fontsize_row = 5,
         fontsize_col = 5,
         main = "Structural Similarity")
dev.off()

DR <- cellchat@netP$similarity$structural$dr[["1-2"]]

df <- as.data.frame(DR)
colnames(df) <- c("UMAP1", "UMAP2")
df$pathway <- rownames(DR)
df$condition <- sub(".*--", "", df$pathway)
df$pathway_clean <- sub("--.*", "", df$pathway)

# Calculate shifts
ctl <- df %>% filter(condition == "adctl") %>% 
  select(pathway = pathway_clean, UMAP1, UMAP2)
dys <- df %>% filter(condition == "addys") %>% 
  select(pathway = pathway_clean, UMAP1, UMAP2)

# Use full_join to keep condition-specific pathways
merged <- full_join(ctl, dys, by = "pathway", suffix = c("_ctl", "_dys")) %>%
  mutate(
    # Calculate shift only for pathways present in both conditions
    shift = ifelse(
      !is.na(UMAP1_ctl) & !is.na(UMAP1_dys),
      sqrt((UMAP1_ctl - UMAP1_dys)^2 + (UMAP2_ctl - UMAP2_dys)^2),
      NA_real_
    )
  ) %>%
  arrange(desc(shift))


# UMAP with labeled pathways
# Pathways specific to one condition will have NA coordinates for the other
topN <- min(500, nrow(merged))  # Safety check: don't exceed available pathways

to_label <- merged$pathway[1:topN]

p1 <- ggplot(df, aes(UMAP1, UMAP2, color = condition)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_text_repel(
    data = df %>% filter(pathway_clean %in% to_label),
    aes(label = pathway_clean),
    size = 3.5, max.overlaps = Inf, box.padding = 0.5
  ) +
  scale_color_manual(values = c("adctl" = "#E41A1C", "addys" = "#377EB8")) +
  theme_classic(base_size = 13) +
  labs(title = sprintf("Structural Similarity - Top %d Shifted Pathways", topN),
       color = "Condition") +
  coord_fixed()

ggsave(file.path(output_dir_structural, "umap_labeled.png"), 
       p1, width = 12, height = 10, dpi = 300)

# Pathway shift arrows
# Filter to only pathways present in both conditions (where shift is not NA)
merged_both <- merged %>% filter(!is.na(shift))

p2 <- ggplot() +
  geom_segment(
    data = merged_both,
    aes(x = UMAP1_ctl, y = UMAP2_ctl, xend = UMAP1_dys, yend = UMAP2_dys),
    arrow = arrow(length = unit(0.12, "inches")),
    alpha = 0.5, color = "grey60"
  ) +
  geom_point(data = df, aes(UMAP1, UMAP2, color = condition), 
             size = 2, alpha = 0.8) +
  geom_text_repel(
    data = merged_both[1:min(15, nrow(merged_both)), ],
    aes(x = UMAP1_dys, y = UMAP2_dys, label = pathway),
    size = 3.5, max.overlaps = Inf
  ) +
  scale_color_manual(values = c("adctl" = "#E41A1C", "addys" = "#377EB8")) +
  theme_classic(base_size = 13) +
  labs(title = "Structural Shift (adctl → addys)", color = "Condition") +
  coord_fixed()

ggsave(file.path(output_dir_structural, "pathway_shifts.png"), 
       p2, width = 12, height = 10, dpi = 300)

p1

p2

# Shift rankings
write.csv(merged, 
          file.path(output_dir_structural, "structural_pathway_shifts_ranked.csv"), 
          row.names = FALSE)

# Similarity matrix
write.csv(M, 
          file.path(output_dir_structural, "structural_similarity_matrix.csv"), 
          row.names = TRUE)

# UMAP coordinates
write.csv(df, 
          file.path(output_dir_structural, "structural_umap_coordinates.csv"), 
          row.names = FALSE)



print("RANK SIMILARITY")

# rankSimilarity(cellchat, type = "functional")
# Compute the distance of signaling networks between datasets 1 2

tryCatch({
  ranked_pathways_functional <- rankSimilarity(cellchat, type = "functional", slot.name = "netP")
  
  options(repr.plot.width = 12, repr.plot.height = 10)
  options(repr.plot.res = 200)
  
  print(ranked_pathways_functional)
  
  png(file.path(output_dir_functional, "ranked_functional_similarity.png"),
      width = 10, height = 8, units = "in", res = 300)
  print(ranked_pathways_functional)
  dev.off()
  
  cat("✓ Saved: ranked_functional_similarity.png\n")
  
}, error = function(e) {
  cat("❌ Error in functional similarity ranking:\n")
  cat("   ", as.character(e), "\n")
  cat("   Skipping ranked_functional_similarity plot...\n")
  # Ensure device is closed if error occurred during plotting
  if (dev.cur() > 1) dev.off()
  if (dev.cur() > 1) dev.off()
})



print("Ranking STRUCTURAL SIMILARITY")

tryCatch({
  ranked_pathways_structural <- rankSimilarity(cellchat, type = "structural", slot.name = "netP")
  # Compute the distance of signaling networks between datasets 1 2
  
  print(ranked_pathways_structural)
  
  png(file.path(output_dir_structural, "ranked_structural_similarity.png"),
      width = 10, height = 8, units = "in", res = 300)
  print(ranked_pathways_structural)
  dev.off()
  
  cat("✓ Saved: ranked_structural_similarity.png\n")
  
}, error = function(e) {
  cat("❌ Error in structural similarity ranking:\n")
  cat("   ", as.character(e), "\n")
  cat("   Skipping ranked_structural_similarity plot...\n")
  # Ensure device is closed if error occurred during plotting
  if (dev.cur() > 1) dev.off()
  # Ensure device is closed if error occurred during plotting
  if (dev.cur() > 1) dev.off()
})



print("Identify altered signaling with distinct interaction strength")

signaling.type = c("Secreted Signaling", "ECM-Receptor", "Cell–Cell Contact", "Non-protein Signaling")

# ═══════════════════════════════════════════════════════════
# DATA VALIDATION BEFORE RANKING
# ═══════════════════════════════════════════════════════════

# Check for problematic values in the data
cat("Checking data quality...\n")

# Check net slot for NaN/Inf values
if (!is.null(cellchat@net)) {
  for (dataset in names(cellchat@net)) {
    net_data <- cellchat@net[[dataset]]
    if (is.matrix(net_data) || is.array(net_data)) {
      na_count <- sum(is.na(net_data))
      inf_count <- sum(is.infinite(net_data))
      if (na_count > 0) cat(sprintf("⚠️  Found %d NA values in cellchat@net$%s\n", na_count, dataset))
      if (inf_count > 0) cat(sprintf("⚠️  Found %d Inf values in cellchat@net$%s\n", inf_count, dataset))
      
      # Replace problematic values
      if (na_count > 0 || inf_count > 0) {
        net_data[is.na(net_data) | is.infinite(net_data)] <- 0
        cellchat@net[[dataset]] <- net_data
        cat(sprintf("   → Replaced with 0 in cellchat@net$%s\n", dataset))
      }
    }
  }
}

# Check netP slot for problematic values (used by signaling role heatmaps)
if (!is.null(cellchat@netP)) {
  for (dataset in names(cellchat@netP)) {
    if (is.list(cellchat@netP[[dataset]])) {
      # Check prob matrix
      if (!is.null(cellchat@netP[[dataset]]$prob)) {
        prob_data <- cellchat@netP[[dataset]]$prob
        if (is.matrix(prob_data)) {
          na_count <- sum(is.na(prob_data))
          inf_count <- sum(is.infinite(prob_data))
          if (na_count > 0) cat(sprintf("⚠️  Found %d NA values in cellchat@netP$%s$prob\n", na_count, dataset))
          if (inf_count > 0) cat(sprintf("⚠️  Found %d Inf values in cellchat@netP$%s$prob\n", inf_count, dataset))
          
          if (na_count > 0 || inf_count > 0) {
            prob_data[is.na(prob_data) | is.infinite(prob_data)] <- 0
            cellchat@netP[[dataset]]$prob <- prob_data
            cat(sprintf("   → Replaced with 0 in cellchat@netP$%s$prob\n", dataset))
          }
        }
      }
      
      # Check centr (centrality) matrices
      if (!is.null(cellchat@netP[[dataset]]$centr)) {
        centr_data <- cellchat@netP[[dataset]]$centr
        if (is.list(centr_data)) {
          for (centr_name in names(centr_data)) {
            if (is.matrix(centr_data[[centr_name]])) {
              na_count <- sum(is.na(centr_data[[centr_name]]))
              inf_count <- sum(is.infinite(centr_data[[centr_name]]))
              if (na_count > 0 || inf_count > 0) {
                cat(sprintf("⚠️  Found problematic values in cellchat@netP$%s$centr$%s\n", dataset, centr_name))
                centr_data[[centr_name]][is.na(centr_data[[centr_name]]) | is.infinite(centr_data[[centr_name]])] <- 0
                cat(sprintf("   → Replaced with 0\n"))
              }
            }
          }
          cellchat@netP[[dataset]]$centr <- centr_data
        }
      }
    }
  }
}

# Also check object.list for problematic values
for (i in seq_along(object.list)) {
  obj_name <- names(object.list)[i]
  
  # Check net slot
  if (!is.null(object.list[[i]]@net)) {
    net_data <- object.list[[i]]@net
    if (is.list(net_data) && !is.null(net_data$prob)) {
      prob_data <- net_data$prob
      if (is.matrix(prob_data)) {
        na_count <- sum(is.na(prob_data))
        inf_count <- sum(is.infinite(prob_data))
        if (na_count > 0 || inf_count > 0) {
          cat(sprintf("⚠️  Found problematic values in object.list$%s@net$prob\n", obj_name))
          prob_data[is.na(prob_data) | is.infinite(prob_data)] <- 0
          object.list[[i]]@net$prob <- prob_data
          cat(sprintf("   → Replaced with 0\n"))
        }
      }
    }
  }
  
  # Check netP slot
  if (!is.null(object.list[[i]]@netP)) {
    if (is.list(object.list[[i]]@netP) && !is.null(object.list[[i]]@netP$prob)) {
      prob_data <- object.list[[i]]@netP$prob
      if (is.matrix(prob_data)) {
        na_count <- sum(is.na(prob_data))
        inf_count <- sum(is.infinite(prob_data))
        if (na_count > 0 || inf_count > 0) {
          cat(sprintf("⚠️  Found problematic values in object.list$%s@netP$prob\n", obj_name))
          prob_data[is.na(prob_data) | is.infinite(prob_data)] <- 0
          object.list[[i]]@netP$prob <- prob_data
          cat(sprintf("   → Replaced with 0\n"))
        }
      }
    }
  }
}

cat("✓ Data validation complete\n\n")

# ═══════════════════════════════════════════════════════════
# RANK NETWORK COMPARISONS WITH ERROR HANDLING
# ═══════════════════════════════════════════════════════════

# 1. Stacked bar plot with statistics
cat("Generating rankNet plots...\n")

tryCatch({
  gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", 
                 sources.use = NULL, targets.use = NULL, 
                 stacked = TRUE, do.stat = TRUE, 
                 signaling.type = signaling.type)
  
  print(gg1)
  
  ggsave("ranknet_weight_stacked.pdf",
         gg1,
         width = 10, 
         height = 8)
  
  cat("✓ Saved: ranknet_weight_stacked.pdf\n")
  
}, error = function(e) {
  cat("❌ Error in rankNet (stacked):\n")
  cat("   ", as.character(e), "\n")
  cat("   This often happens when data contains NA/NaN/Inf values.\n")
  cat("   Skipping stacked rankNet plot...\n\n")
})

# 2. Grouped bar plot with statistics
tryCatch({
  gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", 
                 sources.use = NULL, targets.use = NULL, 
                 stacked = FALSE, do.stat = TRUE, 
                 signaling.type = signaling.type)
  
  print(gg2)
  
  ggsave("ranknet_weight.pdf",
         gg2,
         width = 10, 
         height = 8)
  
  cat("✓ Saved: ranknet_weight.pdf\n")
  
}, error = function(e) {
  cat("❌ Error in rankNet (grouped):\n")
  cat("   ", as.character(e), "\n")
  cat("   Skipping grouped rankNet plot...\n\n")
})



# 3. Network-based ranking (stacked)
tryCatch({
  gg3 <- rankNet(cellchat, slot.name = "net", mode = "comparison", 
                 measure = "weight", stacked = TRUE, do.stat = TRUE, 
                 signaling.type = signaling.type)
  
  # Set plot size for Jupyter notebook display
  options(repr.plot.width = 16, repr.plot.height = 100)
  options(repr.plot.res = 200)  # Resolution (DPI)
  
  print(gg3)
  
  ggsave("ranknet_weight_net_stacked.pdf",
         gg3,
         width = 10, 
         height = 40)
  
  cat("✓ Saved: ranknet_weight_net_stacked.pdf\n")
  
}, error = function(e) {
  cat("❌ Error in rankNet (network stacked):\n")
  cat("   ", as.character(e), "\n")
  cat("   Skipping network stacked rankNet plot...\n\n")
})

# 4. Network-based ranking (grouped)
tryCatch({
  gg4 <- rankNet(cellchat, slot.name = "net", mode = "comparison", 
                 measure = "weight", stacked = FALSE, do.stat = TRUE, 
                 signaling.type = signaling.type)
  
  # Set plot size for Jupyter notebook display
  options(repr.plot.width = 16, repr.plot.height = 70)
  options(repr.plot.res = 200)  # Resolution (DPI)
  
  print(gg4)
  
  ggsave("ranknet_weight_net.pdf",
         gg4,
         width = 10, 
         height = 40)
  
  cat("✓ Saved: ranknet_weight_net.pdf\n")
  
}, error = function(e) {
  cat("❌ Error in rankNet (network grouped):\n")
  cat("   ", as.character(e), "\n")
  cat("   Skipping network grouped rankNet plot...\n\n")
})

cat("\n=== RANKNET PLOTS COMPLETE ===\n\n")

# ═══════════════════════════════════════════════════════════
# ADDITIONAL DATA CLEANING BEFORE HEATMAPS
# ═══════════════════════════════════════════════════════════
# COMMENTED OUT - Not needed since heatmaps are commented out
# cat("Cleaning data for signaling role heatmaps...\n")
# for (i in seq_along(object.list)) {
#   object.list[[i]] <- clean_cellchat_data(object.list[[i]])
#   cat(sprintf("✓ Cleaned data for object.list[[%d]] (%s)\n", i, names(object.list)[i]))
# }
# cat("✓ Data cleaning complete for heatmaps\n\n")

# ═══════════════════════════════════════════════════════════
# HELPER FUNCTION: Clean data in CellChat object
# ═══════════════════════════════════════════════════════════

clean_cellchat_data <- function(obj) {
  # Clean netP slot (used by signaling role heatmaps)
  if (!is.null(obj@netP)) {
    # Clean prob matrix
    if (!is.null(obj@netP$prob) && is.matrix(obj@netP$prob)) {
      obj@netP$prob[is.na(obj@netP$prob) | is.infinite(obj@netP$prob)] <- 0
    }
    
    # Clean centr (centrality) matrices
    if (!is.null(obj@netP$centr) && is.list(obj@netP$centr)) {
      for (centr_name in names(obj@netP$centr)) {
        if (is.matrix(obj@netP$centr[[centr_name]])) {
          obj@netP$centr[[centr_name]][is.na(obj@netP$centr[[centr_name]]) | 
                                       is.infinite(obj@netP$centr[[centr_name]])] <- 0
        }
      }
    }
    
    # Clean outdeg and indeg matrices
    if (!is.null(obj@netP$outdeg) && is.matrix(obj@netP$outdeg)) {
      obj@netP$outdeg[is.na(obj@netP$outdeg) | is.infinite(obj@netP$outdeg)] <- 0
    }
    if (!is.null(obj@netP$indeg) && is.matrix(obj@netP$indeg)) {
      obj@netP$indeg[is.na(obj@netP$indeg) | is.infinite(obj@netP$indeg)] <- 0
    }
  }
  
  # Clean net slot
  if (!is.null(obj@net)) {
    if (is.list(obj@net)) {
      for (net_name in names(obj@net)) {
        if (is.matrix(obj@net[[net_name]])) {
          obj@net[[net_name]][is.na(obj@net[[net_name]]) | 
                             is.infinite(obj@net[[net_name]])] <- 0
        }
      }
    }
  }
  
  return(obj)
}

# ═══════════════════════════════════════════════════════════
# HELPER FUNCTION: Safe heatmap with clustering fallback
# ═══════════════════════════════════════════════════════════

safe_signalingRole_heatmap <- function(object, pattern, signaling, title, 
                                       cluster.rows = TRUE, cluster.cols = TRUE,
                                       ...) {
  # Clean data before attempting heatmap
  object_clean <- clean_cellchat_data(object)
  
  # Try with clustering first
  result <- tryCatch({
    netAnalysis_signalingRole_heatmap(
      object = object_clean,
      pattern = pattern,
      signaling = signaling,
      title = title,
      cluster.rows = cluster.rows,
      cluster.cols = cluster.cols,
      ...
    )
  }, error = function(e) {
    # If clustering fails, try without clustering
    cat(sprintf("⚠️  Clustering failed for %s (%s), trying without clustering...\n", 
                title, pattern))
    cat("   Error:", as.character(e), "\n")
    
    # Fallback: try without clustering
    tryCatch({
      netAnalysis_signalingRole_heatmap(
        object = object_clean,
        pattern = pattern,
        signaling = signaling,
        title = paste0(title, " (no clustering)"),
        cluster.rows = FALSE,
        cluster.cols = FALSE,
        ...
      )
    }, error = function(e2) {
      cat(sprintf("❌ Complete failure for %s (%s): %s\n", 
                  title, pattern, as.character(e2)))
      # Return NULL instead of stopping to allow script to continue
      return(NULL)
    })
  })
  
  return(result)
}

i = 1 
object.list[[i]]
object.list[[i+1]]

i = 1
pathway.union <- union(object.list[[i]]@netP$pathways,
                       object.list[[i+1]]@netP$pathways)

# ═══════════════════════════════════════════════════════════
# SIGNALING ROLE HEATMAPS - COMMENTED OUT DUE TO CLUSTERING ERRORS
# ═══════════════════════════════════════════════════════════
# cat("\n=== GENERATING SIGNALING ROLE HEATMAPS ===\n")
# cat("Pathway union size:", length(pathway.union), "\n")

# OUTGOING PATTERN
# COMMENTED OUT - FAILS WITH: Error in hclust(get_dist(submat, distance), method = method): NA/NaN/Inf in foreign function call (arg 10)
# cat("\n--- Outgoing pattern ---\n")
# tryCatch({
#   ht1 = safe_signalingRole_heatmap(
#     object = object.list[[i]], 
#     pattern = "outgoing",
#     signaling = pathway.union,
#     title = names(object.list)[i],
#     width = 8, height = 18,
#     font.size = 8,
#     font.size.title = 10,
#     cluster.rows = TRUE,
#     cluster.cols = TRUE
#   )
#   
#   ht2 = safe_signalingRole_heatmap(
#     object = object.list[[i+1]], 
#     pattern = "outgoing",
#     signaling = pathway.union,
#     title = names(object.list)[i+1],
#     width = 8, height = 18,
#     font.size = 8,
#     font.size.title = 10,
#     cluster.rows = TRUE,
#     cluster.cols = TRUE
#   )
#   
#   # Check if heatmaps were created successfully
#   if (is.null(ht1)) {
#     cat("⚠️  Skipping ht1 (outgoing) - heatmap creation failed\n")
#   } else {
#     # Draw separately instead of combining
#     options(repr.plot.width = 20, repr.plot.height = 20)
#     draw(ht1)
#     
#     png("azd_ctl_heatmap_netAnalysis_signalingRole_heatmap_outgoing.png",
#         width = 10, height = 20, units = "in", res = 300)
#     draw(ht1)
#     dev.off()
#     cat("✓ Saved: azd_ctl_heatmap_netAnalysis_signalingRole_heatmap_outgoing.png\n")
#   }
#   
#   if (is.null(ht2)) {
#     cat("⚠️  Skipping ht2 (outgoing) - heatmap creation failed\n")
#   } else {
#     options(repr.plot.width = 20, repr.plot.height = 20)
#     draw(ht2)
#     
#     png("azd_dys_heatmap_netAnalysis_signalingRole_heatmap_outgoing.png",
#         width = 10, height = 20, units = "in", res = 300)
#     draw(ht2)
#     dev.off()
#     cat("✓ Saved: azd_dys_heatmap_netAnalysis_signalingRole_heatmap_outgoing.png\n")
#   }
#   
# }, error = function(e) {
#   cat("❌ Error generating outgoing heatmaps:\n")
#   cat("   ", as.character(e), "\n")
#   if (dev.cur() > 1) dev.off()
# })





# INCOMING PATTERN
# COMMENTED OUT - FAILS WITH CLUSTERING ERRORS
# cat("\n--- Incoming pattern ---\n")
# tryCatch({
#   ht1 = safe_signalingRole_heatmap(
#     object = object.list[[i]], 
#     pattern = "incoming", 
#     signaling = pathway.union, 
#     title = names(object.list)[i], 
#     width = 8, height = 18, 
#     color.heatmap = "GnBu",
#     font.size = 8,
#     font.size.title = 10,
#     cluster.rows = TRUE,
#     cluster.cols = TRUE
#   )
#   
#   ht2 = safe_signalingRole_heatmap(
#     object = object.list[[i+1]], 
#     pattern = "incoming", 
#     signaling = pathway.union, 
#     title = names(object.list)[i+1], 
#     width = 8, height = 18, 
#     color.heatmap = "GnBu",
#     font.size = 8,
#     font.size.title = 10,
#     cluster.rows = TRUE,
#     cluster.cols = TRUE
#   )
#   
#   # Check if heatmaps were created successfully
#   if (is.null(ht1)) {
#     cat("⚠️  Skipping ht1 (incoming) - heatmap creation failed\n")
#   } else {
#     options(repr.plot.width = 20, repr.plot.height = 20)
#     draw(ht1)
#     
#     png("azd_ctl_heatmap_netAnalysis_signalingRole_heatmap_incoming.png",
#         width = 10, height = 20, units = "in", res = 300)
#     draw(ht1)
#     dev.off()
#     cat("✓ Saved: azd_ctl_heatmap_netAnalysis_signalingRole_heatmap_incoming.png\n")
#   }
#   
#   if (is.null(ht2)) {
#     cat("⚠️  Skipping ht2 (incoming) - heatmap creation failed\n")
#   } else {
#     options(repr.plot.width = 20, repr.plot.height = 20)
#     draw(ht2)
#     
#     png("azd_dys_heatmap_netAnalysis_signalingRole_heatmap_incoming.png",
#         width = 10, height = 20, units = "in", res = 300)
#     draw(ht2)
#     dev.off()
#     cat("✓ Saved: azd_dys_heatmap_netAnalysis_signalingRole_heatmap_incoming.png\n")
#   }
#   
# }, error = function(e) {
#   cat("❌ Error generating incoming heatmaps:\n")
#   cat("   ", as.character(e), "\n")
#   if (dev.cur() > 1) dev.off()
# })

# ALL PATTERN
# COMMENTED OUT - FAILS WITH CLUSTERING ERRORS
# cat("\n--- All pattern ---\n")
# tryCatch({
#   ht1 = safe_signalingRole_heatmap(
#     object = object.list[[i]], 
#     pattern = "all", 
#     signaling = pathway.union, 
#     title = names(object.list)[i], 
#     width = 8, height = 18, 
#     color.heatmap = "OrRd",
#     font.size = 8,
#     font.size.title = 10,
#     cluster.rows = TRUE,
#     cluster.cols = TRUE
#   )
#   
#   ht2 = safe_signalingRole_heatmap(
#     object = object.list[[i+1]], 
#     pattern = "all", 
#     signaling = pathway.union, 
#     title = names(object.list)[i+1], 
#     width = 8, height = 18, 
#     color.heatmap = "OrRd",
#     font.size = 8,
#     font.size.title = 10,
#     cluster.rows = TRUE,
#     cluster.cols = TRUE
#   )
#   
#   # Check if heatmaps were created successfully
#   if (is.null(ht1)) {
#     cat("⚠️  Skipping ht1 (all) - heatmap creation failed\n")
#   } else {
#     draw(ht1)
#     
#     options(repr.plot.width = 20, repr.plot.height = 20)
#     draw(ht1)
#     
#     png("azd_ctl_heatmap_netAnalysis_signalingRole_heatmap_all.png",
#         width = 10, height = 20, units = "in", res = 300)
#     draw(ht1)       
#     dev.off()
#     cat("✓ Saved: azd_ctl_heatmap_netAnalysis_signalingRole_heatmap_all.png\n")
#   }
#   
#   if (is.null(ht2)) {
#     cat("⚠️  Skipping ht2 (all) - heatmap creation failed\n")
#   } else {
#     draw(ht2)
#     
#     options(repr.plot.width = 20, repr.plot.height = 20)
#     draw(ht2)
#     
#     png("azd_dys_heatmap_netAnalysis_signalingRole_heatmap_all.png",
#         width = 10, height = 20, units = "in", res = 300)
#     draw(ht2)       
#     dev.off()
#     cat("✓ Saved: azd_dys_heatmap_netAnalysis_signalingRole_heatmap_all.png\n")
#   }
#   
# }, error = function(e) {
#   cat("❌ Error generating 'all' pattern heatmaps:\n")
#   cat("   ", as.character(e), "\n")
#   if (dev.cur() > 1) dev.off()
# })

cat("\n=== SIGNALING ROLE HEATMAPS SKIPPED (COMMENTED OUT DUE TO CLUSTERING ERRORS) ===\n\n")

ptm = Sys.time()

# Set figure size for notebook
options(repr.plot.width = 28, repr.plot.height = 60)
options(repr.plot.res = 200)  # Resolution

netVisual_bubble(cellchat, comparison = c(1, 2), angle.x = 45)
#> Comparing communications on a merged object

# Save as PDF instead of PNG for large plots (more efficient, smaller file size)
pdf("netVisual_bubble_comparison.pdf",
    width = 40, height = 40)
netVisual_bubble(cellchat, comparison = c(1, 2), angle.x = 45)
dev.off()
cat("✓ Saved: netVisual_bubble_comparison.pdf (PDF format for large size)\n")

# Comparing communications on a merged object

# L23_IT_Exc 
# Oligo 
# Micro 
# L4_Exc 
# High_mt
# OPc
# L6_Exc
# MGE_Inh
# Astro
# CGE_Inh
# L6_IT_Exc
# L56_Exc
# mesenchymal
# Immune

ptm = Sys.time()

options(repr.plot.width = 40, repr.plot.height = 40)
options(repr.plot.res = 200)  # Resolution

netVisual_bubble(cellchat, 
                 sources.use = c("L23_IT_Exc", "L4_Exc", "L56_Exc", "L6_Exc", "L6_IT_Exc", "Immune"),
                 targets.use = c("L23_IT_Exc", "L4_Exc", "L56_Exc", "L6_Exc", "L6_IT_Exc", "Immune"),
                 comparison = c(1, 2),
                 angle.x = 45)

png("netVisual_bubble_comparison2excit.png",
    width = 40, height = 40, units = "in", res = 200)
netVisual_bubble(cellchat, 
                 sources.use = c("L23_IT_Exc", "L4_Exc", "L56_Exc", "L6_Exc", "L6_IT_Exc", "Immune"),
                 targets.use = c("L23_IT_Exc", "L4_Exc", "L56_Exc", "L6_Exc", "L6_IT_Exc", "Immune"),
                 comparison = c(1, 2),
                 angle.x = 45)
dev.off()

# Comparing communications on a merged object

# L23_IT_Exc 
# Oligo 
# Micro 
# L4_Exc 
# High_mt
# OPc
# L6_Exc
# MGE_Inh
# Astro
# CGE_Inh
# L6_IT_Exc
# L56_Exc
# mesenchymal
# Immune

ptm = Sys.time()

options(repr.plot.width = 20, repr.plot.height = 20)
options(repr.plot.res = 200)  # Resolution

netVisual_bubble(cellchat, 
                 sources.use = c("MGE_Inh", "CGE_Inh", "Immune"),
                 targets.use = c("MGE_Inh", "CGE_Inh", "Immune"),
                 comparison = c(1, 2),
                 angle.x = 45)

png("netVisual_bubble_comparison2inhib.png", 
    width = 40, height = 40, units = "in", res = 200)

netVisual_bubble(cellchat, 
                 sources.use = c("MGE_Inh", "CGE_Inh", "Immune"),
                 targets.use = c("MGE_Inh", "CGE_Inh", "Immune"),
                 comparison = c(1, 2),
                 angle.x = 45)

dev.off()



# Identify the up-regulated (that is, increased) L–R pairs in the second dataset compared with the first dataset

options(repr.plot.width = 30, repr.plot.height = 40)
options(repr.plot.res = 200)  # Resolution

netVisual_bubble(cellchat, 
                 # sources.use = 4, 
                 # targets.use = c(5:11), 
                 comparison = c(1, 2), 
                 max.dataset = 2, 
                 title.name = "Increased signaling in DYS", 
                 angle.x = 45,
                 remove.isolate = T) 

png("netVisual_bubble_comparison3increased_in_DYS.png", 
    width = 40, height = 40, units = "in", res = 200)

netVisual_bubble(cellchat, 
                 # sources.use = 4, 
                 # targets.use = c(5:11), 
                 comparison = c(1, 2), 
                 max.dataset = 2, 
                 title.name = "Increased signaling in DYS", 
                 angle.x = 45,
                 remove.isolate = T) 

dev.off()

# Identify the down-regulated (that is, decreased) L–R pairs in the second dataset compared with the first dataset

options(repr.plot.width = 30, repr.plot.height = 40)
options(repr.plot.res = 200)  # Resolution

netVisual_bubble(cellchat, 
                 # sources.use = 4, 
                 # targets.use = c(5:11), 
                 comparison = c(1, 2), 
                 max.dataset = 1, 
                 title.name = "Decreased signaling in DYS", 
                 angle.x = 45,
                 remove.isolate = T) 

png("netVisual_bubble_comparison3decreased_in_DYS.png", 
    width = 40, height = 40, units = "in", res = 200)

netVisual_bubble(cellchat, 
                 # sources.use = 4, 
                 # targets.use = c(5:11), 
                 comparison = c(1, 2), 
                 max.dataset = 1, 
                 title.name = "Decreased signaling in DYS", 
                 angle.x = 45,
                 remove.isolate = T)

dev.off()

# (cellchat@idents$joint)

print("Identify dysfunctional signaling by using differential expression analysis.")

library(presto)

# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "addys"
# define a char name used for storing the results of differential expression analysis
features.name = paste0(pos.dataset, ".merged")

# features.name: a char name used for storing the over-expressed
#                signaling genes in `object@var.features[[features.name]]`

cellchat <- identifyOverExpressedGenes(cellchat, 
                                       group.dataset = "datasets", 
                                       pos.dataset = pos.dataset, 
                                       features.name = features.name, 
                                       only.pos = FALSE, 
                                       thresh.pc = 0.1, 
                                       thresh.fc = 0.05,
                                       thresh.p = 0.05) 

# glimpse(cellchat)
# glimpse(cellchat@net)
glimpse(cellchat@netP)


# glimpse(cellchat@netP$adctl$centr)
# glimpse(cellchat@netP$addys$centr)

unique((cellchat@var.features$adctl$features))
length((cellchat@var.features$adctl$features))

unique((cellchat@var.features$addys$features))
length((cellchat@var.features$addys$features))

glimpse(cellchat@var.features)

(cellchat@var.features$addys.merged)
length(cellchat@var.features$addys.merged)

glimpse(cellchat@var.features$adctl$features.info)
head(cellchat@var.features$adctl$features.info)
dim(cellchat@var.features$adctl$features.info)

glimpse(cellchat@var.features$addys$features.info)
head(cellchat@var.features$addys$features.info)

glimpse(cellchat@var.features$addys.merged)
cellchat@var.features$addys.merged
length(cellchat@var.features$addys.merged)

glimpse(cellchat@var.features$addys.merged.info)
head(cellchat@var.features$addys.merged.info, 5)
length(cellchat@var.features$addys.merged.info)

features.name

# glimpse(cellchat)

#  ..@ var.features  :List of 4
#  .. ..$ adctl            :List of 2
#  .. ..$ addys            :List of 2
#  .. ..$ addys.merged     : chr [1:1185] "RELN" "LRP2" "LAMA2" "CNTN1" ...
# .. ..$ addys.merged.info:'data.frame':	1185 obs. of  8 variables:





# Use the joint cell labels from the merged CellChat object

# map the results of differential expression analysis onto the inferred cell-cell communications to 
# easily manage/subset the ligand-receptor pairs of interest

net <- netMappingDEG(cellchat, features.name = features.name)

head(net, 4)
tail(net, 4)

colnames(net)
dim(net)

unique(cellchat@meta$datasets)

unique(net$ligand)
length(unique(net$ligand))

unique(net$receptor)
length(unique(net$receptor))

unique(net$interaction_name)
length(unique(net$interaction_name))

unique(net$interaction_name_2)
length(unique(net$interaction_name_2))

unique(net$pathway_name)
length(unique(net$pathway_name))

# to loop over 62 pathways ?

# Save to CSV

write.csv(
  data.frame(net),
  "results.netMappingDEG.ligand.receptor.net.csv",
  row.names = FALSE
)

# Extracts ligand-receptor pairs where BOTH ligand and receptor are upregulated (positive logFC > 0.05) in the addys (AD-Dyslexia) condition
# Higher expression in addys compared to adctl
net.up <- subsetCommunication(cellchat, net = net, datasets = "addys", ligand.logFC = 0.05, receptor.logFC = 0.05)
dim(net.up)

# Extracts ligand-receptor pairs where BOTH ligand and receptor are downregulated in addys (AD-Dyslexia) condition
# i.e., STRONGER in adctl (datasets = "adctl") - this is the correct logic for CellChat's DEG direction handling
# Note: pos.dataset was set to "addys" earlier, so negative logFC means lower in addys vs adctl
net.down <- subsetCommunication(cellchat, net = net, datasets = "adctl", ligand.logFC = -0.05, receptor.logFC = -0.05)
dim(net.down)



#     subsetCommunication(
#       object = NULL,
#       net = NULL,
#       slot.name = "net",
#       sources.use = NULL,
#       targets.use = NULL,
#       signaling = NULL,
#       pairLR.use = NULL,
#      thresh = 0.05,
#      datasets = NULL,
#       ligand.pvalues = NULL,
#       ligand.logFC = NULL,
#       ligand.pct.1 = NULL,
#       ligand.pct.2 = NULL,
#       receptor.pvalues = NULL,
#       receptor.logFC = NULL,
#      receptor.pct.1 = NULL,
#       receptor.pct.2 = NULL\

cat('

slot.name: the slot name of object: slot.name = "net" when extracting
          the inferred communications at the level of
          ligands/receptors; slot.name = "netP" when extracting the
          inferred communications at the level of signaling pathways

sources.use: a vector giving the index or the name of source cell
          groups

targets.use: a vector giving the index or the name of target cell
          groups.

signaling: a character vector giving the name of signaling pathways of
          interest

pairLR.use: a data frame consisting of one column named either
          "interaction_name" or "pathway_name", defining the
          interactions of interest

  thresh: threshold of the p-value for determining significant
          interaction

datasets: select the inferred cell-cell communications from a
          particular `datasets` when inputing a data frame `net`

ligand.pvalues, ligand.logFC, ligand.pct.1, ligand.pct.2: set threshold
          for ligand genes

          ligand.pvalues: threshold for pvalues in the differential
          expression gene analysis (DEG)

          ligand.logFC: threshold for logFoldChange in the DEG
          analysis; When ligand.logFC > 0, keep upgulated genes;
          otherwise, kepp downregulated genes

          ligand.logFC: threshold for logFoldChange in the DEG
          analysis; When ligand.logFC > 0, keep upgulated genes;
          otherwise, kepp downregulated genes

          ligand.pct.1: threshold for the percent of expressed genes in
          the defined >positive< cell group. keep genes with percent
          greater than ligand.pct.1

          ligand.pct.2: threshold for the percent of expressed genes in
          the cells except for the defined >positive< cell group

receptor.pvalues, receptor.logFC, receptor.pct.1, receptor.pct.2: set
          threshold for receptor genes

Value:

     If input object is created from a single dataset, a data frame of
     the inferred cell-cell communications of interest, consisting of
     source, target, interaction_name, pathway_name, prob and other
     information

     If input object is a merged object from multiple datasets, it will
     return a list and each element is a data frame for one dataset

')

head(net.up, 2)
dim(net.up)
tail(net.up, 2)



head(net.down, 2)
dim(net.down)
tail(net.down, 2)



# Save as CSV
write.csv(net.up, "netMappingDEG.net_upregulated_DYS.csv", row.names = FALSE)
write.csv(net.down, "netMappingDEG.net_downregulated_DYS.csv", row.names = FALSE)



gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

gene.up
length(gene.up)

gene.down
length(gene.down)



write.csv(gene.up, "netMappingDEG.net_upregulated_DYS.gene.csv", row.names = FALSE)
write.csv(gene.down, "netMappingDEG.net_downregulated_DYS.gene.csv", row.names = FALSE)




# to subset manually the dataframe NET

head(net)
colnames(net)

#### to extract manually
# thresh.fc = 0.05,
# thresh.p = 0.05,
# thresh.pc = 0.1, 

p_cut   <- 0.05
lfc_cut <- 0.05
pct_cut <- 0.10

# ---- ADDYS (AD vs Dyslexia) ----
addys_LR_bothUp <- subset(
  net,
  datasets == "addys" &
    ligand.pvalues   < p_cut &
    receptor.pvalues < p_cut &
    ligand.logFC     >  lfc_cut &
    receptor.logFC   >  lfc_cut &
    ligand.pct.1     >= pct_cut &
    receptor.pct.1   >= pct_cut
)

addys_LR_bothDown <- subset(
  net,
  datasets == "addys" &
    ligand.pvalues   < p_cut &
    receptor.pvalues < p_cut &
    ligand.logFC     < -lfc_cut &
    receptor.logFC   < -lfc_cut &
    ligand.pct.1     >= pct_cut &
    receptor.pct.1   >= pct_cut
)

dim(addys_LR_bothUp)
dim(addys_LR_bothDown)

# Write full filtered tables
write.csv(addys_LR_bothUp,   file = "addys_LR_bothUp.csv",   row.names = FALSE)
write.csv(addys_LR_bothDown, file = "addys_LR_bothDown.csv", row.names = FALSE)

# Write just interaction identifiers (useful quick lists)
write.table(addys_LR_bothUp$interaction_name,
            file = "addys_LR_bothUp.interactions.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(addys_LR_bothDown$interaction_name,
            file = "addys_LR_bothDown.interactions.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)



head(addys_LR_bothUp)
dim(addys_LR_bothUp)

head(addys_LR_bothDown)
dim(addys_LR_bothDown)

# ---- ADCTL (AD vs Control) ----
adctl_LR_bothUp <- subset(
  net,
  datasets == "adctl" &
    ligand.pvalues   < p_cut &
    receptor.pvalues < p_cut &
    ligand.logFC     >  lfc_cut &
    receptor.logFC   >  lfc_cut &
    ligand.pct.1     >= pct_cut &
    receptor.pct.1   >= pct_cut
)

adctl_LR_bothDown <- subset(
  net,
  datasets == "adctl" &
    ligand.pvalues   < p_cut &
    receptor.pvalues < p_cut &
    ligand.logFC     < -lfc_cut &
    receptor.logFC   < -lfc_cut &
    ligand.pct.1     >= pct_cut &
    receptor.pct.1   >= pct_cut
)

dim(adctl_LR_bothUp)
dim(adctl_LR_bothDown)

write.csv(adctl_LR_bothUp,   file = "adctl_LR_bothUp.csv",   row.names = FALSE)
write.csv(adctl_LR_bothDown, file = "adctl_LR_bothDown.csv", row.names = FALSE)

write.table(adctl_LR_bothUp$interaction_name,
            file = "adctl_LR_bothUp.interactions.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(adctl_LR_bothDown$interaction_name,
            file = "adctl_LR_bothDown.interactions.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)



head(adctl_LR_bothUp)
dim(adctl_LR_bothUp)

head(adctl_LR_bothDown)
dim(adctl_LR_bothDown)





pairLR.use.up = net.up[, "interaction_name", drop = F]

gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, 
             #    sources.use = c("L23_IT_Exc", "L4_Exc", "L56_Exc", "L6_Exc", "L6_IT_Exc", "Immune"),
             #    targets.use = c("L23_IT_Exc", "L4_Exc", "L56_Exc", "L6_Exc", "L6_IT_Exc", "Immune"),
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = T, 
                        title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object

pairLR.use.down = net.down[, "interaction_name", drop = F]

gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, 
              #   sources.use = c("L23_IT_Exc", "L4_Exc", "L56_Exc", "L6_Exc", "L6_IT_Exc", "Immune"),
              #   targets.use = c("L23_IT_Exc", "L4_Exc", "L56_Exc", "L6_Exc", "L6_IT_Exc", "Immune"),
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = T, 
                        title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
gg1 + gg2

# UP-regulated
ggsave("pairLR.use.up.netBubble.png", gg1, 
       width = 12, height = 10, dpi = 300)

# DOWN-regulated
ggsave("pairLR.use.down.netBubble.png", gg2, 
       width = 12, height = 10, dpi = 300)



# Chord diagram
# Visualize the upregulated signaling in the second dataset

# par(mfrow = c(1,2), xpd=TRUE)

# netVisual_chord_gene(object.list[[2]], 
 #                sources.use = c("L23_IT_Exc", "L4_Exc", "L56_Exc", "L6_Exc", "L6_IT_Exc", "Immune"),
 #                targets.use = c("L23_IT_Exc", "L4_Exc", "L56_Exc", "L6_Exc", "L6_IT_Exc", "Immune"),
#                     slot.name = 'net', 
#                     net = net.up, 
#                     lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))

# Visualize the downregulated signaling in the second dataset

# netVisual_chord_gene(object.list[[1]], 
  #               sources.use = c("L23_IT_Exc", "L4_Exc", "L56_Exc", "L6_Exc", "L6_IT_Exc", "Immune"),
  #               targets.use = c("L23_IT_Exc", "L4_Exc", "L56_Exc", "L6_Exc", "L6_IT_Exc", "Immune"), 
#                     slot.name = 'net', 
#                     net = net.down, 
#                     lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))

# You may try the function `netVisual_chord_cell` for visualizing individual signaling pathway



# CHORD DIAGRAMS COMMENTED OUT DUE TO SPACE ERRORS
# Uncomment and adjust parameters if needed

# png(file.path("cellchat_chord_up_addys.selected.excitatory.png"),
#     width = 10, height = 10, units = "in", res = 300)
# 
# par(xpd = TRUE)
# 
# tryCatch({
#   library(circlize)
#   circos.clear()
#   circos.par(gap.after = c(rep(2, length(unique(c("L23_IT_Exc", "L4_Exc", "L56_Exc", "L6_Exc", "L6_IT_Exc", "Immune")))-1), 10))
#   
#   netVisual_chord_gene(
#     object.list[[2]],
#     sources.use = c("L23_IT_Exc", "L4_Exc", "L56_Exc", "L6_Exc", "L6_IT_Exc", "Immune"),
#     targets.use = c("L23_IT_Exc", "L4_Exc", "L56_Exc", "L6_Exc", "L6_IT_Exc", "Immune"),
#     slot.name = "net",
#     net = net.up,
#     lab.cex = 0.6,
#     small.gap = 2,
#     big.gap = 10,
#     title.name = paste0("Up-regulated signaling in ", names(object.list)[2])
#   )
#   
#   circos.clear()
# }, error = function(e) {
#   cat("❌ Error creating chord diagram (up-regulated):", as.character(e), "\n")
#   plot.new()
#   text(0.5, 0.5, paste("Chord diagram failed:\n", as.character(e)), cex = 0.8)
#   circos.clear()
# })
# 
# dev.off()
# 


# png(file.path("cellchat_chord_down_addys.selected.excitatory.png"),
#     width = 10, height = 10, units = "in", res = 300)
# 
# par(xpd = TRUE)
# 
# tryCatch({
#   library(circlize)
#   circos.clear()
#   circos.par(gap.after = c(rep(2, length(unique(c("L23_IT_Exc", "L4_Exc", "L56_Exc", "L6_Exc", "L6_IT_Exc", "Immune")))-1), 10))
#   
#   netVisual_chord_gene(
#     object.list[[1]],
#     sources.use = c("L23_IT_Exc", "L4_Exc", "L56_Exc", "L6_Exc", "L6_IT_Exc", "Immune"),
#     targets.use = c("L23_IT_Exc", "L4_Exc", "L56_Exc", "L6_Exc", "L6_IT_Exc", "Immune"),
#     slot.name = "net",
#     net = net.down,
#     lab.cex = 0.6,
#     small.gap = 2,
#     big.gap = 10,
#     title.name = paste0("Down-regulated signaling in ", names(object.list)[2])
#   )
#   
#   circos.clear()
# }, error = function(e) {
#   cat("❌ Error creating chord diagram (down-regulated):", as.character(e), "\n")
#   plot.new()
#   text(0.5, 0.5, paste("Chord diagram failed:\n", as.character(e)), cex = 0.8)
#   circos.clear()
# })
# 
# dev.off()
# 


# tryCatch({
#   library(circlize)
#   circos.clear()
#   circos.par(gap.after = c(rep(2, length(unique(c("L23_IT_Exc", "L4_Exc", "L56_Exc", "L6_Exc", "L6_IT_Exc", "Immune")))-1), 10))
#   
#   netVisual_chord_gene(
#     object.list[[2]],
#     sources.use = c("L23_IT_Exc", "L4_Exc", "L56_Exc", "L6_Exc", "L6_IT_Exc", "Immune"),
#     targets.use = c("L23_IT_Exc", "L4_Exc", "L56_Exc", "L6_Exc", "L6_IT_Exc", "Immune"),
#     slot.name = "net",
#     net = net.up,
#     lab.cex = 0.6,
#     small.gap = 2,
#     big.gap = 10,
#     title.name = paste0("Up-regulated signaling in ", names(object.list)[2])
#   )
#   
#   circos.clear()
# }, error = function(e) {
#   cat("❌ Error creating chord diagram (up-regulated console):", as.character(e), "\n")
#   circos.clear()
# })
# 
# 
# tryCatch({
#   library(circlize)
#   circos.clear()
#   circos.par(gap.after = c(rep(2, length(unique(c("L23_IT_Exc", "L4_Exc", "L56_Exc", "L6_Exc", "L6_IT_Exc", "Immune")))-1), 10))
#   
#   netVisual_chord_gene(
#     object.list[[1]],
#     sources.use = c("L23_IT_Exc", "L4_Exc", "L56_Exc", "L6_Exc", "L6_IT_Exc", "Immune"),
#     targets.use = c("L23_IT_Exc", "L4_Exc", "L56_Exc", "L6_Exc", "L6_IT_Exc", "Immune"),
#     slot.name = "net",
#     net = net.down,
#     lab.cex = 0.6,
#     small.gap = 2,
#     big.gap = 10,
#     title.name = paste0("Down-regulated signaling in ", names(object.list)[2])
#   )
#   
#   circos.clear()
# }, error = function(e) {
#   cat("❌ Error creating chord diagram (down-regulated console):", as.character(e), "\n")
#   circos.clear()
# })



# visualize the enriched ligands in the second condition
computeEnrichmentScore(net.up, species = 'human', variable.both = TRUE)

png(file.path("cellchat_enrichedScore.net.up.addys.png"),
    width = 7, height = 6, units = "in", res = 300)

computeEnrichmentScore(
  net.up,
  species = "human",
  variable.both = TRUE
)

dev.off()

# visualize the enriched ligands in the first condition
computeEnrichmentScore(net.down, species = 'human')

png(file.path("cellchat_enrichedScore.net.down.addys.png"),
    width = 7, height = 6, units = "in", res = 300)

computeEnrichmentScore(
  net.down,
  species = "human",
  variable.both = TRUE
)

dev.off()



# to select the pathways to show !!



pathways.show <- c("SPP1") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) 

par(mfrow = c(1,2), xpd=TRUE)

for (i in 1:length(object.list)) {

netVisual_aggregate(object.list[[i]], 
                    signaling = pathways.show, 
                    layout = "circle", 
                    edge.weight.max = weight.max[1], 
                    edge.width.max = 10, 
                    signaling.name = paste(pathways.show, names(object.list)[i]))

} 

pathways.show <- "SPP1"

weight.max <- getMaxWeight(
  object.list,
  slot.name = "netP",
  attribute = pathways.show
)

png(file.path("SPP1_aggregate_circle.png"),
    width = 14, height = 7, units = "in", res = 300)

par(
  mfrow = c(1, 2),
  xpd = TRUE,
  mar = c(2, 2, 4, 2)  # good margins for circle plots
)

for (i in seq_along(object.list)) {
  netVisual_aggregate(
    object.list[[i]],
    signaling = pathways.show,
    layout = "circle",
    edge.weight.max = weight.max[1],
    edge.width.max = 10,
    signaling.name = paste(pathways.show, names(object.list)[i])
  )
}

dev.off()


print("Visually compare inferred cell–cell communication networks")

pathways.show <- c("SPP1") 

par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], 
                               signaling = pathways.show, 
                               color.heatmap = "Reds", 
                               title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}

options(repr.plot.width = 6, repr.plot.height = 6)
options(repr.plot.res = 200)  # Resolution

ComplexHeatmap::draw(ht[[1]])

png(file.path("SPP1_heatmap_adctl.png"),
    width = 6, height = 6, units = "in", res = 300)

ComplexHeatmap::draw(ht[[1]])

dev.off()


options(repr.plot.width = 6, repr.plot.height = 6)
options(repr.plot.res = 200)  # Resolution

ComplexHeatmap::draw(ht[[2]])

png(file.path("SPP1_heatmap_addys.png"),
    width = 6, height = 6, units = "in", res = 300)

ComplexHeatmap::draw(ht[[2]])

dev.off()


# Chord diagram

pathways.show <- c("SPP1") 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], 
                      signaling = pathways.show, 
                      layout = "chord", 
                      signaling.name = paste(pathways.show, names(object.list)[i]))
}

pathways.show <- "SPP1"

png(file.path("SPP1_chord.png"),
    width = 14, height = 7, units = "in", res = 300)

par(
  mfrow = c(1, 2),
  xpd = TRUE,
  mar = c(2, 2, 4, 2)  # good margins for chord plots
)

for (i in seq_along(object.list)) {
  netVisual_aggregate(
    object.list[[i]],
    signaling = pathways.show,
    layout = "chord",
    signaling.name = paste(pathways.show, names(object.list)[i])
  )
}

dev.off()





# to select tha pathways to show !!!


cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("adctl", "addys")) # set factor level



plotGeneExpression(cellchat, signaling = "SPP1", split.by = "datasets", colors.ggplot = T, type = "violin")

p <- plotGeneExpression(
  cellchat,
  signaling = "SPP1",
  split.by = "datasets",
  colors.ggplot = TRUE,
  type = "violin"
)

ggsave(
  filename = "SPP1_gene_expression_violin.png",
  plot = p,
  width = 7,
  height = 5,
  units = "in",
  dpi = 300
)


head(cellchat@meta, 2 )

tail(cellchat@meta, 2 )


################################################
# FUNCTION: Generate pathway-specific visualizations
################################################

generate_pathway_plots <- function(pathway_list, output_label, object.list, cellchat, output_dir = ".") {
  
  cat(sprintf("\n=== Generating plots for %s (%d pathways) ===\n", 
              output_label, length(pathway_list)))
  
  # Create subdirectory for this category
  category_dir <- file.path(output_dir, output_label)
  dir.create(category_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Progress tracking
  success_count <- 0
  failed_pathways <- c()
  
  for (pathway in pathway_list) {
    
    cat(sprintf("Processing pathway: %s\n", pathway))
    
    tryCatch({
      
      # 1. CIRCLE/NETWORK PLOTS
      tryCatch({
        weight.max <- getMaxWeight(object.list, slot.name = "netP", attribute = pathway)
        
        png(file.path(category_dir, sprintf("%s_%s_aggregate_circle.png", pathway, output_label)),
            width = 14, height = 7, units = "in", res = 300)
        
        par(mfrow = c(1, 2), xpd = TRUE, mar = c(2, 2, 4, 2))
        
        for (i in seq_along(object.list)) {
          netVisual_aggregate(
            object.list[[i]],
            signaling = pathway,
            layout = "circle",
            edge.weight.max = weight.max[1],
            edge.width.max = 10,
            signaling.name = paste(pathway, names(object.list)[i])
          )
        }
        
        dev.off()
        cat(sprintf("  ✓ Circle plot saved\n"))
        
      }, error = function(e) {
        if (dev.cur() > 1) dev.off()
        cat(sprintf("  ⚠️  Circle plot failed: %s\n", as.character(e)))
      })
      
      # 2. HEATMAPS
      tryCatch({
        ht <- list()
        
        for (i in 1:length(object.list)) {
          ht[[i]] <- netVisual_heatmap(
            object.list[[i]], 
            signaling = pathway, 
            color.heatmap = "Reds", 
            title.name = paste(pathway, "signaling", names(object.list)[i])
          )
        }
        
        # Save heatmap for condition 1
        png(file.path(category_dir, sprintf("%s_%s_heatmap_%s.png", 
                                             pathway, output_label, names(object.list)[1])),
            width = 6, height = 6, units = "in", res = 300)
        ComplexHeatmap::draw(ht[[1]])
        dev.off()
        
        # Save heatmap for condition 2
        png(file.path(category_dir, sprintf("%s_%s_heatmap_%s.png", 
                                             pathway, output_label, names(object.list)[2])),
            width = 6, height = 6, units = "in", res = 300)
        ComplexHeatmap::draw(ht[[2]])
        dev.off()
        
        cat(sprintf("  ✓ Heatmaps saved\n"))
        
      }, error = function(e) {
        if (dev.cur() > 1) dev.off()
        cat(sprintf("  ⚠️  Heatmaps failed: %s\n", as.character(e)))
      })
      
      # 3. CHORD DIAGRAMS
      tryCatch({
        png(file.path(category_dir, sprintf("%s_%s_chord.png", pathway, output_label)),
            width = 14, height = 7, units = "in", res = 300)
        
        par(mfrow = c(1, 2), xpd = TRUE, mar = c(2, 2, 4, 2))
        
        for (i in seq_along(object.list)) {
          netVisual_aggregate(
            object.list[[i]],
            signaling = pathway,
            layout = "chord",
            signaling.name = paste(pathway, names(object.list)[i])
          )
        }
        
        dev.off()
        cat(sprintf("  ✓ Chord diagram saved\n"))
        
      }, error = function(e) {
        if (dev.cur() > 1) dev.off()
        cat(sprintf("  ⚠️  Chord diagram failed: %s\n", as.character(e)))
      })
      
      # 4. GENE EXPRESSION VIOLIN PLOTS
      tryCatch({
        # Set factor levels
        cellchat@meta$datasets <- factor(cellchat@meta$datasets, 
                                          levels = c(names(object.list)[1], names(object.list)[2]))
        
        p <- plotGeneExpression(
          cellchat,
          signaling = pathway,
          split.by = "datasets",
          colors.ggplot = TRUE,
          type = "violin"
        )
        
        ggsave(
          filename = file.path(category_dir, sprintf("%s_%s_gene_expression_violin.png", 
                                                      pathway, output_label)),
          plot = p,
          width = 7,
          height = 5,
          units = "in",
          dpi = 300
        )
        
        cat(sprintf("  ✓ Gene expression plot saved\n"))
        
      }, error = function(e) {
        cat(sprintf("  ⚠️  Gene expression plot failed: %s\n", as.character(e)))
      })
      
      success_count <- success_count + 1
      
    }, error = function(e) {
      cat(sprintf("  ❌ Complete failure for pathway %s: %s\n", pathway, as.character(e)))
      failed_pathways <- c(failed_pathways, pathway)
    })
    
    # Garbage collection every 5 pathways
    if (success_count %% 5 == 0) {
      gc(verbose = FALSE)
    }
  }
  
  # Summary
  cat("\n╔════════════════════════════════════════════╗\n")
  cat(sprintf("║  Summary for %s\n", output_label))
  cat("╚════════════════════════════════════════════╝\n")
  cat(sprintf("Total pathways: %d\n", length(pathway_list)))
  cat(sprintf("✓ Successfully processed: %d\n", success_count))
  cat(sprintf("❌ Failed: %d\n", length(failed_pathways)))
  
  if (length(failed_pathways) > 0) {
    cat("\nFailed pathways:\n")
    cat(paste("  -", failed_pathways, collapse = "\n"))
    cat("\n")
  }
  
  cat(sprintf("\n✓ All plots saved to: %s\n\n", category_dir))
  
  return(invisible(list(
    success_count = success_count,
    failed_pathways = failed_pathways,
    output_dir = category_dir
  )))
}

################################################
# RUN THE FUNCTION FOR ALL PATHWAY CATEGORIES
################################################

cat("\n")
cat("╔═══════════════════════════════════════════════════════════╗\n")
cat("║  GENERATING PATHWAY-SPECIFIC VISUALIZATIONS              ║\n")
cat("╚═══════════════════════════════════════════════════════════╝\n")
cat("\n")

# ═══════════════════════════════════════════════════════════
# ALL CELLCHAT DATABASE PATHWAYS
# ═══════════════════════════════════════════════════════════

cellchatdb_pathways <- c(
  "ACTIVIN", "ADGRE5", "ADIPONECTIN", "AGRN", "AGT", "ALCAM", "AMH", "ANGPT",
  "ANGPTL", "ANNEXIN", "ANXA1", "APELIN", "APJ", "APP", "APRIL", "AVP",
  "BAFF", "BAG", "BMP", "BMP10", "BRADYKININ", "BSP", "BTLA", "CADM",
  "CALCR", "CCK", "CCL", "CD137", "CD22", "CD226", "CD23", "CD30",
  "CD34", "CD39", "CD40", "CD45", "CD46", "CD48", "CD6", "CD70",
  "CD80", "CD86", "CD96", "CD99", "CDH", "CDH1", "CDH5", "CEACAM",
  "CHAD", "CHEMERIN", "CLDN", "CLEC", "CNTN", "COLLAGEN", "COMPLEMENT", "CRH",
  "CSF", "CSF3", "CSPG4", "CX3C", "CXCL", "DESMOSOME", "DMP1", "DSPP",
  "EDA", "EDN", "EGF", "ENHO", "EPGN", "EPHA", "EPHB", "EPO",
  "ESAM", "FASLG", "FGF", "FLT3", "FN1", "FSH", "GALANIN", "GALECTIN",
  "GAS", "GCG", "GDF", "GDNF", "GH", "GHRELIN", "GHRH", "GIPR",
  "GITRL", "GNRH", "GP1BA", "GPR", "GRN", "GUCA", "HCRT", "HGF",
  "HH", "HSPG", "ICAM", "ICOS", "IFN-I", "IFN-II", "IGF", "IL1",
  "IL10", "IL12", "IL16", "IL17", "IL2", "IL4", "IL6", "INSULIN",
  "ITGB2", "JAM", "KISS1", "KIT", "L1CAM", "LAMININ", "LCK", "LEP",
  "LHB", "LIFR", "LIGHT", "LT", "MADCAM", "MAG", "MELANOCORTIN", "MHC-I",
  "MHC-II", "MIF", "MK", "MPZ", "MSTN", "NCAM", "ncWNT", "NECTIN",
  "NEGR", "NGF", "NGL", "NKG2D", "NMU", "NODAL", "NOTCH", "NPFF",
  "NPNT", "NPR1", "NPR2", "NPS", "NPVF", "NPW-B", "NPY", "NRG",
  "NRXN", "NT", "NTS", "OCLN", "OPIOID", "OSM", "OSTN", "OX40",
  "OXT", "PACAP", "PARs", "PDGF", "PD-L1", "PDL2", "PECAM1", "PERIOSTIN",
  "PMCH", "PRL", "PRLH", "PROK", "PROS", "PSAP", "PTH", "PTN",
  "PTPRM", "PVR", "QRFP", "RANKL", "RELN", "RESISTIN", "RLN", "SAA",
  "SCT", "SELE", "SELL", "SELPLG", "SEMA3", "SEMA4", "SEMA5", "SEMA6",
  "SEMA7", "SLURP", "SN", "SOMATOSTATIN", "SPP1", "TAC", "TENASCIN", "TGFb",
  "THBS", "THPO", "THY1", "TIGIT", "TNF", "TRAIL", "TRH", "TSH",
  "TWEAK", "UCN", "UGRP1", "UROTENSIN", "UTS2", "VCAM", "VEGF", "VEGI",
  "VIP", "VISFATIN", "VISTA", "VTN", "VWF", "WNT", "XCR"
)

cat("=== CELLCHAT DATABASE PATHWAYS ===\n")
cat("Total pathways in CellChatDB:", length(cellchatdb_pathways), "\n\n")

# Check which pathways are actually present in your data
if (exists("cellchat")) {
  detected_pathways <- unique(cellchat@netP$pathways)
  cat("Pathways detected in your data:", length(detected_pathways), "\n")
  
  # Find which CellChatDB pathways are in your data
  pathways_in_data <- cellchatdb_pathways[cellchatdb_pathways %in% detected_pathways]
  pathways_not_in_data <- cellchatdb_pathways[!cellchatdb_pathways %in% detected_pathways]
  
  cat("CellChatDB pathways found in your data:", length(pathways_in_data), "\n")
  cat("CellChatDB pathways NOT in your data:", length(pathways_not_in_data), "\n\n")
}

# Main output directory for pathway plots
pathway_plots_dir <- "pathway_specific_plots"
dir.create(pathway_plots_dir, showWarnings = FALSE)

# ═══════════════════════════════════════════════════════════
# GENERATE PLOTS FOR ALL CELLCHATDB PATHWAYS IN YOUR DATA
# ═══════════════════════════════════════════════════════════

cat("\n=== GENERATING PLOTS FOR ALL CELLCHATDB PATHWAYS ===\n")

# Use only the pathways that are actually in your data
if (exists("pathways_in_data") && length(pathways_in_data) > 0) {
  generate_pathway_plots(
    pathway_list = pathways_in_data,
    output_label = "all_cellchatdb_pathways",
    object.list = object.list,
    cellchat = cellchat,
    output_dir = pathway_plots_dir
  )
} else {
  cat("⚠️  No CellChatDB pathways found in data\n")
}

# ═══════════════════════════════════════════════════════════
# PATHWAYS FROM NETMAPPINGDEG RESULTS (using variable 'net')
# ═══════════════════════════════════════════════════════════

cat("\n=== EXTRACTING PATHWAYS FROM NETMAPPINGDEG RESULTS ===\n")

# Check if the 'net' variable exists (created by netMappingDEG function)
if (exists("net") && !is.null(net)) {
  cat("✓ Found 'net' variable with", nrow(net), "rows\n")
  
  # Check if pathway_name column exists
  if ("pathway_name" %in% colnames(net)) {
    # Extract unique pathway names
    pathways_netmapping <- unique(net$pathway_name)
    # Remove NA values if any
    pathways_netmapping <- pathways_netmapping[!is.na(pathways_netmapping)]
    # Remove empty strings
    pathways_netmapping <- pathways_netmapping[pathways_netmapping != ""]
    
    cat("✓ Found", length(pathways_netmapping), "unique pathways in pathway_name column\n")
    
    if (length(pathways_netmapping) > 0) {
      # Filter to only pathways that exist in the CellChat object
      if (exists("cellchat") && !is.null(cellchat@netP$pathways)) {
        detected_pathways <- unique(cellchat@netP$pathways)
        pathways_netmapping_valid <- pathways_netmapping[pathways_netmapping %in% detected_pathways]
        pathways_netmapping_not_in_data <- pathways_netmapping[!pathways_netmapping %in% detected_pathways]
        
        cat("  - Pathways found in CellChat data:", length(pathways_netmapping_valid), "\n")
        if (length(pathways_netmapping_not_in_data) > 0) {
          cat("  - Pathways NOT in CellChat data:", length(pathways_netmapping_not_in_data), "\n")
          cat("    (These will be skipped)\n")
        }
        
        if (length(pathways_netmapping_valid) > 0) {
          cat("\n=== GENERATING PLOTS FOR NETMAPPINGDEG PATHWAYS ===\n")
          
          generate_pathway_plots(
            pathway_list = pathways_netmapping_valid,
            output_label = "netMappingDEG_pathways",
            object.list = object.list,
            cellchat = cellchat,
            output_dir = pathway_plots_dir
          )
        } else {
          cat("⚠️  No valid pathways found in both netMappingDEG results and CellChat data\n")
        }
      } else {
        cat("⚠️  CellChat object not available - cannot validate pathways\n")
        cat("    Attempting to generate plots anyway...\n")
        
        generate_pathway_plots(
          pathway_list = pathways_netmapping,
          output_label = "netMappingDEG_pathways",
          object.list = object.list,
          cellchat = cellchat,
          output_dir = pathway_plots_dir
        )
      }
    } else {
      cat("⚠️  No pathways found in pathway_name column\n")
    }
  } else {
    cat("❌ Error: 'pathway_name' column not found in 'net' variable\n")
    cat("   Available columns:", paste(colnames(net), collapse = ", "), "\n")
  }
} else {
  cat("⚠️  Variable 'net' not found or is NULL\n")
  cat("   Make sure netMappingDEG() has been run before this section\n")
  cat("   Skipping netMappingDEG pathway plots...\n")
}

# ═══════════════════════════════════════════════════════════
# ADDITIONAL PATHWAY CATEGORIES
# ═══════════════════════════════════════════════════════════

cat("\n=== GENERATING PLOTS FOR DIFFERENTIAL PATHWAY CATEGORIES ===\n")

# 1. net.down pathways
if (exists("net.down") && nrow(net.down) > 0) {
  pathways_net_down <- unique(net.down$pathway_name)
  cat(sprintf("Found %d pathways in net.down\n", length(pathways_net_down)))
  
  generate_pathway_plots(
    pathway_list = pathways_net_down,
    output_label = "net.down",
    object.list = object.list,
    cellchat = cellchat,
    output_dir = pathway_plots_dir
  )
} else {
  cat("⚠️  net.down not found or empty\n")
}

# 2. net.up pathways
if (exists("net.up") && nrow(net.up) > 0) {
  pathways_net_up <- unique(net.up$pathway_name)
  cat(sprintf("Found %d pathways in net.up\n", length(pathways_net_up)))
  
  generate_pathway_plots(
    pathway_list = pathways_net_up,
    output_label = "net.up",
    object.list = object.list,
    cellchat = cellchat,
    output_dir = pathway_plots_dir
  )
} else {
  cat("⚠️  net.up not found or empty\n")
}

# 3. addys_LR_bothDown pathways
if (exists("addys_LR_bothDown") && nrow(addys_LR_bothDown) > 0) {
  pathways_addys_down <- unique(addys_LR_bothDown$pathway_name)
  cat(sprintf("Found %d pathways in addys_LR_bothDown\n", length(pathways_addys_down)))
  
  generate_pathway_plots(
    pathway_list = pathways_addys_down,
    output_label = "addys_LR_bothDown",
    object.list = object.list,
    cellchat = cellchat,
    output_dir = pathway_plots_dir
  )
} else {
  cat("⚠️  addys_LR_bothDown not found or empty\n")
}

# 4. addys_LR_bothUp pathways
if (exists("addys_LR_bothUp") && nrow(addys_LR_bothUp) > 0) {
  pathways_addys_up <- unique(addys_LR_bothUp$pathway_name)
  cat(sprintf("Found %d pathways in addys_LR_bothUp\n", length(pathways_addys_up)))
  
  generate_pathway_plots(
    pathway_list = pathways_addys_up,
    output_label = "addys_LR_bothUp",
    object.list = object.list,
    cellchat = cellchat,
    output_dir = pathway_plots_dir
  )
} else {
  cat("⚠️  addys_LR_bothUp not found or empty\n")
}

# 5. adctl_LR_bothDown pathways
if (exists("adctl_LR_bothDown") && nrow(adctl_LR_bothDown) > 0) {
  pathways_adctl_down <- unique(adctl_LR_bothDown$pathway_name)
  cat(sprintf("Found %d pathways in adctl_LR_bothDown\n", length(pathways_adctl_down)))
  
  generate_pathway_plots(
    pathway_list = pathways_adctl_down,
    output_label = "adctl_LR_bothDown",
    object.list = object.list,
    cellchat = cellchat,
    output_dir = pathway_plots_dir
  )
} else {
  cat("⚠️  adctl_LR_bothDown not found or empty\n")
}

# 6. adctl_LR_bothUp pathways
if (exists("adctl_LR_bothUp") && nrow(adctl_LR_bothUp) > 0) {
  pathways_adctl_up <- unique(adctl_LR_bothUp$pathway_name)
  cat(sprintf("Found %d pathways in adctl_LR_bothUp\n", length(pathways_adctl_up)))
  
  generate_pathway_plots(
    pathway_list = pathways_adctl_up,
    output_label = "adctl_LR_bothUp",
    object.list = object.list,
    cellchat = cellchat,
    output_dir = pathway_plots_dir
  )
} else {
  cat("⚠️  adctl_LR_bothUp not found or empty\n")
}

cat("\n")
cat("╔═══════════════════════════════════════════════════════════╗\n")
cat("║  PATHWAY VISUALIZATIONS COMPLETE                          ║\n")
cat("╚═══════════════════════════════════════════════════════════╝\n")
cat("\n")
cat("NOTE: All CellChatDB pathways in your data have been visualized.\n")
cat("      Additional differential pathway categories (net.up, net.down,\n")
cat("      LR_bothUp, LR_bothDown) are also available above.\n")
cat("      You can comment out the differential categories if you only\n")
cat("      want the all_cellchatdb_pathways plots.\n\n")




# ═══════════════════════════════════════════════════════════
# SAVE FINAL RESULTS
# ═══════════════════════════════════════════════════════════
cat("=== SAVING FINAL CELLCHAT OBJECT ===\n")
saveRDS(cellchat, file = "cellchat.adctl.addys.rds")  # Use saveRDS() for .rds files
cat("✓ Saved: cellchat.adctl.addys.rds\n\n")

# ═══════════════════════════════════════════════════════════
# END LOGGING
# ═══════════════════════════════════════════════════════════
cat("\n=================================================\n")
cat("Analysis completed at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Log file:", log_file, "\n")
cat("=================================================\n")

# ============================================================
# CellChat @data.signaling: save + split by condition + preview
# + print & save selected genes (SPP1, RELN) by condition
# ============================================================

suppressPackageStartupMessages({
  library(Matrix)
})

stopifnot(!is.null(cellchat@data.signaling))
stopifnot(inherits(cellchat@data.signaling, "dgCMatrix"))

mat_all <- cellchat@data.signaling

# ----------------------------
# 0) Align per-cell condition
# ----------------------------
stopifnot(!is.null(cellchat@meta))
stopifnot("datasets" %in% colnames(cellchat@meta))

cond <- cellchat@meta$datasets
names(cond) <- rownames(cellchat@meta)  # names are cell IDs

# align condition labels to matrix columns (cells)
cond_ds <- cond[colnames(mat_all)]

cat("Condition labels (table):\n")
print(table(cond_ds, useNA = "ifany"))

cat("\nSanity checks:\n")
cat("All matrix cells in meta? ", all(colnames(mat_all) %in% rownames(cellchat@meta)), "\n")
cat("Missing labels after alignment (NA)? ", sum(is.na(cond_ds)), "\n\n")

# ----------------------------
# 1) Split cells by condition
# ----------------------------
# Adjust these strings if your labels differ
cells_adctl <- colnames(mat_all)[cond_ds == "adctl"]
cells_addys <- colnames(mat_all)[cond_ds == "addys"]  # or "ad_dys" depending on your naming

cat("Cells adctl:", length(cells_adctl), "\n")
cat("Cells addys:", length(cells_addys), "\n\n")

# Guardrails (helpful if labels don't match exactly)
if (length(cells_adctl) == 0) {
  cat("⚠️  No cells matched label 'adctl'. Available labels:\n")
  print(unique(na.omit(cond_ds)))
}
if (length(cells_addys) == 0) {
  cat("⚠️  No cells matched label 'addys'. Available labels:\n")
  print(unique(na.omit(cond_ds)))
}

stopifnot(length(cells_adctl) > 0, length(cells_addys) > 0)

# ----------------------------
# 2) Split sparse matrices
# ----------------------------
mat_adctl <- mat_all[, cells_adctl, drop = FALSE]
mat_addys <- mat_all[, cells_addys, drop = FALSE]

cat("Dims adctl: ", paste(dim(mat_adctl), collapse = " x "), "\n")
cat("Dims addys: ", paste(dim(mat_addys), collapse = " x "), "\n\n")

# ----------------------------
# 3) Safe previews (DON'T print full matrix)
# ----------------------------
cat("Preview [1:10 genes x 1:10 cells] (sparse print):\n")
print(mat_adctl[1:min(10, nrow(mat_adctl)), 1:min(10, ncol(mat_adctl)), drop = FALSE])
print(mat_addys[1:min(10, nrow(mat_addys)), 1:min(10, ncol(mat_addys)), drop = FALSE])

cat("\nDense preview [1:5 x 1:5] (small only):\n")
print(as.matrix(mat_adctl[1:min(5, nrow(mat_adctl)), 1:min(5, ncol(mat_adctl)), drop = FALSE]))
print(as.matrix(mat_addys[1:min(5, nrow(mat_addys)), 1:min(5, ncol(mat_addys)), drop = FALSE]))

# Sparsity stats
cat("\nSparsity stats:\n")
cat("nnzero adctl:", nnzero(mat_adctl), "fraction:", nnzero(mat_adctl)/length(mat_adctl), "\n")
cat("nnzero addys:", nnzero(mat_addys), "fraction:", nnzero(mat_addys)/length(mat_addys), "\n\n")

# ----------------------------
# 4) Optional: gene-level summaries per condition
# ----------------------------
mean_adctl <- Matrix::rowMeans(mat_adctl)
mean_addys <- Matrix::rowMeans(mat_addys)

cat("Top genes by mean (adctl):\n")
print(head(sort(mean_adctl, decreasing = TRUE), 20))

cat("\nTop genes by mean (addys):\n")
print(head(sort(mean_addys, decreasing = TRUE), 20))

# ----------------------------
# 5) Optional: print selected genes safely (SPP1/RELN/SEMA5)
# ----------------------------
genes <- c("SPP1", "RELN", "SEMA5")

genes_present <- intersect(genes, rownames(mat_all))
genes_missing <- setdiff(genes, rownames(mat_all))

cat("\nSelected genes present:", paste(genes_present, collapse = ", "), "\n")
if (length(genes_missing) > 0) cat("Selected genes missing:", paste(genes_missing, collapse = ", "), "\n")

if (length(genes_present) > 0) {
  cat("\nSelected genes preview (first 10 cells):\n")
  print(mat_adctl[genes_present, 1:min(10, ncol(mat_adctl)), drop = FALSE])
  print(mat_addys[genes_present, 1:min(10, ncol(mat_addys)), drop = FALSE])
}

# ============================================================
# Saving options (pick what you need)
# ============================================================

# ----------------------------
# A) BEST: Save sparse matrix as RDS
# ----------------------------
saveRDS(mat_all,   file = "cellchat_data.signaling.rds")
saveRDS(mat_adctl, file = "cellchat_data.signaling_adctl.rds")
saveRDS(mat_addys, file = "cellchat_data.signaling_addys.rds")

# Reload later:
# mat_all <- readRDS("cellchat_data.signaling.rds")

# ----------------------------
# B) Matrix Market .mtx + names (for Python/Scanpy)
# ----------------------------
writeMM(mat_all, file = "cellchat_data.signaling.mtx")
write.table(rownames(mat_all),
            file = "cellchat_data.signaling_genes.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(colnames(mat_all),
            file = "cellchat_data.signaling_cells.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# ----------------------------
# C) OPTIONAL: Save matrix + metadata together (recommended)
# ----------------------------
saveRDS(
  list(
    data.signaling = mat_all,
    genes = rownames(mat_all),
    cells = colnames(mat_all),
    datasets = cellchat@meta$datasets,
    cond_ds = cond_ds
  ),
  file = "cellchat_data.signaling_with_metadata.rds"
)

# ----------------------------
# D) CSV (NOT recommended for 87k+ cells; huge & slow)
# ----------------------------
# write.csv(as.matrix(mat_all), file = "cellchat_data.signaling.csv")

# ============================================================
# Extra: Print & save SPP1/RELN exactly (and previews)
# ============================================================

genes2 <- c("SPP1", "RELN")
genes2_present <- intersect(genes2, rownames(mat_all))
genes2_missing <- setdiff(genes2, rownames(mat_all))

cat("\nSelected genes2 present:", paste(genes2_present, collapse = ", "), "\n")
if (length(genes2_missing) > 0) {
  cat("Selected genes2 missing:", paste(genes2_missing, collapse = ", "), "\n")
}

stopifnot(length(genes2_present) > 0)

# Print exactly as you showed (first 10 cells)
cat("\nSPP1/RELN — adctl (first 10 cells):\n")
print(mat_adctl[genes2_present, 1:min(10, ncol(mat_adctl)), drop = FALSE])

cat("\nSPP1/RELN — addys (first 10 cells):\n")
print(mat_addys[genes2_present, 1:min(10, ncol(mat_addys)), drop = FALSE])

# Save sparse (lossless, best for R)
saveRDS(
  mat_adctl[genes2_present, , drop = FALSE],
  file = "cellchat_data.signaling_genes_SPP1_RELN_adctl.rds"
)
saveRDS(
  mat_addys[genes2_present, , drop = FALSE],
  file = "cellchat_data.signaling_genes_SPP1_RELN_addys.rds"
)

# Save small CSV previews (human-readable; keep small)
write.csv(
  as.matrix(mat_adctl[genes2_present, 1:min(50, ncol(mat_adctl)), drop = FALSE]),
  file = "cellchat_data.signaling_genes_SPP1_RELN_adctl_preview.csv"
)
write.csv(
  as.matrix(mat_addys[genes2_present, 1:min(50, ncol(mat_addys)), drop = FALSE]),
  file = "cellchat_data.signaling_genes_SPP1_RELN_addys_preview.csv"
)

# Optional: long-format table (great for plotting)
make_long_df <- function(mat, condition) {
  df <- as.data.frame(as.matrix(mat))
  df$gene <- rownames(df)
  df_long <- reshape(
    df,
    direction = "long",
    varying = setdiff(colnames(df), "gene"),
    v.names = "expression",
    timevar = "cell",
    times = setdiff(colnames(df), "gene")
  )
  df_long$condition <- condition
  df_long
}

df_long_adctl <- make_long_df(
  mat_adctl[genes2_present, 1:min(50, ncol(mat_adctl)), drop = FALSE],
  "adctl"
)
df_long_addys <- make_long_df(
  mat_addys[genes2_present, 1:min(50, ncol(mat_addys)), drop = FALSE],
  "addys"
)

write.csv(
  rbind(df_long_adctl, df_long_addys),
  file = "cellchat_data.signaling_genes_SPP1_RELN_long_preview.csv",
  row.names = FALSE
)

cat("\n✓ Saved: full matrices (RDS/MTX) + SPP1/RELN sparse & preview files\n")

sink()  # Close the log file


