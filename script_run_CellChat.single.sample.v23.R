##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################

# conda activate cellbender

# Set Python path BEFORE loading any libraries that use reticulate
# Must be set before reticulate initializes (which happens on first import)
python_path <- "~/.virtualenvs/r-reticulate/bin/python"
python_path_expanded <- path.expand(python_path)

if (file.exists(python_path_expanded)) {
  Sys.setenv(RETICULATE_PYTHON = python_path_expanded)
  cat("✓ Python environment set to:", python_path_expanded, "\n")
} else {
  cat("⚠️  Custom Python path not found:", python_path_expanded, "\n")
  cat("   Will use system default Python\n")
}

##########################################################################################################

library(NMF)
library(CellChat)
library(ggalluvial)
library(patchwork)
library(Seurat)
options(stringsAsFactors = FALSE)

library(future)
library(future.apply)
library(progress)
library(BiocParallel)

# SETUP BIOCPARALLEL
n_cores <- parallel::detectCores()
n_workers <- min(4, max(1, n_cores - 2))

cat(sprintf("Setting up BiocParallel: %d cores detected, using %d workers\n", 
            n_cores, n_workers))

if (.Platform$OS.type == "unix" && !grepl("darwin", Sys.info()["sysname"], ignore.case = TRUE)) {
  BPPARAM <- MulticoreParam(workers = n_workers, RNGseed = 1337)
  cat("✓ Using MulticoreParam (fork-based parallel processing)\n")
} else {
  BPPARAM <- SnowParam(workers = n_workers, RNGseed = 1337, type = "SOCK")
  cat("✓ Using SnowParam (socket-based parallel processing)\n")
}

register(BPPARAM)
cat("✓ BiocParallel registered as default parallel backend\n")

options(future.globals.maxSize = 20 * 1024^3)
options(future.seed = TRUE)
future::plan(sequential)

packageVersion("Seurat")
set.seed(1337)

if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(1)
  RhpcBLASctl::omp_set_num_threads(1)
  cat("✓ BLAS/OMP threads set to 1 (safe for parallel processing)\n")
} else {
  cat("⚠️  RhpcBLASctl not available - cannot control BLAS threads\n")
}

options(repr.plot.width = 8, repr.plot.height = 6, repr.plot.res = 160)
options(bitmapType = "cairo")

library(ggplot2)
theme_set(theme_classic(base_size = 14))

update_geom_defaults("point", list(size = 1.2, alpha = 0.8))
update_geom_defaults("text",  list(size = 4))

##########################################################################################################
# MEMORY MONITORING AND GARBAGE COLLECTION FUNCTIONS
##########################################################################################################

check_memory_usage <- function(threshold = 90, force_gc = TRUE) {
  if (force_gc) {
    gc(verbose = FALSE, full = TRUE)
  }
  
  if (file.exists("/proc/meminfo")) {
    meminfo <- readLines("/proc/meminfo")
    mem_total <- as.numeric(sub(".*:\\s+([0-9]+).*", "\\1", 
                                 grep("^MemTotal:", meminfo, value = TRUE)))
    mem_available <- as.numeric(sub(".*:\\s+([0-9]+).*", "\\1", 
                                     grep("^MemAvailable:", meminfo, value = TRUE)))
    mem_used <- mem_total - mem_available
    mem_usage_pct <- (mem_used / mem_total) * 100
    
    cat(sprintf("💾 Memory: %.1f%% used (%.2f GB / %.2f GB)\n", 
                mem_usage_pct, mem_used / 1024^2, mem_total / 1024^2))
    
    if (mem_usage_pct > threshold) {
      stop(sprintf(
        "\n❌ MEMORY THRESHOLD EXCEEDED: %.1f%% > %d%%\n   Used: %.2f GB / %.2f GB\n   Script stopped to prevent system crash.",
        mem_usage_pct, threshold, mem_used / 1024^2, mem_total / 1024^2
      ))
    }
    return(invisible(mem_usage_pct))
  } else {
    cat("⚠️  Memory monitoring not available on this system\n")
    if (force_gc) cat("✓ Garbage collection performed\n")
    return(invisible(NA))
  }
}

do_gc <- function(verbose = TRUE) {
  if (verbose) cat("🧹 Running garbage collection...\n")
  gc(verbose = FALSE, full = TRUE)
  invisible(NULL)
}

safe_dev_off <- function() {
  tryCatch({
    while (grDevices::dev.cur() > 1) {
      grDevices::dev.off()
    }
  }, error = function(e) {
    invisible(NULL)
  })
}

is_gg <- function(x) {
  inherits(x, "ggplot")
}

safe_save_plot <- function(plot_obj, filename, width = 10, height = 8, 
                           units = "in", res = 300, ...) {
  if (is_gg(plot_obj)) {
    ggsave(filename = filename, plot = plot_obj, width = width,
           height = height, units = units, dpi = res, ...)
    cat("✓ ggplot saved to:", filename, "\n")
  } else {
    cat("⚠️  Non-ggplot object - ensure plot is called within png()/dev.off()\n")
  }
  invisible(filename)
}

cat("\n✓ Memory monitoring functions loaded\n")
cat("  - check_memory_usage(threshold = 90): Check and stop if memory > 90%\n")
cat("  - do_gc(): Force garbage collection\n")
cat("  - safe_dev_off(): Safely close all graphics devices\n")
cat("  - is_gg(x): Check if object is ggplot\n")
cat("  - safe_save_plot(): Smart plot saving (auto-detects ggplot vs base)\n\n")

cat("=== Initial Memory Status ===\n")
check_memory_usage(threshold = 90)
cat("\n")

##########################################################################################################

Sys.setenv(TMPDIR = "/mnt/nfs/CX000008_DS1/projects/btanasa/tmp")
tempdir()

library(Matrix)
library(dplyr)
library(ggrepel)
library(reshape2)
library(gridExtra)
library(data.table)
library(vioplot)
library(harmony)
library(cowplot)
library(foreach)
library(doParallel)
library(fs)
library(MAST)
library(DoubletFinder)
library("destiny")
library(scran)
library(scater)
library(ggbeeswarm)

################################################
################################################

setwd("/mnt/nfs/CX000008_DS1/projects/jaeyeon/fastq_file_Dyslexia_r1/zanalysis_bogdan/samples.merged_AG_Harmony_res0.1_anno_091625")
output_dir <- "/mnt/nfs/CX000008_DS1/projects/jaeyeon/fastq_file_Dyslexia_r1/zanalysis_bogdan/samples.merged_AG_Harmony_res0.1_anno_091625"
print(list.files())

# Load the Seurat object
rds_filename <- "merged_AG_Harmony_res0.1_anno_091625.AD_CTL.rds"
data_human <- readRDS(rds_filename)

check_memory_usage(threshold = 90)

sample_name <- sub("\\.rds$", "", rds_filename)
cat("Sample name extracted:", sample_name, "\n")

plot_dir <- paste0(sample_name, '_cellchat')
dir.create(plot_dir, showWarnings = FALSE)

##########################################################################################################
# LOGGING SETUP
##########################################################################################################

# Create log file with timestamp in the plot directory
log_file <- file.path(plot_dir, paste0("cellchat_analysis_", 
                                       format(Sys.time(), "%Y%m%d_%H%M%S"), 
                                       ".log"))

# Record start time for runtime calculation (before logging starts)
start_time <- Sys.time()

# Start logging (both to console and file)
# split=TRUE writes to both console and file simultaneously
sink(log_file, split = TRUE)

cat("╔═══════════════════════════════════════════════════════════╗\n")
cat("║  CELLCHAT ANALYSIS LOG                                    ║\n")
cat("╚═══════════════════════════════════════════════════════════╝\n")
cat("Analysis started:", format(start_time, "%Y-%m-%d %H:%M:%S"), "\n")
cat("Log file:", log_file, "\n")
cat("Sample:", sample_name, "\n")
cat("Output directory:", plot_dir, "\n\n")

# Function to ensure log is closed at the end (even on errors)
on.exit({
  if (sink.number() > 0) {
    cat("\n")
    cat("Analysis ended:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
    total_runtime <- difftime(Sys.time(), start_time, units = "mins")
    cat("Total runtime:", round(total_runtime, 2), "minutes\n")
    sink()  # Close log file
  }
}, add = TRUE)

head(data_human@meta.data, 2)

assay <- "RNA"
DefaultAssay(data_human) <- assay

cat(sprintf("Using assay: %s\n", assay))

genes <- rownames(data_human[[assay]])
is_ensembl <- grepl("^ENS[A-Z]*G\\d+(\\.\\d+)?$", genes)
keep <- genes[!is_ensembl]

cat(sprintf("Removing %d Ensembl gene IDs, keeping %d gene symbols\n", 
            sum(is_ensembl), length(keep)))

data_human <- subset(data_human, features = keep)

cat("Genes after filtering: ")
sum(grepl("^ENS[A-Z]*G\\d+(\\.\\d+)?$", rownames(data_human[[assay]])))
cat("\n")

dim(data_human)

data.input <- LayerData(data_human, assay = assay, layer = "data")
dim(data.input)

head(rownames(data.input))

head(data_human@meta.data, 2)
data_human@meta.data$labels = data_human@meta.data$celltype
head(data_human@meta.data, 2)

meta <- data_human@meta.data
print(dim(meta))
head(meta, 3)

unique(meta$labels)
dim(data.input)

################################################
################################################

sketchData <- function(object, n_cells, do.PCA = TRUE, dimPC = 30) {
    geosketch <- reticulate::import('geosketch')
    object = t(object)
    
    if (do.PCA) {
        X.pcs <- runPCA(object, dimPC = dimPC)
    } else {
        X.pcs <- object
    }

    sketch.size <- as.integer(n_cells)
    cells.all <- rownames(object)

    sketch.index <- geosketch$gs(X.pcs, sketch.size)
    sketch.index <- unlist(sketch.index) + 1
    sketch.cells <- cells.all[sketch.index]
    return(sketch.cells)
}

################################################
################################################

print("Creating CellChat object...")

cellchat <- createCellChat(
    object = data.input,
    meta = meta,
    group.by = 'labels'
)

cellchat

check_memory_usage(threshold = 90)

celltypes <- sort(unique(cellchat@idents))
n_celltypes <- length(celltypes)
cat(sprintf("\n✓ Detected %d cell types: %s\n\n", 
            n_celltypes, paste(celltypes, collapse=", ")))

################################################
# SET THE LIGAND RECEPTOR DATABASE
################################################

CellChatDB <- CellChatDB.human

names(CellChatDB)

head(CellChatDB$interaction, 1)
head(CellChatDB$complex, 1)
head(CellChatDB$cofactor, 1)
head(CellChatDB$geneInfo, 1)

options(repr.plot.width=10, repr.plot.height=8)     
showDatabaseCategory(CellChatDB)

tail(CellChatDB$interaction, 3)
tail(CellChatDB$interaction$annotation, 3)
unique(CellChatDB$interaction$annotation)

CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use

print("Subset the expression data with the genes that appear in the interaction database only.")

cellchat <- subsetData(cellchat)

dim(cellchat@data)
dim(cellchat@data.signaling)
dim(cellchat@DB$interaction)
dim(cellchat@DB$complex)

head(cellchat@DB$interaction, 2)
head(cellchat@DB$complex, 2)

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

head(cellchat@LR$LRsig, 2)

################################################
################################################

cellchat <- projectData(
    object = cellchat,
    adjMatrix = PPI.human
)

head(cellchat@data.project, 2)

cat("Inference of cell-cell communication network\n")

################################################
################################################

cat("\n=== System Diagnostics ===\n")
cat("R version:", as.character(getRversion()), "\n")
cat("CellChat version:", as.character(packageVersion("CellChat")), "\n")
cat("Number of cells:", ncol(cellchat@data), "\n")
cat("Number of genes:", nrow(cellchat@data), "\n")
cat("Number of cell types:", length(unique(cellchat@idents)), "\n")
cat("Cell types:", paste(unique(cellchat@idents), collapse=", "), "\n")
cat("\n")

cat("\n=== Computing Communication Probabilities ===\n")
cat("Using BiocParallel with", bpnworkers(BPPARAM), "workers\n")

nboot_iterations <- 100
cat("Bootstrap iterations (nboot):", nboot_iterations, "\n")
cat("(Tip: Use nboot=50 for exploration, nboot=100 for final analysis)\n")
cat("This may take 10-30 minutes depending on dataset size...\n\n")

do_gc(verbose = TRUE)
check_memory_usage(threshold = 85)

computation_success <- FALSE
attempt <- 1
max_attempts <- 3

while (!computation_success && attempt <= max_attempts) {
  
  cat(sprintf("\n--- Attempt %d/%d ---\n", attempt, max_attempts))
  
  tryCatch({
    
    if (attempt == 1) {
      cat("Using BiocParallel for parallel computation...\n")
      cat("(BiocParallel is optimized for bioinformatics workflows)\n")
      
      cellchat <- computeCommunProb(
        object = cellchat,
        type = 'triMean',
        raw.use = TRUE,
        nboot = nboot_iterations
      )
      
    } else if (attempt == 2) {
      cat("⚠️  BiocParallel failed. Trying doParallel fallback...\n")
      
      cl <- makeCluster(n_workers)
      registerDoParallel(cl)
      
      cellchat <- computeCommunProb(
        object = cellchat,
        type = 'triMean',
        raw.use = TRUE,
        nboot = nboot_iterations
      )
      
      stopCluster(cl)
      
    } else {
      cat("⚠️  Parallel processing failed. Falling back to sequential mode...\n")
      cat("This will be slower but more stable.\n")
      
      old_bpparam <- bpparam()
      register(SerialParam())
      
      cellchat <- computeCommunProb(
        object = cellchat,
        type = 'triMean',
        raw.use = TRUE,
        nboot = nboot_iterations
      )
      
      register(old_bpparam)
    }
    
    computation_success <- TRUE
    cat("\n✓ Communication probability computation successful!\n")
    
  }, error = function(e) {
    cat("\n❌ Error in attempt", attempt, ":\n")
    cat("   ", as.character(e), "\n")
    
    if (attempt < max_attempts) {
      cat("\n🔄 Cleaning up and preparing for fallback...\n")
      do_gc(verbose = TRUE)
      Sys.sleep(2)
    } else {
      cat("\n💥 All attempts failed. Cannot continue.\n")
      stop(sprintf("computeCommunProb failed after %d attempts: %s", 
                   max_attempts, as.character(e)))
    }
  })
  
  attempt <- attempt + 1
}

cellchat@options$parameter
cellchat@net
dim(cellchat@net$prob)
dim(cellchat@net$pval)

check_memory_usage(threshold = 90)

###############################################
# Save RDS after computeCommunProb
###############################################

if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
  cat("✓ Created directory:", plot_dir, "\n")
}

outfile <- file.path(
  plot_dir,
  paste0("cellchat_", sample_name, "_after_computeCommunProb.rds")
)

saveRDS(cellchat, file = outfile)

cat("✓ CellChat object saved after computeCommunProb to:", outfile, "\n")

###############################################
###############################################

cellchat <- filterCommunication(
    cellchat, 
    min.cells = 10
) 

cat("Compute communication probabilty at pathway level\n")

cellchat <- computeCommunProbPathway(
    object=cellchat,
    thresh=0.05
)

cat("Compute the aggregated network\n")

cellchat <- aggregateNet(
    object=cellchat,
    remove.isolate=TRUE,
    thresh=0.05
)

cellchat@net$count
cellchat@net$weight

cellchat@net
cellchat@netP

check_memory_usage(threshold = 90)

################################################
# Identify signaling role of cells 
################################################

cellchat <- netAnalysis_computeCentrality(cellchat)

cat("Visualizations\n")

png(paste0(plot_dir, '/signaling_role_heatmap_incoming.png'), 
    width=8, height=8, units='in', res=300)

netAnalysis_signalingRole_heatmap(
    object=cellchat,
    pattern='incoming'
)

dev.off()

options(repr.plot.width=12, repr.plot.height=12)     
netAnalysis_signalingRole_heatmap(
    object=cellchat,
    pattern='incoming'
)

png(paste0(plot_dir, '/signaling_role_heatmap_outgoing.png'), 
    width=8, height=8, units='in', res=300)
netAnalysis_signalingRole_heatmap(
    object=cellchat,
    pattern='outgoing'
)
dev.off()

options(repr.plot.width=12, repr.plot.height=12)     
netAnalysis_signalingRole_heatmap(
    object=cellchat,
    pattern='outgoing'
)

################################################
################################################

cellchat@netP$pathways

slotNames(cellchat)

cat('Visualize centrality of cell types with scatter plot\n')

options(repr.plot.width=7, repr.plot.height=7)     

netAnalysis_signalingRole_scatter(
    object=cellchat,
    x.measure='outdeg',
    y.measure='indeg',
    xlabel='Outgoing interaction strength',
    ylabel='Incoming interaction strength',
    signaling=cellchat@netP$pathways
)

outfile <- file.path(plot_dir, "signalingRole_scatter_all_pathways.png")

png(outfile, width = 7, height = 7, units = "in", res = 300)

netAnalysis_signalingRole_scatter(
  object = cellchat,
  x.measure = "outdeg",
  y.measure = "indeg",
  xlabel = "Outgoing interaction strength",
  ylabel = "Incoming interaction strength",
  signaling = cellchat@netP$pathways
)

dev.off()

cat("✓ Saved scatter plot to:", outfile, "\n")

################################################
# SCATTER PLOT - SIGNALING ROLE
################################################

options(repr.plot.width=7, repr.plot.height=7)     

p <- netAnalysis_signalingRole_scatter(
    object=cellchat,
    x.measure='outdeg',
    y.measure='indeg',
    xlabel='Outgoing interaction strength',
    ylabel='Incoming interaction strength',
    signaling=cellchat@netP$pathways
)

print(p)

outfile <- file.path(plot_dir, "signalingRole_scatter_all_pathways.png")
ggsave(
  filename = outfile,
  plot = p,
  width = 7,
  height = 7,
  units = "in",
  dpi = 300
)
cat("✓ Saved scatter plot to:", outfile, "\n")

ggsave(
  filename = file.path(plot_dir, "signalingRole_scatter_all_pathways.pdf"),
  plot = p,
  width = 7,
  height = 7,
  units = "in"
)

################################################
################################################

cellchat@netP$pathways

################################################
################################################
# SIGNALING ROLE NETWORK - ALL PATHWAYS WITH PROGRESS BAR
################################################

present = cellchat@netP$pathways

has_signal <- function(p) {
  df <- subsetCommunication(cellchat, slot.name = "netP", signaling = p, thresh = 0.05)
  nrow(df) > 0
}
usable <- Filter(has_signal, present)

cat("\n=== Signaling Role Networks for All Pathways ===\n")
cat("Will save:", paste(usable, collapse = ", "), "\n\n")

pb <- progress_bar$new(
  format = "[:bar] :current/:total (:percent) ETA: :eta | :pathway",
  total = length(usable),
  clear = FALSE,
  width = 80
)

for (p in usable) {
  pb$tick(tokens = list(pathway = sprintf("%-15s", p)))
  
  png(file.path(plot_dir, paste0(p, "_signalingRole_network.png")),
      width = 15, height = 8, units = "in", res = 300)
  ok <- try(
    netAnalysis_signalingRole_network(cellchat, signaling = p,
                                      width = 15, height = 8, font.size = 15),
    silent = TRUE
  )
  dev.off()
  
  if (!inherits(ok, "try-error")) {
    grid::grid.newpage()
    netAnalysis_signalingRole_network(cellchat, signaling = p,
                                      width = 15, height = 8, font.size = 15)
  }
}

cat("\n✓ Signaling role networks complete\n\n")

################################################
################################################
# Identify communication patterns - OUTGOING
################################################

# Run selectK
set.seed(12345)

png(file.path(plot_dir, "selectK_outgoing_patterns.png"),
    width = 14, height = 6, units = "in", res = 300)

set.seed(12345)
select_out <- selectK(cellchat, pattern='outgoing')

dev.off()

cat("✓ Saved selectK plot (outgoing) to:", 
    file.path(plot_dir, "selectK_outgoing_patterns.png"), "\n")

# Try to extract k from select_out if valid
selected_k_out <- NULL

if (!is.null(select_out) && !is.null(select_out$data) && nrow(select_out$data) > 0) {
  # Save the selection data to file
  write.csv(select_out$data, 
            file.path(plot_dir, "selectK_outgoing_scores.csv"),
            row.names = FALSE)
  cat("✓ Saved selectK scores to: selectK_outgoing_scores.csv\n")
  
  # Extract k based on max Cophenetic
  cophenetic_data <- select_out$data[select_out$data$Measure == "Cophenetic", ]
  cophenetic_filtered <- cophenetic_data[cophenetic_data$k > 1, ]
  
  if (nrow(cophenetic_filtered) > 0) {
    selected_k_out <- cophenetic_filtered$k[which.max(cophenetic_filtered$score)]
    cat("Auto-selected k for outgoing:", selected_k_out, 
        "(Cophenetic =", max(cophenetic_filtered$score), ")\n")
  }
}

# Loop through k from 1 to 7
cat("\n=== Running identifyCommunicationPatterns for k = 1 to 7 (OUTGOING) ===\n")

for (k_val in 1:7) {
  cat(sprintf("\nProcessing k = %d...\n", k_val))
  
  tryCatch({
    cellchat <- identifyCommunicationPatterns(
      object = cellchat,
      pattern = 'outgoing',
      k = k_val
    )
    
    # Save this k version
    outfile <- file.path(plot_dir, paste0("cellchat_outgoing_k", k_val, ".rds"))
    saveRDS(cellchat, file = outfile)
    cat("✓ Saved cellchat with k =", k_val, "\n")
    
    # Generate plots
    tryCatch({
      png(file.path(plot_dir, paste0("netAnalysis_river_outgoing_k", k_val, ".png")),
          width = 10, height = 8, units = "in", res = 300)
      netAnalysis_river(object = cellchat, pattern = "outgoing")
      dev.off()
      cat("  ✓ River plot for k =", k_val, "\n")
    }, error = function(e) {
      safe_dev_off()
    })
    
    tryCatch({
      png(file.path(plot_dir, paste0("netAnalysis_dot_outgoing_k", k_val, ".png")),
          width = 14, height = 6, units = "in", res = 300)
      netAnalysis_dot(object = cellchat, pattern = "outgoing", dot.size = c(1, 12))
      dev.off()
      cat("  ✓ Dot plot for k =", k_val, "\n")
    }, error = function(e) {
      safe_dev_off()
    })
    
  }, error = function(e) {
    cat("❌ Failed for k =", k_val, ":", as.character(e), "\n")
  })
  
  do_gc(verbose = FALSE)
}

cat("\n✓ Completed outgoing patterns for k = 1 to 7\n")
if (!is.null(selected_k_out)) {
  cat("★ Recommended k based on Cophenetic:", selected_k_out, "\n")
}

################################################
################################################

################################################
# Identify communication patterns - INCOMING
################################################

# Run selectK
set.seed(12345)

png(file.path(plot_dir, "selectK_incoming_patterns.png"),
    width = 14, height = 6, units = "in", res = 300)

set.seed(12345)
select_in <- selectK(cellchat, pattern='incoming')

dev.off()

cat("✓ Saved selectK plot (incoming) to:", 
    file.path(plot_dir, "selectK_incoming_patterns.png"), "\n")

# Try to extract k from select_in if valid
selected_k_in <- NULL

if (!is.null(select_in) && !is.null(select_in$data) && nrow(select_in$data) > 0) {
  # Save the selection data to file
  write.csv(select_in$data, 
            file.path(plot_dir, "selectK_incoming_scores.csv"),
            row.names = FALSE)
  cat("✓ Saved selectK scores to: selectK_incoming_scores.csv\n")
  
  # Extract k based on max Cophenetic
  cophenetic_data <- select_in$data[select_in$data$Measure == "Cophenetic", ]
  cophenetic_filtered <- cophenetic_data[cophenetic_data$k > 1, ]
  
  if (nrow(cophenetic_filtered) > 0) {
    selected_k_in <- cophenetic_filtered$k[which.max(cophenetic_filtered$score)]
    cat("Auto-selected k for incoming:", selected_k_in, 
        "(Cophenetic =", max(cophenetic_filtered$score), ")\n")
  }
}

# Loop through k from 1 to 7
cat("\n=== Running identifyCommunicationPatterns for k = 1 to 7 (INCOMING) ===\n")

for (k_val in 1:7) {
  cat(sprintf("\nProcessing k = %d...\n", k_val))
  
  tryCatch({
    cellchat <- identifyCommunicationPatterns(
      object = cellchat,
      pattern = 'incoming',
      k = k_val
    )
    
    # Save this k version
    outfile <- file.path(plot_dir, paste0("cellchat_incoming_k", k_val, ".rds"))
    saveRDS(cellchat, file = outfile)
    cat("✓ Saved cellchat with k =", k_val, "\n")
    
    # Generate plots
    tryCatch({
      png(file.path(plot_dir, paste0("netAnalysis_river_incoming_k", k_val, ".png")),
          width = 10, height = 8, units = "in", res = 300)
      netAnalysis_river(object = cellchat, pattern = "incoming")
      dev.off()
      cat("  ✓ River plot for k =", k_val, "\n")
    }, error = function(e) {
      safe_dev_off()
    })
    
    tryCatch({
      png(file.path(plot_dir, paste0("netAnalysis_dot_incoming_k", k_val, ".png")),
          width = 10, height = 6, units = "in", res = 300)
      netAnalysis_dot(object = cellchat, pattern = "incoming", dot.size = c(1, 12))
      dev.off()
      cat("  ✓ Dot plot for k =", k_val, "\n")
    }, error = function(e) {
      safe_dev_off()
    })
    
  }, error = function(e) {
    cat("❌ Failed for k =", k_val, ":", as.character(e), "\n")
  })
  
  do_gc(verbose = FALSE)
}

cat("\n✓ Completed incoming patterns for k = 1 to 7\n")
if (!is.null(selected_k_in)) {
  cat("★ Recommended k based on Cophenetic:", selected_k_in, "\n")
}

################################################
################################################

print("Compute the functional similarity")

cellchat <- computeNetSimilarity(cellchat, type='functional')

if (!is.null(cellchat@netP$similarity$functional$matrix)) {
  cat("✓ Functional similarity computed successfully\n")
} else {
  cat("⚠️  Functional similarity matrix is NULL\n")
}

library(reticulate)

py_config()

cat("\n=== Functional Similarity Embedding (UMAP) ===\n")

umap_available <- FALSE
umap_install_attempted <- FALSE

tryCatch({
  umap_module <- reticulate::import("umap")
  umap_available <- TRUE
  cat("✓ umap-learn is already installed\n")
}, error = function(e) {
  cat("⚠️  umap-learn not found. Attempting installation...\n")
  umap_install_attempted <- TRUE
  
  tryCatch({
    reticulate::py_install("umap-learn", pip = TRUE)
    Sys.sleep(2)
    umap_module <- reticulate::import("umap")
    umap_available <- TRUE
    cat("✓ umap-learn installed and verified successfully\n")
  }, error = function(e2) {
    cat("❌ Failed to install umap-learn automatically.\n")
    cat("   Error:", as.character(e2), "\n")
    cat("   Please install manually:\n")
    cat("   reticulate::py_install('umap-learn', pip = TRUE)\n")
    cat("   OR: ~/.virtualenvs/r-reticulate/bin/pip install umap-learn\n")
    cat("   Skipping embedding step - continuing with rest of analysis...\n")
    umap_available <<- FALSE
  })
})

if (umap_available) {
  tryCatch({
    cat("Attempting functional similarity embedding...\n")
    cellchat <- netEmbedding(cellchat, type ='functional')
    cat("✓ Functional embedding completed\n")
    
    tryCatch({
      cellchat <- netClustering(cellchat, type='functional', do.parallel = FALSE)
      cat("✓ Functional clustering completed\n")
    }, error = function(e) {
      cat("⚠️  Functional clustering failed:", as.character(e), "\n")
      cat("   Continuing without clustering results...\n")
    })
    
    tryCatch({
      options(repr.plot.width=7, repr.plot.height=6)     
      netVisual_embedding(
          object=cellchat,
          type='functional',
          label.size=3.5,
          dot.size=c(5,15)
      )
      
      png(file.path(plot_dir, "netVisual_embedding_functional.png"),
          width = 7, height = 6, units = "in", res = 300)
      
      netVisual_embedding(
        object = cellchat,
        type = "functional",
        label.size = 3.5,
        dot.size = c(5, 15)
      )
      
      dev.off()
      
      cat("✓ Functional embedding plot saved to:",
          file.path(plot_dir, "netVisual_embedding_functional.png"), "\n")
    }, error = function(e) {
      cat("⚠️  Failed to create functional embedding plot:", as.character(e), "\n")
      safe_dev_off()
    })
    
  }, error = function(e) {
    cat("❌ Functional embedding failed:", as.character(e), "\n")
    cat("   This may be due to UMAP issues or insufficient data.\n")
    cat("   Continuing with rest of analysis...\n")
  })
} else {
  cat("⚠️  Skipping functional embedding (UMAP not available)\n")
  cat("   The rest of the analysis will continue normally.\n")
}

print("Compute the structural similarity")

cellchat <- computeNetSimilarity(cellchat, type='structural')

if (!is.null(cellchat@netP$similarity$structural$matrix)) {
  cat("✓ Structural similarity computed successfully\n")
} else {
  cat("⚠️  Structural similarity matrix is NULL\n")
}

cat("\n=== Structural Similarity Embedding (UMAP) ===\n")

structural_umap_available <- FALSE
tryCatch({
  umap_module <- reticulate::import("umap")
  structural_umap_available <- TRUE
  cat("✓ umap-learn is available for structural embedding\n")
}, error = function(e) {
  cat("⚠️  umap-learn not available. Attempting installation...\n")
  tryCatch({
    reticulate::py_install("umap-learn", pip = TRUE)
    Sys.sleep(2)
    umap_module <- reticulate::import("umap")
    structural_umap_available <<- TRUE
    cat("✓ umap-learn installed and verified\n")
  }, error = function(e2) {
    cat("❌ Failed to install umap-learn:", as.character(e2), "\n")
    cat("   Skipping structural embedding...\n")
    structural_umap_available <<- FALSE
  })
})

if (structural_umap_available) {
  tryCatch({
    cat("Attempting structural similarity embedding...\n")
    cellchat <- netEmbedding(cellchat, type='structural')
    cat("✓ Structural embedding completed\n")
    
    tryCatch({
      cellchat <- netClustering(cellchat, type='structural', do.parallel = FALSE)
      cat("✓ Structural clustering completed\n")
    }, error = function(e) {
      cat("⚠️  Structural clustering failed:", as.character(e), "\n")
      cat("   Continuing without clustering results...\n")
    })
    
    tryCatch({
      options(repr.plot.width=7, repr.plot.height=6)     
      netVisual_embedding(
          object=cellchat,
          type='structural',
          label.size=3.5,
          dot.size=c(5,10)
      )
      
      png(file.path(plot_dir, "netVisual_embedding_structural.png"),
          width = 7, height = 6, units = "in", res = 300)
      
      netVisual_embedding(
        object = cellchat,
        type = "structural",
        label.size = 3.5,
        dot.size = c(5, 10)
      )
      
      dev.off()
      
      cat("✓ Saved structural embedding to:",
          file.path(plot_dir, "netVisual_embedding_structural.png"), "\n")
    }, error = function(e) {
      cat("⚠️  Failed to create structural embedding plot:", as.character(e), "\n")
      safe_dev_off()
    })
    
  }, error = function(e) {
    cat("❌ Structural embedding failed:", as.character(e), "\n")
    cat("   This may be due to UMAP issues or insufficient data.\n")
    cat("   Continuing with rest of analysis...\n")
  })
} else {
  cat("⚠️  Skipping structural embedding (UMAP not available)\n")
  cat("   The rest of the analysis will continue normally.\n")
}

################################################
################################################

data.frame(
    celltype=unique(cellchat@idents),
    index=as.numeric(unique(cellchat@idents))
)

cellchat@netP$pathways

writeLines(
  cellchat@netP$pathways, 
  con = file.path(plot_dir, "active_pathways.txt")
)

cat("✓ Pathways saved to:", 
    file.path(plot_dir, "active_pathways.txt"), "\n")

##########################################################################################################
# Adding plots in ggplot2 format
##########################################################################################################

################################################
# RIVER PLOT - OUTGOING PATTERNS (GGPLOT2)
################################################

options(repr.plot.width=10, repr.plot.height=8)

p_river <- netAnalysis_river(cellchat, pattern='outgoing')
print(p_river)

outfile <- file.path(plot_dir, "netAnalysis_river_outgoing.png")

ggsave(
  filename = outfile,
  plot = p_river,
  width = 10,
  height = 8,
  units = "in",
  dpi = 300,
  bg = "white"
)

cat("✓ Saved river plot (outgoing) to:", outfile, "\n")

################################################
# DOT PLOT - OUTGOING PATTERNS (GGPLOT2)
################################################

options(repr.plot.width=14, repr.plot.height=6)

p_dot <- netAnalysis_dot(cellchat, pattern='outgoing', dot.size=c(1,12))
print(p_dot)

outfile <- file.path(plot_dir, "netAnalysis_dot_outgoing.png")

ggsave(
  filename = outfile,
  plot = p_dot,
  width = 14,
  height = 6,
  units = "in",
  dpi = 300,
  bg = "white"
)

cat("✓ Saved dot plot (outgoing) to:", outfile, "\n")

##########################################################################################################

################################################
# RIVER PLOT - INCOMING PATTERNS (GGPLOT2)
################################################
options(repr.plot.width=10, repr.plot.height=8)

p_river <- netAnalysis_river(cellchat, pattern='incoming')
print(p_river)

outfile <- file.path(plot_dir, "netAnalysis_river_incoming.png")

ggsave(
  filename = outfile,
  plot = p_river,
  width = 10,
  height = 8,
  units = "in",
  dpi = 300,
  bg = "white"
)

cat("✓ Saved river plot (incoming) to:", outfile, "\n")

################################################
# DOT PLOT - INCOMING PATTERNS (GGPLOT2)
################################################
options(repr.plot.width=10, repr.plot.height=6)

p_dot <- netAnalysis_dot(cellchat, pattern='incoming', dot.size=c(1,12))
print(p_dot)

outfile <- file.path(plot_dir, "netAnalysis_dot_incoming.png")

ggsave(
  filename = outfile,
  plot = p_dot,
  width = 10,
  height = 6,
  units = "in",
  dpi = 300,
  bg = "white"
)

cat("✓ Saved dot plot (incoming) to:", outfile, "\n")

##########################################################################################################

################################################
# FUNCTIONAL EMBEDDING PLOT (GGPLOT2)
################################################
cat("\n=== Creating ggplot2 functional embedding plot ===\n")

has_functional_embedding <- tryCatch({
  !is.null(cellchat@net$similarity$functional) || 
  !is.null(cellchat@netP$similarity$functional)
}, error = function(e) FALSE)

if (has_functional_embedding) {
  tryCatch({
    options(repr.plot.width=7, repr.plot.height=6)
    
    p_embed_func <- netVisual_embedding(
      object = cellchat,
      type = "functional",
      label.size = 3.5,
      dot.size = c(5, 15)
    )
    print(p_embed_func)
    
    outfile <- file.path(plot_dir, "netVisual_embedding_functional.png")
    
    ggsave(
      filename = outfile,
      plot = p_embed_func,
      width = 7,
      height = 6,
      units = "in",
      dpi = 300,
      bg = "white"
    )
    
    cat("✓ Functional embedding plot saved to:", outfile, "\n")
  }, error = function(e) {
    cat("⚠️  Failed to create functional embedding plot:", as.character(e), "\n")
    cat("   This may be because embedding was not computed successfully.\n")
    cat("   Skipping this plot...\n")
  })
} else {
  cat("⚠️  Skipping functional embedding plot (embedding not available)\n")
}

################################################
# STRUCTURAL EMBEDDING PLOT (GGPLOT2)
################################################
cat("\n=== Creating ggplot2 structural embedding plot ===\n")

has_structural_embedding <- tryCatch({
  !is.null(cellchat@net$similarity$structural) || 
  !is.null(cellchat@netP$similarity$structural)
}, error = function(e) FALSE)

if (has_structural_embedding) {
  tryCatch({
    options(repr.plot.width=7, repr.plot.height=6)
    
    p_embed_struct <- netVisual_embedding(
      object = cellchat,
      type = "structural",
      label.size = 3.5,
      dot.size = c(5, 10)
    )
    print(p_embed_struct)
    
    outfile <- file.path(plot_dir, "netVisual_embedding_structural.png")
    
    ggsave(
      filename = outfile,
      plot = p_embed_struct,
      width = 7,
      height = 6,
      units = "in",
      dpi = 300,
      bg = "white"
    )
    
    cat("✓ Saved structural embedding to:", outfile, "\n")
  }, error = function(e) {
    cat("⚠️  Failed to create structural embedding plot:", as.character(e), "\n")
    cat("   This may be because embedding was not computed successfully.\n")
    cat("   Skipping this plot...\n")
  })
} else {
  cat("⚠️  Skipping structural embedding plot (embedding not available)\n")
}

################################################
# BUBBLE PLOT - SAVE AS PDF (RECOMMENDED)
################################################

options(repr.plot.width=40, repr.plot.height=100)

p_bubble <- netVisual_bubble(object = cellchat)
print(p_bubble)

outfile_pdf <- file.path(plot_dir, "netVisual_bubble_all.pdf")

if (is_gg(p_bubble)) {
  ggsave(
    filename = outfile_pdf,
    plot = p_bubble,
    width = 40,
    height = 100,
    units = "in",
    limitsize = FALSE
  )
  cat("✓ Bubble plot PDF saved (full resolution, vector) to:", outfile_pdf, "\n")
  
  ggsave(
    filename = file.path(plot_dir, "netVisual_bubble_all.png"),
    plot = p_bubble,
    width = 16,
    height = 40,
    units = "in",
    dpi = 150,
    bg = "white"
  )
  cat("✓ Bubble plot PNG saved (preview version, 16x40 in @ 150 dpi)\n")
  cat("   (Use PDF for publications - vector format, infinite zoom)\n")
} else {
  cat("⚠️  Bubble plot is not ggplot - skipping ggsave\n")
}

##########################################################################################################
# HEATMAPS FOR ALL PATHWAYS WITH PROGRESS BAR
##########################################################################################################

library(ComplexHeatmap)

present = cellchat@netP$pathways

has_signal <- function(p) {
  df <- subsetCommunication(cellchat, slot.name = "netP", signaling = p, thresh = 0.05)
  nrow(df) > 0
}
usable <- Filter(has_signal, present)

cat("\n=== Heatmaps for All Pathways ===\n")
cat("Saving heatmaps for", length(usable), "pathways\n\n")

pb <- progress_bar$new(
  format = "[:bar] :current/:total (:percent) ETA: :eta | :pathway",
  total = length(usable),
  clear = FALSE,
  width = 80
)

for (p in usable) {
  pb$tick(tokens = list(pathway = sprintf("%-15s", p)))
  
  tryCatch({
    ht <- netVisual_heatmap(
      object = cellchat,
      signaling = c(p),
      color.heatmap = "Reds"
    )
    
    png(file.path(plot_dir, paste0(p, "_netVisual_heatmap.png")),
        width = 12, height = 8, units = "in", res = 300, type = "cairo")
    print(draw(ht))
    invisible(safe_dev_off())
    
    pdf(file.path(plot_dir, paste0(p, "_netVisual_heatmap.pdf")),
        width = 12, height = 8)
    print(draw(ht))
    invisible(safe_dev_off())
    
  }, error = function(e) {
    invisible(dev.off())
  })
}

cat("\n✓ Heatmap generation complete\n\n")

##########################################################################################################
# HEATMAPS + GENE EXPRESSION FOR ALL PATHWAYS WITH PROGRESS BAR
##########################################################################################################

cat("\n=== Heatmaps + Gene Expression for All Pathways ===\n")
cat("Saving heatmaps and gene expression plots for", length(usable), "pathways\n\n")

success_heatmap <- 0
success_geneexpr <- 0
failed_pathways <- c()

pb <- progress_bar$new(
  format = "[:bar] :current/:total (:percent) ETA: :eta | :pathway",
  total = length(usable),
  clear = FALSE,
  width = 80
)

for (i in seq_along(usable)) {
  p <- usable[i]
  pb$tick(tokens = list(pathway = sprintf("%-15s", p)))
  
  # ===== 1. HEATMAP =====
  tryCatch({
    ht <- netVisual_heatmap(
      object = cellchat,
      signaling = c(p),
      color.heatmap = "Reds"
    )
    
    png(file.path(plot_dir, paste0(p, "_netVisual_heatmap.png")),
        width = 12, height = 8, units = "in", res = 300, type = "cairo")
    print(draw(ht))
    invisible(dev.off())
    
    pdf(file.path(plot_dir, paste0(p, "_netVisual_heatmap.pdf")),
        width = 12, height = 8)
    print(draw(ht))
    invisible(dev.off())
    
    success_heatmap <- success_heatmap + 1
    
  }, error = function(e) {
    invisible(safe_dev_off())
  })
  
  # ===== 2. GENE EXPRESSION =====
  ok <- try({
    p_expr <- plotGeneExpression(
      object = cellchat,
      signaling = p
    )
    
    outfile <- file.path(plot_dir, paste0(p, "_plotGeneExpression.png"))
    
    ggsave(
      filename = outfile,
      plot = p_expr,
      width = 12,
      height = 8,
      units = "in",
      dpi = 300,
      bg = "white"
    )
    
    success_geneexpr <- success_geneexpr + 1
    
  }, silent = TRUE)
  
  if (inherits(ok, "try-error")) {
    if (success_heatmap != i) {
      failed_pathways <- c(failed_pathways, p)
    }
  }
}

cat("\n╔════════════════════════════════════════════╗\n")
cat("║           SUMMARY REPORT                  ║\n")
cat("╚════════════════════════════════════════════╝\n")
cat(sprintf("Total pathways: %d\n", length(usable)))
cat(sprintf("✓ Heatmaps: %d/%d\n", success_heatmap, length(usable)))
cat(sprintf("✓ Gene expression: %d/%d\n", success_geneexpr, length(usable)))

if (length(failed_pathways) > 0) {
  cat(sprintf("\n⚠️  Completely failed: %d pathways\n", length(failed_pathways)))
  cat("   ", paste(failed_pathways, collapse = ", "), "\n")
}
cat("\n✓ All plots saved to:", plot_dir, "\n\n")

##########################################################################################################
# CHORD DIAGRAMS FOR ALL PATHWAYS WITH PROGRESS BAR
##########################################################################################################

library(circlize)

cat("\n=== Chord Diagrams for All Pathways ===\n")
cat("Generating chord diagrams for", length(usable), "pathways\n\n")

pb <- progress_bar$new(
  format = "[:bar] :current/:total (:percent) ETA: :eta | :pathway",
  total = length(usable),
  clear = FALSE,
  width = 80
)

for (i in seq_along(usable)) {
  p <- usable[i]
  pb$tick(tokens = list(pathway = sprintf("%-15s", p)))
  
  if (i %% 10 == 0) {
    do_gc(verbose = FALSE)
    check_memory_usage(threshold = 90)
  }
  
  circos.clear()
  outfile <- file.path(plot_dir, paste0(p, "_netVisual_chord_gene.png"))
  
  png(outfile, width = 20, height = 20, units = "in", res = 300)
  ok <- try({
    circos.clear()
    circos.par(gap.degree = 1, track.margin = c(0.005, 0.005))
    netVisual_chord_gene(cellchat, sources.use = 1:n_celltypes, 
                         targets.use = 1:n_celltypes, signaling = c(p))
  }, silent = TRUE)
  dev.off()
  circos.clear()
}

cat("\n✓ Chord diagrams complete\n\n")

##########################################################################################################
# CONTRIBUTION PLOTS FOR ALL PATHWAYS WITH PROGRESS BAR
##########################################################################################################

cat("\n=== Contribution Plots for All Pathways ===\n")
cat("Generating contribution plots for", length(usable), "pathways\n\n")

pb <- progress_bar$new(
  format = "[:bar] :current/:total (:percent) ETA: :eta | :pathway",
  total = length(usable),
  clear = FALSE,
  width = 80
)

for (i in seq_along(usable)) {
  p <- usable[i]
  pb$tick(tokens = list(pathway = sprintf("%-15s", p)))
  
  if (i %% 10 == 0) {
    do_gc(verbose = FALSE)
    check_memory_usage(threshold = 90)
  }
  
  ok <- try({
    p_contrib <- netAnalysis_contribution(cellchat, signaling = p)
    
    if (is_gg(p_contrib)) {
      ggsave(
        file.path(plot_dir, paste0(p, "_netAnalysis_contribution.png")),
        plot = p_contrib,
        width = 6, height = 4, dpi = 300, bg = "white"
      )
    }
  }, silent = TRUE)
}

cat("\n✓ Contribution plots complete\n\n")

##########################################################################################################
# SIGNALING ROLE ANALYSIS - COMBINED HEATMAP
##########################################################################################################

ht1 <- netAnalysis_signalingRole_heatmap(
  cellchat, 
  pattern = "outgoing",
  font.size = 4,
  font.size.title = 10
)

ht2 <- netAnalysis_signalingRole_heatmap(
  cellchat, 
  pattern = "incoming",
  font.size = 4,
  font.size.title = 10
)

ht_combined <- ht1 + ht2

print(ht_combined)

outfile_png <- file.path(plot_dir, paste0(sample_name, "_signalingRole_heatmap_combined.png"))

png(outfile_png, width = 16, height = 16, units = "in", res = 300, type = "cairo")
print(ht_combined)
dev.off()

cat("✓ Combined heatmap (PNG) saved to:", outfile_png, "\n")

outfile_pdf <- file.path(plot_dir, paste0(sample_name, "_signalingRole_heatmap_combined.pdf"))

pdf(outfile_pdf, width = 16, height = 16)
print(ht_combined)
dev.off()

cat("✓ Combined heatmap (PDF) saved to:", outfile_pdf, "\n")

##########################################################################################################
# SAVE FINAL RDS
##########################################################################################################

outfile <- file.path(
  plot_dir,
  paste0("cellchat_", sample_name, "_at_the_end.rds")
)

saveRDS(cellchat, file = outfile)

cat("✓ CellChat object saved at the very end:", outfile, "\n")

##########################################################################################################
# SAVE EVERYTHING - COMPLETE EXPORT
##########################################################################################################

cat("\n=== EXPORTING ALL CELLCHAT DATA ===\n\n")

write.csv(
  cellchat@LR$LRsig,
  file.path(plot_dir, paste0(sample_name, "_01_all_LR_pairs.csv")),
  row.names = FALSE
)

all_lr_comm <- subsetCommunication(cellchat, slot.name = "net", thresh = 0.05)
write.csv(
  all_lr_comm,
  file.path(plot_dir, paste0(sample_name, "_02_all_LR_communications.csv")),
  row.names = FALSE
)

all_pathway_comm <- subsetCommunication(cellchat, slot.name = "netP", thresh = 0.05)
write.csv(
  all_pathway_comm,
  file.path(plot_dir, paste0(sample_name, "_03_all_pathway_communications.csv")),
  row.names = FALSE
)

write.csv(
  cellchat@net$prob[,,1],
  file.path(plot_dir, paste0(sample_name, "_04_communication_probability_matrix.csv"))
)

write.csv(
  cellchat@net$count,
  file.path(plot_dir, paste0(sample_name, "_05_communication_count_matrix.csv"))
)

write.csv(
  cellchat@net$weight,
  file.path(plot_dir, paste0(sample_name, "_06_communication_weight_matrix.csv"))
)

writeLines(
  cellchat@netP$pathways,
  file.path(plot_dir, paste0(sample_name, "_07_active_pathways.txt"))
)

cat("\n=== Exporting per-pathway interactions ===\n")

pb <- progress_bar$new(
  format = "[:bar] :current/:total (:percent) | :pathway",
  total = length(cellchat@netP$pathways),
  clear = FALSE,
  width = 80
)

for (i in seq_along(cellchat@netP$pathways)) {
  pathway <- cellchat@netP$pathways[i]
  pb$tick(tokens = list(pathway = sprintf("%-15s", pathway)))
  
  if (i %% 20 == 0) {
    do_gc(verbose = FALSE)
    check_memory_usage(threshold = 90)
  }
  
  pw_data <- subsetCommunication(cellchat, signaling = pathway, thresh = 0.05)
  if (nrow(pw_data) > 0) {
    write.csv(
      pw_data,
      file.path(plot_dir, paste0(sample_name, "_", pathway, "_pathway_interactions.csv")),
      row.names = FALSE
    )
  }
}

cat("\n--- Export Summary ---\n")
cat("L-R pairs:", nrow(cellchat@LR$LRsig), "\n")
cat("Total communications:", nrow(all_lr_comm), "\n")
cat("Active pathways:", length(cellchat@netP$pathways), "\n")
cat("\n✓ All data exported to:", plot_dir, "\n\n")

do_gc(verbose = TRUE)
check_memory_usage(threshold = 90)

cat("\n")
cat("╔═══════════════════════════════════════════════════════════╗\n")
cat("║                                                           ║\n")
cat("║   ✓✓✓  CELLCHAT ANALYSIS COMPLETE  ✓✓✓                  ║\n")
cat("║                                                           ║\n")
cat("╚═══════════════════════════════════════════════════════════╝\n")
cat("\n")

# Calculate and log total runtime
end_time <- Sys.time()
total_runtime <- difftime(end_time, start_time, units = "mins")
cat("Total analysis time:", round(total_runtime, 2), "minutes\n")
cat("Analysis completed:", format(end_time, "%Y-%m-%d %H:%M:%S"), "\n")

# Close log file (on.exit will also handle this, but explicit is good)
sink()

cat("✓ Log file saved to:", log_file, "\n")
cat("  You can review the complete analysis log in:", log_file, "\n\n")

##########################################################################################################
