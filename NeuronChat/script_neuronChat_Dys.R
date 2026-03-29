# SINGLE-SAMPLE ANALYSIS — Dyslexia
# Cell type column : subcelltype (28 cell types)

library(uwot)
library(NeuronChat)
library(Seurat)
library(future)
library(Matrix)
options(stringsAsFactors = FALSE)

plan("sequential")

# ── Paths ─────────────────────────────────────────────────────────
base_dir <- "/mnt/nfs/CX000008_DS1/projects/jaeyeon/fastq_file_Dyslexia_r1/zanalysis_bogdan/samples.merged_AG_Harmony_res0.x_anno_031826"

rds_seurat <- file.path(base_dir,
			  "18_merged_AG_031826.28celltypes.Dyslexia.rds")

out_dir <- file.path(base_dir,
		       "18_merged_AG_031826.28celltypes.Dyslexia_neuronchat")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

tmpdir <- "/mnt/nfs/CX000008_DS1/projects/btanasa/tmp"
Sys.setenv(TMPDIR = tmpdir)
options(future.globals.maxSize = 200 * 1024^3)

# ── Logging ───────────────────────────────────────────────────────
log_file   <- file.path(out_dir,
			  paste0("neuronchat_Dyslexia_28celltypes_",
				          format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
start_time <- Sys.time()
sink(log_file, split = TRUE)

cat("╔═══════════════════════════════════════════════════════════╗\n")
cat("║  NeuronChat Single-Sample: Dyslexia (subcelltype / 28ct)  ║\n")
cat("╚═══════════════════════════════════════════════════════════╝\n")
cat("Started:", format(start_time), "\n\n")

on.exit({
	  if (sink.number() > 0) {
		      cat("\nEnded  :", format(Sys.time()), "\n")
	      cat("Runtime:", round(difftime(Sys.time(), start_time,
					                                        units = "mins"), 2), "min\n")
	          sink()
	        }
}, add = TRUE)

# ── Memory check ──────────────────────────────────────────────────
check_memory <- function(threshold = 0.90) {
	  meminfo <- tryCatch(readLines("/proc/meminfo", n = 20),
			                            error = function(e) NULL)
  if (is.null(meminfo)) return(invisible(NULL))
    total <- as.numeric(gsub("[^0-9]", "", meminfo[1]))
    avail <- as.numeric(gsub("[^0-9]", "", meminfo[3]))
      pct   <- (total - avail) / total
      cat(sprintf("💾 Memory: %.1f%% used (%.1f / %.1f GB)\n",
		                pct * 100,
				              (total - avail) / 1024^2,
				              total / 1024^2))
        if (pct > threshold) {
		    cat("❌ Memory threshold exceeded — stopping\n")
	    sink(); stop("Memory threshold exceeded")
	      }
        invisible(pct)
}

# ── Load Seurat object ────────────────────────────────────────────
cat("=== Loading Seurat Object ===\n")
cat("File:", rds_seurat, "\n")
seurat_obj <- readRDS(rds_seurat)
cat("✓ Loaded\n")
cat(sprintf("Cells    : %d\n", ncol(seurat_obj)))
cat(sprintf("Genes    : %d\n", nrow(seurat_obj)))
cat(sprintf("Assays   : %s\n", paste(names(seurat_obj@assays), collapse = ", ")))

# ── Verify subcelltype column ────────────────────────────────────
cat("\n=== Checking subcelltype column ===\n")
if (!"subcelltype" %in% colnames(seurat_obj@meta.data)) {
	  cat("❌ 'subcelltype' column NOT found in meta.data\n")
  cat("Available columns:\n")
    print(colnames(seurat_obj@meta.data))
    sink(); stop("subcelltype column missing")
}

celltypes <- sort(unique(seurat_obj@meta.data$subcelltype))
cat(sprintf("✓ subcelltype found: %d cell types\n", length(celltypes)))
print(celltypes)
print(table(seurat_obj@meta.data$subcelltype))

# ── Convert Assay5 → Assay if needed (Seurat v5) ─────────────────
cat("\n=== Checking Assay version ===\n")
DefaultAssay(seurat_obj) <- "RNA"

if (inherits(seurat_obj[["RNA"]], "Assay5")) {
	  cat("Seurat v5 Assay5 detected — converting to Assay...\n")
  seurat_obj[["RNA"]] <- as(seurat_obj[["RNA"]], "Assay")
    cat("✓ Converted to Assay\n")
} else {
	  cat("✓ Assay type OK\n")
}

# ── Extract count matrix ──────────────────────────────────────────
cat("\n=== Extracting Count Matrix ===\n")

# Try counts layer first, fall back to data layer
count_mat <- tryCatch({
	  m <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
	    if (is.null(m) || nrow(m) == 0) stop("counts layer empty")
	    cat("✓ Extracted from 'counts' layer\n")
	      m
}, error = function(e) {
	  cat("⚠️  counts layer failed:", as.character(e), "\n")
  cat("   Trying 'data' layer...\n")
    m <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")
    cat("✓ Extracted from 'data' layer\n")
      m
})
cat(sprintf("Matrix: %d genes × %d cells\n", nrow(count_mat), ncol(count_mat)))

# ── Set identity to subcelltype ───────────────────────────────────
Idents(seurat_obj) <- "subcelltype"
group_info <- Idents(seurat_obj)
cat(sprintf("✓ Identity set to subcelltype (%d levels)\n",
	                length(levels(group_info))))

check_memory()

# ── Create NeuronChat object ──────────────────────────────────────
cat("\n=== Creating NeuronChat Object ===\n")
nc_obj <- createNeuronChat(
			     object   = count_mat,
			       DB       = "human",
			       group.by = as.character(group_info))

cat("✓ NeuronChat object created\n")
cat(sprintf("Pathways in DB: %d\n", length(nc_obj@DB)))

# ── Run NeuronChat ────────────────────────────────────────────────
cat("\n=== Running NeuronChat (M=100 permutations) ===\n")
cat("This will take ~30-60 minutes for 28 cell types...\n")

check_memory()

nc_obj <- run_NeuronChat(
			   object = nc_obj,
			     M      = 100)

cat("✓ NeuronChat run complete\n")
cat(sprintf("Inferred pathways: %d\n", length(nc_obj@net)))

check_memory()

# ── Pattern analysis ──────────────────────────────────────────────
cat("\n=== Pattern Analysis (k=4) ===\n")
tryCatch({
	  nc_obj <- identifyCommunicationPatterns_Neuron(
						      nc_obj,
						          slot.name  = "net",
						          pattern    = "outgoing",
							      k          = 4,
							      height     = 18)
	    cat("✓ Outgoing patterns done\n")
}, error = function(e) {
	  cat("⚠️  Outgoing patterns failed:", as.character(e), "\n")
})

tryCatch({
	  nc_obj <- identifyCommunicationPatterns_Neuron(
						      nc_obj,
						          slot.name  = "net",
						          pattern    = "incoming",
							      k          = 4,
							      height     = 18)
	    cat("✓ Incoming patterns done\n")
}, error = function(e) {
	  cat("⚠️  Incoming patterns failed:", as.character(e), "\n")
})

check_memory()

# ── Manifold learning ─────────────────────────────────────────────
cat("\n=== Manifold Learning ===\n")
tryCatch({
	  nc_obj <- computeNetSimilarity_Neuron(nc_obj,
					                                   slot.name = "net",
									                                     type      = "functional")
	    nc_obj <- netEmbedding(nc_obj,
				                             slot.name  = "net_analysis",
							                               type       = "functional",
							                               umap.method = "uwot")
	    nc_obj <- netClustering(nc_obj,
				                               slot.name = "net_analysis",
							                                  type      = "functional",
							                                  k         = 5)
	      cat("✓ Manifold learning done\n")
}, error = function(e) {
	  cat("⚠️  Manifold learning failed:", as.character(e), "\n")
})

check_memory()

# ── Save RDS ──────────────────────────────────────────────────────
cat("\n=== Saving RDS ===\n")
rds_out <- file.path(out_dir,
		       "neuronchat_18_merged_AG_031826.28celltypes.Dyslexia_at_the_end.rds")

saveRDS(nc_obj, file = rds_out)
cat("✓ Saved:", rds_out, "\n")

# ── Verify saved object cell types ───────────────────────────────
cat("\n=== Verifying saved object ===\n")
cat(sprintf("Cell types in saved object: %d\n",
	                length(levels(nc_obj@idents))))
cat("Cell types:\n")
print(sort(levels(nc_obj@idents)))

cat("\n")
cat("╔═══════════════════════════════════════════════════════════╗\n")
cat("║  ✓✓✓  NeuronChat Dyslexia COMPLETE  ✓✓✓                   ║\n")
cat("╚═══════════════════════════════════════════════════════════╝\n")
cat("Output:", rds_out, "\n")

sink()
