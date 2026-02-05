# lsf.str("package:NeuronChat")

# Load libraries
library(NeuronChat)
library(CellChat)
library(Seurat)
library(ggplot2)
library(patchwork)
library(ggalluvial)
library(ComplexHeatmap)
library(circlize)
library(uwot)

# Check the built-in database
data(package = "NeuronChat")

ls("package:NeuronChat")

# Define working directory where the files are located
work_dir <- "/mnt/nfs/CX000008_DS1/projects/jaeyeon/fastq_file_Dyslexia_r1/zanalysis_bogdan/samples.merged_AG_Harmony_res0.1_anno_091625"
setwd(work_dir)

list.files(pattern = "\\.rds$")

# Load Seurat object
rds_file <- "merged_AG_Harmony_res0.1_anno_091625.AD_CTL.rds"
data_human <- readRDS(rds_file)

# Set assay and extract normalized data
assay <- "RNA"
DefaultAssay(data_human) <- assay

sample_name <- tools::file_path_sans_ext(basename(rds_file))
sample_name

# Remove Ensembl IDs, keep gene symbols
genes <- rownames(data_human[[assay]])
is_ensembl <- grepl("^ENS[A-Z]*G\\d+(\\.\\d+)?$", genes)
keep <- genes[!is_ensembl]

if (sum(is_ensembl) > 0) {
  cat("Removing", sum(is_ensembl), "Ensembl gene IDs, keeping", length(keep), "gene symbols\n")
  data_human <- subset(data_human, features = keep)
}

length(genes)
length(keep)

# Extract normalized data (genes by cells)
data_matrix <- LayerData(data_human, assay = assay, layer = "data")

# unique(rownames(data_matrix))

head(data_human@meta.data,3)
tail(data_human@meta.data,3)
cat("the number of cells is :")
dim(data_human@meta.data)
cat("the number of cell types is :")
unique(data_human@meta.data$celltype)

table(is.na(rownames(data_matrix)))

unique(data_human@meta.data$celltype)

# Get cell groups (adjust column name as needed)
if ("celltype" %in% colnames(data_human@meta.data)) {
  cell_groups <- data_human@meta.data$celltype
} else {
  stop("Cannot find cell type annotation. Please ensure 'celltype' or 'labels' column exists in meta.data")
}

head(cell_groups)
tail(cell_groups)
length(cell_groups)
unique(cell_groups)

# Set species
species <- 'human'

# Create group structure for plotting
# cell_groups <- cell_groups[colnames(data_matrix)]
group <- structure(cell_groups, names = colnames(data_matrix))

# Diagnostic checks
cat("\n=== Diagnostic Checks ===\n")

# Check column names
cat("Column names check:\n")
table(is.na(colnames(data_matrix)))

# Check cell_groups
cat("\nCell groups check:\n")
table(is.na(cell_groups))
cat("Cell groups class:", class(cell_groups), "\n")
cat("Cell groups levels (if factor):", levels(cell_groups), "\n")

# Check if names match
cat("\nDimensions match check:\n")
cat("data_matrix columns:", ncol(data_matrix), "\n")
cat("cell_groups length:", length(cell_groups), "\n")
cat("Names alignment:", identical(names(cell_groups), colnames(data_matrix)), "\n")

# Check for empty strings
cat("\nEmpty string check:\n")
cat("Empty cell groups:", sum(cell_groups == ""), "\n")

# Check the group structure
cat("\nGroup structure check:\n")
cat("Group length:", length(group), "\n")
cat("Group names present:", !any(is.na(names(group))), "\n")
table(is.na(names(group)))

head(group)
length(group)

cat("Data loaded:\n")
cat("  Dimensions:", dim(data_matrix), "\n")
cat("  Cell groups:", length(unique(cell_groups)), "\n")
cat("  Species:", species, "\n")

# ==============================================================================
# CREATE NEURONCHAT OBJECT
# ==============================================================================

x <- createNeuronChat(
  object = data_matrix,
  DB = species,
  group.by = cell_groups
)

# Loads the Neurotransmitter/Neuropeptide Database

# Loads a built-in database of neuronal signaling interactions
# For humans: includes ligand-receptor pairs like:

# Neurotransmitters: GABA, Glutamate, Dopamine, Serotonin, etc.
# Neuropeptides: Oxytocin, Vasopressin, NPY, etc.
# Their receptors: GABBR1, GABBR2, DRD1, DRD2, HTR1A, etc.

# Filters Expression Data

# Keeps only genes that are part of the neuronal signaling database
# Example: If your data has 20,000 genes, it might keep only ~500 that are relevant ligands/receptors

# Organizes Cell Groups

# Associates each cell with its cell type label
# Creates internal structures to analyze communication between cell types

glimpse(x)

# Data slots

# x@data.raw          # Original expression data
# x@data              # Processed expression data
# x@data.signaling    # Expression of signaling molecules only

# Annotation slots

# x@idents            # Cell type labels
# x@meta              # Cell metadata
# x@DB                # The neuronal signaling database

# Analysis slots (filled later by run_NeuronChat)

# x@net               # Communication networks (empty initially)
# x@pvalue            # Statistical significance (empty initially)
# x@net_analysis      # Pattern analysis (empty initially)

# x@data.signaling — signaling-only expression

# What it contains
# Subset of x@data
# Only genes involved in signaling, i.e.:
# ligands
# receptors
# downstream target genes

# Why this matters
# NeuronChat does not infer communication from whole transcriptomes.

# Annotation slots

# x@idents — cell identity labels

# What it is
# The group.by vector
# One label per cell
# Defines
# Nodes in the communication network
# Sender/receiver populations
# Change this → you change the network topology.

# x@meta — cell-level metadata

# What it includes
# Sample IDs
# Region labels

# Class / subclass info (if provided)
# Any metadata passed during construction
# Used for
# Plot labeling

# Subsetting
# Group mapping
# Not directly used in inference math, but critical for interpretation.



# ==============================================================================
# RUN NEURONCHAT ANALYSIS
# ==============================================================================

x <- run_NeuronChat(x, M = 100)

glimpse(x, max.level=2)

slotNames(x)

# x@data
# x@data.signaling
# x@net
# x@pvalue
# x@idents
# x@DB[["GABA_GABBR1"]]

### data.raw
cat("x@data.raw\n")
dim(x@data.raw)
head(x@data.raw)

### data
cat("\nx@data\n")
dim(x@data)
x@data[1:5, 1:5]

### data.signaling
cat("\nx@data.signaling\n")
dim(x@data.signaling)
head(x@data.signaling[, 1:5])

### net0
cat("\nx@net0\n")
length(x@net0)
names(x@net0)[1:5]
str(x@net0[[1]], max.level = 2)

### pvalue
cat("\nx@pvalue\n")
length(x@pvalue)

str(x@pvalue[[1]], max.level = 2)

head(x@pvalue)
tail(x@pvalue)

# names(x@pvalue)[1:5]

dim(x@pvalue)
sum(!is.na(x@pvalue))
range(x@pvalue, na.rm=TRUE)
length(x@pvalue)

### net
cat("\nx@net\n")
length(x@net)
names(x@net)[1:5]

str(x@net[[1]], max.level = 2)

head(x@net)
tail(x@net)

# Your x@net is a named list of 190 pathways, and each pathway is a 14×14 matrix (sender groups × receiver groups). 
# The names you’re seeing like 5HT_HTR1A, 5HT_HTR1B, etc. are your “pathways / interaction pairs”.
# Now you can get “significant pathways” in two ways:

# A) If NeuronChat already zeroed-out non-significant edges
# B) If you want to use p-values explicitly

# pathways with any inferred communication (non-zero)

# Code : chatGPT

sig_paths <- names(x@net)[sapply(x@net, function(m) any(m > 0, na.rm=TRUE))]

length(sig_paths)
head(sig_paths, 30)

# If NeuronChat already zeroed-out non-significant edges
# Then “significant” = any pathway matrix has any non-zero edge.

# Explanation:

# 1. x@net is a LIST of 190 matrices

# 190 pathways/interactions (ligand-receptor pairs)
# Each pathway name follows the pattern: Ligand_Receptor

# Example: 5HT_HTR1A = Serotonin (5-HT) → Serotonin receptor 1A (HTR1A)

# 2. Each matrix is 14 × 14

# 14 rows = 14 sender cell types
# 14 columns = 14 receiver cell types
# Values = communication weight/probability between cell types

# 3. The matrix shows cell type → cell type communication
# For pathway 5HT_HTR1A:

# Row (Sender): Which cell type is sending 5-HT (serotonin)
# Column (Receiver): Which cell type is receiving 5-HT via HTR1A receptor
# Value: Strength/probability of this communication

# 4. All zeros = No significant communication detected
# The 5HT_HTR1A pathway shows all zeros, meaning:

# Either the cells don't express enough 5-HT or HTR1A
# Or the communication wasn't statistically significant
# Or it was filtered out during the analysis

# To rank them by total flow (sum of weights):

sig_strength <- sort(sapply(x@net[sig_paths], function(m) sum(m, na.rm=TRUE)), decreasing=TRUE)
head(sig_strength, 30)

# Code Claude :

# Quick check - which pathways have communication?
pathway_max <- sapply(x@net, function(mat) max(mat, na.rm = TRUE))

# Count non-zero pathways
sum(pathway_max > 0)

# Show pathways with communication
active_pathways <- names(x@net)[pathway_max > 0]
cat("Active pathways:\n")
print(active_pathways)
length(active_pathways)

# Show the values
cat("\nMax communication weights:\n")
print(sort(pathway_max[pathway_max > 0], decreasing = TRUE))

# Returns TRUE if they contain the same elements (order doesn't matter)
setequal(sig_paths, active_pathways)
# Check if sorted versions are identical
identical(sort(sig_paths), sort(active_pathways))

# x <- run_NeuronChat(x, M = 100)
# This only fills: x@net
# x@pvalue (in your case as a vector)
# x@net_analysis is 0 (length 0, empty list) because the communication pattern analysis hasn't been run yet!

### net_analysis
cat("\nx@net_analysis\n")
length(x@net_analysis)
str(x@net_analysis, max.level = 2)

head(x@net_analysis)

### fc
cat("\nx@fc\n")
length(x@fc)
head(x@fc)

# Interpretation (important for NeuronChat)
# x@fc is typically functional contribution / information flow score per pathway : chatGPT

# Zero values mean:
# pathway inferred but contributes nothing at the global level
# or filtered out during aggregation

# Non-zero values = pathways that actually drive communication patterns

# Count how many fc values are non-zero
sum(x@fc != 0, na.rm = TRUE)
# If you want the indices of non-zero entries
which(x@fc != 0)
table(x@fc == 0)

cat("

| Use case                        | Use this slot |
| ------------------------------- | ------------- |
| Show *which pathways exist*     | `x@net`       |
| Show *statistical support*      | `x@pvalue`    |
| Show *drivers of communication* | **`x@fc`** ✅  |
| Rank pathways                   | `x@fc`        |
| Pattern analysis / river plots  | `x@net`       |

")

# ✅ Where “significant pathways” usually come from
# 1️⃣ Active / inferred pathways

# Slot: x@net
active_paths <- names(x@net)[
  sapply(x@net, function(m) any(m > 0, na.rm = TRUE))
]
         
# A pathway is active if any edge is non-zero:
# ➡️ These are pathways with inferred communication
# ➡️ This is what produced your “Active pathways:” list

# 2️⃣ Statistically significant pathways
# Slot: x@pvalue

# If you want explicit statistical filtering

sig_paths_net <- names(x@net)[
  sapply(x@net, function(m) any(m > 0, na.rm = TRUE))
]

# x@pvalue is pathway-level
# Often used implicitly (edges already zeroed in x@net)
# Many workflows do not explicitly threshold this

# 3️⃣ Functionally important pathways (drivers)
# Slot: x@fc
# These are pathways that actually shape the global communication structure:

sig_paths_fc <- names(x@net)[x@fc != 0 & !is.na(x@fc)]

# ➡️ This is what your fc_df contains
# ➡️ These are high-impact pathways, not just present ones

# This is the strongest biological definition of “significant” in NeuronChat.

# ✅ Final answer (one sentence)
# NeuronChat does not store “significant pathways” in a single slot — 
# significance is inferred from x@net (active edges), 
# optionally filtered by x@pvalue, and biologically prioritized using x@fc.

# Slot meanings involved:

# names(x@net)
# → the canonical list of all pathways (≈190)

# x@fc
# → pathway-level functional contribution / information flow

# x@fc != 0
# → pathway contributes non-trivially to the global communication structure

# !is.na(x@fc)
# → removes undefined pathways

# So sig_paths_fc =
# pathways that actually drive the inferred communication network

active_paths
length(active_paths)

sig_paths_net
length(sig_paths_net)

sig_paths_fc
length(sig_paths_fc)

# To map the non-zero x@fc values back to pathway names (which are the names of x@net), and then rank them.

# Pathway names (should be length 190)
pw <- names(x@net)

stopifnot(length(pw) == length(x@fc))

# Indices of non-zero fc
idx <- which(x@fc != 0 & !is.na(x@fc))

# Table: pathway + fc
fc_df <- data.frame(
  pathway = pw[idx],
  fc      = x@fc[idx],
  stringsAsFactors = FALSE
)

# fc_df
# dim(fc_df)

fc_df_sorted <- fc_df[order(fc_df$fc, decreasing = TRUE), ]
fc_df_sorted$rank <- seq_len(nrow(fc_df_sorted))

# fc_df_sorted

# Show top 48 pathways by fc
fc_df_sorted[1:50, ]



# What this table is
# Each row in your fc_df is a NeuronChat signaling pathway (interaction) that has a non-zero functional contribution (fc).
# In plain language:
# These are the signaling pathways that actually drive the global neural communication structure in your dataset.

# x@info is the per-pathway metadata table that NeuronChat uses internally to link each interaction (pathway) to its statistical 
# support, global contribution, and bookkeeping indices ?

# These numbers are internal indices that NeuronChat uses to:

# Track which pathways are active
# Link pathway names to their data
# Support global communication calculations
# Enable cross-referencing between different data slots

# It’s the “ledger” behind x@net, x@pvalue, and x@fc.

### info
cat("\nx@info\n")
length(x@info)
head(x@info)
tail(x@info)

sum(x@info != 0, na.rm = TRUE)

# If it's a list
length(x@DB)

# View structure first
# str(x@DB)
class(x@DB)

# View the data
# head(x@DB)

# If it's a list
length(x@LR)

# View structure first
# str(x@DB)
class(x@LR)

# View the data
# head(x@DB)



# Your object 'x' should contain the database used
# x@DB
# x@LR

dim(x@data.signaling)
head(x@data.signaling)
tail(x@data.signaling)

length(x@info)    # ~190 pathways
length(x@info)    # ~10–15 metrics

### ligand.abundance
cat("\nx@ligand.abundance\n")
dim(x@ligand.abundance)
head(x@ligand.abundance[, 1:5])

library(ggplot2)
library(reshape2)

# Convert to long format for ggplot2
ligand_df <- as.data.frame(x@ligand.abundance)
ligand_long <- melt(ligand_df, variable.name = "Ligand", value.name = "Abundance")

# Log-transformed histogram (add pseudocount to avoid log(0))
hist(x@ligand.abundance, 
     main = "Distribution of Ligand Abundance (log scale)",
     xlab = "Log(Abundance + 1)",
     col = "coral",
     breaks = 50)

# Density plot
plot(density(as.vector(x@ligand.abundance)), 
     main = "Density of Ligand Abundance",
     xlab = "Abundance",
     lwd = 2,
     col = "darkblue")

### target.abundance
cat("\nx@target.abundance\n")
dim(x@target.abundance)
head(x@target.abundance[, 1:5])

library(ggplot2)
library(reshape2)

# Convert to long format for ggplot2
target_df <- as.data.frame(x@target.abundance)
target_long <- melt(target_df, variable.name = "Target", value.name = "Abundance")

# Optional: ggplot2 histogram
ggplot(target_long, aes(x = Abundance)) +
  geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
  labs(title = "Distribution of Target Abundance",
       x = "Abundance", y = "Count") +
  theme_minimal()

# Density plot
plot(density(as.vector(x@target.abundance)), 
     main = "Density of Target Abundance",
     xlab = "Abundance",
     lwd = 2,
     col = "darkblue")

### meta
cat("\nx@meta\n")
dim(x@meta)
head(x@meta)

# What x@idents is :
# x@idents is the cell identity vector used by NeuronChat to define who can talk to whom.

# Concretely:
# It is a factor of length = number of cells
# Each entry assigns a cell to a cell group

### idents
cat("\nx@idents\n")
length(x@idents)
levels(x@idents)
head(x@idents)

### DB
cat("\nx@DB\n")
length(x@DB)
names(x@DB)[1:5]
str(x@DB[[1]], max.level = 2)

### LR
cat("\nx@LR\n")
length(x@LR)
head(x@LR)

### dr
cat("\nx@dr\n")
length(x@dr)
str(x@dr, max.level = 2)

### options
cat("\nx@options\n")
length(x@options)
x@options

# str(x@pvalue)
length(x@pvalue)

#str(x@net)
length(x@net)

glimpse(x@net_analysis)
# print(x@net_analysis)

# What x@idents is (key concept)

# In NeuronChat:
# x@idents = the cell identity vector
# One entry per cell
# Length = number of cells
# Values = cell type labels
# Used to define:
# sender populations
# receiver populations
# network nodes in x@net

# In your output:
# 41,976 cells
# 14 cell types

# ?net_aggregation
# Aggregation of communication networks over all ligand-target pairs
# net_list	
# a list containing communication strength matrices for all ligand-target pairs, e.g., 
# 'net' slot of NeuronChat object after run 'run_NeuronChat'

# Aggregation of communication networks over all ligand-target pairs
# ?net_aggregation

# Aggregate networks
net_aggregated_x <- net_aggregation(x@net, method = 'weight')

net_aggregated_x

# Get interaction names
interaction_names <- names(x@net)
n_celltypes <- length(unique(x@idents))

interaction_names
n_celltypes
length(interaction_names)
length(n_celltypes)

# Check for NaN/NA values
cat("Checking for NaN/NA values...\n")
cat("NaN in net_aggregated_x:", sum(is.nan(net_aggregated_x)), "\n")
cat("NA in net_aggregated_x:", sum(is.na(net_aggregated_x)), "\n")

length(group)
any(is.na(group))

options(
  repr.plot.width = 10,
  repr.plot.height = 8,
  repr.plot.res = 150
)

# optional: output directory
fig_dir <- "results_NeuronChat_AD_CTL_figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# Circle plot – aggregated network
p <- netVisual_circle_neuron(
  net_aggregated_x,
  # group = group,
  vertex.label.cex = 1
)
print(p)

png(
  file.path(fig_dir, paste0("netVisual_circle_aggregated_", sample_name, ".png")),
  width = 12, height = 10, units = "in", res = 300
)
print(p)
dev.off()

# Chord plot – aggregated network
p <- netVisual_chord_neuron(
  x,
  method = "weight",
  # group = group,
  lab.cex = 1
)
print(p)

png(
  file.path(fig_dir, paste0("netVisual_chord_aggregated_", sample_name, ".png")),
  width = 12, height = 10, units = "in", res = 300
)
print(p)
dev.off()


# Heatmap – aggregated
p <- heatmap_aggregated(
  x,
  method = "weight",
 # group = group
)
print(p)

png(
  file.path(fig_dir, paste0("heatmap_aggregated_", sample_name, ".png")),
  width = 12, height = 10, units = "in", res = 300
)
print(p)
dev.off()


# See what's in the net slot
# str(x@net)

# See the names of elements
names(x@net)

# Get all pathway names
pathways <- names(x@net)
pathways

length(pathways)

mat <- x@net[["VIP_VIPR2"]]

sum(mat, na.rm = TRUE)
sum(mat > 0, na.rm = TRUE)
range(mat, na.rm = TRUE)

active_idx <- which(pathway_max > 0)
active_idx
length(active_idx)

cat("the active indexes are :")
length(active_idx)

names(x@net)[active_idx]

active_idx <- match(active_pathways, names(x@net))
active_idx

active_pathways
active_pathways[1]
active_pathways[2]
active_pathways[44]
active_pathways[45]

par(mfrow=c(1,2))
# Visualization for the single interaction pair, circle plot  
netVisual_circle_neuron(x@net$"CO_GUCY1A1",
                        # group=group,
                        vertex.label.cex = 1)
# Visualization for the single interaction pair, chord diagram 
netVisual_chord_neuron(x,
                      interaction_use='CO_GUCY1A1',
                      # group=group,
                      lab.cex = 1)

heatmap_aggregated(
x,
method = "count",
interaction_use = "CO_GUCY1A1",
  # group = group,         # optional: for grouping cell types into classes
  # sender.names = ...,     # optional: subset senders
  # receiver.names = ...    # optional: subset receivers
)

heatmap_aggregated(
x,
method = "weight",
interaction_use = "CO_GUCY1A1",
  # group = group,         # optional: for grouping cell types into classes
  # sender.names = ...,     # optional: subset senders
  # receiver.names = ...    # optional: subset receivers
)

# heatmap_single(x, interaction_name = "NRXN3_NLGN1")
heatmap_single(x, 
               interaction_name = "CO_GUCY1A1")

# Select pathway

# pw <- "NRXN3_NLGN1"

# Extract its matrix
# mat <- x@net[[pw]]

# Plot as heatmap
# ComplexHeatmap::Heatmap(
#  mat,
#  name = pw,
#  cluster_rows = FALSE,
#  cluster_columns = FALSE
# )

# p <- heatmap_aggregated(
#  x,
#  method = "weight",
#  interaction_use = "NRXN3_NLGN1",
  # group = group,         # optional: for grouping cell types into classes
  # sender.names = ...,     # optional: subset senders
  # receiver.names = ...    # optional: subset receivers
# )

# print(p)

# Select pathway : another pathway
# pw <- "NRXN3_NLGN1"

# Extract its matrix
# mat <- x@net[[pw]]

# Plot as heatmap
# ComplexHeatmap::Heatmap(
#  mat,
#  name = pw,
#  cluster_rows = FALSE,
#  cluster_columns = FALSE
# )



# group

# Yes, group is supposed to contain all cells
# (one label per cell, named by cell barcode).
# NeuronChat then aggregates cells → cell types → coarse groups internally when plotting.

group_coarse <- as.character(group)

# Glia
group_coarse[group %in% c("Astro", "Oligo", "OPc", "Micro")] <- "Glia"

# Microglia / immune-like
group_coarse[group %in% c("Immune")] <- "Immune"

# Neurons: Excitatory (any Exc* or IT_Exc*)
group_coarse[grepl("Exc", group)] <- "Neuron_Exc"

# Neurons: Inhibitory (any *Inh*)
group_coarse[grepl("Inh", group)] <- "Neuron_Inh"

# QC / other / unknown
group_coarse[group %in% c("High_mt")] <- "QC_high_mt"
group_coarse[group %in% c("mesenchymal")] <- "mesenchymal"

# Anything still unmapped -> Other
group_coarse[is.na(group_coarse) | group_coarse == ""] <- "Other"

group_coarse <- factor(group_coarse, levels = c("Neuron_Exc","Neuron_Inh","Glia",
                                                "Immune","Other","mesenchymal", "QC_high_mt"))

# Sanity check
table(group, group_coarse)
table(group_coarse)


unique(group_coarse)
unique(group_coarse)

# ?lig_tar_heatmap
# A set of plots illustrating the communication network as well as violin plots for related genes for single ligand-target pair

options(
  repr.plot.width = 30,
  repr.plot.height = 8,
  repr.plot.res = 150
)

# Heatmap for a specific pathway
lig_tar_heatmap(
  x,
  interaction_name='CO_GUCY1A1'
)

ls("package:NeuronChat")

options(
  repr.plot.width = 10,
  repr.plot.height = 10,
  repr.plot.res = 150
)


# ==============================================================================
# Loop over active pathways: circle + chord + heatmaps + ligand–target heatmap
# (Filnames start with safe_name(pw))
# ==============================================================================

# fig_dir <- "AD_CTL_NeuronChat_active_pathways"
# if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

active_pathways <- as.character(active_pathways)

safe_name <- function(x) gsub("[^A-Za-z0-9_\\-]+", "_", x)

for (pw in active_pathways) {
  message("Plotting pathway: ", pw)
  pw_safe <- safe_name(pw)

  # --------------------------------------------------------------------------
  # 1) Circle + Chord (side-by-side)
  # --------------------------------------------------------------------------
  png(
    file.path(fig_dir, paste0(pw_safe, "_circle_chord_", sample_name, ".png")),
    width = 16, height = 8, units = "in", res = 300
  )

  par(mfrow = c(1, 2), mar = c(3, 3, 3, 3))

  try({
    netVisual_circle_neuron(
      x@net[[pw]],
      vertex.label.cex = 1
    )
    title(main = paste0("Circle: ", pw))
  }, silent = TRUE)

  try({
    netVisual_chord_neuron(
      x,
      interaction_use = pw,
      lab.cex = 1
    )
    title(main = paste0("Chord: ", pw))
  }, silent = TRUE)

  dev.off()

  # --------------------------------------------------------------------------
  # 2) Aggregated heatmap – COUNT
  # --------------------------------------------------------------------------
  png(
    file.path(fig_dir, paste0(pw_safe, "_heatmap_count_", sample_name, ".png")),
    width = 10, height = 9, units = "in", res = 300
  )

  try({
    heatmap_aggregated(
      x,
      method = "count",
      interaction_use = pw
    )
    title(main = paste0("Aggregated count: ", pw))
  }, silent = TRUE)

  dev.off()

  # --------------------------------------------------------------------------
  # 3) Aggregated heatmap – WEIGHT
  # --------------------------------------------------------------------------
  png(
    file.path(fig_dir, paste0(pw_safe, "_heatmap_weight_", sample_name, ".png")),
    width = 10, height = 9, units = "in", res = 300
  )

  try({
    heatmap_aggregated(
      x,
      method = "weight",
      interaction_use = pw
    )
    title(main = paste0("Aggregated weight: ", pw))
  }, silent = TRUE)

  dev.off()

  # --------------------------------------------------------------------------
  # 4) Single-interaction heatmap
  # --------------------------------------------------------------------------
  png(
    file.path(fig_dir, paste0(pw_safe, "_heatmap_single_", sample_name, ".png")),
    width = 10, height = 9, units = "in", res = 300
  )

  try({
    heatmap_single(
      x,
      interaction_name = pw
    )
    title(main = paste0("Single interaction: ", pw))
  }, silent = TRUE)

  dev.off()

  # --------------------------------------------------------------------------
  # 5) Ligand–target heatmap (wide figure; also sets Jupyter repr options)
  # --------------------------------------------------------------------------
  options(
    repr.plot.width  = 30,
    repr.plot.height = 8,
    repr.plot.res    = 150
  )

  png(
    file.path(fig_dir, paste0(pw_safe, "_ligand_target_heatmap_", sample_name, ".png")),
    width = 30, height = 8, units = "in", res = 300
  )

  try({
    lig_tar_heatmap(
      x,
      interaction_name = pw
    )
    title(main = paste0("Ligand–target heatmap: ", pw))
  }, silent = TRUE)

  dev.off()
}

message("All pathway plots saved to: ", normalizePath(fig_dir))




options(
  repr.plot.width = 30,
  repr.plot.height = 20,
  repr.plot.res = 150
)

# cutoff.pvalue	
# the cutoff of pvalue when doing Wilcoxon test; Default = 0.05

rankNet_Neuron(x, slot.name = "net", measure = "count", mode = "single")

# ------------------------------------------------------------------------------
# Save rankNet_Neuron plot (count, single)
# ------------------------------------------------------------------------------

png(
  file.path(fig_dir, paste0("rankNet_single_count_", sample_name, ".png")),
  width = 12, height = 8, units = "in", res = 300
)

rankNet_Neuron(
  x,
  slot.name = "net",
  measure = "count",
  mode = "single"
)

dev.off()


# Rank signaling networks based on the information flow or the number of interactions; adapted from CellChat https://github.com/sqjin/CellChat

options(
  repr.plot.width = 30,
  repr.plot.height = 20,
  repr.plot.res = 150
)

# cutoff.pvalue	
# the cutoff of pvalue when doing Wilcoxon test; Default = 0.05
rankNet_Neuron(x, slot.name = "net", measure = "weight", mode = "single")


png(
  file.path(fig_dir, paste0("rankNet_single_weight_", sample_name, ".png")),
  width = 12, height = 8, units = "in", res = 300
)

rankNet_Neuron(
  x,
  slot.name = "net",
  measure = "weight",
  mode = "single"
)

dev.off()


png(
  file.path(fig_dir, paste0("rankNet_single_count_vs_weight_", sample_name, ".png")),
  width = 16, height = 8, units = "in", res = 300
)

par(mfrow = c(1, 2))

rankNet_Neuron(x, slot.name = "net", measure = "count",  mode = "single")
rankNet_Neuron(x, slot.name = "net", measure = "weight", mode = "single")

dev.off()


options(
  repr.plot.width = 20,
  repr.plot.height = 20,
  repr.plot.res = 150
)

# ==============================================================================
# COMMUNICATION PATTERNS
# ==============================================================================

# For pattern identification (outgoing)
# Open PNG device first, then call function which both modifies object and creates plot
png(
  file.path(fig_dir, paste0("pattern_outgoing_", sample_name, ".png")),
  width = 7, height = 6, units = "in", res = 300
)
x <- identifyCommunicationPatterns_Neuron(
  x, 
  slot.name = "net",      # ✓ This exists
  pattern = "outgoing", 
  k = 4,
  height = 18
)
dev.off()

# For pattern identification (incoming)
# Open PNG device first, then call function which both modifies object and creates plot
png(
  file.path(fig_dir, paste0("pattern_incoming_", sample_name, ".png")),
  width = 7, height = 6, units = "in", res = 300
)
x <- identifyCommunicationPatterns_Neuron(
  x, 
  slot.name = "net",      # ✓ This exists
  pattern = "incoming", 
  k = 4, 
  height = 18
)
dev.off()



# x

# 1️⃣ $pattern$incoming / $pattern$outgoing

# This is a cell group × signaling (or LR / pathway) matrix.
# Example from your screenshot:

# Rows   = cell groups (Oligo, L23_IT_Exc, MGE_Inh, …)
# Cols   = signaling entities (e.g. 5HT_HTR1A, 5HT_HTR1B, …)
# Values = contribution / importance (often 0 after filtering)

# 2️⃣ $pattern (data.frame, 56 × 3)
# CellGroup | Pattern   | Contribution
# ------------------------------------
# Oligo     | Pattern 1 | 3.3e-10
# L23_IT_Exc| Pattern 1 | 7.2e-01
# ...

# This tells you:
# “How much each cell group contributes to each communication pattern”

# 3️⃣ $signaling (data.frame, 176 × 3)
# Pattern   | Signaling     | Contribution
# ----------------------------------------
# Pattern 1 | CO_GUCY1A1   | 5.9e-05
# ...

# This tells you:
# “Which signaling molecules / LR pairs / pathways define each pattern”

# 🚫 Why you see mostly zeros
# This is expected and normal:
# Pattern analysis uses NMF / SVD-like decomposition
# Only top contributors survive
# Most entries are near-zero → displayed as 0

# This does not mean:
# no interactions exist
# signaling is absent

x@net_analysis$pattern

x@net_analysis$pattern

x@net_analysis$pattern$outgoing
x@net_analysis$pattern$incoming

x@net_analysis$pattern$incoming$data
x@net_analysis$pattern$incoming$pattern$cell
x@net_analysis$pattern$incoming$pattern$signaling

x@net_analysis$pattern$outgoing$data
x@net_analysis$pattern$outgoing$pattern$cell
x@net_analysis$pattern$outgoing$pattern$signaling



options(
  repr.plot.width = 50,
  repr.plot.height = 10,
  repr.plot.res = 150
)

# Libraries
library(pheatmap)
library(dplyr)
library(tidyr)
library(tibble)

# ---- output folder ----
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# ---- robust getters (handles both nested and top-level structures) ----
get_incoming_data <- function(x) {
  if (!is.null(x@net_analysis$pattern$incoming$data)) return(x@net_analysis$pattern$incoming$data)
  stop("Could not find x@net_analysis$pattern$incoming$data")
}

get_incoming_cell_df <- function(x) {
  # preferred (as you requested)
  if (!is.null(x@net_analysis$pattern$incoming$pattern$cell)) return(x@net_analysis$pattern$incoming$pattern$cell)
  # fallback (common structure)
  if (!is.null(x@net_analysis$cell)) return(x@net_analysis$cell)
  stop("Could not find incoming cell data.frame (nested or top-level).")
}

get_incoming_signaling_df <- function(x) {
  # preferred (as you requested)
  if (!is.null(x@net_analysis$pattern$incoming$pattern$signaling)) return(x@net_analysis$pattern$incoming$pattern$signaling)
  # fallback (common structure)
  if (!is.null(x@net_analysis$signaling)) return(x@net_analysis$signaling)
  stop("Could not find incoming signaling data.frame (nested or top-level).")
}

# =========================
# 1) incoming$data  (matrix)
# =========================
mat_in <- get_incoming_data(x)

# Write numbers to file (tab-separated)
write.table(
  mat_in,
  file = file.path(fig_dir, paste0("net_analysis_incoming_data_matrix_", sample_name, ".txt")),
  sep = "\t", quote = FALSE, col.names = NA
)

# Display in Jupyter (simple pheatmap, with numbers)
# pheatmap(
#  mat_in,
#  cluster_rows = FALSE,
#  cluster_cols = FALSE,
#  display_numbers = TRUE,
#  number_format = "%.2g",
#  fontsize_number = 7,
#  fontsize_row = 9,
#  fontsize_col = 9,
#  main = "Incoming: pattern$data"
# )

# Save the same heatmap
png(file.path(fig_dir, paste0("net_analysis_incoming_data_matrix_heatmap_", sample_name, ".png")),
    width = 12, height = 8, units = "in", res = 300)
pheatmap(
  mat_in,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  number_format = "%.2g",
  fontsize_number = 7,
  fontsize_row = 9,
  fontsize_col = 9,
  main = "Incoming: pattern$data"
)
dev.off()


# ==========================================
# 2) incoming$pattern$cell (data.frame long)
# ==========================================
cell_df_in <- get_incoming_cell_df(x)

# Save table to file
write.table(
  cell_df_in,
  file = file.path(fig_dir, paste0("net_analysis_incoming_cell_table_", sample_name, ".txt")),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Convert to matrix: rows=CellGroup, cols=Pattern, values=Contribution
cell_mat_in <- cell_df_in %>%
  mutate(CellGroup = as.character(CellGroup),
         Pattern   = as.character(Pattern)) %>%
  select(CellGroup, Pattern, Contribution) %>%
  pivot_wider(names_from = Pattern, values_from = Contribution, values_fill = 0) %>%
  column_to_rownames("CellGroup") %>%
  as.matrix()

# Also save the matrix numbers
write.table(
  cell_mat_in,
  file = file.path(fig_dir, paste0("net_analysis_incoming_cell_matrix_", sample_name, ".txt")),
  sep = "\t", quote = FALSE, col.names = NA
)

# Display in Jupyter
# pheatmap(
#  cell_mat_in,
#  cluster_rows = FALSE,
# cluster_cols = FALSE,
# display_numbers = TRUE,
#  number_format = "%.2g",
#  fontsize_number = 7,
#  fontsize_row = 9,
#  fontsize_col = 10,
#  main = "Incoming: cell contributions"
# )

# Save PNG
png(file.path(fig_dir, paste0("net_analysis_incoming_cell_heatmap_", sample_name, ".png")),
    width = 7, height = 6, units = "in", res = 300)
pheatmap(
  cell_mat_in,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  number_format = "%.2g",
  fontsize_number = 7,
  fontsize_row = 9,
  fontsize_col = 10,
  main = "Incoming: cell contributions"
)
dev.off()


# ===============================================
# 3) incoming$pattern$signaling (data.frame long)
# ===============================================
sig_df_in <- get_incoming_signaling_df(x)

# Save table to file
write.table(
  sig_df_in,
  file = file.path(fig_dir, paste0("net_analysis_incoming_signaling_table_", sample_name, ".txt")),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Convert to matrix: rows=Signaling, cols=Pattern, values=Contribution
sig_mat_in <- sig_df_in %>%
  mutate(Signaling = as.character(Signaling),
         Pattern   = as.character(Pattern)) %>%
  select(Signaling, Pattern, Contribution) %>%
  pivot_wider(names_from = Pattern, values_from = Contribution, values_fill = 0) %>%
  column_to_rownames("Signaling") %>%
  as.matrix()

# Save matrix numbers
write.table(
  sig_mat_in,
  file = file.path(fig_dir, paste0("net_analysis_incoming_signaling_matrix_", sample_name, ".txt")),
  sep = "\t", quote = FALSE, col.names = NA
)

# Display in Jupyter (often huge; show top 40 rows by max contribution for readability)
topN <- 100
ord <- order(apply(sig_mat_in, 1, max), decreasing = TRUE)
sig_mat_in_top <- sig_mat_in[ord[seq_len(min(topN, nrow(sig_mat_in)))], , drop = FALSE]

# pheatmap(
#  sig_mat_in_top,
#  cluster_rows = FALSE,
#  cluster_cols = FALSE,
#  display_numbers = TRUE,
#  number_format = "%.2g",
#  fontsize_number = 7,
#  fontsize_row = 7,
#  fontsize_col = 10,
#  main = paste0("Incoming: signaling contributions (top ", nrow(sig_mat_in_top), ")")
# )

# Save PNG
png(file.path(fig_dir, paste0("net_analysis_incoming_signaling_heatmap_top_", sample_name, ".png")),
    width = 10, height = 10, units = "in", res = 300)
pheatmap(
  sig_mat_in_top,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  number_format = "%.2g",
  fontsize_number = 7,
  fontsize_row = 7,
  fontsize_col = 10,
  main = paste0("Incoming: signaling contributions (top ", nrow(sig_mat_in_top), ")")
)
dev.off()


options(
  repr.plot.width = 10,
  repr.plot.height = 50,
  repr.plot.res = 150
)

# Display in Jupyter (simple pheatmap, with numbers)
pheatmap(
  t(mat_in),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  number_format = "%.2g",
  fontsize_number = 7,
  fontsize_row = 9,
  fontsize_col = 9,
  main = "Incoming: pattern$data"
)


png(
  file.path(fig_dir, paste0("incoming_pattern_data_heatmap_", sample_name, ".png")),
  width = 12,
  height = 8,
  units = "in",
  res = 300
)

pheatmap(
  t(mat_in),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  number_format = "%.2g",
  fontsize_number = 7,
  fontsize_row = 9,
  fontsize_col = 9,
  main = "Incoming: pattern$data"
)

dev.off()


options(
  repr.plot.width = 10,
  repr.plot.height = 20,
  repr.plot.res = 150
)

# Keep only rows with at least one non-zero value
rows_keep <- rowSums(mat_in != 0) > 0

# Keep only columns with at least one non-zero value
cols_keep <- colSums(mat_in != 0) > 0

mat_in_nz <- mat_in[rows_keep, cols_keep, drop = FALSE]

dim(mat_in)       # original
dim(mat_in_nz)    # filtered

pheatmap(
  t(mat_in_nz),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  number_format = "%.2g",
  fontsize_number = 7,
  fontsize_row = 19,
  fontsize_col = 19,
  main = "Incoming: pattern$data (non-zero rows/columns only)"
)

png(
  file.path(fig_dir, paste0("incoming_pattern_data_nonzero_heatmap_", sample_name, ".png")),
  width = 10,
  height = 50,
  units = "in",
  res = 150
)

pheatmap(
  t(mat_in_nz),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  number_format = "%.2g",
  fontsize_number = 7,
  fontsize_row = 19,
  fontsize_col = 19,
  main = "Incoming: pattern$data (non-zero rows/columns only)"
)

dev.off()


options(
  repr.plot.width = 15,
  repr.plot.height = 20,
  repr.plot.res = 150
)

pheatmap(
  t(mat_in_nz),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  display_numbers = TRUE,
  number_format = "%.2g",
  fontsize_number = 7,
  fontsize_row = 19,
  fontsize_col = 19,
  main = "Incoming: pattern$data (clustered)"
)

png(
  file.path(fig_dir, paste0("incoming_pattern_data_nonzero_clustered_heatmap_", sample_name, ".png")),
  width = 50,
  height = 15,
  units = "in",
  res = 150
)

pheatmap(
  t(mat_in_nz),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  display_numbers = TRUE,
  number_format = "%.2g",
  fontsize_number = 7,
  fontsize_row = 19,
  fontsize_col = 19,
  main = "Incoming: pattern$data (clustered)"
)

dev.off()




options(
  repr.plot.width = 10,
  repr.plot.height = 10,
  repr.plot.res = 150
)

pheatmap(
  cell_mat_in,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  number_format = "%.2g",
  fontsize_number = 7,
  fontsize_row = 9,
  fontsize_col = 10,
  main = "Incoming: cell contributions"
)

png(
  file.path(fig_dir, paste0("incoming_cell_contributions_heatmap_", sample_name, ".png")),
  width = 10,
  height = 10,
  units = "in",
  res = 150
)

pheatmap(
  cell_mat_in,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  number_format = "%.2g",
  fontsize_number = 7,
  fontsize_row = 9,
  fontsize_col = 10,
  main = "Incoming: cell contributions"
)

dev.off()




options(
  repr.plot.width = 8,
  repr.plot.height = 10,
  repr.plot.res = 150
)

pheatmap(
  sig_mat_in_top,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  number_format = "%.2g",
  fontsize_number = 7,
  fontsize_row = 7,
  fontsize_col = 10,
  main = paste0("Incoming: signaling contributions (top ", nrow(sig_mat_in_top), ")")

)

png(
  file.path(fig_dir, paste0("incoming_signaling_contributions_top_heatmap_", sample_name, ".png")),
  width = 8,
  height = 10,
  units = "in",
  res = 150
)

pheatmap(
  sig_mat_in_top,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  number_format = "%.2g",
  fontsize_number = 7,
  fontsize_row = 7,
  fontsize_col = 10,
  main = paste0(
    "Incoming: signaling contributions (top ",
    nrow(sig_mat_in_top),
    ")"
  )
)

dev.off()




options(
  repr.plot.width = 50,
  repr.plot.height = 10,
  repr.plot.res = 150
)

# =========================
# OUTGOING patterns (same as incoming)
# =========================

# ---- robust getters (outgoing) ----
get_outgoing_data <- function(x) {
  if (!is.null(x@net_analysis$pattern$outgoing$data)) return(x@net_analysis$pattern$outgoing$data)
  stop("Could not find x@net_analysis$pattern$outgoing$data")
}

get_outgoing_cell_df <- function(x) {
  # preferred (if present)
  if (!is.null(x@net_analysis$pattern$outgoing$pattern$cell)) return(x@net_analysis$pattern$outgoing$pattern$cell)
  # fallback (most common: stored at top-level net_analysis)
  if (!is.null(x@net_analysis$cell)) return(x@net_analysis$cell)
  stop("Could not find outgoing cell data.frame (nested or top-level).")
}

get_outgoing_signaling_df <- function(x) {
  # preferred (if present)
  if (!is.null(x@net_analysis$pattern$outgoing$pattern$signaling)) return(x@net_analysis$pattern$outgoing$pattern$signaling)
  # fallback (most common: stored at top-level net_analysis)
  if (!is.null(x@net_analysis$signaling)) return(x@net_analysis$signaling)
  stop("Could not find outgoing signaling data.frame (nested or top-level).")
}

# =========================
# 1) outgoing$data  (matrix)
# =========================
mat_out <- get_outgoing_data(x)

# Write numbers to file (tab-separated)
write.table(
  mat_out,
  file = file.path(fig_dir, paste0("net_analysis_outgoing_data_matrix_", sample_name, ".txt")),
  sep = "\t", quote = FALSE, col.names = NA
)

# Save heatmap PNG
png(file.path(fig_dir, paste0("net_analysis_outgoing_data_matrix_heatmap_", sample_name, ".png")),
    width = 12, height = 8, units = "in", res = 300)
pheatmap(
  mat_out,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  number_format = "%.2g",
  fontsize_number = 7,
  fontsize_row = 9,
  fontsize_col = 9,
  main = "Outgoing: pattern$data"
)
dev.off()

# ==========================================
# 2) outgoing$pattern$cell (data.frame long)
# ==========================================
cell_df_out <- get_outgoing_cell_df(x)

# Save table to file
write.table(
  cell_df_out,
  file = file.path(fig_dir, paste0("net_analysis_outgoing_cell_table_", sample_name, ".txt")),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Convert to matrix: rows=CellGroup, cols=Pattern, values=Contribution
cell_mat_out <- cell_df_out %>%
  mutate(CellGroup = as.character(CellGroup),
         Pattern   = as.character(Pattern)) %>%
  select(CellGroup, Pattern, Contribution) %>%
  pivot_wider(names_from = Pattern, values_from = Contribution, values_fill = 0) %>%
  column_to_rownames("CellGroup") %>%
  as.matrix()

# Save matrix numbers
write.table(
  cell_mat_out,
  file = file.path(fig_dir, paste0("net_analysis_outgoing_cell_matrix_", sample_name, ".txt")),
  sep = "\t", quote = FALSE, col.names = NA
)

# Save PNG
png(file.path(fig_dir, paste0("net_analysis_outgoing_cell_heatmap_", sample_name, ".png")),
    width = 7, height = 6, units = "in", res = 300)
pheatmap(
  cell_mat_out,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  number_format = "%.2g",
  fontsize_number = 7,
  fontsize_row = 9,
  fontsize_col = 10,
  main = "Outgoing: cell contributions"
)
dev.off()

# ===============================================
# 3) outgoing$pattern$signaling (data.frame long)
# ===============================================
sig_df_out <- get_outgoing_signaling_df(x)

# Save table to file
write.table(
  sig_df_out,
  file = file.path(fig_dir, paste0("net_analysis_outgoing_signaling_table_", sample_name, ".txt")),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Convert to matrix: rows=Signaling, cols=Pattern, values=Contribution
sig_mat_out <- sig_df_out %>%
  mutate(Signaling = as.character(Signaling),
         Pattern   = as.character(Pattern)) %>%
  select(Signaling, Pattern, Contribution) %>%
  pivot_wider(names_from = Pattern, values_from = Contribution, values_fill = 0) %>%
  column_to_rownames("Signaling") %>%
  as.matrix()

# Save matrix numbers
write.table(
  sig_mat_out,
  file = file.path(fig_dir, paste0("net_analysis_outgoing_signaling_matrix_", sample_name, ".txt")),
  sep = "\t", quote = FALSE, col.names = NA
)

# Show top N for readability
topN <- 100
ord <- order(apply(sig_mat_out, 1, max), decreasing = TRUE)
sig_mat_out_top <- sig_mat_out[ord[seq_len(min(topN, nrow(sig_mat_out)))], , drop = FALSE]

# Save PNG
png(file.path(fig_dir, paste0("net_analysis_outgoing_signaling_heatmap_top_", sample_name, ".png")),
    width = 10, height = 10, units = "in", res = 300)
pheatmap(
  sig_mat_out_top,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  number_format = "%.2g",
  fontsize_number = 7,
  fontsize_row = 7,
  fontsize_col = 10,
  main = paste0("Outgoing: signaling contributions (top ", nrow(sig_mat_out_top), ")")
)
dev.off()


options(
  repr.plot.width = 10,
  repr.plot.height = 30,
  repr.plot.res = 150
)


# Display in Jupyter (simple pheatmap, with numbers)
pheatmap(
  t(mat_out),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  number_format = "%.2g",
  fontsize_number = 7,
  fontsize_row = 9,
  fontsize_col = 9,
  main = "Outgoing: pattern$data"
)

png(
  file.path(fig_dir, paste0("outgoing_pattern_data_heatmap_", sample_name, ".png")),
  width = 12,
  height = 8,
  units = "in",
  res = 300
)

pheatmap(
  t(mat_out),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  number_format = "%.2g",
  fontsize_number = 7,
  fontsize_row = 9,
  fontsize_col = 9,
  main = "Outgoing: pattern$data"
)

dev.off()


options(
  repr.plot.width = 20,
  repr.plot.height = 30,
  repr.plot.res = 150
)


# Keep only rows with at least one non-zero value
rows_keep <- rowSums(mat_out != 0) > 0

# Keep only columns with at least one non-zero value
cols_keep <- colSums(mat_out != 0) > 0

mat_out_nz <- mat_out[rows_keep, cols_keep, drop = FALSE]

dim(mat_out)       # original
dim(mat_out_nz)    # filtered


pheatmap(
  t(mat_out_nz),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  number_format = "%.2g",
  fontsize_number = 7,
  fontsize_row = 19,
  fontsize_col = 19,
  main = "Outgoing: pattern$data (non-zero rows/columns only)"
)

png(
  file.path(fig_dir, paste0("outgoing_pattern_data_nonzero_heatmap_", sample_name, ".png")),
  width = 50,
  height = 10,
  units = "in",
  res = 150
)

pheatmap(
  t(mat_out_nz),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  number_format = "%.2g",
  fontsize_number = 7,
  fontsize_row = 19,
  fontsize_col = 19,
  main = "Outgoing: pattern$data (non-zero rows/columns only)"
)

dev.off()



options(
  repr.plot.width = 20,
  repr.plot.height = 20,
  repr.plot.res = 150
)

pheatmap(
  t(mat_out_nz),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  display_numbers = TRUE,
  number_format = "%.2g",
  fontsize_number = 7,
  fontsize_row = 19,
  fontsize_col = 19,
  main = "Outgoing: pattern$data (clustered)"
)

png(
  file.path(fig_dir, paste0("outgoing_pattern_data_nonzero_clustered_heatmap_", sample_name, ".png")),
  width = 50,
  height = 15,
  units = "in",
  res = 150
)

pheatmap(
  t(mat_out_nz),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  display_numbers = TRUE,
  number_format = "%.2g",
  fontsize_number = 7,
  fontsize_row = 19,
  fontsize_col = 19,
  main = "Outgoing: pattern$data (clustered)"
)

dev.off()


options(
  repr.plot.width = 10,
  repr.plot.height = 10,
  repr.plot.res = 150
)

pheatmap(
  cell_mat_out,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  number_format = "%.2g",
  fontsize_number = 7,
  fontsize_row = 9,
  fontsize_col = 10,
  main = "Outgoing: cell contributions"
)

png(
  file.path(fig_dir, paste0("outgoing_cell_contributions_heatmap_", sample_name, ".png")),
  width = 10,
  height = 10,
  units = "in",
  res = 150
)

pheatmap(
  cell_mat_out,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  number_format = "%.2g",
  fontsize_number = 7,
  fontsize_row = 9,
  fontsize_col = 10,
  main = "Outgoing: cell contributions"
)

dev.off()


options(
  repr.plot.width = 8,
  repr.plot.height = 10,
  repr.plot.res = 150
)

pheatmap(
  sig_mat_out_top,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  number_format = "%.2g",
  fontsize_number = 7,
  fontsize_row = 7,
  fontsize_col = 10,
  main = paste0("Outgoing: signaling contributions (top ", nrow(sig_mat_out_top), ")")
)

png(
  file.path(fig_dir, paste0("outgoing_signaling_contributions_top_heatmap_", sample_name, ".png")),
  width = 8,
  height = 10,
  units = "in",
  res = 150
)

pheatmap(
  sig_mat_out_top,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  number_format = "%.2g",
  fontsize_number = 7,
  fontsize_row = 7,
  fontsize_col = 10,
  main = paste0(
    "Outgoing: signaling contributions (top ",
    nrow(sig_mat_out_top),
    ")"
  )
)

dev.off()






options(
  repr.plot.width = 16,
  repr.plot.height = 10,
  repr.plot.res = 150
)

# show in Jupyter
netAnalysis_river_Neuron(
  x,
  slot.name = "net",
  pattern = "outgoing",
  font.size = 2.5,
  cutoff.1 = 0.5,
  cutoff.2 = 0.5
)

# save PNG
png(
  file.path(fig_dir, paste0("netAnalysis_river_outgoing_", sample_name, ".png")),
  width = 10, height = 8, units = "in", res = 300
)
netAnalysis_river_Neuron(
  x,
  slot.name = "net",
  pattern = "outgoing",
  font.size = 2.5,
  cutoff.1 = 0.5,
  cutoff.2 = 0.5
)
dev.off()


options(
  repr.plot.width = 16,
  repr.plot.height = 10,
  repr.plot.res = 150
)

# show in Jupyter
netAnalysis_river_Neuron(
  x,
  slot.name = "net",
  pattern = "incoming",
  font.size = 2.5,
  cutoff.1 = 0.5,
  cutoff.2 = 0.5
)

# save PNG
png(
  file.path(fig_dir, paste0("netAnalysis_river_incoming_", sample_name, ".png")),
  width = 10, height = 8, units = "in", res = 300
)
netAnalysis_river_Neuron(
  x,
  slot.name = "net",
  pattern = "incoming",
  font.size = 2.5,
  cutoff.1 = 0.5,
  cutoff.2 = 0.5
)
dev.off()




# Alias (for clarity / safety)
neta <- x@net_analysis

# 1) pattern matrices (cell group × signaling)
out_mat <- neta$pattern$outgoing$data
in_mat  <- neta$pattern$incoming$data

write.table(
  out_mat,
  file = file.path(fig_dir, paste0("net_analysis_pattern_outgoing_matrix_", sample_name, ".txt")),
  sep = "\t", quote = FALSE, col.names = NA
)

write.table(
  in_mat,
  file = file.path(fig_dir, paste0("net_analysis_pattern_incoming_matrix_", sample_name, ".txt")),
  sep = "\t", quote = FALSE, col.names = NA
)

# 2) cell contributions (CellGroup × Pattern)
# Note: Access from incoming or outgoing pattern (they should be the same)
# Try incoming first, fallback to outgoing if needed
if (!is.null(neta$pattern$incoming$pattern$cell)) {
  write.table(
    neta$pattern$incoming$pattern$cell,
    file = file.path(fig_dir, paste0("net_analysis_cell_contribution_", sample_name, ".txt")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
} else if (!is.null(neta$pattern$outgoing$pattern$cell)) {
  write.table(
    neta$pattern$outgoing$pattern$cell,
    file = file.path(fig_dir, paste0("net_analysis_cell_contribution_", sample_name, ".txt")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
} else if (!is.null(neta$cell)) {
  write.table(
    neta$cell,
    file = file.path(fig_dir, paste0("net_analysis_cell_contribution_", sample_name, ".txt")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
} else {
  warning("Could not find cell contribution data")
}

# 3) signaling contributions (Pattern × Signaling)
# Try incoming first, fallback to outgoing if needed
if (!is.null(neta$pattern$incoming$pattern$signaling)) {
  write.table(
    neta$pattern$incoming$pattern$signaling,
    file = file.path(fig_dir, paste0("net_analysis_signaling_contribution_", sample_name, ".txt")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
} else if (!is.null(neta$pattern$outgoing$pattern$signaling)) {
  write.table(
    neta$pattern$outgoing$pattern$signaling,
    file = file.path(fig_dir, paste0("net_analysis_signaling_contribution_", sample_name, ".txt")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
} else if (!is.null(neta$signaling)) {
  write.table(
    neta$signaling,
    file = file.path(fig_dir, paste0("net_analysis_signaling_contribution_", sample_name, ".txt")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
} else {
  warning("Could not find signaling contribution data")
}





# ==============================================================================
# FUNCTIONAL SIMILARITY AND EMBEDDING (using uwot)
# ==============================================================================

# Compute functional similarity
x <- computeNetSimilarity_Neuron(x, type = 'functional')

# Embedding using uwot (instead of umap-learn)
x <- netEmbedding(x, slot.name = "net_analysis", type = "functional", umap.method = "uwot")

# Clustering
x <- netClustering(x, type = 'functional', slot.name = "net_analysis", k = 5)

netVisual_embedding_Neuron(x, type = "functional", label.size = 3.5, pathway.remove.show = F)

# Visualization
png(file.path(fig_dir, paste0("netVisual_embedding_functional_", sample_name, ".png")), width = 7, height = 6, units = "in", res = 300)
netVisual_embedding_Neuron(x, type = "functional", label.size = 3.5, pathway.remove.show = F)
dev.off()

# Visualization
png(file.path(fig_dir, paste0("netVisual_embedding_functional_zoom_", sample_name, ".png")), width = 7, height = 6, units = "in", res = 300)
netVisual_embeddingZoomIn_Neuron(x, type = "functional", nCol = 2,label.size = 3)
dev.off()



# ==============================================================================
# STRUCTURAL SIMILARITY AND EMBEDDING (using uwot)
# ==============================================================================

# Compute structural similarity
x <- computeNetSimilarity_Neuron(x, type = 'structural')

# Embedding using uwot
x <- netEmbedding(x, slot.name = "net_analysis", type = "structural", umap.method = "uwot")

# Clustering
x <- netClustering(x, type = 'structural', slot.name = "net_analysis", k = 5)

netVisual_embedding_Neuron(x, type = "structural", label.size = 3.5, pathway.remove.show = F)

netVisual_embedding_Neuron(x, type = "structural", label.size = 3.5, pathway.remove.show = F)

# Visualization
png(file.path(fig_dir, paste0("netVisual_embedding_structural_", sample_name, ".png")), width = 7, height = 6, units = "in", res = 300)
netVisual_embedding_Neuron(x, type = "structural", label.size = 3.5, pathway.remove.show = F)
dev.off()

# Visualization
png(file.path(fig_dir, paste0("netVisual_embedding_structural_zoom_", sample_name, ".png")), width = 7, height = 6, units = "in", res = 300)
netVisual_embeddingZoomIn_Neuron(x, type = "structural", nCol = 2,label.size = 3)
dev.off()

options(
  repr.plot.width = 16,
  repr.plot.height = 16,
  repr.plot.res = 150
)
netVisual_embeddingZoomIn_Neuron(x, type = "structural", nCol = 2,label.size = 3)



interaction_names[1]
interaction_names[190]



# x@net_analysis$similarity

# x@net_analysis$similarity$functional
# x@net_analysis$similarity$functional$matrix
# x@net_analysis$similarity$functional$matrix$single

# x@net_analysis$similarity

# x@net_analysis$similarity$structural
# x@net_analysis$similarity$structural$matrix
# x@net_analysis$similarity$structural$matrix$single

# x@net_analysis
# x@net_analysis$similarity$functional$matrix$single
# x@net_analysis$similarity$structural$matrix$single

# Inline plotting
options(
  repr.plot.width = 12,
  repr.plot.height = 12,
  repr.plot.res = 150
)


netVisual_embeddingZoomIn_Neuron(x, type = "functional", nCol = 2,label.size = 3)

png(
  file.path(fig_dir, paste0("netVisual_embeddingZoomIn_functional_", sample_name, ".png")),
  width = 12, height = 12, units = "in", res = 300
)

netVisual_embeddingZoomIn_Neuron(
  x,
  type = "functional",
  nCol = 2,
  label.size = 3
)

dev.off()

# Inline plotting
options(
  repr.plot.width = 12,
  repr.plot.height = 12,
  repr.plot.res = 150
)


netVisual_embeddingZoomIn_Neuron(x, 
                                 type = "structural", 
                                 nCol = 2,label.size = 3)

png(
  file.path(fig_dir, paste0("netVisual_embeddingZoomIn_structural_", sample_name, ".png")),
  width = 12, height = 12, units = "in", res = 300
)

netVisual_embeddingZoomIn_Neuron(
  x,
  type = "structural",
  nCol = 2,
  label.size = 3
)

dev.off()




# ==============================================================================
# HEATMAPS FOR INDIVIDUAL INTERACTIONS (SAVE ALL)
# ==============================================================================

# Optional: where to save (set to "." if you want current folder)
# fig_dir <- "AD_CTL_interaction_heatmaps"
# if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# Helper to make safe filenames
safe_name <- function(x) gsub("[^A-Za-z0-9_\\-\\.]+", "_", x)

if (length(interaction_names) > 0) {
  for (interaction_example in interaction_names) {

    # Get the interaction matrix (skip if missing)
    if (!interaction_example %in% names(x@net)) {
      cat("Skipping", interaction_example, "- not found in x@net\n")
      next
    }
    interaction_matrix <- x@net[[interaction_example]]

    # Check if there's variation in the data
    vals <- as.vector(interaction_matrix)
    vals <- vals[is.finite(vals)]  # drop NA/Inf
    if (length(vals) == 0) {
      cat("Skipping", interaction_example, "- empty / all NA\n")
      next
    }

    if (length(unique(vals)) > 1) {
      out_file <- file.path(fig_dir, paste0(safe_name(interaction_example), "_lig_tar_heatmap_", sample_name, ".png"))

      png(out_file, width = 12, height = 8, units = "in", res = 300)
      tryCatch({
        lig_tar_heatmap(
          x,
          interaction_name = interaction_example,
          width.vector = c(0.38, 0.35, 0.27)
        )
      }, error = function(e) {
        cat("Failed for", interaction_example, ":", conditionMessage(e), "\n")
      })
      dev.off()

      cat("Created heatmap for:", interaction_example, "->", out_file, "\n")
    } else {
      cat("Skipping", interaction_example, "- no variation in values\n")
    }
  }
} else {
  cat("No interaction_names found.\n")
}




active_paths



# ==============================================================================
# SAVE RESULTS
# ==============================================================================

# Save RDS
saveRDS(x, file = paste0(sample_name, "_neuronchat_object_final.rds"))

# Export interaction names
writeLines(interaction_names, paste0("active_interactions_", sample_name, ".txt"))

# Export aggregated network
write.csv(net_aggregated_x, paste0("aggregated_network_", sample_name, ".csv"))

cat("Analysis complete!\n")


