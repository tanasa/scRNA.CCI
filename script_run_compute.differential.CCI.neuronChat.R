# ============================================================================
# NeuronChat Comparison Script - Version 17
# ============================================================================
# 
# PURPOSE: Compare cell-cell communication networks between two conditions
# (AD_CTL vs AD_Dyslexia) using NeuronChat analysis
#

# Key sections:
# 1. Part I: Individual pathway visualization
# 2. Part II: Merge NeuronChat objects for comparison
# 3. Part III: Compare overall interaction counts and weights
# 4. Part IV: Functional similarity analysis (UMAP embedding)
# 5. Part V: Differential pathway analysis (fold-changes, scatter plots)


setwd("/mnt/nfs/CX000008_DS1/projects/jaeyeon/fastq_file_Dyslexia_r1/zanalysis_bogdan/samples.merged_AG_Harmony_res0.1_anno_091625")

# Make sure reticulate can find Python
library(reticulate)
py_config()  # Check your Python configuration
library("uwot")

# Load libraries
library(NeuronChat)
library(CellChat)

library(grid)
library(pheatmap)
library(ComplexHeatmap)
library(RColorBrewer)
library(gplots)  # For heatmap.2 (better legend visibility)

set.seed(1234)



# Part I: Load NeuronChat objects

# Load pre-computed NeuronChat objects
adctl <- readRDS("merged_AG_Harmony_res0.1_anno_091625.AD_CTL_neuronchat_object_final.rds")
addys <- readRDS("merged_AG_Harmony_res0.1_anno_091625.AD_Dyslexia_neuronchat_object_final.rds")

# Create list with named objects
cortex_list <- list(AD_CTL = adctl, AD_Dyslexia = addys)

head(adctl@data.signaling, 2)

head(addys@data.signaling, 2)

# cortex_list
# str(glimpse(cortex_list))

# Individually plotting the communication networks for AD_Dyslexia and AD_CTL  

# Rows = sender cell types
# Cells that produce / release the signal (e.g. serotonin-producing neurons)

# Columns = receiver cell types
# Cells that express the receptor (e.g. Htr1a-expressing neurons)

# Values = communication weight

# A continuous score inferred from:

# ligand contributors (e.g. synthesis + vesicular transport)
# receptor expression
# stoichiometry rules
# cell-type–specific expression

# length(names(cortex_list[[1]]@net))
# cortex_list[[1]]@net

# length(cortex_list[[2]]@net)
# cortex_list[[2]]@net

# The matrix represents cell–cell communication strength for the serotonin (5HT) → Htr1a receptor interaction.

# Rows = sender cell types
# Cells that produce / release the signal (e.g. serotonin-producing neurons)

# Columns = receiver cell types
# Cells that express the receptor (e.g. Htr1a-expressing neurons)

# Values = communication weight

# A continuous score inferred from:
# ligand contributors (e.g. synthesis + vesicular transport)
# receptor expression
# stoichiometry rules
# cell-type–specific expression

identical(names(cortex_list[[1]]@net), names(cortex_list[[2]]@net))

all_pathways <- union(names(cortex_list[[1]]@net), names(cortex_list[[2]]@net))
length(all_pathways)
print(all_pathways)



# ============================================================================
# NeuronChat Comparison Script - Version 13
# ============================================================================

# Create output directory for all results
output_dir <- "neuronchat_diff_comp"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
cat("Created output directory:", output_dir, "\n")

# ============================================================================
# PART I: INDIVIDUAL PATHWAY VISUALIZATION
# ============================================================================
# Create output directory for all results - stores all plots and data files
output_dir <- "neuronchat_diff_comp"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
cat("Created output directory:", output_dir, "\n")

# Generate side-by-side comparison plots for each pathway
# Creates visual network diagrams showing cell-cell communication patterns
for(pathway in all_pathways) {
  
  # Create filename starting with pathway name
  filename <- file.path(output_dir, paste0(pathway, "_comparison.png"))
  
  # Open PNG device - high resolution for publication quality
  png(filename, width = 1400, height = 700, res = 100)
  
  # Set up side-by-side plotting - left panel AD_CTL, right panel AD_Dyslexia
  par(mfrow=c(1,2))
  
  # Loop through both datasets (AD_CTL and AD_Dyslexia)
  for(j in c(1,2)) {
    
    # Check if pathway exists in this condition - some pathways are condition-specific
    if(pathway %in% names(cortex_list[[j]]@net)) {
      
      # Plot the network - circular diagram showing sender->receiver communication
      netVisual_circle_neuron(cortex_list[[j]]@net[[pathway]], 
                              title.name = paste(pathway, '--', names(cortex_list)[j]), 
                              arrow.size = 0.5, 
                              margin = 0.3, 
                              edge.width.max = 8)
    } else {
      # If pathway doesn't exist in this region, create empty plot with message
      plot.new()
      title(main = paste(pathway, '--', names(cortex_list)[j], '\n(Not detected)'))
    }
  }
  
  # Close the PNG device - saves the plot to file
  dev.off()
  
  # Optional: print progress
  cat("Saved:", filename, "\n")
}

# Reset plotting parameters
par(mfrow=c(1,1))

cat("All plots saved! Total files:", length(all_pathways), "\n")

# An example :

par(mfrow=c(1,2))

for(j in c(1,2)){
  netVisual_circle_neuron(cortex_list[[j]]@net[["NRXN1_NLGN1"]], 
                          title.name = paste(names(cortex_list)[j]), 
                          arrow.size = 0.5, 
                          margin=0.3, 
                          edge.width.max=8
                         )
}

for(j in c(1,2)){
  netVisual_circle_neuron(cortex_list[[j]]@net[["NRXN1_NLGN3"]], 
                          title.name = paste(names(cortex_list)[j]), 
                          arrow.size = 0.5, 
                          margin=0.3, 
                          edge.width.max=8
                         )
}



# ============================================================================
# NETWORK AGGREGATION: Combine all pathways into single communication matrix
# ============================================================================
# Aggregates all individual pathways (neurotransmitters, neuropeptides, etc.)
# into one matrix showing total communication strength between cell types

par(mfrow=c(1,2))

for(j in c(1,2)){
  # Aggregate by counting number of interactions (method='count')
  # Shows how many different pathways connect each cell type pair
  net_aggregated_x <- net_aggregation(cortex_list[[j]]@net, 
                                      method='count')
    
  netVisual_circle_neuron(net_aggregated_x, 
                          title.name = names(cortex_list)[j], 
                          arrow.size = 0.5, 
                          margin=0.3, 
                          edge.width.max=8)
}

# Save the aggregated network comparison (count method)
png(file.path(output_dir, "fig_aggregated_network_count_comparison.png"), width = 1400, height = 700, res = 100)

par(mfrow=c(1,2))
for(j in c(1,2)){
  net_aggregated_x <- net_aggregation(cortex_list[[j]]@net, 
                                      method='count')
  netVisual_circle_neuron(net_aggregated_x, 
                          title.name = names(cortex_list)[j], 
                          arrow.size = 0.5, 
                          margin=0.3, 
                          edge.width.max=8)
}

dev.off()

# Reset plotting parameters
par(mfrow=c(1,1))

cat("Saved: aggregated_network_count_comparison.png\n")

par(mfrow=c(1,2))
for(j in c(1,2)){
  net_aggregated_x <- net_aggregation(cortex_list[[j]]@net, 
                                      method='weight')
  netVisual_circle_neuron(net_aggregated_x, 
                          title.name = names(cortex_list)[j], 
                          arrow.size = 0.5, 
                          margin=0.3, 
                          edge.width.max=8)
}

# Save the aggregated network comparison (weight method)
png(file.path(output_dir, "fig_aggregated_network_weight_comparison.png"), width = 1400, height = 700, res = 100)

par(mfrow=c(1,2))
for(j in c(1,2)){
  net_aggregated_x <- net_aggregation(cortex_list[[j]]@net, 
                                      method='weight')
  netVisual_circle_neuron(net_aggregated_x, 
                          title.name = names(cortex_list)[j], 
                          arrow.size = 0.5, 
                          margin=0.3, 
                          edge.width.max=8)
}

dev.off()

# Reset plotting parameters
par(mfrow=c(1,1))

cat("Saved: aggregated_network_weight_comparison.png\n")




# Save both aggregation methods
methods <- c('count', 'weight')

for(method in methods) {
  
  filename <- file.path(output_dir, paste0("fig_aggregated_network_", method, "_comparison.png"))
  
  png(filename, width = 1400, height = 700, res = 100)
  
  par(mfrow=c(1,2))
  for(j in c(1,2)){
    net_aggregated_x <- net_aggregation(cortex_list[[j]]@net, 
                                        method=method)
      
    netVisual_circle_neuron(net_aggregated_x, 
                            title.name = names(cortex_list)[j], 
                            arrow.size = 0.5, 
                            margin=0.3, 
                            edge.width.max=8)
  }
  
  dev.off()
  
  cat("Saved:", filename, "\n")
}


# ?net_aggregation 
# Aggregation of communication networks over all ligand-target pairs 

# How to read the matrix

# Rows = sender cell types
# Cells that send signal
# Columns = receiver cell types
# Cells that receive signals

# Values = total communication strength
# A scalar summarizing:
# neurotransmitter signaling
# neuropeptide signaling
# electrical synapses (gap junctions)
# other supported NeuronChat programs

# net_aggregated_x is an aggregated cell–cell communication matrix, where all individual signaling programs 
# (neurotransmitters, neuropeptides, gap junctions, etc.) have been summed together.

# In one sentence:

# It represents the total inferred communication strength from each sender cell type to each receiver cell type, 
# across all ligand–receptor programs.

# ?mergeNeuronChat # Merge NeuronChat objects Merge NeuronChat objects; adapted from CellChat https://github.com/sqjin/CellChat 

names(cortex_list)

head(cortex_list[[1]]@meta, 2)

head(cortex_list[[2]]@meta, 2)




# Part II: Merge the NeuronChat object list

# Fix: Ensure unique row names in meta slots before merging
# This prevents "invalid 'row.names' length" error

# cortex_list[[1]]@meta

head(colnames(cortex_list[[1]]@data))

# Fix: Populate the meta slot for both objects before merging
# The meta slot needs to contain cell-level metadata as a data.frame

cat("Populating meta slot for both datasets...\n")

# For AD_CTL
if(is.null(cortex_list[[1]]@meta) || nrow(cortex_list[[1]]@meta) == 0) {
  cat("Creating meta for AD_CTL...\n")
  
  # Get cell barcodes from @data slot
  cell_barcodes <- colnames(cortex_list[[1]]@data)
  
  # Create a basic metadata data.frame
  cortex_list[[1]]@meta <- data.frame(
    cell_barcode = cell_barcodes,
    dataset = "AD_CTL",
    sample = "AD_CTL",
    row.names = cell_barcodes,
    stringsAsFactors = FALSE
  )
  
  cat("  AD_CTL meta: ", nrow(cortex_list[[1]]@meta), " cells\n")
}

# cortex_list[[2]]@meta

head(colnames(cortex_list[[2]]@data))

# Part II: Merge the NeuronChat object list

# Fix: Ensure unique row names in meta slots before merging
# This prevents "invalid 'row.names' length" error

# For AD_Dyslexia
if(is.null(cortex_list[[2]]@meta) || nrow(cortex_list[[2]]@meta) == 0) {
  cat("Creating meta for AD_Dyslexia...\n")
  
  # Get cell barcodes from @data slot
  cell_barcodes <- colnames(cortex_list[[2]]@data)
  
  # Create a basic metadata data.frame
  cortex_list[[2]]@meta <- data.frame(
    cell_barcode = cell_barcodes,
    dataset = "AD_Dyslexia",
    sample = "AD_Dyslexia",
    row.names = cell_barcodes,
    stringsAsFactors = FALSE
  )
  
  cat("  AD_Dyslexia meta: ", nrow(cortex_list[[2]]@meta), " cells\n")
}




# Part II: Merge the NeuronChat object list

neuronchat_list <- mergeNeuronChat(cortex_list, 
                                   add.names = names(cortex_list))

# neuronchat_list
glimpse(neuronchat_list, max.level = 2)

# neuronchat_list@net$AD_CTL
# neuronchat_list@net$AD_Dyslexia

# str(cortex_list[1]$AD_CTL)
# str(cortex_list[1]$AD_CTL@data)
# str(cortex_list[1]$AD_CTL@data.signaling)

# str(cortex_list[1]$AD_DYS)
# str(cortex_list[1]$AD_CTL@data)
# str(cortex_list[1]$AD_CTL@data.signaling)

# neuronchat_list@net$AD_CTL
# neuronchat_list@net$AD_Dyslexia

glimpse(neuronchat_list@meta)
dim(neuronchat_list@meta)


head(neuronchat_list@meta, 2)
tail(neuronchat_list@meta, 2)





# Part III: Barplots to compare link count and weight of interaction pairs between AD_CTL and AD_Dyslexia

# Overall communication

p1 <- compareInteractions_Neuron(neuronchat_list, measure = c("count"), comparison = c(1,2), group=c(1,2), show.legend = F)
p2 <- compareInteractions_Neuron(neuronchat_list, measure = c("weight"), comparison = c(1,2), group=c(1,2), show.legend = F)

p1 + p2

# Create the plots
p1 <- compareInteractions_Neuron(neuronchat_list, 
                                 measure = c("count"), 
                                 comparison = c(1,2), 
                                 group=c(1,2), 
                                 show.legend = F)
p2 <- compareInteractions_Neuron(neuronchat_list, 
                                 measure = c("weight"), 
                                 comparison = c(1,2), 
                                 group=c(1,2), 
                                 show.legend = F)

# Combine and save
library(ggplot2)  # Make sure patchwork or cowplot is loaded for the + operator

combined_plot <- p1 + p2

ggsave(file.path(output_dir, "fig_interaction_comparison_count_weight.png"), 
       plot = combined_plot, 
       width = 10, 
       height = 5, 
       dpi = 300)

cat("Saved: interaction_comparison_count_weight.png\n")



# Communication for individual interaction pairs

# ?rankNet_Neuron 
# Rank signaling networks based on the information flow or the number of interactions; 
# adapted from CellChat https://github.com/sqjin/CellChat

# INFORMATION FLOW

# Information flow is a summary score of how strong a signaling pathway is across the whole network.
# In plain terms:
# Information flow = total communication strength carried by a signaling pathway across all sender–receiver cell-type pairs.

# For a given signaling pathway (e.g. Glu_Grin3a, Wnt, Notch):
# NeuronChat / CellChat builds a cell-type × cell-type communication matrix

# rows = sender cell types
# columns = receiver cell types
# values = communication strength (probability / weight)

# Information flow = sum of all edge weights in that matrix:

# | Metric                     | What it measures        | Interpretation                          |
# | -------------------------- | ----------------------- | --------------------------------------- |
# | **Number of interactions** | Count of non-zero edges | How widespread a pathway is             |
# | **Information flow**       | Sum of edge weights     | How dominant / influential a pathway is |
# | **High count, low flow**   | Many weak signals       | Broad but weak                          |
# | **Low count, high flow**   | Few strong signals      | Focused but powerful                    |



options(repr.plot.width = 14, repr.plot.height = 16)

g1 <- rankNet_Neuron(neuronchat_list, mode='comparison', 
                     measure = c("count"), 
                     comparison = 1:2, do.stat = F, tol = 0.1, stacked = F, font.size = 11)
g2 <- rankNet_Neuron(neuronchat_list, mode='comparison', 
                     measure = c("weight"), 
                     comparison = 1:2, do.stat = F, tol = 0.1, stacked = F, font.size = 11) 
g1+g2

# Combine and save
combined_plot <- g1 + g2

ggsave(file.path(output_dir, "fig_pathway_ranking_count_weight.png"), 
       plot = combined_plot, 
       width = 14, 
       height = 16, 
       dpi = 300)

cat("Saved: pathway_ranking_count_weight.png\n")


# the differences in the Information Flow : 

# NRXN1_NLGN4Y
# NRXN3_NLGN4Y
# NRXN3_NLGN1
# NRXN1_NLGN1
# NRXN1_NLGN4X
# NRXN3_NLGN4X


# ?computeNetSimilarityPairwise_Neuron
# Compute signaling network similarity for any pair of datasets
# Description
# Compute signaling network similarity for any pair of datasets; adapted from CellChat https://github.com/sqjin/CellChat 

# rankNet_Neuron accesses:
# neuronchat_list@net

# Structure after merging:
# neuronchat_list@net$AD_CTL    # Contains all pathway matrices for AD_CTL
# neuronchat_list@net$AD_Dyslexia   # Contains all pathway matrices for AD_Dyslexia


# Each pathway is an 8×8 matrix, for example:

# neuronchat_list@net[["AD_CTL"]][["5HT_HTR1A"]]
# neuronchat_list@net[["AD_Dyslexia"]][["5HT_HTR1A"]]

neuronchat_list@net[["AD_CTL"]][["NRXN3_NLGN1"]]
neuronchat_list@net[["AD_Dyslexia"]][["NRXN3_NLGN1"]]

glimpse(neuronchat_list@net, max.level =1)



# Part IV: Shared and specific interaction patterns across AD_CTL and AD_DYS

# Big picture

# When NeuronChat compares two communication networks, it can ask two different questions:
# Do they use the same signaling pathways? → functional similarity      : FUNCTIONAL
# Do they connect cell types in the same way? → structural similarity   : STRUCTURAL

# Old notes - neuronChat tutorial :

cat("

1️⃣ Functional similarity

Functional similarity compares networks based on signaling identity.

In plain terms:

Are the same signaling pathways (ligand–receptor programs) active in both networks, regardless of which cell types talk to each other?

What is compared

Presence/strength of pathways (e.g., Glu, GABA, Wnt, Notch, cytokines)

Aggregated signaling activity per pathway

What is ignored

Exact sender → receiver cell-type wiring

Interpretation

High functional similarity means:

The two regions/conditions rely on similar biological signaling programs

Even if different cell types are responsible

Example

ALM and VISp both show strong glutamatergic + calcium signaling

But ALM uses L2/3 IT → L5 IT

VISp uses L4 IT → L2/3 IT

➡ High functional similarity, low structural similarity

")

# Old notes - neuronChat tutorial :

cat("

2️⃣ Structural similarity

Structural similarity compares the wiring diagram of the network.

In plain terms:

Do the same sender cell types communicate with the same receiver cell types, with similar strengths?

What is compared

Sender × receiver adjacency matrices

Directionality and weight of interactions

What is ignored

Identity of signaling pathways

Interpretation

High structural similarity means:

Cell-type communication architecture is preserved

Even if signaling molecules differ

Example

ALM and VISp both have:

L2/3 IT → L5 IT

L5 IT → interneurons

But one uses Glu signaling, the other uses neuromodulators

➡ High structural similarity, low functional similarity

")




# neuronchat_list@net$AD_CTL

# neuronchat_list@net$AD_Dyslexia

# Manifold learning

# Functional similarity

# Compute functional similarity 
neuronchat_list <- computeNetSimilarityPairwise_Neuron(neuronchat_list, 
                                                       slot.name = "net", 
                                                       type = "functional", 
                                                       comparison = c(1,2))

glimpse(neuronchat_list@net_analysis)
glimpse(neuronchat_list@net_analysis$similarity$functional$matrix)
glimpse(neuronchat_list@net_analysis$similarity$functional$matrix$"1-2")
# glimpse(neuronchat_list)

# reticulate::py_install("umap-learn", pip = TRUE)
neuronchat_list <- netEmbedding(neuronchat_list, 
                                slot.name = "net_analysis", 
                                type = "functional",
                                umap.method = "uwot",
                                comparison = c(1,2))

# Clustering on interactions 
neuronchat_list <- netClustering(neuronchat_list, 
                                 slot.name = "net_analysis", 
                                 type = "functional", 
                                 comparison = c(1,2), k = 5)

# Visualization
netVisual_embeddingPairwise_Neuron(neuronchat_list, 
                                   slot.name = "net_analysis", 
                                   type = "functional", 
                                   label.size = 3.5, 
                                   comparison=c(1,2), 
                                   pathway.remove.show = FALSE, 
                                   pathway.labeled = F)

# Visualization zoom in
netVisual_embeddingPairwiseZoomIn_Neuron(neuronchat_list, 
                                         slot.name = "net_analysis", 
                                         type = "functional", 
                                         label.size = 5, 
                                         comparison=c(1,2), 
                                         nCol=3)



# If the ○ and △ for a pathway are close together → similar behavior
# If they are far apart → region-specific regulation

# UMAP distance reflects:

# differences in overall pathway activity
# differences in information flow
# differences in sender–receiver usage aggregated at pathway level

glimpse(neuronchat_list@net_analysis$similarity)
glimpse(neuronchat_list@net_analysis$similarity$functional$matrix)

# neuronchat_list@net_analysis$similarity$functional$matrix

matrix_data <- neuronchat_list@net_analysis$similarity$functional$matrix$`1-2`
head(matrix_data)
dim(matrix_data)


# Diagonal values = 1.0 (each condition is perfectly similar to itself)
# Off-diagonal values = similarity score between AD_CTL and AD_Dyslexia (range 0-1)

# Values closer to 1 = high functional similarity
# Values closer to 0 = low functional similarity

# Save the FULL matrix without filtering (no sparsity removal)

# Check dimensions
cat("Matrix dimensions:", dim(matrix_data), "\n")

# Plot the full matrix
options(repr.plot.width = 60, repr.plot.height = 60, repr.plot.res = 100)
pheatmap(matrix_data,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(c("white", "red"))(100),
         border_color = NA,
         fontsize_row = 12,
         fontsize_col = 12,
         main = "Functional Similarity: AD_CTL vs AD_Dyslexia (Full Matrix)",
         show_colnames = TRUE,
         show_rownames = TRUE)

# Save the full heatmap
png(file.path(output_dir, "fig_similarity_heatmap_full.png"), width = 6000, height = 6000, res = 100)
pheatmap(matrix_data,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(c("white", "red"))(100),
         border_color = NA,
         fontsize_row = 12,
         fontsize_col = 12,
         main = "Functional Similarity: AD_CTL vs AD_Dyslexia (Full Matrix)",
         show_colnames = TRUE,
         show_rownames = TRUE)
dev.off()

cat("Saved: similarity_heatmap_full.png\n")


# To remove the sparsity: 
# Keep only rows/columns with MORE than just the diagonal
# (i.e., rows with at least 2 non-zero values)
rows_to_keep <- rowSums(matrix_data != 0) > 1
cols_to_keep <- colSums(matrix_data != 0) > 1

# Subset the matrix
matrix_filtered <- matrix_data[rows_to_keep, cols_to_keep]

# Check dimensions
cat("Original matrix:", dim(matrix_data), "\n")
cat("Filtered matrix:", dim(matrix_filtered), "\n")
cat("Rows removed:", sum(!rows_to_keep), "\n")
cat("Cols removed:", sum(!cols_to_keep), "\n")

# Plot the filtered matrix
options(repr.plot.width = 60, repr.plot.height = 60, repr.plot.res = 100)
pheatmap(matrix_filtered,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(c("white", "red"))(100),
         border_color = NA,
         fontsize_row = 12,
         fontsize_col = 12,
         main = "Functional Similarity: AD_CTL vs AD_Dyslexia (Pathways with cross-similarities)",
         show_colnames = TRUE,
         show_rownames = TRUE)

# Save the filtered heatmap
png(file.path(output_dir, "fig_similarity_heatmap_filtered.png"), width = 3000, height = 3000, res = 100)
pheatmap(matrix_filtered,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(c("white", "red"))(100),
         border_color = NA,
         fontsize_row = 12,
         fontsize_col = 12,
         main = "Functional Similarity: AD_CTL vs AD_Dyslexia (Pathways with cross-similarities)",
         show_colnames = TRUE,
         show_rownames = TRUE)
dev.off()

cat("Saved: similarity_heatmap_filtered.png\n")





# What the clustering shows:

# This matrix represents functional similarity between cell type pairs across two groups (AD_CTL and AD_Dyslexia).
# Here's what it shows:

# Matrix Structure:

# Rows: Cell type pairs from AD_CTL (e.g., "5HT_Htr1a--AD_CTL", "5HT_Htr1b--AD_CTL")
# Columns: Cell type pairs from both AD_CTL and AD_Dyslexia
# Values: Binary (0 or 1) indicating whether the functional interactions are similar

# What the values mean:

# 1 = The two cell type pairs have similar functional interactions (similar ligand-receptor communication patterns)
# 0 = The cell type pairs do not have similar functional interactions

# Interpretation:

# This is a cross-group comparison matrix that shows which cell-cell communication patterns in AD_CTL are functionally similar to patterns in AD_Dyslexia. For example:

# The diagonal pattern of 1s suggests that each cell type pair is similar to itself
# Off-diagonal 1s would indicate that different cell type pairs share similar functional interaction profiles

# ?computeNetSimilarityPairwise_Neuron

# Row clustering (vertical dendrogram):

# Groups pathways from ALM that have similar functional profiles
# Pathways that cluster together have correlated communication patterns across cell types

# Column clustering (horizontal dendrogram):

# Groups pathways from both ALM and VISp based on their similarity patterns
# Reveals which pathways behave similarly within and across regions

# Biological interpretations:

# 1. Conserved signaling modules (diagonal blocks of red)

# Pathways that are functionally similar in both ALM and VISp
# Indicates core communication programs that are preserved across cortical regions
# Example: If glutamate receptors cluster together, they share similar cell-type expression patterns in both regions

# 2. Region-specific modules (off-diagonal patterns)

# ALM pathways that are more similar to certain VISp pathways than to their own ALM counterparts
# Suggests functional rewiring or specialization between regions





# Part II: Merge the NeuronChat object list

# Fix: Ensure unique row names in meta slots before merging
# This prevents "invalid 'row.names' length" error

# cortex_list[[1]]@meta

head(colnames(cortex_list[[1]]@data))

# Fix: Populate the meta slot for both objects before merging
# The meta slot needs to contain cell-level metadata as a data.frame

cat("Populating meta slot for both datasets...\n")

# For AD_CTL
if(is.null(cortex_list[[1]]@meta) || nrow(cortex_list[[1]]@meta) == 0) {
  cat("Creating meta for AD_CTL...\n")
  
  # Get cell barcodes from @data slot
  cell_barcodes <- colnames(cortex_list[[1]]@data)
  
  # Create a basic metadata data.frame
  cortex_list[[1]]@meta <- data.frame(
    cell_barcode = cell_barcodes,
    dataset = "AD_CTL",
    sample = "AD_CTL",
    row.names = cell_barcodes,
    stringsAsFactors = FALSE
  )
  
  cat("  AD_CTL meta: ", nrow(cortex_list[[1]]@meta), " cells\n")
}



# Part II: Merge the NeuronChat object list

# Fix: Ensure unique row names in meta slots before merging
# This prevents "invalid 'row.names' length" error

# For AD_Dyslexia
if(is.null(cortex_list[[2]]@meta) || nrow(cortex_list[[2]]@meta) == 0) {
  cat("Creating meta for AD_Dyslexia...\n")
  
  # Get cell barcodes from @data slot
  cell_barcodes <- colnames(cortex_list[[2]]@data)
  
  # Create a basic metadata data.frame
  cortex_list[[2]]@meta <- data.frame(
    cell_barcode = cell_barcodes,
    dataset = "AD_Dyslexia",
    sample = "AD_Dyslexia",
    row.names = cell_barcodes,
    stringsAsFactors = FALSE
  )
  
  cat("  AD_Dyslexia meta: ", nrow(cortex_list[[2]]@meta), " cells\n")
}


head(cortex_list[[1]]@meta)

head(cortex_list[[2]]@meta)

# Part II: Merge the NeuronChat object list

neuronchat_list <- mergeNeuronChat(cortex_list, 
                                   add.names = names(cortex_list))

# neuronchat_list
glimpse(neuronchat_list, max.level = 2)

# neuronchat_list@net$AD_CTL
# neuronchat_list@net$AD_Dyslexia

# str(cortex_list[1]$AD_CTL)
# str(cortex_list[1]$AD_CTL@data)
# str(cortex_list[1]$AD_CTL@data.signaling)

# str(cortex_list[1]$AD_DYS)
# str(cortex_list[1]$AD_CTL@data)
# str(cortex_list[1]$AD_CTL@data.signaling)

# neuronchat_list@net$AD_CTL
# neuronchat_list@net$AD_Dyslexia

glimpse(neuronchat_list@meta)
dim(neuronchat_list@meta)


head(neuronchat_list@meta, 2)
tail(neuronchat_list@meta, 2)




# Part III: Barplots to compare link count and weight of interaction pairs between AD_CTL and AD_Dyslexia

# Overall communication

p1 <- compareInteractions_Neuron(neuronchat_list, measure = c("count"), comparison = c(1,2), group=c(1,2), show.legend = F)
p2 <- compareInteractions_Neuron(neuronchat_list, measure = c("weight"), comparison = c(1,2), group=c(1,2), show.legend = F)

p1 + p2

# Create the plots
p1 <- compareInteractions_Neuron(neuronchat_list, 
                                 measure = c("count"), 
                                 comparison = c(1,2), 
                                 group=c(1,2), 
                                 show.legend = F)
p2 <- compareInteractions_Neuron(neuronchat_list, 
                                 measure = c("weight"), 
                                 comparison = c(1,2), 
                                 group=c(1,2), 
                                 show.legend = F)

# Combine and save
library(ggplot2)  # Make sure patchwork or cowplot is loaded for the + operator

combined_plot <- p1 + p2

ggsave(file.path(output_dir, "fig_interaction_comparison_count_weight.png"), 
       plot = combined_plot, 
       width = 10, 
       height = 5, 
       dpi = 300)

cat("Saved: interaction_comparison_count_weight.png\n")









# Communication for individual interaction pairs

# ?rankNet_Neuron 
# Rank signaling networks based on the information flow or the number of interactions; 
# adapted from CellChat https://github.com/sqjin/CellChat

# INFORMATION FLOW

# Information flow is a summary score of how strong a signaling pathway is across the whole network.
# In plain terms:
# Information flow = total communication strength carried by a signaling pathway across all sender–receiver cell-type pairs.

# For a given signaling pathway (e.g. Glu_Grin3a, Wnt, Notch):
# NeuronChat / CellChat builds a cell-type × cell-type communication matrix

# rows = sender cell types
# columns = receiver cell types
# values = communication strength (probability / weight)

# Information flow = sum of all edge weights in that matrix:

# | Metric                     | What it measures        | Interpretation                          |
# | -------------------------- | ----------------------- | --------------------------------------- |
# | **Number of interactions** | Count of non-zero edges | How widespread a pathway is             |
# | **Information flow**       | Sum of edge weights     | How dominant / influential a pathway is |
# | **High count, low flow**   | Many weak signals       | Broad but weak                          |
# | **Low count, high flow**   | Few strong signals      | Focused but powerful                    |




options(repr.plot.width = 14, repr.plot.height = 16)

g1 <- rankNet_Neuron(neuronchat_list, mode='comparison', 
                     measure = c("count"), 
                     comparison = 1:2, do.stat = F, tol = 0.1, stacked = F, font.size = 11)
g2 <- rankNet_Neuron(neuronchat_list, mode='comparison', 
                     measure = c("weight"), 
                     comparison = 1:2, do.stat = F, tol = 0.1, stacked = F, font.size = 11) 
g1+g2

# Combine and save
combined_plot <- g1 + g2

ggsave(file.path(output_dir, "fig_pathway_ranking_count_weight.png"), 
       plot = combined_plot, 
       width = 14, 
       height = 16, 
       dpi = 300)

cat("Saved: pathway_ranking_count_weight.png\n")

# Some pathwways to look into !!!

# NRXN1_NLGN4Y
# NRXN3_NLGN4Y
# NRXN3_NLGN1
# NRXN1_NLGN1
# NRXN1_NLGN4X
# NRXN3_NLGN4X



# ?computeNetSimilarityPairwise_Neuron
# Compute signaling network similarity for any pair of datasets
# Description
# Compute signaling network similarity for any pair of datasets; adapted from CellChat https://github.com/sqjin/CellChat 

# rankNet_Neuron accesses:
# neuronchat_list@net

# Structure after merging:
# neuronchat_list@net$AD_CTL    # Contains all pathway matrices for AD_CTL
# neuronchat_list@net$AD_Dyslexia   # Contains all pathway matrices for AD_Dyslexia




# Each pathway is an 8×8 matrix, for example:

neuronchat_list@net[["AD_CTL"]][["NRXN3_NLGN1"]]
neuronchat_list@net[["AD_Dyslexia"]][["NRXN3_NLGN1"]]

glimpse(neuronchat_list@net, max.level =1)

# neuronchat_list@net$AD_Ctl
# neuronchat_list@net$AD_Dyslexia






# Part IV: Shared and specific interaction patterns across AD_CTL and AD_DYS

# Big picture

# When NeuronChat compares two communication networks, it can ask two different questions:
# Do they use the same signaling pathways? → functional similarity      : FUNCTIONAL
# Do they connect cell types in the same way? → structural similarity   : STRUCTURAL

# Old notes - neuronChat tutorial :

cat("

1️⃣ Functional similarity

Functional similarity compares networks based on signaling identity.

In plain terms:

Are the same signaling pathways (ligand–receptor programs) active in both networks, regardless of which cell types talk to each other?

What is compared

Presence/strength of pathways (e.g., Glu, GABA, Wnt, Notch, cytokines)

Aggregated signaling activity per pathway

What is ignored

Exact sender → receiver cell-type wiring

Interpretation

High functional similarity means:

The two regions/conditions rely on similar biological signaling programs

Even if different cell types are responsible

Example

ALM and VISp both show strong glutamatergic + calcium signaling

But ALM uses L2/3 IT → L5 IT

VISp uses L4 IT → L2/3 IT

➡ High functional similarity, low structural similarity

")

# Old notes - neuronChat tutorial :

cat("

2️⃣ Structural similarity

Structural similarity compares the wiring diagram of the network.

In plain terms:

Do the same sender cell types communicate with the same receiver cell types, with similar strengths?

What is compared

Sender × receiver adjacency matrices

Directionality and weight of interactions

What is ignored

Identity of signaling pathways

Interpretation

High structural similarity means:

Cell-type communication architecture is preserved

Even if signaling molecules differ

Example

ALM and VISp both have:

L2/3 IT → L5 IT

L5 IT → interneurons

But one uses Glu signaling, the other uses neuromodulators

➡ High structural similarity, low functional similarity

")




# Manifold learning

# Functional similarity

# Compute functional similarity 
neuronchat_list <- computeNetSimilarityPairwise_Neuron(neuronchat_list, 
                                                       slot.name = "net", 
                                                       type = "functional", 
                                                       comparison = c(1,2))

glimpse(neuronchat_list@net_analysis)
glimpse(neuronchat_list@net_analysis$similarity$functional$matrix)
glimpse(neuronchat_list@net_analysis$similarity$functional$matrix$"1-2")
# glimpse(neuronchat_list)

# reticulate::py_install("umap-learn", pip = TRUE)
neuronchat_list <- netEmbedding(neuronchat_list, 
                                slot.name = "net_analysis", 
                                type = "functional",
                                umap.method = "uwot",
                                comparison = c(1,2))

# Clustering on interactions 
neuronchat_list <- netClustering(neuronchat_list, 
                                 slot.name = "net_analysis", 
                                 type = "functional", 
                                 comparison = c(1,2), k = 5)

# Visualization
netVisual_embeddingPairwise_Neuron(neuronchat_list, 
                                   slot.name = "net_analysis", 
                                   type = "functional", 
                                   label.size = 3.5, 
                                   comparison=c(1,2), 
                                   pathway.remove.show = FALSE, 
                                   pathway.labeled = F)

# Visualization zoom in
netVisual_embeddingPairwiseZoomIn_Neuron(neuronchat_list, 
                                         slot.name = "net_analysis", 
                                         type = "functional", 
                                         label.size = 5, 
                                         comparison=c(1,2), 
                                         nCol=3)

# If the ○ and △ for a pathway are close together → similar behavior
# If they are far apart → region-specific regulation

# UMAP distance reflects:

# differences in overall pathway activity
# differences in information flow
# differences in sender–receiver usage aggregated at pathway level



glimpse(neuronchat_list@net_analysis$similarity)

glimpse(neuronchat_list@net_analysis$similarity$functional$matrix)

neuronchat_list@net_analysis$similarity$functional$matrix

matrix_data <- neuronchat_list@net_analysis$similarity$functional$matrix$`1-2`
# head(matrix_data)
dim(matrix_data)



# Diagonal values = 1.0 (each condition is perfectly similar to itself)
# Off-diagonal values = similarity score between AD_CTL and AD_Dyslexia (range 0-1)

# Values closer to 1 = high functional similarity
# Values closer to 0 = low functional similarity

# Save the FULL matrix without filtering (no sparsity removal)

# Check dimensions
cat("Matrix dimensions:", dim(matrix_data), "\n")

# Plot the full matrix
options(repr.plot.width = 60, repr.plot.height = 60, repr.plot.res = 100)
pheatmap(matrix_data,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(c("white", "red"))(100),
         border_color = NA,
         fontsize_row = 12,
         fontsize_col = 12,
         main = "Functional Similarity: AD_CTL vs AD_Dyslexia (Full Matrix)",
         show_colnames = TRUE,
         show_rownames = TRUE)

# Save the full heatmap
png(file.path(output_dir, "fig_similarity_heatmap_full.png"), width = 6000, height = 6000, res = 100)
pheatmap(matrix_data,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(c("white", "red"))(100),
         border_color = NA,
         fontsize_row = 12,
         fontsize_col = 12,
         main = "Functional Similarity: AD_CTL vs AD_Dyslexia (Full Matrix)",
         show_colnames = TRUE,
         show_rownames = TRUE)
dev.off()

cat("Saved: similarity_heatmap_full.png\n")






# To remove the sparsity: 
# Keep only rows/columns with MORE than just the diagonal
# (i.e., rows with at least 2 non-zero values)
rows_to_keep <- rowSums(matrix_data != 0) > 1
cols_to_keep <- colSums(matrix_data != 0) > 1

# Subset the matrix
matrix_filtered <- matrix_data[rows_to_keep, cols_to_keep]

# Check dimensions
cat("Original matrix:", dim(matrix_data), "\n")
cat("Filtered matrix:", dim(matrix_filtered), "\n")
cat("Rows removed:", sum(!rows_to_keep), "\n")
cat("Cols removed:", sum(!cols_to_keep), "\n")

# Plot the filtered matrix
options(repr.plot.width = 60, repr.plot.height = 60, repr.plot.res = 100)
pheatmap(matrix_filtered,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(c("white", "red"))(100),
         border_color = NA,
         fontsize_row = 12,
         fontsize_col = 12,
         main = "Functional Similarity: AD_CTL vs AD_Dyslexia (Pathways with cross-similarities)",
         show_colnames = TRUE,
         show_rownames = TRUE)

# Save the filtered heatmap
png(file.path(output_dir, "fig_similarity_heatmap_filtered.png"), width = 3000, height = 3000, res = 100)
pheatmap(matrix_filtered,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(c("white", "red"))(100),
         border_color = NA,
         fontsize_row = 12,
         fontsize_col = 12,
         main = "Functional Similarity: AD_CTL vs AD_Dyslexia (Pathways with cross-similarities)",
         show_colnames = TRUE,
         show_rownames = TRUE)
dev.off()

cat("Saved: similarity_heatmap_filtered.png\n")





# What the clustering shows:

# This matrix represents functional similarity between cell type pairs across two groups (AD_CTL and AD_Dyslexia).
# Here's what it shows:

# Matrix Structure:

# Rows: Cell type pairs from AD_CTL (e.g., "5HT_Htr1a--AD_CTL", "5HT_Htr1b--AD_CTL")
# Columns: Cell type pairs from both AD_CTL and AD_Dyslexia
# Values: Binary (0 or 1) indicating whether the functional interactions are similar

# What the values mean:

# 1 = The two cell type pairs have similar functional interactions (similar ligand-receptor communication patterns)
# 0 = The cell type pairs do not have similar functional interactions

# Interpretation:

# This is a cross-group comparison matrix that shows which cell-cell communication patterns in AD_CTL are functionally similar to patterns in AD_Dyslexia. For example:

# The diagonal pattern of 1s suggests that each cell type pair is similar to itself
# Off-diagonal 1s would indicate that different cell type pairs share similar functional interaction profiles

# ?computeNetSimilarityPairwise_Neuron

# Row clustering (vertical dendrogram):

# Groups pathways from ALM that have similar functional profiles
# Pathways that cluster together have correlated communication patterns across cell types

# Column clustering (horizontal dendrogram):

# Groups pathways from both ALM and VISp based on their similarity patterns
# Reveals which pathways behave similarly within and across regions

# Biological interpretations:

# 1. Conserved signaling modules (diagonal blocks of red)

# Pathways that are functionally similar in both ALM and VISp
# Indicates core communication programs that are preserved across cortical regions
# Example: If glutamate receptors cluster together, they share similar cell-type expression patterns in both regions

# 2. Region-specific modules (off-diagonal patterns)

# ALM pathways that are more similar to certain VISp pathways than to their own ALM counterparts
# Suggests functional rewiring or specialization between regions




length(all_pathways)

embedding_data <- neuronchat_list@net_analysis$similarity$functional$dr$`1-2`
head(embedding_data)

cluster_info <- neuronchat_list@net_analysis$similarity$functional$group$`1-2`
head(cluster_info)

# Extract the embedding coordinates and cluster information
embedding_data <- neuronchat_list@net_analysis$similarity$functional$dr$`1-2`
cluster_info <- neuronchat_list@net_analysis$similarity$functional$group$`1-2`

# Get pathway names and their regions
pathway_names <- rownames(embedding_data)

# head(embedding_data)
# head(cluster_info)
# tail(cluster_info)
# pathway_names
length(pathway_names)





# Code written by Claude to find the SHIFTED PATHWAYS :

# Separate AD_CTL and AD_Dyslexia pathways
adctl_pathways  <- grep("--AD_CTL$",  pathway_names, value = TRUE)
addys_pathways <- grep("--AD_Dyslexia$", pathway_names, value = TRUE)

# Extract base pathway names (without dataset suffix)
base_names_adctl  <- gsub("--AD_CTL$",  "", adctl_pathways)
base_names_addys <- gsub("--AD_Dyslexia$", "", addys_pathways)

# Use UNION of pathways (all pathways present in either dataset)
all_pathways <- union(base_names_adctl, base_names_addys)
adctl_only        <- setdiff(base_names_adctl,  base_names_addys)
addys_only       <- setdiff(base_names_addys, base_names_adctl)

# Summary counts
cat("Total AD_CTL pathways:        ", length(base_names_adctl),  "\n")
cat("Total AD_Dyslexia pathways:       ", length(base_names_addys), "\n")
cat("Union pathways (all_pathways):           ", length(all_pathways), "\n")
cat("AD_CTL-only pathways:         ", length(adctl_only),        "\n")
cat("AD_Dyslexia-only pathways:        ", length(addys_only),       "\n")


cat("AD_CTL-only pathways:         ", adctl_only,        "\n")
cat("AD_Dyslexia-only pathways:        ", addys_only,       "\n")


# Initialize distance table
distances <- data.frame(
  pathway = character(),
  distance = numeric(),
  adctl_cluster = integer(),
  addys_cluster = integer(),
  stringsAsFactors = FALSE
)

# all_pathways
length(all_pathways)



for (pathway in all_pathways) {
  adctl_name <- paste0(pathway, "--AD_CTL")
  addys_name <- paste0(pathway, "--AD_Dyslexia")

  # Get coordinates (handle missing pathways with NA)
  adctl_coords <- if(adctl_name %in% rownames(embedding_data)) {
    embedding_data[adctl_name, ]
  } else {
    NA
  }
  addys_coords <- if(addys_name %in% rownames(embedding_data)) {
    embedding_data[addys_name, ]
  } else {
    NA
  }
  
  # Calculate Euclidean distance (only if both exist)
  if(!any(is.na(adctl_coords)) && !any(is.na(addys_coords))) {
    dist <- sqrt(sum((adctl_coords - addys_coords)^2))
  } else {
    dist <- NA  # Pathway missing in one condition
  }
  
  # Get cluster assignments
  adctl_clust <- if(adctl_name %in% names(cluster_info)) {
    cluster_info[adctl_name]
  } else {
    NA
  }
  addys_clust <- if(addys_name %in% names(cluster_info)) {
    cluster_info[addys_name]
  } else {
    NA
  }
  
  distances <- rbind(distances, data.frame(
    pathway = pathway,
    distance = dist,
    adctl_cluster = adctl_clust,
    addys_cluster = addys_clust
  ))
}




# distances

# Modified code to work with ALL pathways instead of top_n_pathways
# Replace the corresponding section in neuronchat.diff.comp.v2_5.r

# distances

# Sort by distance (largest change first)
distances <- distances[order(distances$distance, decreasing = TRUE), ]

# Display all pathways with largest location changes
print("All pathways with UMAP location changes:")
print(distances)

# Visualize all pathways
library(ggplot2)

# Use all pathways (no top_n_pathways filtering)
all_pathways <- distances$pathway

ggplot(distances, aes(x = reorder(pathway, distance), y = distance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "All Pathways with UMAP Distance Between AD_CTL and AD_Dyslexia",
       x = "Pathway",
       y = "Euclidean Distance in UMAP Space") +
  theme_minimal()

# Also identify pathways that switched clusters
cluster_switches <- distances[distances$adctl_cluster != distances$addys_cluster, ]
print("\nPathways that switched clusters between AD_CTL and AD_Dyslexia:")
print(cluster_switches[order(cluster_switches$distance, decreasing = TRUE), ])

# Create a summary table for ALL pathways
summary_table <- data.frame(
  Pathway = distances$pathway,
  UMAP_Distance = round(distances$distance, 3),
  AD_CTL_Cluster = distances$adctl_cluster,
  AD_Dyslexia_Cluster = distances$addys_cluster,
  Cluster_Change = distances$adctl_cluster != distances$addys_cluster
)

print("\nSummary of all pathways:")
print(summary_table)

p <- ggplot(distances,
            aes(x = reorder(pathway, distance), y = distance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "All Pathways with UMAP Distance Between AD_CTL and AD_Dyslexia",
       x = "Pathway",
       y = "Euclidean Distance in UMAP Space") +
  theme_minimal()

ggsave(
  filename = file.path(output_dir, "fig_all_pathways_umap_distance_AD_CTL_vs_AD_Dyslexia.png"),
  plot = p,
  width = 7,
  height = max(6, nrow(distances) * 0.1),  # Adjust height based on number of pathways
  dpi = 300, 
  bg = "white"
)

print(p)










# Long bar = Large UMAP distance = Pathway has very different functional characteristics between conditions
# Short bar = Small UMAP distance = Pathway is functionally similar/conserved between conditions



# Interpretation:

# AD_Dyslexia → AD_CTL Similarities shows that certain cell-cell communication patterns (ligand-receptor interactions) 
# that occur in AD_Dyslexia are functionally similar to communication patterns in AD_CTL.
# Specific Examples:
# Row 1: Cort_Sstr2--AD_Dyslexia is similar to Cort_Sstr2--AD_CTL

# Meaning: The communication between Cort and Sstr2 cells in AD_Dyslexia uses similar ligand-receptor pathways as 
# the communication between Cort and Sstr2 cells in AD_CTL
# Interpretation: This same cell type pair has conserved functional communication across both conditions

# Row 2: Glu_Grik1--AD_Dyslexia is similar to Glu_Grin3a--AD_CTL

# Meaning: The communication between Glu and Grik1 cells in AD_Dyslexia is functionally similar to the communication between Glu and Grin3a cells 
# in AD_CTL
# Interpretation: Even though the target cells are different (Grik1 vs Grin3a), they use similar ligand-receptor communication mechanisms

# Biological Significance:

# This tells you:

# Conserved Communication Patterns: These interactions are maintained across different conditions (AD_CTL and AD_Dyslexia), 
# suggesting they're important for brain function
# Functional Equivalence: Even when different cell types are involved (like Grik1 in AD_Dyslexia vs Grin3a in AD_CTL), 
# they can perform similar communicative roles

# Condition Specificity vs Conservation:

# Same cell types (like Cort_Sstr2 in both conditions) = conserved specific pathway
# Different cell types (like Grik1 vs Grin3a) = functionally equivalent pathways using different molecular players

# Measures how far apart the same pathway is in UMAP space between regions

# Large distance = pathway has very different functional characteristics between AD_CTL and AD_Dyslexia
# Small distance = pathway is functionally conserved

# Bar plot visualization
# Top 50 pathways with largest UMAP distances
# Longer bars = more region-specific functional profiles

# cluster_switches <- distances[distances$adctl_cluster != distances$addys_cluster, ]

# Identifies pathways that moved to a different functional module between regions
# Example: A pathway might cluster with glutamate signaling in AD_CTL but with GABA signaling in AD_Dyslexia

# Biological interpretation:

# High distance + same cluster = Quantitative change (same function, different strength)
# High distance + cluster switch = Qualitative change (different functional role in each region)
# Low distance = Conserved pathway across regions

# Next, the script calculates:

# Average distance of high-change pathways
# Median and maximum distances
# How many switched clusters
# What percentage switched clusters




# Filter pathways with distance > 1
high_distance_pathways <- distances[distances$distance > 1, ]

# Remove duplicate rows (you have duplicates in your output)
high_distance_pathways <- high_distance_pathways[!duplicated(high_distance_pathways$pathway), ]

print("Pathways with UMAP distance > 1:")
print(high_distance_pathways)

# Count them
print(paste("\nTotal number of pathways with distance > 1:", nrow(high_distance_pathways)))

# is this OK to add this code here ? 

  adctl_name <- paste0(pathway, "--AD_CTL")
  addys_name <- paste0(pathway, "--AD_Dyslexia")

  alm_name = adctl_name
  visp_name = addys_name
    
  # Get coordinates
#   alm_coords <- embedding_data[alm_name, ]
#  visp_coords <- embedding_data[visp_name, ]
  
  # Calculate Euclidean distance
#  dist <- sqrt(sum((alm_coords - visp_coords)^2))
  
  # Get cluster assignments
#  alm_clust <- cluster_info[alm_name]
#  visp_clust <- cluster_info[visp_name]






# Visualize these high-distance pathways
library(ggplot2)

ggplot(high_distance_pathways, aes(x = reorder(pathway, distance), y = distance, 
                                    fill = factor(adctl_cluster != addys_cluster))) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "coral", "FALSE" = "steelblue"),
                    labels = c("Same Cluster", "Cluster Switch"),
                    name = "Cluster Change") +
  labs(title = "Pathways with UMAP Distance > 1",
       subtitle = "Between AD_CTL and AD_Dyslexia",
       x = "Pathway",
       y = "Euclidean Distance in UMAP Space") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12))

# Create a detailed summary table
high_distance_summary <- data.frame(
  Pathway = high_distance_pathways$pathway,
  Distance = round(high_distance_pathways$distance, 3),
  AD_CTL_Cluster = high_distance_pathways$adctl_cluster,
  AD_Dyslexia_Cluster = high_distance_pathways$addys_cluster,
  Cluster_Switch = ifelse(high_distance_pathways$adctl_cluster != high_distance_pathways$addys_cluster, 
                          "Yes", "No"),
  Distance_Category = cut(high_distance_pathways$distance, 
                          breaks = c(1, 2, 3, 4, Inf),
                          labels = c("1-2", "2-3", "3-4", ">4"))
)

print("\nDetailed summary of high-distance pathways:")
print(high_distance_summary)



# Get the embedding data and extract coordinates for high-distance pathways

embedding_coords <- as.data.frame(neuronchat_list@net_analysis$similarity$functional$dr$`1-2`)
colnames(embedding_coords) <- c("UMAP1", "UMAP2")
embedding_coords$pathway_full <- rownames(embedding_coords)

# Create arrow data for high-distance pathways
arrow_data <- data.frame()

for(i in 1:nrow(high_distance_pathways)) {
  pathway <- high_distance_pathways$pathway[i]
  adctl_name <- paste0(pathway, "--AD_CTL")
  addys_name <- paste0(pathway, "--AD_Dyslexia")
  
  # Get coordinates (handle missing pathways)
  adctl_coord <- if(adctl_name %in% rownames(embedding_coords)) {
    embedding_coords[adctl_name, c("UMAP1", "UMAP2")]
  } else {
    data.frame(UMAP1 = NA, UMAP2 = NA)
  }
  addys_coord <- if(addys_name %in% rownames(embedding_coords)) {
    embedding_coords[addys_name, c("UMAP1", "UMAP2")]
  } else {
    data.frame(UMAP1 = NA, UMAP2 = NA)
  }
  
  # Create arrow data (only if both exist)
  if(!any(is.na(adctl_coord)) && !any(is.na(addys_coord))) {
    arrow_data <- rbind(arrow_data, data.frame(
      pathway = pathway,
      x_start = adctl_coord$UMAP1,
      y_start = adctl_coord$UMAP2,
      x_end = addys_coord$UMAP1,
      y_end = addys_coord$UMAP2,
      x_mid = (adctl_coord$UMAP1 + addys_coord$UMAP1) / 2,
      y_mid = (adctl_coord$UMAP2 + addys_coord$UMAP2) / 2,
      distance = high_distance_pathways$distance[i],
      adctl_cluster = high_distance_pathways$adctl_cluster[i],
      addys_cluster = high_distance_pathways$addys_cluster[i],
      cluster_switch = high_distance_pathways$adctl_cluster[i] != high_distance_pathways$addys_cluster[i]
    ))
  }
}

library("ggrepel")

# Add region and cluster info to embedding coords
embedding_coords$region <- ifelse(grepl("--AD_CTL$", embedding_coords$pathway_full), "AD_CTL", "AD_Dyslexia")
embedding_coords$pathway <- gsub("--AD_CTL$|--AD_Dyslexia$", "", embedding_coords$pathway_full)
embedding_coords$cluster <- sapply(embedding_coords$pathway_full, function(x) cluster_info[x])

# Filter to only show high-distance pathways in the base plot
high_dist_pathway_names <- high_distance_pathways$pathway
embedding_coords_filtered <- embedding_coords[embedding_coords$pathway %in% high_dist_pathway_names, ]

# Create the UMAP plot with arrows and labels
p <- ggplot() +
  # Plot all points (faded in background)
  geom_point(data = embedding_coords, 
             aes(x = UMAP1, y = UMAP2, color = factor(cluster), shape = region),
             alpha = 0.2, size = 4) +
  # Highlight high-distance pathways
  geom_point(data = embedding_coords_filtered, 
             aes(x = UMAP1, y = UMAP2, color = factor(cluster), shape = region),
             size = 4, alpha = 0.8) +
  # Add arrows
  geom_segment(data = arrow_data,
               aes(x = x_start, y = y_start, xend = x_end, yend = y_end, 
                   linewidth = distance, alpha = distance),
               arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
               color = "black") +
  # Add pathway labels at arrow midpoints
  geom_text_repel(data = arrow_data,
                  aes(x = x_mid, y = y_mid, label = pathway),
                  size = 3, fontface = "bold",
                  box.padding = 0.5,
                  point.padding = 0.5,
                  segment.color = "grey50",
                  segment.size = 0.2,
                  max.overlaps = 20) +
  scale_linewidth_continuous(range = c(0.5, 2), name = "Distance") +
  scale_alpha_continuous(range = c(0.3, 0.8), guide = "none") +
  scale_shape_manual(values = c("AD_CTL" = 16, "AD_Dyslexia" = 17), name = "Dataset") +
  scale_color_brewer(palette = "Set1", name = "Cluster") +
  labs(title = "UMAP with Arrows Showing High-Distance Pathways (Distance > 1)",
       subtitle = "Arrows point from AD_CTL to AD_Dyslexia with pathway labels",
       x = "UMAP1", y = "UMAP2") +
  theme_minimal() +
  theme(legend.position = "right")

print(p)

# Alternative: Show only pathways that switch clusters with labels
arrow_data_switch <- arrow_data[arrow_data$cluster_switch == TRUE, ]

p2 <- ggplot() +
  geom_point(data = embedding_coords, 
             aes(x = UMAP1, y = UMAP2, color = factor(cluster), shape = region),
             alpha = 0.15, size = 4) +
  geom_point(data = embedding_coords[embedding_coords$pathway %in% arrow_data_switch$pathway, ], 
             aes(x = UMAP1, y = UMAP2, color = factor(cluster), shape = region),
             size = 4, alpha = 0.8) +
  geom_segment(data = arrow_data_switch,
               aes(x = x_start, y = y_start, xend = x_end, yend = y_end, 
                   linewidth = distance),
               arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
               color = "red", alpha = 0.6) +
  geom_text_repel(data = arrow_data_switch,
                  aes(x = x_mid, y = y_mid, label = pathway),
                  size = 4, fontface = "bold", color = "darkred",
                  box.padding = 0.5,
                  point.padding = 0.5,
                  segment.color = "grey50",
                  max.overlaps = 25) +
  scale_linewidth_continuous(range = c(0.5, 2), name = "Distance") +
  scale_shape_manual(values = c("AD_CTL" = 16, "AD_Dyslexia" = 17), name = "Dataset") +
  scale_color_brewer(palette = "Set1", name = "Cluster") +
  labs(title = "UMAP: Pathways that Switch Clusters (Distance > 1)",
       subtitle = "Red arrows show AD_CTL → AD_Dyslexia transitions with pathway labels",
       x = "UMAP1", y = "UMAP2") +
  theme_minimal() +
  theme(legend.position = "right")

print(p2)




# Option 3: Top 10 most divergent with larger labels
top_divergent <- arrow_data[order(arrow_data$distance, decreasing = TRUE), ][1:20, ]

p3 <- ggplot() +
  geom_point(data = embedding_coords, 
             aes(x = UMAP1, y = UMAP2, color = factor(cluster), shape = region),
             alpha = 0.2, size = 4) +
  geom_segment(data = top_divergent,
               aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
               arrow = arrow(length = unit(0.4, "cm"), type = "closed"),
               color = "darkred", linewidth = 1.2, alpha = 0.7) +
  geom_text_repel(data = top_divergent,
                  aes(x = x_mid, y = y_mid, label = pathway),
                  size = 4, fontface = "bold", color = "darkred",
                  box.padding = 1,
                  point.padding = 0.5,
                  segment.color = "grey50",
                  max.overlaps = 15) +
  scale_shape_manual(values = c("AD_CTL" = 16, "AD_Dyslexia" = 17), name = "Dataset") +
  scale_color_brewer(palette = "Set1", name = "Cluster") +
  labs(title = "Top 10 Most Divergent Pathways Between AD_CTL and AD_Dyslexia",
       subtitle = "Arrows show direction from AD_CTL to AD_Dyslexia with pathway labels",
       x = "UMAP1", y = "UMAP2") +
  theme_minimal() +
  theme(legend.position = "right")

print(p3)








# Option 4: Simpler version without ggrepel (direct labels)
p4 <- ggplot() +
  geom_point(data = embedding_coords, 
             aes(x = UMAP1, y = UMAP2, color = factor(cluster), shape = region),
             alpha = 0.2, size = 4) +
  geom_segment(data = arrow_data,
               aes(x = x_start, y = y_start, xend = x_end, yend = y_end, 
                   linewidth = distance),
               arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
               color = "black", alpha = 0.5) +
  # Direct labels (may overlap)
  geom_text(data = arrow_data,
            aes(x = x_mid, y = y_mid, label = pathway),
            size = 3, fontface = "bold", hjust = 0.5, vjust = -0.5) +
  scale_linewidth_continuous(range = c(0.5, 2), name = "Distance") +
  scale_shape_manual(values = c("AD_CTL" = 16, "AD_Dyslexia" = 17), name = "Dataset") +
  scale_color_brewer(palette = "Set1", name = "Cluster") +
  labs(title = "UMAP with All High-Distance Pathways Labeled",
       subtitle = "Direct labels at arrow midpoints",
       x = "UMAP1", y = "UMAP2") +
  theme_minimal() +
  theme(legend.position = "right")

print(p4)

# Save plots with larger dimensions to accommodate labels
ggsave(file.path(output_dir, "fig_umap_with_arrows_labeled_all.png"), plot = p, width = 14, height = 10, dpi = 300,  bg = "white")
ggsave(file.path(output_dir, "fig_umap_with_arrows_labeled_cluster_switch.png"), plot = p2, width = 14, height = 10, dpi = 300,  bg = "white")
ggsave(file.path(output_dir, "fig_umap_with_arrows_labeled_simple.png"), plot = p4, width = 14, height = 10, dpi = 300,  bg = "white")

print("Plots with labels saved!")


# UMAP_all_divergent_pathways.csv - Complete list sorted by distance
# UMAP_top10_most_divergent_pathways.csv - Quick reference for top 10
# UMAP_analysis_summary.txt - Text file with BOTH top 10 section AND complete list
# UMAP_complete_analysis.xlsx - Excel file with 6 sheets including both

# 1. Save arrow data (pathways with distance > 1)
arrow_data_export <- arrow_data
arrow_data_export <- arrow_data_export[order(arrow_data_export$distance, decreasing = TRUE), ]

write.csv(arrow_data_export, 
          file.path(output_dir, "res_UMAP_high_distance_pathways_arrow_data.csv"), 
          row.names = FALSE)

# 2. Save all UMAP coordinates with metadata
embedding_export <- data.frame(
  Pathway_Full = embedding_coords$pathway_full,
  Pathway = embedding_coords$pathway,
  Region = embedding_coords$region,
  Cluster = embedding_coords$cluster,
  UMAP1 = embedding_coords$UMAP1,
  UMAP2 = embedding_coords$UMAP2
)

write.csv(embedding_export, 
          file.path(output_dir, "res_UMAP_all_pathways_coordinates.csv"), 
          row.names = FALSE)

# 3. Save high-distance pathways with their coordinates
high_distance_coords <- data.frame()

for(i in 1:nrow(high_distance_pathways)) {
  pathway <- high_distance_pathways$pathway[i]
  adctl_name <- paste0(pathway, "--AD_CTL")
  addys_name <- paste0(pathway, "--AD_Dyslexia")
  
  adctl_data <- embedding_coords[embedding_coords$pathway_full == adctl_name, ]
  addys_data <- embedding_coords[embedding_coords$pathway_full == addys_name, ]
  
  # Only add if AD_CTL data exists
  if(nrow(adctl_data) > 0) {
    high_distance_coords <- rbind(high_distance_coords, data.frame(
      Pathway = pathway,
      Region = "AD_CTL",
      UMAP1 = adctl_data$UMAP1,
      UMAP2 = adctl_data$UMAP2,
      Cluster = adctl_data$cluster,
      Distance = high_distance_pathways$distance[i],
      AD_CTL_Cluster = high_distance_pathways$adctl_cluster[i],
      AD_Dyslexia_Cluster = high_distance_pathways$addys_cluster[i],
      Cluster_Switch = high_distance_pathways$adctl_cluster[i] != high_distance_pathways$addys_cluster[i]
    ))
  }
  
  # Only add if AD_Dyslexia data exists
  if(nrow(addys_data) > 0) {
    high_distance_coords <- rbind(high_distance_coords, data.frame(
      Pathway = pathway,
      Region = "AD_Dyslexia",
      UMAP1 = addys_data$UMAP1,
      UMAP2 = addys_data$UMAP2,
      Cluster = addys_data$cluster,
      Distance = high_distance_pathways$distance[i],
      AD_CTL_Cluster = high_distance_pathways$adctl_cluster[i],
      AD_Dyslexia_Cluster = high_distance_pathways$addys_cluster[i],
      Cluster_Switch = high_distance_pathways$adctl_cluster[i] != high_distance_pathways$addys_cluster[i]
    ))
  }
}

write.csv(high_distance_coords, 
          file.path(output_dir, "res_UMAP_high_distance_pathways_coordinates.csv"), 
          row.names = FALSE)

# Save all divergent pathways (sorted by distance)
all_divergent <- arrow_data[order(arrow_data$distance, decreasing = TRUE), ]

# Save top 10 divergent pathways
top_10_divergent <- arrow_data[order(arrow_data$distance, decreasing = TRUE), ][1:10, ]

# Optionally: Create summary report and open sink
sink(file.path(output_dir, "UMAP_analysis_summary.txt"))

cat("Cluster Distribution:\n")
cluster_dist <- table(embedding_coords$cluster)
for(i in 1:length(cluster_dist)) {
  cat("  Cluster", names(cluster_dist)[i], ":", cluster_dist[i], "pathways\n")
}

cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("TOP 10 MOST DIVERGENT PATHWAYS\n")
cat(paste(rep("=", 60), collapse=""), "\n\n")
for(i in 1:10) {
  cat(i, ". ", top_10_divergent$pathway[i], "\n", sep="")
  cat("   Distance: ", round(top_10_divergent$distance[i], 3), "\n")
  cat("   AD_CTL Cluster: ", top_10_divergent$adctl_cluster[i], 
      " -> AD_Dyslexia Cluster: ", top_10_divergent$addys_cluster[i], "\n")
  cat("   Cluster Switch: ", ifelse(top_10_divergent$cluster_switch[i], "YES", "NO"), "\n\n")
}

cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("ALL DIVERGENT PATHWAYS (sorted by distance)\n")
cat(paste(rep("=", 60), collapse=""), "\n\n")
for(i in 1:nrow(all_divergent)) {
  cat(i, ". ", all_divergent$pathway[i], "\n", sep="")
  cat("   Distance: ", round(all_divergent$distance[i], 3), "\n")
  cat("   AD_CTL Cluster: ", all_divergent$adctl_cluster[i], 
      " -> AD_Dyslexia Cluster: ", all_divergent$addys_cluster[i], "\n")
  cat("   Cluster Switch: ", ifelse(all_divergent$cluster_switch[i], "YES", "NO"), "\n\n")
}

sink()




# 8. Create an Excel workbook with multiple sheets (if openxlsx is available)
if(require(openxlsx)) {
  wb <- createWorkbook()
  
  # Sheet 1: Arrow data
  addWorksheet(wb, "High_Distance_Arrows")
  writeData(wb, "High_Distance_Arrows", arrow_data_export)
  
  # Sheet 2: All coordinates
  addWorksheet(wb, "All_UMAP_Coordinates")
  writeData(wb, "All_UMAP_Coordinates", embedding_export)
  
  # Sheet 3: High distance coordinates
  addWorksheet(wb, "High_Distance_Coords")
  writeData(wb, "High_Distance_Coords", high_distance_coords)
  
  # Sheet 4: Cluster switches
  addWorksheet(wb, "Cluster_Switches")
  writeData(wb, "Cluster_Switches", cluster_switch_data)
  
  # Sheet 5: All divergent pathways
  addWorksheet(wb, "All_Divergent")
  writeData(wb, "All_Divergent", all_divergent)
  
  # Sheet 6: Top 10 divergent pathways
  addWorksheet(wb, "Top_10_Divergent")
  writeData(wb, "Top_10_Divergent", top_10_divergent)
  
  # Save workbook
  saveWorkbook(wb, "UMAP_complete_analysis.xlsx", overwrite = TRUE)
  
  print("Excel workbook created: UMAP_complete_analysis.xlsx")
} else {
  print("Package 'openxlsx' not installed. Individual CSV files created instead.")
}

# Print summary of saved files
cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("FILES SAVED:\n")
cat(paste(rep("=", 60), collapse=""), "\n")
cat("1. UMAP_high_distance_pathways_arrow_data.csv - Arrow coordinates and distances\n")
cat("2. UMAP_all_pathways_coordinates.csv - All UMAP coordinates\n")
cat("3. UMAP_high_distance_pathways_coordinates.csv - High-distance pathway details\n")
cat("4. UMAP_cluster_switching_pathways.csv - Pathways that switch clusters\n")
cat("5. UMAP_all_divergent_pathways.csv - ALL divergent pathways (sorted by distance)\n")
cat("6. UMAP_top10_most_divergent_pathways.csv - Top 10 divergent pathways\n")
cat("7. UMAP_analysis_summary.txt - Text summary report with TOP 10 and ALL pathways\n")
if(require(openxlsx)) {
  cat("8. UMAP_complete_analysis.xlsx - All data in one Excel file (6 sheets)\n")
}
cat(paste(rep("=", 60), collapse=""), "\n\n")





# Heatmap for each interaction pattern

net12 <- neuronchat_list@net[c(1,2)]

net1 <- net12[[1]]
names(net1) <- paste(names(net1),'--AD_CTL',sep='')

net2 <- net12[[2]]
names(net2) <- paste(names(net2),'--AD_Dyslexia',sep='')

net12_list <- append(net1,net2)

# net1
# net2

# str(net1)
# str(net2)

length(net1)
length(net2)

interaction_group <- neuronchat_list@net_analysis$similarity$functional
head(interaction_group, 2)

length(interaction_group)

# First, check the structure of interaction_group
str(interaction_group)

# Check what names/elements are available
names(interaction_group)

# Or if it's a list without names:
length(interaction_group)

# Look at the first few elements
head(interaction_group)

M <- interaction_group$matrix[["1-2"]]
head(M, 2)
dim(M)

DR <- interaction_group$dr[["1-2"]]
head(DR, 2)
dim(DR)

interaction_group <- neuronchat_list@net_analysis$similarity$functional$group$`1-2`

# interaction_group 
table(interaction_group)








cat("

What a “pattern” is in NeuronChat

NeuronChat does this in three steps:

Each interaction program (e.g. Sst_Sstr2, Glu_Grm2, Vip_Vipr1) produces a cell-type × cell-type communication matrix.

NeuronChat computes functional similarity between all these matrices (based on sender→receiver structure).

It clusters interactions with similar matrices into pattern clusters.

Each cluster is labeled:

pattern cluster 1

pattern cluster 2

pattern cluster 3

")

cat("

What “pattern cluster 1” specifically represents

Pattern cluster 1 = one canonical communication motif, characterized by:

A specific subset of sender cell types

A specific subset of receiver cell types

A characteristic directional flow (who talks to whom)

In your heatmap, pattern 1 is showing:

Rows = sender cell types

Columns = receiver cell types

Color = aggregated communication weight

Aggregation = sum of all interaction programs assigned to pattern 1

")


pattern1_pathways <- names(interaction_group)[interaction_group == 1]
# pattern1_pathways

pattern2_pathways <- names(interaction_group)[interaction_group == 2]
# pattern2_pathways

pattern3_pathways <- names(interaction_group)[interaction_group == 3]
# pattern3_pathways

pattern4_pathways <- names(interaction_group)[interaction_group == 4]
# pattern4_pathways

pattern5_pathways <- names(interaction_group)[interaction_group == 5]
# pattern5_pathways

# net1

# net2

# Save all pathway matrices for AD_CTL (net1) and AD_Dyslexia (net2)

# Save all pathway matrices for AD_CTL (net1) and AD_Dyslexia (net2)


# Get pathway names
pathways_adctl <- names(net1)
pathways_addys <- names(net2)

# Get union of all pathways
all_pathways <- union(pathways_adctl, pathways_addys)

cat("Total pathways in AD_CTL:", length(pathways_adctl), "\n")
cat("Total pathways in AD_Dyslexia:", length(pathways_addys), "\n")
cat("Total unique pathways:", length(all_pathways), "\n\n")

pattern1_pathways 

pattern2_pathways 

pattern3_pathways 

pattern4_pathways 

pattern5_pathways 






# ============================================================================
# Save each pathway matrix as a separate CSV file in current directory
# Filename format: PathwayName_DatasetName.csv
# ============================================================================

# Save AD_CTL matrices
cat("Saving AD_CTL pathway matrices...\n")
for(pathway in pathways_adctl) {
  filename <- file.path(output_dir, paste0(pathway, "_AD_CTL.csv"))
  write.csv(net1[[pathway]], filename, row.names = TRUE)
}

# Save AD_Dyslexia matrices
cat("Saving AD_Dyslexia pathway matrices...\n")
for(pathway in pathways_addys) {
  filename <- file.path(output_dir, paste0(pathway, "_AD_Dyslexia.csv"))
  write.csv(net2[[pathway]], filename, row.names = TRUE)
}

cat("\n=============================================================\n")
cat("SUMMARY OF SAVED FILES\n")
cat("=============================================================\n")
cat("Total files saved in current directory: ", length(pathways_adctl) + length(pathways_addys), "\n", sep="")
cat("  - AD_CTL matrices: ", length(pathways_adctl), " files\n", sep="")
cat("  - AD_Dyslexia matrices: ", length(pathways_addys), " files\n", sep="")
cat("\nFilename format: PathwayName_DatasetName.csv\n")
cat("Examples:\n")
cat("  - 5HT_Htr1a--AD_CTL_AD_CTL.csv\n")
cat("  - 5HT_Htr1a--AD_Dyslexia_AD_Dyslexia.csv\n")
cat("  - Glu_Grin3a--AD_CTL_AD_CTL.csv\n")
cat("  - Glu_Grin3a--AD_Dyslexia_AD_Dyslexia.csv\n")
cat("=============================================================\n")





# Part IV: Shared and specific interaction patterns across AD_CTL and AD_Dyslexia

# Compute structural similarity 
# neuronchat_list <- computeNetSimilarityPairwise_Neuron(neuronchat_list, 
#                                                       slot.name = "net", 
#                                                       type = "structural", 
#                                                       comparison = c(1,2))

# Why Structural Similarity Isn't Available in NeuronChat
# This is likely because:

# NeuronChat is adapted from CellChat but doesn't implement all features
# Neuron-specific focus: The package may prioritize comparing signaling mechanisms (functional) 
# over network topology (structural) for neuron communication
# Simplified implementation: The developers may have chosen to focus only on functional similarity 
# for neuron-neuron communication analysis

# packageVersion("NeuronChat")
# args(computeNetSimilarityPairwise_Neuron)
# or
# ?computeNetSimilarityPairwise_Neuron


# Why are most off-diagonal values zero?
# Because NeuronChat is very conservative.
# Zero means:
# These two interaction programs do not produce similar communication patterns across cell types.
# Example:
# 5HT_Htr1a--AD_CTL vs Ach_Chrm1--AD_CTL
# → completely different signaling modes → similarity = 0








# Measures how far apart the same pathway is in UMAP space between regions

# Large distance = pathway has very different functional characteristics between AD_CTL and AD_Dyslexia
# Small distance = pathway is functionally conserved

# Bar plot visualization
# Top 50 pathways with largest UMAP distances
# Longer bars = more region-specific functional profiles

# cluster_switches <- distances[distances$adctl_cluster != distances$addys_cluster, ]

# Identifies pathways that moved to a different functional module between regions
# Example: A pathway might cluster with glutamate signaling in AD_CTL but with GABA signaling in AD_Dyslexia

# Biological interpretation:

# High distance + same cluster = Quantitative change (same function, different strength)
# High distance + cluster switch = Qualitative change (different functional role in each region)
# Low distance = Conserved pathway across regions

# Next, the script calculates:

# Average distance of high-change pathways
# Median and maximum distances
# How many switched clusters
# What percentage switched clusters




# Part IV: Shared and specific interaction patterns across ALM and VISp

# Compute structural similarity 
# neuronchat_list <- computeNetSimilarityPairwise_Neuron(neuronchat_list, 
#                                                       slot.name = "net", 
#                                                       type = "structural", 
#                                                       comparison = c(1,2))

# Why Structural Similarity Isn't Available in NeuronChat
# This is likely because:

# NeuronChat is adapted from CellChat but doesn't implement all features
# Neuron-specific focus: The package may prioritize comparing signaling mechanisms (functional) 
# over network topology (structural) for neuron communication
# Simplified implementation: The developers may have chosen to focus only on functional similarity 
# for neuron-neuron communication analysis

# packageVersion("NeuronChat")
# args(computeNetSimilarityPairwise_Neuron)
# or
# ?computeNetSimilarityPairwise_Neuron










# Heatmap for each interaction pattern

net12 <- neuronchat_list@net[c(1,2)]

net1 <- net12[[1]]
names(net1) <- paste(names(net1),'--AD_CTL',sep='')

net2 <- net12[[2]]
names(net2) <- paste(names(net2),'--AD_Dyslexia',sep='')

net12_list <- append(net1,net2)

# net1
# net2

# str(net1)
# str(net2)

length(net1)
length(net2)

interaction_group <- neuronchat_list@net_analysis$similarity$functional
interaction_group

length(interaction_group)

# First, check the structure of interaction_group
str(interaction_group)

# Check what names/elements are available
names(interaction_group)

# Or if it's a list without names:
length(interaction_group)

# Look at the first few elements
head(interaction_group)

M <- interaction_group$matrix[["1-2"]]
# head(M)
dim(M)

DR <- interaction_group$dr[["1-2"]]
# head(DR)
dim(DR)




interaction_group <- neuronchat_list@net_analysis$similarity$functional$group$`1-2`

interaction_group 
table(interaction_group)

cat("

What a “pattern” is in NeuronChat

NeuronChat does this in three steps:

Each interaction program (e.g. Sst_Sstr2, Glu_Grm2, Vip_Vipr1) produces a cell-type × cell-type communication matrix.

NeuronChat computes functional similarity between all these matrices (based on sender→receiver structure).

It clusters interactions with similar matrices into pattern clusters.

Each cluster is labeled:

pattern cluster 1

pattern cluster 2

pattern cluster 3

")

cat("

What “pattern cluster 1” specifically represents

Pattern cluster 1 = one canonical communication motif, characterized by:

A specific subset of sender cell types

A specific subset of receiver cell types

A characteristic directional flow (who talks to whom)

In your heatmap, pattern 1 is showing:

Rows = sender cell types

Columns = receiver cell types

Color = aggregated communication weight

Aggregation = sum of all interaction programs assigned to pattern 1

")




pattern1_pathways <- names(interaction_group)[interaction_group == 1]
# pattern1_pathways

pattern2_pathways <- names(interaction_group)[interaction_group == 2]
# pattern2_pathways

pattern3_pathways <- names(interaction_group)[interaction_group == 3]
# pattern3_pathways

pattern4_pathways <- names(interaction_group)[interaction_group == 4]
# pattern4_pathways

pattern5_pathways <- names(interaction_group)[interaction_group == 5]
# pattern5_pathways

# net1

# net2

# Save all pathway matrices for AD_CTL (net1) and AD_Dyslexia (net2)

# Save all pathway matrices for AD_CTL (net1) and AD_Dyslexia (net2)

# Get pathway names
pathways_adctl <- names(net1)
pathways_addys <- names(net2)

# Get union of all pathways
all_pathways <- union(pathways_adctl, pathways_addys)

cat("Total pathways in AD_CTL:", length(pathways_adctl), "\n")
cat("Total pathways in AD_Dyslexia:", length(pathways_addys), "\n")
cat("Total unique pathways:", length(all_pathways), "\n\n")




# ============================================================================
# Save each pathway matrix as a separate CSV file in current directory
# Filename format: PathwayName_DatasetName.csv
# ============================================================================

# Save AD_CTL matrices
cat("Saving AD_CTL pathway matrices...\n")
for(pathway in pathways_adctl) {
  filename <- file.path(output_dir, paste0(pathway, "_AD_CTL.csv"))
  write.csv(net1[[pathway]], filename, row.names = TRUE)
}

# Save AD_Dyslexia matrices
cat("Saving AD_Dyslexia pathway matrices...\n")
for(pathway in pathways_addys) {
  filename <- file.path(output_dir, paste0(pathway, "_AD_Dyslexia.csv"))
  write.csv(net2[[pathway]], filename, row.names = TRUE)
}

cat("\n=============================================================\n")
cat("SUMMARY OF SAVED FILES\n")
cat("=============================================================\n")
cat("Total files saved in current directory: ", length(pathways_adctl) + length(pathways_addys), "\n", sep="")
cat("  - AD_CTL matrices: ", length(pathways_adctl), " files\n", sep="")
cat("  - AD_Dyslexia matrices: ", length(pathways_addys), " files\n", sep="")
cat("\nFilename format: PathwayName_DatasetName.csv\n")
cat("Examples:\n")
cat("  - 5HT_Htr1a--AD_CTL_AD_CTL.csv\n")
cat("  - 5HT_Htr1a--AD_Dyslexia_AD_Dyslexia.csv\n")
cat("  - Glu_Grin3a--AD_CTL_AD_CTL.csv\n")
cat("  - Glu_Grin3a--AD_Dyslexia_AD_Dyslexia.csv\n")
cat("=============================================================\n")

# interaction_group



# heatmap for each interaction pattern : this is the code from the tutorial

net12 <- neuronchat_list@net[c(1,2)]

net1 <- net12[[1]];names(net1) <- paste(names(net1),'--AD_CTL',sep='')
net2 <- net12[[2]];names(net2) <- paste(names(net2),'--AD_Dyslexia',sep='')

net12_list <- append(net1,net2)

interaction_group <- neuronchat_list@net_analysis$similarity$functional$group$`1-2`

hlist <- list();gb_heatmap <- list()
library(grid); library(ComplexHeatmap);grid.newpage(); x_seq <- c(0,0.2,0.4,0.6,0.8)

#> ========================================
#> ComplexHeatmap version 2.8.0
#> Bioconductor page: http://bioconductor.org/packages/ComplexHeatmap/
#> Github page: https://github.com/jokergoo/ComplexHeatmap
#> Documentation: http://jokergoo.github.io/ComplexHeatmap-reference
#> 
#> If you use it in published research, please cite:
#> Gu, Z. Complex heatmaps reveal patterns and correlations in multidimensional 
#>   genomic data. Bioinformatics 2016.
#> 
#> The new InteractiveComplexHeatmap package can directly export static 
#> complex heatmaps into an interactive Shiny app with zero effort. Have a try!
#> 
#> This message can be suppressed by:
#>   suppressPackageStartupMessages(library(ComplexHeatmap))
#> ========================================

for(j in 1:length(sort(unique(interaction_group),decreasing = F))){
  net_aggregated_group2 <- net_aggregation(net12_list[names(interaction_group[interaction_group==j])],
                                           method = 'weight')
  library(RColorBrewer);
  col_map = brewer.pal(8,"YlOrBr");
  h <- Heatmap(net_aggregated_group2, name = "Weight",
                        col = col_map,
                        cluster_rows = FALSE,
                        cluster_columns=FALSE,
                        row_names_side='left',
                        column_names_side='bottom',
                        row_title='Sender',
                        row_title_side='left',
                        row_title_gp = gpar(fontsize = 16),
                        column_title='Receiver',
               column_title_side = "bottom",
               column_title_gp = gpar(fontsize = 16),
               column_names_rot = 60)
  
  gb_heatmap[[j]] = grid.grabExpr(draw(h,column_title=paste('pattern cluster',j), padding = unit(c(2, 2, 2, 2), "mm")) )
  pushViewport(viewport(x = x_seq[j], 
                        y = 0.30, 
                        width = 0.30, 
                        height = 0.5, 
                        just = c("left", "top"), 
                        xscale = c(0, 1), 
                        yscale = c(0, 1)));
    grid.draw(gb_heatmap[[j]]);
    popViewport()
}





# saving the figures based on the code provided in the tutorial 

library(ComplexHeatmap)
library(RColorBrewer)
library(grid)

# Initialize list to store heatmaps
gb_heatmap <- list()

# Calculate number of clusters
n_clusters <- length(sort(unique(interaction_group), decreasing = FALSE))

# Calculate x positions for horizontal layout
x_seq <- seq(0.05, 0.95 - 0.30, length.out = n_clusters)

# Save as PNG
png(file.path(output_dir, "fig_pattern_clusters_horizontal_layout.png"), 
    width = 1000 * n_clusters, 
    height = 1200, 
    res = 150, 
    bg = "white")

grid.newpage()

for(j in 1:length(sort(unique(interaction_group), decreasing = FALSE))) {
  net_aggregated_group2 <- net_aggregation(
    net12_list[names(interaction_group[interaction_group == j])],
    method = 'weight'
  )
  
  library(RColorBrewer)
  col_map <- brewer.pal(8, "YlOrBr")
  
  h <- Heatmap(
    net_aggregated_group2, 
    name = "Weight",
    col = col_map,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_names_side = 'left',
    column_names_side = 'bottom',
    row_title = 'Sender',
    row_title_side = 'left',
    row_title_gp = gpar(fontsize = 16),
    column_title = 'Receiver',
    column_title_side = "bottom",
    column_title_gp = gpar(fontsize = 16),
    column_names_rot = 60
  )
  
  gb_heatmap[[j]] <- grid.grabExpr(
    draw(h, 
         column_title = paste('pattern cluster', j), 
         padding = unit(c(2, 2, 2, 2), "mm"))
  )
  
  pushViewport(viewport(
    x = x_seq[j], 
    y = 0.30, 
    width = 0.30, 
    height = 0.5, 
    just = c("left", "top"), 
    xscale = c(0, 1), 
    yscale = c(0, 1)
  ))
  
  grid.draw(gb_heatmap[[j]])
  popViewport()
}

dev.off()

cat("Saved: pattern_clusters_horizontal_layout_format_tutorial.png\n")

# Save as PDF (publication quality)
pdf(file.path(output_dir, "fig_pattern_clusters_horizontal_layout.pdf"), 
    width = 5 * n_clusters, 
    height = 6, 
    bg = "white")

grid.newpage()


for(j in 1:length(sort(unique(interaction_group), decreasing = FALSE))) {
  net_aggregated_group2 <- net_aggregation(
    net12_list[names(interaction_group[interaction_group == j])],
    method = 'weight'
  )
  
  col_map <- brewer.pal(8, "YlOrBr")
  
  h <- Heatmap(
    net_aggregated_group2, 
    name = "Weight",
    col = col_map,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_names_side = 'left',
    column_names_side = 'bottom',
    row_title = 'Sender',
    row_title_side = 'left',
    row_title_gp = gpar(fontsize = 16),
    column_title = 'Receiver',
    column_title_side = "bottom",
    column_title_gp = gpar(fontsize = 16),
    column_names_rot = 60
  )
  
  gb_heatmap[[j]] <- grid.grabExpr(
    draw(h, 
         column_title = paste('pattern cluster', j), 
         padding = unit(c(2, 2, 2, 2), "mm"))
  )
  
  pushViewport(viewport(
    x = x_seq[j], 
    y = 0.30, 
    width = 0.30, 
    height = 0.5, 
    just = c("left", "top"), 
    xscale = c(0, 1), 
    yscale = c(0, 1)
  ))
  
  grid.draw(gb_heatmap[[j]])
  popViewport()
}

dev.off()

cat("Saved: pattern_clusters_horizontal_layout_format_tutorial.pdf\n")

# str(net12)
glimpse(net12, max.level = 1)





# This code is creating and displaying heatmaps for different pattern clusters of signaling pathways. Let me break it down:

# Step 1: Loop through each pattern cluster
# for(j in 1:length(sort(unique(interaction_group), decreasing = F)))

# Loops through each unique cluster ID in interaction_group
# interaction_group contains cluster assignments for pathways (e.g., pathway X is in cluster 1, pathway Y is in cluster 2)

# Step 2: Aggregate networks within each cluster
# net_aggregated_group2 <- net_aggregation(net12_list[names(interaction_group[interaction_group==j])],
#                                         method = 'weight')

# Selects all pathways that belong to cluster j
# net12_list contains all pathway communication matrices
# net_aggregation() combines them by summing their weights
# Creates a single 8×8 aggregated matrix representing the total communication for all pathways in this cluster

# Heatmap settings:

# No clustering (rows/columns stay in original order)
# Rows = Sender cell types (labeled on left)
# Columns = Receiver cell types (labeled on bottom, rotated 60°)
# Colors represent communication weight strength

# Example interpretation:
# If you have 3 pattern clusters:

# Pattern cluster 1: Pathways with similar neurotransmitter signaling profiles
# Pattern cluster 2: Pathways with similar neuropeptide signaling profiles
# Pattern cluster 3: Pathways with similar gap junction profiles

# net1

# net2





# Set plot size
options(repr.plot.width = 12, repr.plot.height = 16, repr.plot.res = 120)

library(ComplexHeatmap)
library(RColorBrewer)

# Set plot size
options(repr.plot.width = 12, repr.plot.height = 16)

n_clusters <- length(sort(unique(interaction_group), decreasing = FALSE))

# Create a list to store heatmaps
heatmap_list <- list()

for(j in 1:n_clusters) {
    
  net_aggregated_group2 <- net_aggregation(
    net12_list[names(interaction_group[interaction_group == j])], 
    method = 'weight'
  )
  
  col_map <- brewer.pal(8, "YlOrBr")
  
  h <- Heatmap(
    net_aggregated_group2, 
    name = paste("Weight", j),
    col = col_map,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    column_title = paste('Pattern Cluster', j),
    column_title_gp = gpar(fontsize = 16, fontface = "bold"),
    row_title = 'Sender',
    column_names_rot = 45,
    border = TRUE
  )
  
  heatmap_list[[j]] <- h
}

# Draw all heatmaps vertically stacked
ht_list <- NULL
for(i in 1:length(heatmap_list)) {
  if(is.null(ht_list)) {
    ht_list <- heatmap_list[[i]]
  } else {
    ht_list <- ht_list %v% heatmap_list[[i]]  # Stack vertically
  }
}

# Save as PNG
png(file.path(output_dir, "fig_pattern_clusters_heatmaps.png"), width = 2400, height = 3200, res = 150)
draw(ht_list)
dev.off()

cat("Saved: pattern_clusters_heatmaps.png\n")


# Optional: Save as PDF (publication quality)
pdf(file.path(output_dir, "fig_pattern_clusters_heatmaps.pdf"), width = 12, height = 16)
draw(ht_list)
dev.off()

cat("Saved: pattern_clusters_heatmaps.pdf\n")

# another visualization

library(ComplexHeatmap)
library(RColorBrewer)
library(grid)



# Save as PNG
png(file.path(output_dir, "fig_pattern_clusters_grid_layout.png"), width = 3000, height = 4000, res = 200, bg = "white")

grid.newpage()

n_clusters <- length(sort(unique(interaction_group), decreasing = FALSE))

# More spacing between heatmaps
spacing <- 0.06
height_available <- 0.92  # 92% of page height
height_each <- (height_available - spacing * (n_clusters - 1)) / n_clusters

for(j in 1:n_clusters) {
  net_aggregated_group2 <- net_aggregation(
    net12_list[names(interaction_group[interaction_group == j])], 
    method = 'weight'
  )
  
  col_map <- brewer.pal(8, "YlOrBr")
  
  h <- Heatmap(
    net_aggregated_group2, 
    name = "Weight",
    col = col_map,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    
    # Good dimensions for PDF
    width = unit(14, "cm"),
    height = unit(7, "cm"),
    
    row_names_side = 'left',
    row_names_gp = gpar(fontsize = 13),
    column_names_side = 'bottom',
    column_names_gp = gpar(fontsize = 13),
    column_names_rot = 45,
    
    row_title = 'Sender',
    row_title_side = 'left',
    row_title_gp = gpar(fontsize = 16, fontface = "bold"),
    row_title_rot = 90,
    
    column_title = 'Receiver',
    column_title_side = "bottom",
    column_title_gp = gpar(fontsize = 16, fontface = "bold"),
    
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 14, fontface = "bold"),
      labels_gp = gpar(fontsize = 12)
    ),
    
    border = TRUE,
    rect_gp = gpar(col = "white", lwd = 1)
  )
  
  # Calculate position from top
  y_pos <- 0.96 - (j - 1) * (height_each + spacing)
  
  pushViewport(viewport(
    x = 0.5,
    y = y_pos,
    width = 0.88,
    height = height_each,
    just = c("center", "top")
  ))
  
  draw(h, 
       column_title = paste('Pattern Cluster', j),
       column_title_gp = gpar(fontsize = 18, fontface = "bold"),
       padding = unit(c(8, 8, 8, 8), "mm"),
       newpage = FALSE)
  
  popViewport()
}

dev.off()

cat("Saved: pattern_clusters_grid_layout.png\n")

# Also save as PDF (publication quality)
pdf(file.path(output_dir, "fig_pattern_clusters_grid_layout.pdf"), width = 15, height = 20, bg = "white")


grid.newpage()

n_clusters <- length(sort(unique(interaction_group), decreasing = FALSE))

spacing <- 0.06
height_available <- 0.92
height_each <- (height_available - spacing * (n_clusters - 1)) / n_clusters

for(j in 1:n_clusters) {
  net_aggregated_group2 <- net_aggregation(
    net12_list[names(interaction_group[interaction_group == j])], 
    method = 'weight'
  )
  
  col_map <- brewer.pal(8, "YlOrBr")
  
  h <- Heatmap(
    net_aggregated_group2, 
    name = "Weight",
    col = col_map,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    
    width = unit(14, "cm"),
    height = unit(7, "cm"),
    
    row_names_side = 'left',
    row_names_gp = gpar(fontsize = 13),
    column_names_side = 'bottom',
    column_names_gp = gpar(fontsize = 13),
    column_names_rot = 45,
    
    row_title = 'Sender',
    row_title_side = 'left',
    row_title_gp = gpar(fontsize = 16, fontface = "bold"),
    row_title_rot = 90,
    
    column_title = 'Receiver',
    column_title_side = "bottom",
    column_title_gp = gpar(fontsize = 16, fontface = "bold"),
    
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 14, fontface = "bold"),
      labels_gp = gpar(fontsize = 12)
    ),
    
    border = TRUE,
    rect_gp = gpar(col = "white", lwd = 1)
  )
  
  y_pos <- 0.96 - (j - 1) * (height_each + spacing)
  
  pushViewport(viewport(
    x = 0.5,
    y = y_pos,
    width = 0.88,
    height = height_each,
    just = c("center", "top")
  ))
  
  draw(h, 
       column_title = paste('Pattern Cluster', j),
       column_title_gp = gpar(fontsize = 18, fontface = "bold"),
       padding = unit(c(8, 8, 8, 8), "mm"),
       newpage = FALSE)
  
  popViewport()
}

dev.off()

cat("Saved: pattern_clusters_grid_layout.pdf\n")




cat("

What the numerical values in the heatmap are

Each number in the heatmap is a communication weight between a sender cell type (row) and a receiver cell type (column), 
aggregated across all interaction programs in a given pattern cluster.

In symbols:

Where the numbers come from (step by step)

1️⃣ Individual pathway networks (before aggregation)

For each pathway / interaction program (e.g. Glu_Grin3a, Vip_Vipr1, Nrxn1_Nlgn3), NeuronChat computes a sender × receiver matrix:

Rows = sender cell types

Columns = receiver cell types

Values = inferred communication strength
(based on expression of ligands, receptors, cofactors, etc.)

Example for one pathway:

          L2/3 IT   L4/5 IT   L5 IT
L2/3 IT     0.10      0.03     0.01
L4/5 IT     0.05      0.20     0.08
L5 IT       0.02      0.06     0.15

2️⃣ Pathways are grouped into clusters

Your earlier steps (netClustering) grouped pathways into functional pattern clusters:

Cluster 1: mostly glutamatergic synaptic signaling

Cluster 2: neuromodulatory peptides

etc.

Each cluster contains multiple pathway-specific matrices.

3️⃣ net_aggregation(..., method = weight)

This is the key line:

net_aggregated_group2 <- net_aggregation(..., method = weight)


What it does:

Takes all pathway matrices in that cluster

Sums (or equivalently aggregates) the edge weights

Produces one matrix per cluster

So the value you see is:

Total signaling strength from sender → receiver carried by all pathways in that cluster

How to interpret a single cell in the heatmap

Example:

Sender:   L4/5 IT

Receiver: L5 IT

Value:    0.8

This means:

Across all pathways in this pattern cluster, L4/5 IT neurons send strong cumulative signals to L5 IT neurons.

It does not mean:

probability of a synapse

number of cells

gene expression level

It does mean:

inferred cell–cell communication intensity

summed over multiple ligand–receptor programs

What the magnitude means

Absolute scale

Values are unitless

They depend on:

expression levels

number of pathways in the cluster

how NeuronChat normalizes internally

So:

0.8 vs 0.4 → roughly 2× stronger signaling

0 → no inferred signaling

Relative interpretation (the correct way)

The heatmap is best read comparatively:

Within a panel:

Which sender–receiver pairs dominate?

Across panels:

Which pattern cluster emphasizes which layers/cell types?

Directionality matters

Rows → columns means:

Sender (row)  →  Receiver (column)


So the matrix is directed, not symmetric.

L2/3 → L5 can be strong

L5 → L2/3 can be weak (or vice versa.

")



cat("

How These Values Are Calculated

Step 1: Individual Pathway Communication

For each ligand-receptor pair (e.g., Glu_Grin1--AD_CTL), NeuronChat calculates communication probability for every sender-receiver pair:

Example: Glu_Grin1--AD_CTL pathway
                L2/3 IT  L4/5 IT  L5 IT
L2/3 IT          0.1      0.05     0.02
L4/5 IT          0.03     0.2      0.08
L5 IT            0.01     0.1      0.15

Step 2: Aggregation Across Cluster

The code aggregates all pathways in cluster j:
rnet_aggregated_group2 <- net_aggregation(
  net12_list[names(interaction_group[interaction_group==j])], 
  method = weight
)

If cluster 1 contains 3 pathways: Glu_Grin1--AD_CTL, Glu_Grin2--AD_CTL, Glu_Grik1--AD_CTL

**The values are summed:**

Cluster 1 aggregated = Glu_Grin1 + Glu_Grin2 + Glu_Grik1

                L2/3 IT  L4/5 IT  L5 IT
L2/3 IT          0.5      0.2      0.1    (sum of 3 pathways)
L4/5 IT          0.3      0.8      0.4
L5 IT            0.1      0.5      0.6

## Biological Interpretation

### High Values (dark brown/orange)
- **Strong total signaling** through pathways in this cluster
- These sender-receiver pairs are **key communication hubs** for this pattern
- Example: [0.8] means L4/5 IT cells have strong autocrine signaling via these pathways

### Low Values (light yellow)
- **Weak or absent signaling** through these pathways
- Less important for this communication pattern
- Example: [0.1] means minimal L2/3 IT → L5 IT communication

### Diagonal Values (sender = receiver)
- **Autocrine signaling** - cells signaling to themselves
- Important for self-regulation

### Off-diagonal Values (sender ≠ receiver)
- **Paracrine signaling** - cells signaling to other cell types
- Important for cell-cell communication

## Color Scale
The YlOrBr (Yellow-Orange-Brown) palette maps values:
- **0 or very low** → Light yellow (almost white)
- **Medium** → Orange
- **High** → Dark brown

## Real Example from Your Data

If you see in Pattern Cluster 3:

                L5/6 NP CTX
L5 IT CTX       [dark brown]  → High value (e.g., 0.8+)

This means: L5 IT CTX neurons send very strong signals to L5/6 NP CTX neurons 
through the ligand-receptor pathways grouped in cluster 3.

")


# an example :

# pathway_differences[["Glu_Grin3a"]] = list(
#  AD_CTL = 8x8 matrix,           # Communication matrix for AD_CTL
#  AD_Dyslexia = 8x8 matrix,      # Communication matrix for AD_Dyslexia
#  Difference = 8x8 matrix,       # AD_Dyslexia - AD_CTL (element-wise)
#  AD_CTL_total = 12.5,           # sum(AD_CTL matrix)
#  AD_Dyslexia_total = 18.3,      # sum(AD_Dyslexia matrix)
#  Diff_total = 8.2,              # sum(abs(Difference matrix))
#  Status = "Both"                # "Both", "AD_CTL only", or "AD_Dyslexia only"
# )

# What Each Column Represents

# Pathway
# names(pathway_differences)
# → Pathway name (e.g., "Glu_Grin3a")

# AD_CTL_strength
# sum(AD_CTL matrix)
# → Total communication strength in AD_CTL across all cell type pairs

# AD_Dyslexia_strength
# sum(AD_Dyslexia matrix)
# → Total communication strength in AD_Dyslexia across all cell type pairs

# Total_diff
# sum(abs(Difference matrix))
# → Total absolute change (how much the pathway changed overall)

# Status
# Metadata
# → Indicates whether the pathway exists in:

# "Both"
# "AD_CTL only"
# "AD_Dyslexia only"

# ============================================================================
# 1. EXTRACT UNIQUE PATHWAY NAMES
# ============================================================================

# Get pathway names from net1 and net2 (which have suffixes --AD_CTL and --AD_Dyslexia)
pathways_net1 <- names(net1)
pathways_net2 <- names(net2)

# Remove suffixes to get base pathway names
base_names_net1 <- gsub("--AD_CTL$", "", pathways_net1)
base_names_net2 <- gsub("--AD_Dyslexia$", "", pathways_net2)

# Get union of all unique pathway names (without suffixes)
unique_pathways <- unique(c(base_names_net1, base_names_net2))

cat("Total unique pathways to compare:", length(unique_pathways), "\n\n")


# ============================================================================
# 2. EXTRACT AND COMPARE MATRICES FOR EACH UNIQUE PATHWAY
# ============================================================================

# Create comparison for each pathway
pathway_differences <- list()

for(pathway in unique_pathways) {
  
  # Construct the full names with suffixes
  adctl_name <- paste0(pathway, "--AD_CTL")
  addys_name <- paste0(pathway, "--AD_Dyslexia")
  
  # Check if both exist
  adctl_exists <- adctl_name %in% names(net1)
  addys_exists <- addys_name %in% names(net2)
  
  cat("Processing:", pathway, "\n")
  cat("  AD_CTL name:", adctl_name, "- Exists:", adctl_exists, "\n")
  cat("  AD_Dyslexia name:", addys_name, "- Exists:", addys_exists, "\n")
  
  # Extract matrices
  if(adctl_exists && addys_exists) {
    # Both exist - extract and compare
    adctl_matrix <- net1[[adctl_name]]
    addys_matrix <- net2[[addys_name]]
    
    # Calculate difference (AD_Dyslexia - AD_CTL)
    # diff_matrix <- adctl_matrix - addys_matrix  # Old version (AD_CTL - AD_Dyslexia)

    diff_matrix <- addys_matrix - adctl_matrix 
    
    # Store
    pathway_differences[[pathway]] <- list(
      AD_CTL = adctl_matrix,
      AD_Dyslexia = addys_matrix,
      Difference = diff_matrix,
      AD_CTL_total = sum(adctl_matrix),
      AD_Dyslexia_total = sum(addys_matrix),
      Diff_total = sum(abs(diff_matrix)),
      Status = "Both"
    )
    
    cat("  ✓ Both matrices extracted and subtracted\n\n")
    
  } else if(adctl_exists && !addys_exists) {
    # Only AD_CTL exists
    adctl_matrix <- net1[[adctl_name]]
    
    pathway_differences[[pathway]] <- list(
      AD_CTL = adctl_matrix,
      AD_Dyslexia = NULL,
      Difference = adctl_matrix,  # All difference is from AD_CTL
      AD_CTL_total = sum(adctl_matrix),
      AD_Dyslexia_total = 0,
      Diff_total = sum(abs(adctl_matrix)),
      Status = "AD_CTL only"
    )
    
    cat("  ⚠ Only AD_CTL exists\n\n")
    
  } else if(!adctl_exists && addys_exists) {
    # Only AD_Dyslexia exists
    addys_matrix <- net2[[addys_name]]
    
    pathway_differences[[pathway]] <- list(
      AD_CTL = NULL,
      AD_Dyslexia = addys_matrix,
      Difference = -addys_matrix,  # All difference is from AD_Dyslexia (negative)
      AD_CTL_total = 0,
      AD_Dyslexia_total = sum(addys_matrix),
      Diff_total = sum(abs(addys_matrix)),
      Status = "AD_Dyslexia only"
    )
    
    cat("  ⚠ Only AD_Dyslexia exists\n\n")
    
  } else {
    # Neither exists (shouldn't happen if unique_pathways is correct)
    cat("  ✗ ERROR: Neither exists!\n\n")
    
    pathway_differences[[pathway]] <- list(
      AD_CTL = NULL,
      AD_Dyslexia = NULL,
      Difference = NULL,
      AD_CTL_total = 0,
      AD_Dyslexia_total = 0,
      Diff_total = 0,
      Status = "Neither"
    )
  }
}

# pathway_differences
# str(pathway_differences)
# head(pathway_differences)
# tail(pathway_differences)


pathway_scores <- data.frame(
  Pathway = names(pathway_differences),
  AD_CTL_strength = sapply(pathway_differences, function(x) x$AD_CTL_total),
  AD_Dyslexia_strength = sapply(pathway_differences, function(x) x$AD_Dyslexia_total),
  Total_diff = sapply(pathway_differences, function(x) x$Diff_total),
  Status = sapply(pathway_differences, function(x) x$Status)
) %>%
  arrange(desc(Total_diff))

                      # View top differences
cat("\n=== Top 30 Pathways by Difference ===\n")
print(head(pathway_scores, 30))

# Save summary
write.csv(pathway_scores, file.path(output_dir, "res_pathway_comparison_summary.csv"), row.names = FALSE)

cat("\nSaved: pathway_comparison_summary.csv\n")

# ============================================================================
# CREATE pathway_comparison DATA FRAME FOR SCATTER PLOTS
# ============================================================================
# This data frame contains pathway-level metrics for creating scatter plots
# comparing AD_CTL vs AD_Dyslexia pathway strengths

cat("\n=== Creating pathway_comparison data frame ===\n")

# Initialize empty data frame with proper column structure
pathway_comparison <- data.frame(
  pathway = character(),
  AD_CTL_strength = numeric(),
  AD_Dyslexia_strength = numeric(),
  difference = numeric(),
  fold_change = numeric(),
  log2FC = numeric(),
  stringsAsFactors = FALSE
)

# Loop through all unique pathways and calculate metrics
for(pathway in unique_pathways) {
  
  # Construct the full names with suffixes
  adctl_name <- paste0(pathway, "--AD_CTL")
  addys_name <- paste0(pathway, "--AD_Dyslexia")
  
  # Get pathway matrices (or zero if missing)
  if(adctl_name %in% names(net1)) {
    adctl_matrix <- net1[[adctl_name]]
    adctl_total <- sum(adctl_matrix, na.rm = TRUE)
  } else {
    adctl_total <- 0
  }
  
  if(addys_name %in% names(net2)) {
    addys_matrix <- net2[[addys_name]]
    addys_total <- sum(addys_matrix, na.rm = TRUE)
  } else {
    addys_total <- 0
  }
  
  # Calculate metrics
  difference <- addys_total - adctl_total
  fold_change <- (addys_total + 1e-9) / (adctl_total + 1e-9)  # Avoid division by zero
  log2FC <- log2(fold_change)
  
  # Add to data frame
  pathway_comparison <- rbind(pathway_comparison, data.frame(
    pathway = pathway,
    AD_CTL_strength = adctl_total,
    AD_Dyslexia_strength = addys_total,
    difference = difference,
    fold_change = fold_change,
    log2FC = log2FC,
    stringsAsFactors = FALSE
  ))
}

# Sort by absolute difference
pathway_comparison <- pathway_comparison[order(abs(pathway_comparison$difference), decreasing = TRUE), ]

cat("  Created pathway_comparison with", nrow(pathway_comparison), "pathways\n")
cat("  Top 5 pathways by difference:\n")
print(head(pathway_comparison, 5))

# ============================================================================
# 4. SUMMARY STATISTICS
# ============================================================================

cat("\n=== PATHWAY STATUS SUMMARY ===\n")
status_table <- table(pathway_scores$Status)
print(status_table)

cat("\nPathways with differences > 0:", sum(pathway_scores$Total_diff > 0), "\n")
cat("Pathways with no differences:", sum(pathway_scores$Total_diff == 0), "\n")

head(pathway_scores, 50)

cat("\nPathways with differences < 0:", sum(pathway_scores$Total_diff < 0), "\n")
cat("Pathways with no differences:", sum(pathway_scores$Total_diff == 0), "\n")

tail(pathway_scores, 50)

sum(pathway_scores$Total_diff < 0)



# ============================================================================
# 5. RANK PATHWAYS BY DIFFERENCE
# ============================================================================

# Create ranking based on Total_diff (already sorted in descending order)
pathway_scores_ranked <- pathway_scores %>%
  filter(Total_diff > 0) %>%
  mutate(Rank = row_number()) %>%
  arrange(Rank)

cat("\n=== PATHWAY RANKING BY DIFFERENCE ===\n")
cat("Total pathways with differences > 0:", nrow(pathway_scores_ranked), "\n\n")

# Display top 20 ranked pathways
cat("Top 50 ranked pathways:\n")
print(head(pathway_scores_ranked[, c("Rank", "Pathway", "Total_diff", "AD_CTL_strength", "AD_Dyslexia_strength")], 50))

# Save ranked pathways
write.csv(pathway_scores_ranked, file.path(output_dir, "res_pathway_ranking_by_difference.csv"), row.names = FALSE)

cat("\nSaved: pathway_ranking_by_difference.csv\n")








# ============================================================================
# 6. CREATE HEATMAPS FOR RANKED PATHWAYS
# ============================================================================

cat("\n=== Creating heatmaps for", nrow(pathway_scores_ranked), "ranked pathways ===\n\n")

# Create directory for heatmaps inside main output directory
heatmap_dir <- file.path(output_dir, "pathway_difference_heatmaps")
dir.create(heatmap_dir, showWarnings = FALSE, recursive = TRUE)

# Loop through ranked pathways
files_created <- 0
files_skipped <- 0

for(i in 1:nrow(pathway_scores_ranked)) {
  
  rank <- pathway_scores_ranked$Rank[i]
  pathway <- pathway_scores_ranked$Pathway[i]
  total_diff <- pathway_scores_ranked$Total_diff[i]
  
  cat("Processing Rank", rank, ":", pathway, "(Total diff =", round(total_diff, 3), ")\n")
  
  # Check if pathway exists in pathway_differences
  if(!pathway %in% names(pathway_differences)) {
    cat("  ⚠ Skipping - pathway not found in pathway_differences\n\n")
    files_skipped <- files_skipped + 1
    next
  }
  
  # Get difference matrix
  diff_mat <- pathway_differences[[pathway]]$Difference
  
  # Skip if difference matrix is NULL
  if(is.null(diff_mat)) {
    cat("  ⚠ Skipping - no difference matrix\n\n")
    files_skipped <- files_skipped + 1
    next
  }
  
  # Check if matrix is valid
  if(!is.matrix(diff_mat) || nrow(diff_mat) == 0 || ncol(diff_mat) == 0) {
    cat("  ⚠ Skipping - invalid matrix dimensions\n\n")
    files_skipped <- files_skipped + 1
    next
  }
  
  # Get max absolute difference for symmetric color scale
  max_diff <- max(abs(diff_mat), na.rm = TRUE)
  
  # Skip if max_diff is 0 or invalid
  if(is.na(max_diff) || max_diff == 0 || !is.finite(max_diff)) {
    cat("  ⚠ Skipping - max_diff is 0 or invalid (", max_diff, ")\n\n")
    files_skipped <- files_skipped + 1
    next
  }
  
  # Create filename: rank_number_pathway_name.png
  rank_padded <- sprintf("%03d", rank)  # Pads to 3 digits (001, 002, etc.)
  # Sanitize pathway name for filename (replace special characters with underscores)
  pathway_clean <- gsub("[^A-Za-z0-9_-]", "_", pathway)
  filename <- file.path(heatmap_dir, paste0("rank_", rank_padded, "_", pathway_clean, ".png"))
  
  # Try to create heatmap with error handling
  tryCatch({
    # Save as PNG with white background - wider to accommodate legend
    png(filename, width = 1800, height = 1400, res = 150, bg = "white")
    
    # Create color palette: Blue (negative) -> White (zero) -> Red (positive)
    col_palette <- colorRampPalette(c("blue", "white", "red"))(100)
    
    # Create breaks for color scale
    breaks <- seq(-max_diff, max_diff, length.out = 101)
    
    # Create heatmap using heatmap.2 (optimized for visible legend)
    heatmap.2(diff_mat,
              main = paste("Rank", rank, "-", pathway, "\nDifference (AD_Dyslexia - AD_CTL) | Total Diff =", round(total_diff, 3)),
              col = col_palette,
              breaks = breaks,
              dendrogram = "none",  # No clustering
              Rowv = FALSE,
              Colv = FALSE,
              trace = "none",  # No trace lines
              key = TRUE,  # Show color key/legend - CRITICAL FOR VISIBILITY
              keysize = 1.5,  # Larger legend size
              key.title = "Difference\n(AD_Dyslexia - AD_CTL)",
              key.xlab = "",
              key.ylab = "",
              key.par = list(cex = 0.8),  # Legend text size
              density.info = "none",  # No density plot in legend
              margins = c(12, 12),  # Bottom and right margins (increased)
              cexRow = 0.9,  # Row label size
              cexCol = 0.9,  # Column label size
              srtCol = 45,  # Column label angle
              lmat = rbind(c(0,3,0), c(2,1,0), c(0,4,0)),  # Layout: legend on right (3 columns)
              lhei = c(0.5, 4, 0.5),  # Row heights
              lwid = c(0.5, 4, 1.2))  # Column widths (wider for legend)
    
    dev.off()
    
    cat("  ✓ Saved:", filename, "\n\n")
    files_created <- files_created + 1
    
  }, error = function(e) {
    # Close device if still open
    if(dev.cur() != 1) dev.off()
    cat("  ✗ ERROR creating heatmap:", e$message, "\n\n")
    files_skipped <- files_skipped + 1
  })
}

cat("\n=============================================================\n")
cat("HEATMAP CREATION SUMMARY\n")
cat("=============================================================\n")
cat("Files created:", files_created, "\n")
cat("Files skipped:", files_skipped, "\n")
cat("Total pathways processed:", nrow(pathway_scores_ranked), "\n")
cat("=============================================================\n")

cat("\n=============================================================\n")
cat("SUMMARY\n")
cat("=============================================================\n")
cat("Total ranked pathways:", nrow(pathway_scores_ranked), "\n")
cat("Files saved to:", heatmap_dir, "\n")
cat("Filename format: rank_001_PathwayName.png, rank_002_PathwayName.png, ...\n")
cat("Ranking CSV: pathway_ranking_by_difference.csv\n")
cat("=============================================================\n")




pathway_scores$Delta_strength <- pathway_scores$AD_Dyslexia_strength - pathway_scores$AD_CTL_strength

head(pathway_scores)
tail(pathway_scores)



# 🔺 Up in AD_Dyslexia

up_in_dyslexia <- pathway_scores %>%
  filter(Delta_strength > 0) %>%
  arrange(desc(Delta_strength))

head(up_in_dyslexia, 20)

up_in_ctl <- pathway_scores %>%
  filter(Delta_strength < 0) %>%
  arrange(Delta_strength)   # most negative first

head(up_in_ctl, 20)


# ============================================================================
# CREATE HEATMAPS BASED ON DELTA_STRENGTH (USING BOTH pheatmap AND heatmap.2)
# ============================================================================

cat("\n=== Creating heatmaps based on Delta_strength ===\n\n")

# Create ranking based on absolute Delta_strength
pathway_scores_delta_ranked <- pathway_scores %>%
  filter(abs(Delta_strength) > 0) %>%
  arrange(desc(abs(Delta_strength))) %>%
  mutate(Rank = row_number()) %>%
  arrange(Rank)

cat("Total pathways with |Delta_strength| > 0:", nrow(pathway_scores_delta_ranked), "\n")
cat("Top 10 pathways by |Delta_strength|:\n")
print(head(pathway_scores_delta_ranked[, c("Rank", "Pathway", "Delta_strength", "AD_CTL_strength", "AD_Dyslexia_strength")], 10))

# Save ranked pathways by Delta_strength
write.csv(pathway_scores_delta_ranked, file.path(output_dir, "res_pathway_ranking_by_delta_strength.csv"), row.names = FALSE)
cat("\nSaved: pathway_ranking_by_delta_strength.csv\n\n")

# Create directory for Delta_strength heatmaps
delta_heatmap_dir <- file.path(output_dir, "pathway_delta_strength_heatmaps")
dir.create(delta_heatmap_dir, showWarnings = FALSE, recursive = TRUE)

# Initialize counters
files_created_pheatmap_delta <- 0
files_created_heatmap2_delta <- 0
files_skipped_delta <- 0

# Loop through ranked pathways by Delta_strength
for(i in 1:nrow(pathway_scores_delta_ranked)) {
  
  rank <- pathway_scores_delta_ranked$Rank[i]
  pathway <- pathway_scores_delta_ranked$Pathway[i]
  delta_strength <- pathway_scores_delta_ranked$Delta_strength[i]
  
  cat("Processing Rank", rank, ":", pathway, "(Delta_strength =", round(delta_strength, 3), ")\n")
  
  # Check if pathway exists in pathway_differences
  if(!pathway %in% names(pathway_differences)) {
    cat("  ⚠ Skipping - pathway not found in pathway_differences\n\n")
    files_skipped_delta <- files_skipped_delta + 1
    next
  }
  
  # Get difference matrix
  diff_mat <- pathway_differences[[pathway]]$Difference
  
  # Skip if difference matrix is NULL
  if(is.null(diff_mat)) {
    cat("  ⚠ Skipping - no difference matrix\n\n")
    files_skipped_delta <- files_skipped_delta + 1
    next
  }
  
  # Check if matrix is valid
  if(!is.matrix(diff_mat) || nrow(diff_mat) == 0 || ncol(diff_mat) == 0) {
    cat("  ⚠ Skipping - invalid matrix dimensions\n\n")
    files_skipped_delta <- files_skipped_delta + 1
    next
  }
  
  # Get max absolute difference for symmetric color scale
  max_diff <- max(abs(diff_mat), na.rm = TRUE)
  
  # Skip if max_diff is 0 or invalid
  if(is.na(max_diff) || max_diff == 0 || !is.finite(max_diff)) {
    cat("  ⚠ Skipping - max_diff is 0 or invalid (", max_diff, ")\n\n")
    files_skipped_delta <- files_skipped_delta + 1
    next
  }
  
  # Create filename base
  rank_padded <- sprintf("%03d", rank)
  pathway_clean <- gsub("[^A-Za-z0-9_-]", "_", pathway)
  
  # ========================================================================
  # VERSION 1: Using pheatmap
  # ========================================================================
  filename_pheatmap <- file.path(delta_heatmap_dir, paste0("rank_", rank_padded, "_", pathway_clean, "_pheatmap.png"))
  
  tryCatch({
    # Save as PNG with white background - wider to accommodate legend
    png(filename_pheatmap, width = 1800, height = 1400, res = 150, bg = "white")
    
    # Create heatmap using pheatmap
    pheatmap::pheatmap(diff_mat,
             main = paste("Rank", rank, "-", pathway, "\nDifference (AD_Dyslexia - AD_CTL) | Delta_strength =", round(delta_strength, 3)),
             color = colorRampPalette(c("blue", "white", "red"))(100),
             breaks = seq(-max_diff, max_diff, length.out = 101),
             cluster_rows = FALSE,
             cluster_cols = FALSE,
             fontsize = 12,
             fontsize_row = 10,
             fontsize_col = 10,
             border_color = "grey80",
             angle_col = 45,
             display_numbers = FALSE,
             cellwidth = 40,
             cellheight = 40,
             legend = TRUE,
             legend_breaks = c(-max_diff, -max_diff/2, 0, max_diff/2, max_diff),
             legend_labels = c(paste0(round(-max_diff, 3)), 
                              paste0(round(-max_diff/2, 3)), 
                              "0", 
                              paste0(round(max_diff/2, 3)), 
                              paste0(round(max_diff, 3))))
    
    dev.off()
    cat("  ✓ Saved pheatmap version:", basename(filename_pheatmap), "\n")
    files_created_pheatmap_delta <- files_created_pheatmap_delta + 1
    
  }, error = function(e) {
    if(dev.cur() != 1) dev.off()
    cat("  ✗ ERROR creating pheatmap:", e$message, "\n")
  })
  
  # ========================================================================
  # VERSION 2: Using heatmap.2
  # ========================================================================
  filename_heatmap2 <- file.path(delta_heatmap_dir, paste0("rank_", rank_padded, "_", pathway_clean, "_heatmap2.png"))
  
  tryCatch({
    # Save as PNG with white background - wider to accommodate legend
    png(filename_heatmap2, width = 1800, height = 1400, res = 150, bg = "white")
    
    # Create color palette: Blue (negative) -> White (zero) -> Red (positive)
    col_palette <- colorRampPalette(c("blue", "white", "red"))(100)
    
    # Create breaks for color scale
    breaks <- seq(-max_diff, max_diff, length.out = 101)
    
    # Create heatmap using heatmap.2 (optimized for visible legend)
    heatmap.2(diff_mat,
              main = paste("Rank", rank, "-", pathway, "\nDifference (AD_Dyslexia - AD_CTL) | Delta_strength =", round(delta_strength, 3)),
              col = col_palette,
              breaks = breaks,
              dendrogram = "none",  # No clustering
              Rowv = FALSE,
              Colv = FALSE,
              trace = "none",  # No trace lines
              key = TRUE,  # Show color key/legend - CRITICAL FOR VISIBILITY
              keysize = 1.5,  # Larger legend size
              key.title = "Difference\n(AD_Dyslexia - AD_CTL)",
              key.xlab = "",
              key.ylab = "",
              key.par = list(cex = 0.8),  # Legend text size
              density.info = "none",  # No density plot in legend
              margins = c(12, 12),  # Bottom and right margins (increased)
              cexRow = 0.9,  # Row label size
              cexCol = 0.9,  # Column label size
              srtCol = 45,  # Column label angle
              lmat = rbind(c(0,3,0), c(2,1,0), c(0,4,0)),  # Layout: legend on right (3 columns)
              lhei = c(0.5, 4, 0.5),  # Row heights
              lwid = c(0.5, 4, 1.2))  # Column widths (wider for legend)
    
    dev.off()
    cat("  ✓ Saved heatmap.2 version:", basename(filename_heatmap2), "\n\n")
    files_created_heatmap2_delta <- files_created_heatmap2_delta + 1
    
  }, error = function(e) {
    if(dev.cur() != 1) dev.off()
    cat("  ✗ ERROR creating heatmap.2:", e$message, "\n\n")
  })
}

cat("\n=============================================================\n")
cat("DELTA_STRENGTH HEATMAP CREATION SUMMARY\n")
cat("=============================================================\n")
cat("Files created (pheatmap):", files_created_pheatmap_delta, "\n")
cat("Files created (heatmap.2):", files_created_heatmap2_delta, "\n")
cat("Files skipped:", files_skipped_delta, "\n")
cat("Total pathways processed:", nrow(pathway_scores_delta_ranked), "\n")
cat("Output directory:", delta_heatmap_dir, "\n")
cat("=============================================================\n\n")


# ============================================================================
# CREATE HEATMAPS FOR ALL PATHWAYS WITH DIFFERENCES > 0
# ============================================================================

# Filter pathways with differences
pathways_with_diff <- pathway_scores$Pathway[pathway_scores$Total_diff > 0]

cat("\n=== Creating heatmaps for", length(pathways_with_diff), "pathways with differences ===\n\n")

# Create directory for heatmaps inside output directory (reuse existing)
# heatmap_dir already created above

files_created_alt <- 0
files_skipped_alt <- 0

for(pathway in pathways_with_diff) {
  
  cat("Creating heatmap for:", pathway, "\n")
  
  # Check if pathway exists in pathway_differences
  if(!pathway %in% names(pathway_differences)) {
    cat("  ⚠ Skipping - pathway not found in pathway_differences\n")
    files_skipped_alt <- files_skipped_alt + 1
    next
  }
  
  # Get difference matrix
  diff_mat <- pathway_differences[[pathway]]$Difference
  
  # Skip if difference matrix is NULL
  if(is.null(diff_mat)) {
    cat("  ⚠ Skipping - no difference matrix\n")
    files_skipped_alt <- files_skipped_alt + 1
    next
  }
  
  # Check if matrix is valid
  if(!is.matrix(diff_mat) || nrow(diff_mat) == 0 || ncol(diff_mat) == 0) {
    cat("  ⚠ Skipping - invalid matrix dimensions\n")
    files_skipped_alt <- files_skipped_alt + 1
    next
  }
  
  # Get max absolute difference for symmetric color scale
  max_diff <- max(abs(diff_mat), na.rm = TRUE)
  
  # Skip if max_diff is 0 or invalid
  if(is.na(max_diff) || max_diff == 0 || !is.finite(max_diff)) {
    cat("  ⚠ Skipping - max_diff is 0 or invalid (", max_diff, ")\n")
    files_skipped_alt <- files_skipped_alt + 1
    next
  }
  
  # Create filename (sanitize pathway name)
  filename <- file.path(heatmap_dir, 
                     paste0(gsub("[^A-Za-z0-9]", "_", pathway), 
                     "_difference.png"))
  
  # Try to create heatmap with error handling
  tryCatch({
    # Save as PNG with white background - wider to accommodate legend
    png(filename, width = 1800, height = 1400, res = 150, bg = "white")
    
    # Create color palette: Blue (negative) -> White (zero) -> Red (positive)
    col_palette <- colorRampPalette(c("blue", "white", "red"))(100)
    
    # Create breaks for color scale
    breaks <- seq(-max_diff, max_diff, length.out = 101)
    
    # Create heatmap using heatmap.2 (optimized for visible legend)
    heatmap.2(diff_mat,
              main = paste("Difference (AD_Dyslexia - AD_CTL):", pathway),
              col = col_palette,
              breaks = breaks,
              dendrogram = "none",  # No clustering
              Rowv = FALSE,
              Colv = FALSE,
              trace = "none",  # No trace lines
              key = TRUE,  # Show color key/legend - CRITICAL FOR VISIBILITY
              keysize = 1.5,  # Larger legend size
              key.title = "Difference\n(AD_Dyslexia - AD_CTL)",
              key.xlab = "",
              key.ylab = "",
              key.par = list(cex = 0.8),  # Legend text size
              density.info = "none",  # No density plot in legend
              margins = c(12, 12),  # Bottom and right margins (increased)
              cexRow = 0.9,  # Row label size
              cexCol = 0.9,  # Column label size
              srtCol = 45,  # Column label angle
              lmat = rbind(c(0,3,0), c(2,1,0), c(0,4,0)),  # Layout: legend on right (3 columns)
              lhei = c(0.5, 4, 0.5),  # Row heights
              lwid = c(0.5, 4, 1.2))  # Column widths (wider for legend)
    
    dev.off()
    
    cat("  ✓ Saved:", filename, "\n")
    files_created_alt <- files_created_alt + 1
    
  }, error = function(e) {
    # Close device if still open
    if(dev.cur() != 1) dev.off()
    cat("  ✗ ERROR creating heatmap:", e$message, "\n")
    files_skipped_alt <- files_skipped_alt + 1
  })
}

cat("\n=== All heatmaps saved to:", heatmap_dir, "===\n")
cat("Files created:", files_created_alt, "\n")
cat("Files skipped:", files_skipped_alt, "\n")
cat("Total pathways processed:", length(pathways_with_diff), "\n")

cat("\n=== All heatmaps saved to:", heatmap_dir, "===\n")
cat("Total files created:", length(pathways_with_diff), "\n")

getwd()





# net12[[1]] 
# net12[[2]] 



# This code is creating a detailed comparison table that shows how every cell-to-cell communication changes between AD_CTL and AD_Dyslexia for all pathways.



# Claude AI

# For each of the 183 languages/pathways:

# a) Get the conversation tables:

# ALM table (who talks to whom in ALM)
# VISp table (who talks to whom in VISp)

# b) Focus on your selected cells (6×6 = 36 conversations):

# Cell 1 → Cell 1
# Cell 1 → Cell 2
# Cell 1 → Cell 3
# ...
# Cell 6 → Cell 6

# c) Write it out as a list:
# Instead of a 6×6 table, make a list:
# SpeakerListenerLanguageALM strengthVISp strengthL5 PT CTXL2/3 IT CTXGlu_Grin3a0.23.5L5 PT CTXL4/5 IT CTXGlu_Grin3a1.11.3...............

# d) Do this for all 183 pathways and stack them:

# 183 pathways × 36 conversations each = 6,588 rows total

# Calculate the differences (lines 47-51)
# For each row, add:

# delta = AD_Dyslexia - AD_CTL (How much stronger in AD_Dyslexia?)
# mean_w = Average of both
# log2FC = Fold change (ratio)

# Find the biggest changes (line 54)
# Sort by delta and show top 20:
# "Which 20 specific cell-to-cell conversations (using which language) changed the most between AD_CTL and AD_Dyslexia?"



# all_pathways





library(dplyr)
library(tidyr)
library(ggplot2)

library(dplyr)
library(tidyr)
library(ggplot2)

# ---- Clean pathway names and get union ----
# Remove the --AD_CTL and --AD_Dyslexia suffixes from pathway names
all_pathways_clean <- unique(gsub("--AD_CTL|--AD_Dyslexia", "", all_pathways))

cat("Total unique pathways (after cleaning and union):", length(all_pathways_clean), "\n")
cat("First few pathways:", paste(head(all_pathways_clean, 5), collapse = ", "), "\n")




# If pathway exists: extracts the actual communication matrix
# If pathway missing: creates a zero matrix (assumes no communication)
# This allows union analysis - missing pathways get zeros instead of being excluded

# METRICS

# delta: Raw difference (positive = increased in Dyslexia)
# mean_w: Average strength across conditions (for weighting)
# log2FC: Log2 fold-change (symmetrical metric, +1 = doubled, -1 = halved)
# 1e-6: Small constant to prevent division by zero

cat("

delta: Raw difference (positive = increased in Dyslexia)

* mean_w: Average strength across conditions (for weighting)

* log2FC: Log2 fold-change (symmetrical metric, +1 = doubled, -1 = halved)

* 1e-6: Small constant to prevent division by zero, give a precise example !!!

")

cat("

Communication Example

L23_IT_Exc → MGE_Inh via GABA_GABRA1 pathway

From your matrix, we can extract the cell types and create a concrete example.

Hypothetical Values

AD_CTL: L23_IT_Exc → MGE_Inh = 0.1129170

AD_Dyslexia: L23_IT_Exc → MGE_Inh = 0.4516680 (increased)

1️⃣ Delta (Raw Difference)

delta = AD_Dyslexia - AD_CTL

Calculation
delta = 0.4516680 - 0.1129170
      = 0.3387510

Interpretation

+0.3387510 means communication from L23_IT_Exc → MGE_Inh increased by 0.34 units in AD_Dyslexia.

This tells you the absolute change in communication strength.

Real-world meaning

Excitatory neurons in Layer 2/3 are sending more GABA-related signals to inhibitory interneurons from the medial ganglionic eminence in the dyslexic brain.

2️⃣ Mean_w (Average Strength)

mean_w = (AD_Dyslexia + AD_CTL) / 2

Calculation
mean_w = (0.4516680 + 0.1129170) / 2
       = 0.5645850 / 2
       = 0.2822925

Interpretation

The average communication strength across both conditions is 0.28.

This is used as a weight in downstream analysis.

Why Mean_w Matters

Compare two scenarios:

Connection	AD_CTL	AD_Dyslexia	delta	mean_w	Importance
L23_IT_Exc → MGE_Inh	0.113	0.452	+0.339	0.282	More important ✓
Astro → Micro	0.000001	0.339001	+0.339	0.170	Less important

Both have the same delta (+0.339), but:

The first has higher mean_w (0.282) → strong baseline connection that got stronger.

The second has lower mean_w (0.170) → appeared from almost nothing.

Biological interpretation

Changes in already-strong connections (high mean_w) are often more functionally relevant than weak connections appearing from zero.

3️⃣ Log2FC (Log2 Fold-Change)

log2FC = log2((AD_Dyslexia + 1e-6) / (AD_CTL + 1e-6))

Calculation
log2FC = log2((0.4516680 + 0.000001) / (0.1129170 + 0.000001))
       = log2(0.4516690 / 0.1129180)
       = log2(4.000)
       = 2.0

Interpretation

log2FC = 2.0 means communication is 2² = 4× stronger in AD_Dyslexia.

L23_IT_Exc → MGE_Inh communication quadrupled in dyslexia.

Why log2FC Is Better Than Raw Ratios
Scenario	AD_CTL	AD_Dyslexia	Raw Ratio	log2FC	Interpretation
Quadrupled	0.113	0.452	4.0	+2.0	4× stronger
Doubled	0.113	0.226	2.0	+1.0	2× stronger
No change	0.113	0.113	1.0	0	Same
Halved	0.226	0.113	0.5	-1.0	2× weaker
Quartered	0.452	0.113	0.25	-2.0	4× weaker

Notice the symmetry:
+1 and -1 represent the same magnitude change (doubling vs halving).

")






library(dplyr)
library(tidyr)
library(ggplot2)

# ---- Clean pathway names and get union ----
# Remove the --AD_CTL and --AD_Dyslexia suffixes from pathway names
all_pathways_clean <- unique(gsub("--AD_CTL|--AD_Dyslexia", "", all_pathways))

cat("Total unique pathways (after cleaning and union):", length(all_pathways_clean), "\n")
cat("First few pathways:", paste(head(all_pathways_clean, 5), collapse = ", "), "\n")

# ---- Extract networks from cortex_list ----
net_AD_CTL      <- cortex_list[[1]]@net  # AD_CTL
net_AD_Dyslexia <- cortex_list[[2]]@net  # AD_Dyslexia

# ---- Determine cell type ordering (rows/cols of the matrices) ----
# Use any one pathway matrix as a template:
any_path <- all_pathways_clean[1]

# Find which condition has this pathway to get cell type names
if (any_path %in% names(net_AD_CTL)) {
  celltypes <- rownames(net_AD_CTL[[any_path]])
} else {
  celltypes <- rownames(net_AD_Dyslexia[[any_path]])
}

# Use ALL cell types as both sources and targets
sources.use <- 1:length(celltypes)
targets.use <- 1:length(celltypes)

source_names <- celltypes[sources.use]
target_names <- celltypes[targets.use]

cat("Using ALL cell types:", length(celltypes), "\n")
cat("Cell types:", paste(celltypes, collapse = ", "), "\n")

# ---- Build long table of weights for ALL source/target pairs ----
cmp_df <- lapply(all_pathways_clean, function(p) {
  
  # Check if pathway exists in each condition
  has_AD_CTL <- p %in% names(net_AD_CTL)
  has_AD_Dyslexia <- p %in% names(net_AD_Dyslexia)
  
  # Get matrices (or create zero matrix if missing)
  if (has_AD_CTL) {
    A <- net_AD_CTL[[p]]
  } else {
    # Create zero matrix with correct dimensions
    A <- matrix(0, nrow = length(celltypes), ncol = length(celltypes))
    rownames(A) <- celltypes
    colnames(A) <- celltypes
  }
  
  if (has_AD_Dyslexia) {
    B <- net_AD_Dyslexia[[p]]
  } else {
    # Create zero matrix with correct dimensions
    B <- matrix(0, nrow = length(celltypes), ncol = length(celltypes))
    rownames(B) <- celltypes
    colnames(B) <- celltypes
  }
  
  # Align just in case (should usually match)
  common_ct <- intersect(rownames(A), rownames(B))
  A <- A[common_ct, common_ct, drop = FALSE]
  B <- B[common_ct, common_ct, drop = FALSE]
  
  # Subset to chosen sources/targets by name (robust)
  A_sub <- A[source_names, target_names, drop = FALSE]
  B_sub <- B[source_names, target_names, drop = FALSE]
  
  # Long format
  expand.grid(sender = rownames(A_sub), receiver = colnames(A_sub), stringsAsFactors = FALSE) %>%
    mutate(
      pathway = p,
      AD_CTL      = as.vector(A_sub),
      AD_Dyslexia = as.vector(B_sub)
    )
}) %>%
  bind_rows() %>%
  mutate(
    delta  = AD_Dyslexia - AD_CTL,                           # AD_Dyslexia minus AD_CTL (positive = AD_Dyslexia increased)
    mean_w = (AD_Dyslexia + AD_CTL)/2,
    log2FC = log2((AD_Dyslexia + 1e-6) / (AD_CTL + 1e-6))    # ratio (stabilized)
  )

# ---- Summary ----
cat("\nTotal rows in comparison table:", nrow(cmp_df), "\n")
cat("  (", length(celltypes), "senders ×", length(celltypes), "receivers ×", length(all_pathways_clean), "pathways )\n")


cmp_df


# This is the NeuronChat analogue of gg1$data
cat("\nTop 30 edges with largest increases (AD_Dyslexia > AD_CTL):\n")
cmp_df %>% arrange(desc(delta)) %>% head(30) %>% print()

cat("\nTop 30 edges with largest decreases (AD_Dyslexia < AD_CTL):\n")
cmp_df %>% arrange(delta) %>% head(30) %>% print()

# Save results
write.csv(cmp_df, file.path(output_dir, "res_cell_to_cell_comparison_AD_CTL_AD_Dyslexia_all_cells_union_pathways.csv"), row.names = FALSE)

cat("\n✓ Saved: res_cell_to_cell_comparison_AD_CTL_AD_Dyslexia_all_cells_union_pathways.csv\n")

head(cmp_df)

# As a reminder :

# Components of the weight:

# Ligand/neurotransmitter production - Expression of genes encoding:

# Neurotransmitter synthesis enzymes
# Vesicular transporters
# Release machinery

# Receptor expression - Expression level of the receptor in target cells
# Stoichiometry - Proper pairing of ligand-receptor components
# Statistical significance - The M=100 permutations test whether the communication is significant

# The resulting values are:

# Continuous weights (not discrete counts)
# Higher values = stronger inferred communication
# Zero = no detectable communication (either no expression or not significant)
# NOT normalized to sum to 1 or any fixed scale
# Each pathway has its own scale based on expression levels

dim(cmp_df)






# Set plot width and height (in inches for notebooks)
options(repr.plot.width = 22, repr.plot.height = 82)

# size_by <- "mean_w"

# You choose which column controls bubble size.
# Here:
# mean_w = (AD_Dyslexia + AD_CTL) / 2

# Bigger dots = stronger overall signaling (regardless of direction)
# you could also use:
# abs(delta) → magnitude of change
# AD_Dyslexia or AD_CTL → condition-specific strength



# Choose what you want bubbles sized by:
# - mean_w  (overall strength)
# - abs(delta) (magnitude of change)

head(cmp_df)
tail(cmp_df)




# Overall Purpose
# Creates one plot per pathway (e.g., GABA_GABRA1, 5HT_HTR1A, etc.) showing:

# Which sender-receiver pairs communicate via that pathway
# How much the communication changed (color)
# How strong the average communication is (bubble size)

library(ggplot2)
library(dplyr)

# output directory (inside main output directory)
plot_dir <- file.path(output_dir, "NeuronChat_pathway_plots")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# shared color scaling (global, consistent across pathways)
max_abs <- quantile(abs(cmp_df$delta), 0.95, na.rm = TRUE)

for (pw in unique(cmp_df$pathway)) {

  df_pw <- cmp_df %>% filter(pathway == pw)

  p_pw <- ggplot(df_pw, aes(x = receiver, y = sender)) +
    geom_point(
      aes(size = mean_w, color = delta),
      alpha = 0.9
    ) +
    scale_color_gradient2(
      low = "#2166ac",
      mid = "white",
      high = "#b2182b",
      midpoint = 0,
      limits = c(-max_abs, max_abs),
      oob = scales::squish,
      name = "AD_Dyslexia − AD_CTL"
    ) +
    scale_size_continuous(name = "Mean weight", range = c(2, 6)) +
    theme_bw(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      title = paste("NeuronChat:", pw, "(AD_Dyslexia − AD_CTL)"),
      x = "Receiver",
      y = "Sender"
    )

  png(
    file.path(plot_dir, paste0(pw, "_AD_Dyslexia_vs_AD_CTL.png")),
    width = 7, height = 6, units = "in", res = 300
  )
  print(p_pw)
  dev.off()
}





cat("

Overall Purpose

Aggregates all cell-to-cell connections within each pathway to answer:

Which pathways show the strongest overall change in AD_Dyslexia?

Step-by-Step Explanation

1️⃣ Filter and Group Data (Lines 2–5)

pathway_rank <- cmp_df %>%
filter(mean_w > 0.01) %>%
group_by(pathway)

Filter: mean_w > 0.01

Removes very weak connections (noise)

Keeps only edges with meaningful average communication strength

Example:
Excludes connections like:

Astro → High_mt mean_w = 0.00001

Why filter?

Weak connections are unreliable and may be artifacts

You want pathway-level conclusions based on biologically meaningful edges

Group by pathway

Groups all sender-receiver pairs belonging to the same pathway together

Example:
All GABA_GABRA1 edges (e.g., Oligo→L23_IT_Exc, L23_IT_Exc→MGE_Inh, etc.) are grouped together.

2️⃣ Calculate Pathway-Level Metrics (Lines 6–11)

summarise(
n_edges = n(),
log2FC_wmean = weighted.mean(log2FC, w = mean_w, na.rm = TRUE),
delta_sum = sum(delta, na.rm = TRUE),
.groups = drop
)

For each pathway, three metrics are calculated:

n_edges — Number of Edges

n_edges = n()

Counts how many sender-receiver pairs use this pathway

Example:
GABA_GABRA1 might have 120 edges
(e.g., 14 cell types × 14 cell types, excluding zeros)

log2FC_wmean — Weighted Mean Log2 Fold-Change ⭐ MOST IMPORTANT

log2FC_wmean = weighted.mean(log2FC, w = mean_w, na.rm = TRUE)

What this does

Calculates the average log2FC across all edges in the pathway

Weights each edge by its mean_w (average strength)

Strong connections contribute more to the pathway-level score

Example Calculation for GABA_GABRA1

Sender Receiver log2FC mean_w Contribution
L23_IT_Exc MGE_Inh +2.0 0.30 2.0 × 0.30 = 0.60
Oligo L4_Exc +0.5 0.10 0.5 × 0.10 = 0.05
L4_Exc CGE_Inh -1.0 0.20 -1.0 × 0.20 = -0.20
Astro Micro +5.0 0.001 5.0 × 0.001 = 0.005

Total weight:

0.30 + 0.10 + 0.20 + 0.001 = 0.601

Weighted mean:

(0.60 + 0.05 - 0.20 + 0.005) / 0.601 = +0.75

Interpretation

GABA_GABRA1 pathway is overall:

2^0.75 ≈ 1.68× stronger in AD_Dyslexia

Why Weighted?

Strong edge (L23_IT_Exc → MGE_Inh, mean_w=0.30) strongly influences result

Weak edge (Astro → Micro, mean_w=0.001) barely influences result

This prevents rare/weak edges from skewing pathway-level conclusions.

delta_sum — Sum of Raw Differences

delta_sum = sum(delta, na.rm = TRUE)

Sums all raw differences (AD_Dyslexia − AD_CTL)

Example:
If a pathway has 100 edges with deltas from -0.2 to +0.5, this adds them all together.

Why less important than log2FC_wmean?

Does not account for proportional changes

Does not normalize for connection strength

Can be biased by a few large raw values

")



# log2FC is computed per sender→receiver edge, not per pathway. 
# So to get “top up/down pathways”, you need to aggregate log2FC across edges within each pathway (and then rank).

# Also: the top rows have AD_CTL = 0, AD_Dyslexia > 0 → log2FC becomes huge (~18) because of the tiny pseudocount (1e-6). 
# That’s not “wrong”, but it will dominate rankings unless you filter.

# log2FC_wmean = the pathway-level, edge-strength-weighted average log2 fold-change between AD_Dyslexia and AD_CTL.
# It tells you:
# "On average, across all sender→receiver interactions in this pathway, is signaling stronger in AD_Dyslexia or in AD_CTL — 
# and by how much, giving more weight to strong interactions?”

# log2FC_wmean 

# Each edge has
# log2FCᵢ = log2((VISpᵢ + ε)/(ALMᵢ + ε))
# mean_wᵢ = (VISpᵢ + ALMᵢ)/2

library(dplyr)

pathway_rank <- cmp_df %>%
  # optional but recommended: ignore tiny/near-zero interactions
  filter(mean_w > 0.01) %>%
  group_by(pathway) %>%
  summarise(
    n_edges = n(),
    # weighted mean log2FC across edges in that pathway
    log2FC_wmean = weighted.mean(log2FC, w = mean_w, na.rm = TRUE),
    delta_sum = sum(delta, na.rm = TRUE),
    .groups = "drop"
  )

top_up   <- pathway_rank %>% arrange(desc(log2FC_wmean)) %>% head(30)
top_down <- pathway_rank %>% arrange(log2FC_wmean) %>% head(30)

top_up
top_down

# log2FC_wmean	Interpretation
# > 0	Pathway overall up-regulated in AD_Dyslexia
# < 0	Pathway overall down-regulated in AD_Dyslexia (AD_CTL-biased)
# ≈ 0	Balanced / mixed effects
# Large magnitude	Strong, consistent bias

display_tbl <- bind_rows(
  top_up   %>% mutate(direction = "UP (AD_Dyslexia > AD_CTL)"),
  top_down %>% mutate(direction = "DOWN (AD_Dyslexia < AD_CTL)")
) %>%
  arrange(direction, desc(abs(log2FC_wmean)))

display_tbl


write.csv(
  display_tbl,
  file.path(output_dir, "res_NeuronChat_top_pathways_AD_Dyslexia_vs_AD_CTL_considering_log2FC_mean_w.csv"),
  row.names = FALSE
)





library(ggplot2)

ggplot(display_tbl,
       aes(x = reorder(pathway, log2FC_wmean),
           y = log2FC_wmean,
           fill = direction)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(
    values = c(
      "UP (AD_Dyslexia > AD_CTL)"   = "#b2182b",
      "DOWN (AD_Dyslexia < AD_CTL)" = "#2166ac"
    )
  ) +
  theme_bw(base_size = 14) +
  labs(
    title = "Top differential NeuronChat pathways (AD_Dyslexia vs AD_CTL)",
    x = "Pathway",
    y = "Weighted mean log2FC"
  )


library(ggplot2)
library(dplyr)

# Symmetric limits around zero (prevents visual bias)
max_abs <- max(abs(display_tbl$log2FC_wmean), na.rm = TRUE)

p_pathways <- ggplot(
  display_tbl,
  aes(
    x = reorder(pathway, log2FC_wmean),
    y = log2FC_wmean,
    fill = direction
  )
) +
  geom_col(width = 0.7, color = "grey30", linewidth = 0.2) +
  coord_flip() +

  # Zero reference line
  geom_hline(yintercept = 0, linewidth = 0.4, color = "black") +

  # Symmetric scale
  scale_y_continuous(
    limits = c(-max_abs * 1.05, max_abs * 1.05),
    expand = expansion(mult = c(0, 0.02))
  ) +

  scale_fill_manual(
    values = c(
      "UP (AD_Dyslexia > AD_CTL)"   = "#b2182b",
      "DOWN (AD_Dyslexia < AD_CTL)" = "#2166ac"
    )
  ) +

  labs(
    title = "Differential NeuronChat pathways (AD_Dyslexia vs AD_CTL)",
    subtitle = "Pathway-level weighted mean log2 fold-change",
    x = NULL,
    y = "Weighted mean log2FC (AD_Dyslexia / AD_CTL)",
    fill = NULL
  ) +

  theme_bw(base_size = 14) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "grey85"),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 11),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 11),
    legend.position = "top",
    legend.justification = "left"
  )

p_pathways

png(
  file = file.path(output_dir, "fig_NeuronChat_pathway_log2FC_AD_Dyslexia_vs_AD_CTL.png"),
  width = 7,
  height = 8,
  units = "in",
  res = 300
)

print(p_pathways)

dev.off()








library(ggplot2)
library(ggrepel)
library(dplyr)
library(IRdisplay)

# Set plot size for Jupyter output
options(repr.plot.width = 12, repr.plot.height = 10)

# Create labeled data and assign colors
pathway_comparison_labeled <- pathway_comparison %>%
  mutate(
    label_text = ifelse(AD_CTL_strength > 1 & AD_Dyslexia_strength > 1, pathway, ""),
    color_group = case_when(
      log2FC > 0 ~ "Increased in AD_Dyslexia",
      log2FC < 0 ~ "Decreased in AD_Dyslexia",
      TRUE ~ "No change"
    )
  )

# Create the plot
pathway_plot <- ggplot(
  pathway_comparison_labeled,
  aes(x = AD_CTL_strength, y = AD_Dyslexia_strength)
) +
  geom_point(
    aes(color = color_group),
    size = 6,              # Large dots
    alpha = 0.9            # Nearly full opacity
  ) +
  scale_color_manual(
    values = c(
      "Increased in AD_Dyslexia" = "#ca0020",    # Red
      "Decreased in AD_Dyslexia" = "#0571b0",    # Blue
      "No change" = "grey70"                      # Grey
    ),
    name = "Change"
  ) +
  geom_text_repel(
    aes(label = label_text),
    size = 4,
    max.overlaps = Inf,
    box.padding = 0.5,
    point.padding = 0.4,
    segment.color = "grey70",
    segment.size = 0.3
  ) +
  geom_abline(
    slope = 1, intercept = 0,
    linetype = "dashed",
    color = "gray50",
    linewidth = 1
  ) +
  labs(
    title = "Pathway Strength Comparison: AD_CTL vs AD_Dyslexia",
    subtitle = "Labels shown for pathways with strength > 1 in both datasets",
    x = "Communication strength (AD_CTL)",
    y = "Communication strength (AD_Dyslexia)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

# Display the plot
pathway_plot

# Save the plot
ggsave(
  filename = file.path(output_dir, "fig_pathway_strength_comparison_red_blue.png"),
  plot = pathway_plot,
  width = 12,
  height = 10,
  dpi = 300,
  bg = "white"
)

cat("✓ Saved: fig_pathway_strength_comparison_red_blue.png\n")






