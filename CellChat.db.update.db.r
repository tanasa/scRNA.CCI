# Load required libraries
library(NeuronChat)
library(CellChat)
library(dplyr)

# After loading the package
library(NeuronChat)
ls("package:NeuronChat")

# ?computeNetSimilarity_Neuron
# ?compareInteractions_Neuron
# ?computeNetSimilarityPairwise_Neuron
# ?netVisual_circle_compare
# ?netVisual_circle_neuron
# ?netVisual_embedding_Neuron

interactionDB_human

interactionDB_human$VIP_VIPR1

interactionDB_human$VIP_VIPR2

# how to update NeuronChat
# https://htmlpreview.github.io/?https://github.com/Wei-BioMath/NeuronChat/blob/main/vignettes/Update_NeuronChat_database.html

# load example data: cortex data from Yao 2021, processed

data(list='cortex_data')
# subset the data by choosing the region VISp

region_name <- 'VISp'; cell_class <- names(table(meta$class_label))
cell_idx <- which(meta$region_label %in% region_name & meta$class_label %in% cell_class & !(meta$subclass_label %in%c('Car3','CR','DG','L2/3 IT PPP','L5/6 IT TPE-ENT')))
target_df_single  <- target_df[cell_idx,] # a data frame: row  cell, column gene (the last column is the cell subclass)

# creat NeuronChat object; choose the database 'mouse' for mouse data; 'human' for human data
# note that the first parameter should be a numeric matrix with row gene and column cell

x <- createNeuronChat(t(as.matrix(target_df_single[,1:(dim(target_df_single)[2]-1)])),DB='mouse',group.by = target_df_single$cell_subclass)
# Create a NeuronChat object from a data matrix


# Access and update the current NeuronChat interaction database
# database is stored in the `DB` slot of NeuronChat object 
DB <- x@DB
# view the number of entries in current database
length(DB)

head(DB)

# another way to access the database;  replace 'mouse' to 'human' for human data 
data(list='interactionDB_mouse'); 
DB <- eval(parse(text = 'interactionDB_mouse'))
length(DB)
head(DB)

# add one entry of user defined interaction pair to DB via function `update_interactionDB`

# user needs to provide  : 
# interaction_name (character), 
# lig_contributor (a vector of genes), 
# receptor_subunit (a vector of genes), 

# We show this functionality by applying it to the interaction between glutamate  & the AMPA receptor 
# (e.g., subunit composition is Gria1:Gria2:Gria3 = 1:2:1) 

# Contributing genes of glutamate are Gls for synthesis and Slc17a6, Slc17a7, Slc17a8 for vesicular transport
# Receptor subunits of AMPA are Gria1, Gria2, Gria3 

cat("TUTORIAL to UPDATE NeuronChat db")

# code Claude AI to flatten the database

collapse_or_na <- function(x) {
  if (is.null(x) || length(x) == 0) return(NA_character_)
  paste(x, collapse = ";")
}

# Helper to safely extract single values
get_single <- function(x) {
  if (is.null(x) || length(x) == 0) return(NA_character_)
  as.character(x[1])
}

db_df <- do.call(rbind, lapply(names(interactionDB_human), function(nm) {
  x <- interactionDB_human[[nm]]
  data.frame(
    interaction_name       = get_single(x$interaction_name),
    lig_contributor        = collapse_or_na(x$lig_contributor),
    receptor_subunit       = collapse_or_na(x$receptor_subunit),
    lig_contributor_group  = collapse_or_na(x$lig_contributor_group),
    lig_contributor_coeff  = collapse_or_na(x$lig_contributor_coeff),
    receptor_subunit_group = collapse_or_na(x$receptor_subunit_group),
    receptor_subunit_coeff = collapse_or_na(x$receptor_subunit_coeff),
    activator              = collapse_or_na(x$activator),
    inhibitor              = collapse_or_na(x$inhibitor),
    targets_up             = collapse_or_na(x$targets_up),
    targets_down           = collapse_or_na(x$targets_down),
    interaction_type       = get_single(x$interaction_type),
    ligand_type            = get_single(x$ligand_type),
    stringsAsFactors = FALSE
  )
}))

head(db_df)
dim(db_df)

write.csv(
  db_df,
  file = "NeuronChat_interactionDB_flat.csv",
  row.names = FALSE
)

split_genes <- function(x) {
  if (is.null(x) || length(x) == 0) return(character(0))
  x <- x[!is.na(x)]
  if (length(x) == 0) return(character(0))
  unique(trimws(unlist(strsplit(x, ";", fixed = TRUE))))
}

lig_genes <- split_genes(db_df$lig_contributor)
rec_genes <- split_genes(db_df$receptor_subunit)

all_lr_genes <- sort(unique(c(lig_genes, rec_genes)))

cat("all_lr_genes")

length(all_lr_genes)
head(all_lr_genes, 30)
tail(all_lr_genes, 30)

writeLines(all_lr_genes, "NeuronChat_LR_genes_list.txt")

# add one entry of user defined interaction pair to DB via function `update_interactionDB`
# user needs to provide  interaction_name (character), lig_contributor (a vector of genes), receptor_subunit (a vector of genes), 
# We show this functionality by applying it to the interaction between glutamate  & 
# the AMPA receptor (e.g., subunit composition is Gria1:Gria2:Gria3 = 1:2:1) 
# Contributing genes of glutamate are Gls for synthesis and Slc17a6, Slc17a7, Slc17a8 for vesicular transport
# Receptor subunits of AMPA are Gria1, Gria2, Gria3 

interaction_name <- 'Glu_AMPA'
lig_contributor <- c('Gls','Slc17a6','Slc17a7','Slc17a8')
receptor_subunit <- c('Gria1', 'Gria2', 'Gria3')

lig_contributor_group <- c(1,2,2,2) 

# 'Gls' corresponds to group 1; 
# 'Slc17a6','Slc17a7','Slc17a8' correspond to group 2 (redundant for vesicular transport)

lig_contributor_coeff <- c(1,1) # stoichiometry coefficients (or weights) for groups in uniuqe(lig_contributor_group)
receptor_subunit_group <- c(1,2,3) # 1 for Gria1, 2 for Gria2, and 3 for Gria3
receptor_subunit_coeff <- c(1,2,1) # subunit composition is Gria1:Gria2:Gria3 = 1:2:1

DB_new <- update_interactionDB(DB,
                               interaction_name=interaction_name,
                               lig_contributor=lig_contributor,
                               receptor_subunit,
                               interaction_type='user_defined',
                               ligand_type='user_defined',
                               lig_contributor_group=lig_contributor_group,
                               lig_contributor_coeff=lig_contributor_coeff,
                               receptor_subunit_group=receptor_subunit_group,
                               receptor_subunit_coeff=receptor_subunit_coeff)
length(DB_new)
tail(DB_new)

dplyr::glimpse(DB_new$Glu_AMPA) # view the newly added entries

interactionDB_human[["VIP_VIPR1"]]

for (nm in names(interactionDB_human)) {
  cat("\n==============================\n")
  cat("Interaction:", nm, "\n")
  print(interactionDB_human[[nm]])
}

cat("Create NeuronChat object using the updated database")


y <- createNeuronChat(t(as.matrix(target_df_single[,1:(dim(target_df_single)[2]-1)])), 
                      DB='mouse',
                      group.by = target_df_single$cell_subclass)

# Create a NeuronChat object from a data matrix

# dplyr::glimpse(y@DB$Glu_AMPA) # check the newly added entries 

#  List of 14
#>  $ interaction_name      : chr "Glu_AMPA"
#>  $ lig_contributor       : chr [1:4] "Gls" "Slc17a6" "Slc17a7" "Slc17a8"
#>  $ receptor_subunit      : chr [1:3] "Gria1" "Gria2" "Gria3"
#>  $ lig_contributor_group : num [1:4] 1 2 2 2
#>  $ lig_contributor_coeff : num [1:2] 1 1
#>  $ receptor_subunit_group: num [1:3] 1 2 3
#>  $ receptor_subunit_coeff: num [1:3] 1 2 1
#>  $ targets_up            : chr(0) 
#>  $ targets_down          : chr(0) 
#>  $ activator             : chr(0) 
#>  $ inhibitor             : chr(0) 
#>  $ interactors           : chr(0) 
#>  $ interaction_type      : chr "user_defined"
#>  $ ligand_type           : chr "user_defined"

cat("

OpenAI : Big picture: how NeuronChat thinks about an interaction

NeuronChat does not treat a ligand–receptor interaction as “one gene meets one gene”.

Instead, it models an interaction as:

A signal exists only if all required functional components are present,
and each component can itself be made of one or more genes, possibly with weights / stoichiometry.

That’s why you see gene lists + groups + coefficients.

1️⃣ What lig_contributor means (ligand side)

lig_contributor <- c(Gls,Slc17a6,Slc17a7,Slc17a8)


This does not mean “all these genes are ligands”.

It means:

These genes collectively enable the ligand-side signal.

For glutamatergic signaling:

Gls → synthesizes glutamate inside the neuron

Slc17a6 / Slc17a7 / Slc17a8 → package glutamate into synaptic vesicles (VGLUTs)

So biologically:

No synthesis → no glutamate

No vesicular transporter → no release

NeuronChat is trying to encode this logic.

2️⃣ What lig_contributor_group means

lig_contributor_group <- c(1,2,2,2)

This assigns each gene to a functional group.

Interpretation :

Gene	Group	Meaning
Gls	1	Glutamate synthesis
Slc17a6	2	Vesicular transport
Slc17a7	2	Vesicular transport
Slc17a8	2	Vesicular transport

Key idea:
Genes with the same group number are alternatives (OR logic).

So for Group 2:

Slc17a6 or Slc17a7 or Slc17a8 is sufficient
(they are redundant VGLUT isoforms)

Across groups, the logic is AND:

Group 1 AND Group 2 must be present

3️⃣ What lig_contributor_coeff means

lig_contributor_coeff <- c(1,1)

This gives weights / stoichiometry per group, not per gene.

Here:

Group 1 weight = 1

Group 2 weight = 1

Biologically:

You need one unit of synthesis capacity

And one unit of vesicular transport capacity

In practice:

NeuronChat uses these coefficients when aggregating expression to compute a ligand-side signal strength

4️⃣ What receptor_subunit means

receptor_subunit <- c(Gria1,Gria2,Gria3)

This represents a multimeric receptor, not independent receptors.

AMPA receptors are heterotetramers composed of GluA subunits.

NeuronChat models this explicitly instead of pretending there is a single “AMPA gene”.

5️⃣ What receptor_subunit_group means

receptor_subunit_group <- c(1,2,3)

Each receptor subunit is assigned to its own group.

Gene	Group	Meaning
Gria1	1	distinct required subunit
Gria2	2	distinct required subunit
Gria3	3	distinct required subunit

This encodes AND logic across all three:

Gria1 AND

Gria2 AND

Gria3

If one is missing, the receptor complex is incomplete.

6️⃣ What receptor_subunit_coeff means :

receptor_subunit_coeff <- c(1,2,1)

This encodes stoichiometry:

Gria1 : Gria2 : Gria3 = 1 : 2 : 1

Interpretation:

The functional AMPA receptor contains:

one Gria1

two Gria2

one Gria3

NeuronChat uses this when computing receptor availability:

Expression of Gria2 contributes twice as much as Gria1 or Gria3

Limiting subunits constrain the final receptor signal

7️⃣ The logical structure (this is the core idea)

What NeuronChat is really encoding is:

Ligand side
(Gls) AND (Slc17a6 OR Slc17a7 OR Slc17a8)

Receptor side
(1×Gria1) AND (2×Gria2) AND (1×Gria3)

Only when both sides are satisfied does NeuronChat assign a strong Glu→AMPA communication score.

8️⃣ Why this matters (and why NeuronChat is different)

This is much richer than typical CCI tools:

CellChat: usually ligand gene × receptor gene

NeuronChat: biochemical assembly logic

That’s why it’s especially well suited for:

Neurotransmitter systems

Ion channels

Multimeric receptors

Metabolic–signaling coupling

")

cat(

"Summary : 

Adding Custom Interactions

Demonstrates adding a new interaction: Glutamate-AMPA receptor signaling
Defines ligand contributors (genes involved in glutamate synthesis and transport)
Defines receptor subunits (AMPA receptor components: Gria1, Gria2, Gria3)
Specifies stoichiometry: the ratio of receptor subunits (1:2:1)
Uses grouping to indicate that certain genes are redundant (multiple transporters)

")



# Now, let's look at CellChat database

cat("TUTORIAL on updating CellChat database ")



# CellChatDB database

CellChatDB <- CellChatDB.human

CellChatDB

str(CellChatDB)

head(CellChatDB$interaction)

head(CellChatDB$complex)

head(CellChatDB$cofactor)

head(CellChatDB$geneInfo)

# help(updateCellChatDB)

# https://github.com/qingchen36/ligand-receptor
# other databases transformed into CellChat format

# for example to use : https://github.com/ZJUFanLab/CellTalkDB

colnames(CellChatDB.human$interaction)
dim(CellChatDB.human$interaction)

getwd()

# The first input db is a data frame with as least two columns named as ligand and receptor. 

# We highly suggest users to provide a column of pathway information named ‘pathway_name’ associated with each L-R pair.
# Other optional columns include ‘interaction_name’ and ‘interaction_name_2’. 
# The default columns of CellChatDB can be checked via ‘colnames(CellChatDB.human$interaction)’.

# The second input gene_info is a data frame with at least one column named as ‘Symbol’. 

# “When setting gene_info = NULL, the input ‘species_target’ should be provided: either ‘human’ or ‘mouse’.

# The third input other_info is a list consisting of other information including a dataframe named as ‘complex’ 
# and a dataframe named as ‘cofactor’. 

# This additional information is not necessary. 
# If other_info is provided, the ‘complex’ and ‘cofactor’ are dataframes with defined rownames.

cat("Step 1: Load the customized ligand-receptor pairs and the gene information")

# Load the dataframe consisting of customized ligand-receptor pairs
# db.user <- readRDS("./CellTalkDB_rds/CellTalkDB/human_lr_pair.rds")
# Load the dataframe consisting of gene information (optional)
# gene_info <- readRDS("./CellTalkDB_rds/CellTalkDB/human_gene_info.rds")

cat("Step 2: Formulate the input files to be compatible with CellChatDB")

# getwd()
# list.files()

# Load the dataframe consisting of customized ligand-receptor pairs
db.user <- readRDS("/home/tanasa/Desktop/NeuronChat/CellTalkDB/database/human_lr_pair.rds")
# Load the dataframe consisting of gene information (optional)
gene_info <- readRDS("/home/tanasa/Desktop/NeuronChat/CellTalkDB/data/human_gene_info.rds")

head(db.user)

head(gene_info)

# Users need to make sure that user-defined complex/cofactor names are the same in CellChatDB$interaction.

# Modify the colnames because of incompatible with `colnames(CellChatDB.human$interaction)`
colnames(db.user) <- plyr::mapvalues(colnames(db.user), 
                                     from = c("ligand_gene_symbol","receptor_gene_symbol","lr_pair"), 
                                     to = c("ligand","receptor","interaction_name"), 
                                     warn_missing = TRUE)

cat("Step 3: Update CellChatDB by running the function updateCellChatDB")

# Create a new database by typing one of the following commands
# Use user-provided gene information

db.new <- updateCellChatDB(db = db.user, 
                           gene_info = gene_info)
# Use built-in gene information of CellChatDB

db.new <- updateCellChatDB(db = db.user, 
                           gene_info = NULL, 
                           species_target = "human")

# Alternatively, users can integrate the customized L-R pairs into the built-in CellChatDB 
db.new <- updateCellChatDB(db = db.user, 
                           merged = TRUE, 
                           species_target = "human")

head(db.new$interaction)
dim(db.new$interaction)

# Warning: Because no pathway information of each L-R pair is provided in db.user here, 
# all pathway-level analysis from CellChat cannot be used, such as computeCommunProb! 

cat("Step 4: Use the new database in CellChat analysis or re-build CellChat package")
# Re-build CellChat package by updating the database as follows

# Users can now use this new database in CellChat analysis 
# cellchat@DB <- db.new

# Users can save the new database for future use
# save(db.new, file = "CellChatDB.human_user.rda")

# Users can also re-build CellChatDB in CellChat package

# If working on a human dataset, do following:
# CellChatDB.human <- CellChatDB
# usethis::use_data(CellChatDB.human, overwrite = TRUE)

cat("Using L-R interactions from CellPhoneDB")

# Since cellphonedb v5, cellphonedb also introduces signalling directionality (ligand is partner A, receptor partner B) 
# and classification of signaling pathways, which makes it easiler to make full use of CellChat’s versatile analysis 
# by taking CellPhoneDB as an input. 

# cellphonedb
# https://github.com/ventolab/cellphonedb-data

# Subunit architecture is included for both ligands and receptors, representing heteromeric complexes accurately. 
# This is crucial, as cell-cell communication relies on multi-subunit protein complexes that go beyond 
# the binary representation used in most databases and studies.

# Includes interactions involving non-proteins (i.e., molecules not encoded by a gene/trasncript) acting as ligands. 
# Examples of these include steroid hormones (e.g., estrogen), small molecules (e.g. histamine) or neurotransmitters. 

# To do so, we have reconstructed the biosynthetic pathways and used the last representative enzyme as a proxy of ligand abundance. 
# We retrieve this information by manually reviewing and curating relevant literature and peer-reviewed pathway resources such as REACTOME. 
# Cellphonedb-data includes more than 200 interactions involving non-peptidic ligands! 
# All non-protein ligands are coded as complexes in CellphoneDB. 
# This is to allow the inclusion of all the key enzymes in their biosynthetic pathway (even if the number of enzymes n=1)

library(CellChat)
options(stringsAsFactors = FALSE)

## Step 1: Load the database files
## load the database files

interaction_input <- read.csv(file = './cellphonedb-data/data/interaction_input.csv')
complex_input <- read.csv(file = './cellphonedb-data/data/complex_input.csv', row.names = 1)
geneInfo <- read.csv(file = './cellphonedb-data/data/gene_input.csv')
geneInfo$Symbol <- geneInfo$hgnc_symbol
geneInfo <- select(geneInfo, -c("ensembl"))
geneInfo <- unique(geneInfo)

head(interaction_input)

head(complex_input)

cat("Step 2: Formulate the input files to be compatible with CellChatDB")

## prepare interaction_input file
## get the ligand information

idx_partnerA <- match(interaction_input$partner_a, geneInfo$uniprot)
idx.use <- !is.na(idx_partnerA)
interaction_input$ligand <- interaction_input$partner_a
interaction_input$ligand[idx.use] <- geneInfo$hgnc_symbol[idx_partnerA[idx.use]]

# get the receptor information

idx_partnerB <- match(interaction_input$partner_b, geneInfo$uniprot)
idx.use <- !is.na(idx_partnerB)
interaction_input$receptor <- interaction_input$partner_b
interaction_input$receptor[idx.use] <- geneInfo$hgnc_symbol[idx_partnerB[idx.use]]

# get other information
interaction_input$interaction_name <- interaction_input$interactors
interaction_input$interaction_name_2 <- interaction_input$interaction_name
interaction_input$pathway_name <- interaction_input$classification
interaction_input$pathway_name <- gsub(".*by ", "", interaction_input$pathway_name)
interaction_input$annotation <- interaction_input$directionality
interaction_input <- select(interaction_input, -c("partner_a","partner_b","protein_name_a","protein_name_b","interactors","classification","directionality"))

## prepare complex_input file
complexsubunits <- dplyr::select(complex_input, starts_with("uniprot"))
for (i in 1:ncol(complexsubunits)) {
  idx_complex <- match(complex_input[,paste0("uniprot_",i)], geneInfo$uniprot)
  idx.use <- !is.na(idx_complex)
  complex_input[idx.use,paste0("uniprot_",i)] <- geneInfo$hgnc_symbol[idx_complex[idx.use]]
}
colnames(complex_input)[1:ncol(complexsubunits)] <- paste0("subunit_",seq_len(ncol(complexsubunits)))
complex_input <- dplyr::select(complex_input, starts_with("subunit"))

## prepare other information 
other_info <- list(complex = complex_input)

head(interaction_input) 

cat("Step 3: Update CellChatDB by running the function updateCellChatDB")

db.new <- updateCellChatDB(db = interaction_input, gene_info = geneInfo, other_info = other_info, trim.pathway = T)

cat("Step 4: Use the new database in CellChat analysis or re-build CellChat package")

# Users can now use this new database in CellChat analysis 
# cellchat@DB <- db.new

# Users can save the new database for future use
# save(db.new, file = "CellChatDB.human_user.rda")

# Users can also re-build CellChatDB in CellChat package

# setwd("/Users/$USERS/Downloads/CellChat-master") # This is the folder of CellChat package downloaded from Github
# CellChatDB.human <- CellChatDB
# usethis::use_data(CellChatDB.human, overwrite = TRUE)



cat("Part II: update CellChatDB by manually modifying the required files")

cat("Step 1: Access the ligand-receptor interaction information in CellChatDB")

# Extract the database information in CellChatDB and then save them in a local computer, 
# including four files: 
# ‘geneInfo.csv’, 
# ‘interaction_input_CellChat.csv’, 
# ‘complex_input_CellChat.csv’, 
# ‘and cofactor_input_CellChat.csv’.

CellChatDB <- CellChatDB.mouse # set CellChatDB <- CellChatDB.human if working on the human dataset

interaction_input <- CellChatDB$interaction
complex_input <- CellChatDB$complex
cofactor_input <- CellChatDB$cofactor
geneInfo <- CellChatDB$geneInfo

# write.csv(interaction_input, file = "interaction_input_CellChatDB.csv")
# write.csv(complex_input, file = "complex_input_CellChatDB.csv")
# write.csv(cofactor_input, file = "cofactor_input_CellChatDB.csv")
# write.csv(geneInfo, file = "geneInfo_input_CellChatDB.csv")

cat("Step 2 :Update the required files by adding users’ curated ligand-receptor pairs")

# Update the above four .csv files by adding users’ curated ligand-receptor pairs.

# The main file is ‘interaction_input_CellChatDB.csv’. Users can first update the ligands, receptors and co-factors 
# in the corresponding columns in ‘interaction_input_CellChatDB.csv’.

# Users can then update “complex_input_CellChatDB.csv” and “cofactor_input_CellChatDB.csv” if any ligand complex, receptor complex and cofactors are updated. 
# Users need to make sure that user-defined complex/cofactor names are the same in ‘interaction_input_CellChatDB.csv’ 
# and “complex_input_CellChatDB.csv”, ‘interaction_input_CellChatDB.csv’ and ” cofactor_input_CellChatDB.csv”.

# “geneInfo_input_CellChatDB.csv” contains all gene information in mouse and it should have a column named ‘Symbol’, which does not need to be changed when updating CellChatDB.


cat("Step 3: Update CellChatDB")

# options(stringsAsFactors = FALSE)

# interaction_input <- read.csv(file = 'interaction_input_CellChatDB.csv', row.names = 1)
# complex_input <- read.csv(file = 'complex_input_CellChatDB.csv', row.names = 1)
# cofactor_input <- read.csv(file = 'cofactor_input_CellChatDB.csv', row.names = 1)
# geneInfo <- read.csv(file = ' geneInfo_input_CellChatDB.csv', row.names = 1)

# CellChatDB <- list()
# CellChatDB$interaction <- interaction_input
# CellChatDB$complex <- complex_input
# CellChatDB$cofactor <- cofactor_input
# CellChatDB$geneInfo <- geneInfo

cat("Step 4: Re-build CellChat package (optional)")

# CellChatDB.mouse <- CellChatDB
# usethis::use_data(CellChatDB.mouse, overwrite = TRUE)

# If working on a human dataset, do following:
# CellChatDB.human <- CellChatDB
# usethis::use_data(CellChatDB.human, overwrite = TRUE)


