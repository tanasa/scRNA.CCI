# use the code from: 
# https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/Update-CellChatDB.html

library(CellChat)
options(stringsAsFactors = FALSE)

getwd()

# CellChatDB.human

dim(CellChatDB.human$interaction)

write.csv(CellChatDB.human$interaction, file = "CellChat_Interactions_before_merging_with_CellPhoneDB.csv", row.names = FALSE)

list.files()



cat("Step 1: Load the database files")

interaction_input <- read.csv(file = './cellphonedb-data/data/interaction_input.csv')
complex_input <- read.csv(file = './cellphonedb-data/data/complex_input.csv', row.names = 1)
geneInfo <- read.csv(file = './cellphonedb-data/data/gene_input.csv')
geneInfo$Symbol <- geneInfo$hgnc_symbol
geneInfo <- select(geneInfo, -c("ensembl"))
geneInfo <- unique(geneInfo)

CellChatDB <- CellChatDB.human 

geneInfo$Symbol <- geneInfo$hgnc_symbol
# geneInfo <- select(geneInfo, -c("ensembl"))
geneInfo <- unique(geneInfo)

head(geneInfo)

cat("Step 2: Formulate the input files to be compatible with CellChatDB")

## prepare interaction_input file

# get the ligand information
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

cat("Step 3: Update CellChatDB by running the function updateCellChatDB")

db.new <- updateCellChatDB(db = interaction_input, 
                           gene_info = geneInfo, 
                           other_info = other_info, 
                           trim.pathway = T)

dim(db.new$interaction)

head(db.new$interaction)

# ?updateCellChatDB
# Update CellChatDB by integrating new L-R pairs from other resources or adding more information 

write.csv(db.new$interaction, 
          file = "CellChat_Interactions_is_it_from_CellPhoneDB.csv", 
          row.names = FALSE)

getwd()

save(db.new, file = "CellPhoneDB.to.Cell.Chat.human.rda")


