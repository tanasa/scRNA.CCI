library(CellChat)
library(patchwork)
library(dplyr)



# Define working directory where the files are located
work_dir <- "/mnt/nfs/CX000008_DS1/projects/jaeyeon/fastq_file_Dyslexia_r1/zanalysis_bogdan/samples.merged_AG_Harmony_res0.1_anno_091625/merged_AG_Harmony_res0.1_anno_091625.to.compare.results_cellchat"
setwd(work_dir)
list.files()

ad_ctl = readRDS("cellchat_merged_AG_Harmony_res0.1_anno_091625.AD_CTL_at_the_end.rds")

# glimpse(ad_ctl)
unique(ad_ctl@meta$orig.ident)
unique(ad_ctl@meta$group)
unique(ad_ctl@meta$celltype)

ad_dys = readRDS("cellchat_merged_AG_Harmony_res0.1_anno_091625.AD_Dyslexia_at_the_end.rds")

# glimpse(ad_dys)
unique(ad_dys@meta$orig.ident)
unique(ad_dys@meta$group)
unique(ad_dys@meta$celltype)

ctl_ctl = readRDS("cellchat_merged_AG_Harmony_res0.1_anno_091625.Ctr_Ctr_at_the_end.rds")

glimpse(ctl_ctl)

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

glimpse(ctl_ctl@net)

ctl_ctl@net

# glimpse(ctl_ctl@netP)
# See what type each component is
lapply(ad_ctl@netP, class)

str(ctl_ctl@netP, max.level = 1)

ctl_ctl@netP$pathways

ctl_ctl@netP$prob

# ctl_ctl@netP$centr
ctl_ctl@netP$centr[["LAMININ"]]

ctl_ctl@netP$pattern

str(ctl_ctl@netP$pattern, max.level = 1)

ctl_ctl@netP$pattern$incoming

ctl_ctl@netP$pattern$outgoing



ctl_ctl@netP$similarity

ctl_ctl@netP$similarity$functional

ctl_ctl@netP$similarity$structural

colnames(ctl_ctl@meta)

head(ctl_ctl@meta)

dim(ctl_ctl@meta)

unique(ctl_ctl@meta$orig.ident)

unique(ctl_ctl@meta$celltype)

unique(ctl_ctl@meta$group)

unique(ctl_ctl@meta$group)
unique(ctl_ctl@meta$celltype)

# unique(ctl_ctl@idents)

# unique(ctl_ctl@meta$orig.ident)
# unique(ctl_ctl@meta$group)
# unique(ctl_ctl@meta$celltype)

head(ctl_ctl@LR$LRsig, 2)

# head(ctl_ctl@DB, 2)
head(ctl_ctl@DB$interaction, 2)
head(ctl_ctl@DB$complex, 2)

# (ctl_ctl@var.features)
head(ctl_ctl@var.features$features, 5)
head(ctl_ctl@var.features$features.info, 5)



# To continue with : Create comparison list
# cellchat_comparison <- list(Control = cellchat_azd_ctr_, Dyslexia = cellchat_azd_dys)
# cellchat_merged <- mergeCellChat(cellchat_comparison, add.names = names(cellchat_comparison))


