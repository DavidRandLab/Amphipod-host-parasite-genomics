library(edgeR)  
library(EBSeq)
library(DESeq2)
library(goseq)
library(seqsetvis)
library(tidyverse)
library(EnhancedVolcano)
library(ggpubr)
library(googledrive)
library(readxl)


options(
  gargle_oauth_cache = "~/.secrets",
  gargle_oauth_email = TRUE
)

source('Functions/getFly.R')
source('Functions/GOamnat.R')
source('~/Documents/R-Projects/R_GlobalFunction/Bioconductor_Functions.R')
source('~/Documents/R-Projects/Amphipod-Analysis/Functions/reactome.amnat.R')
source('~/Documents/R-Projects/Amphipod-Analysis/Functions/GOamnat.R')
source('~/Documents/R-Projects/Amphipod-Analysis/Functions/getFly.R')
source('~/Documents/R-Projects/R_GlobalFunction/GoogleDrive_Functions.R')

rmarkdown::render("DESeq2/DESeq2.Rmd")
rmarkdown::render("EBSeq/EBSeq.Rmd") 
rmarkdown::render("edgeR/edgeR.Rmd")


geneinfo<- read.delim('/Volumes/GoogleDrive/My Drive/2019 Amphipod Analyses/DEanalysis/Current Analysis/Counts_Annotations/Ogril1_RNAevidence_MASTER_ANNOTATION_FILE.txt', stringsAsFactors=FALSE)


#### PCA ######
# 
# vsd<- vst(dds,blind = T)
# 
# pcData<- plotPCA(vsd, intgroup = c("colnames.Ogril.cnts.", "infection.status"),returnData = T)
# percentVar <- round(100 * attr(pcData, "percentVar"))
# 
# ggplot(pcData, aes(PC1, PC2, color = infection.status, label = colnames.Ogril.cnts.)) +
#   geom_point(size=6,alpha = .8) +
#   geom_text(size=2, nudge_y = 6)+
#   xlab(paste0("PC1: ",percentVar[1],"% variance")) +
#   ylab(paste0("PC2: ",percentVar[2],"% variance")) +
#   theme_pubr()
# 
# ggsave('PCA.png', device = 'png', dpi = 300)
# upload_file('PCA.png',
#             path = as_id('https://drive.google.com/open?id=1uf3B1DOlM5h3Jsxk_rBRM2eC1reAxzIv'))

##### Infection Status ####

# read in saved tables 
drive_download(as_id('1gvrQ2u_p04_E3KmV0aXImJRiENlmX3IkNVZjju--Uok'), overwrite = T) # download count table
DEtable_InfectionStatus_Genome_DESeq <- read_excel('DEtable_InfectionStatus_Genome_DESeq.xlsx')
file.remove('DEtable_InfectionStatus_Genome_DESeq.xlsx')


drive_download(as_id('1N07lSDQIeD9v78lV5cXeW-ev27-WhmJLWiQSPsKafdc'), overwrite = T) # download count table
DEtable_InfectionStatus_Genome_edgeR <- read_excel('DEtable_InfectionStatus_Genome_edgeR.xlsx')
file.remove('DEtable_InfectionStatus_Genome_edgeR.xlsx')

drive_download(as_id('107GctPvQYTx5kFVxtqwA6DxXiQqnzr41SPASsPDW27c'), overwrite = T) # download count table
DEtable_InfectionStatus_Genome_EBSeq <- read_excel('DEtable_InfectionStatus_Genome_EBSeq.xlsx')
file.remove('DEtable_InfectionStatus_Genome_EBSeq.xlsx')

overlap.list.InfectionStatus = list(pull(DEtable_InfectionStatus_Genome_DESeq,'PredictedGene'),
                                    pull(DEtable_InfectionStatus_Genome_edgeR,'PredictedGene'),
                                    pull(DEtable_InfectionStatus_Genome_EBSeq,'PredictedGene'))


names(overlap.list.InfectionStatus) = c('DESeq2','edgeR','EBSeq')
overlap.InfectionStatus <- ssvMakeMembTable(overlap.list.InfectionStatus)


ssvFeatureVenn(overlap.InfectionStatus,circle_colors =c("red3",'goldenrod','darkgreen'),
               fill_alpha=.4,
               line_alpha=0)

ggsave('InfectionStatus_MethodsOverlap.png', device = 'png', dpi = 300)
upload_file('InfectionStatus_MethodsOverlap.png',
            path = as_id('https://drive.google.com/open?id=1uf3B1DOlM5h3Jsxk_rBRM2eC1reAxzIv'))




