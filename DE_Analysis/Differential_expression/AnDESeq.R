BiocManager::install("DESeq2")
library(DESeq2)

# read in JCBN count table
# generated from genome 
Ogril.cnts<- read.csv("https://www.dropbox.com/s/z3lxxn4clc4q6lm/Ogril_Read_count.txt?dl=1",sep="\t")
pdata<- read.csv("https://www.dropbox.com/s/sfaxuk2xw019loz/Amphipod_SNP_cluster_codes.csv?dl=1")
#head(pdata)
#head(Ogril.cnts)

# rename columns to match cluster file names 
colnames(Ogril.cnts)<- sub("Og","",colnames(Ogril.cnts))

#Create a data frame where row names match the col names in the count table
# order is important 
pheno<- data.frame(row.names = colnames(Ogril.cnts))
for(SampleID in row.names(pheno)){
  indx<- which(pdata$ind==SampleID)
  pheno[SampleID,1]=pdata[indx,2]
}

# Double check everything is indexed in the correct order 
pheno$Test<- colnames(Ogril.cnts)

####
# filter rows with counts lower than 2/million
rowSumskeep <- rowSums(cpm(Ogril.cnts)>1) >=24
Ogril.cnts<- Ogril.cnts[rowSumskeep,]
table(rowSumskeep)
cnts<- as.matrix(Ogril.cnts)


pheno$Test<- sub(".*uninf","Uf",pheno$Test)
pheno$Test<- sub(".*inf","If",pheno$Test)
colnames(pheno)<- c("Cluster","InfectionStatus")
pheno$Condition<- paste(pheno$InfectionStatus,pheno$Cluster,sep="_")
coldata<- pheno

coldata$Cluster<- factor(coldata$Cluster)
coldata$InfectionStatus<- factor(coldata$InfectionStatus)


#des<- model.matrix(~ 0+InfectionStatus,coldata)
dds <- DESeqDataSetFromMatrix(countData = cnts,
                              colData = coldata,
                              design = ~ InfectionStatus + Cluster )
## If status
dds<- DESeq(dds)
resIf<- results(dds,alpha=0.05,contrast = c("InfectionStatus","If","Uf"))
summary(resIf)

DE.If<- row.names(resIf[which(resIf$padj<0.05),]) 
ncnts<- assay(vsd)
row.names(ncnts)
DE.If

ncnts[which(row.names(ncnts) %in% DE.If),]
my_palette <- colorRampPalette(c("red", "white", "blue"))(n = 299)
heatmap.2(ncnts[which(row.names(ncnts) %in% DE.If),order(coldata$Condition)],col=my_palette,
          density.info="none",  # turns off density plot inside color legend
          trace="none",
          key = F,
          Colv = NA)

plotMA(results)
### Cluster
dds<- DESeq(dds)
results(dds)
resClus<- results(dds,alpha=0.05,lfcThreshold = 1,contrast = c("Cluster","1","2"))
summary(resClus)

DE.Clus<- row.names(resClus[which(resClus$padj<0.05),]) 
vsd<- vst(dds,blind = F)
ncnts<- assay(vsd)
row.names(ncnts)
DE.Clus

# plot a venn diagram to determine overlap

A1= length(DE.Clus)
A2= length(DE.If)
CA= length(intersect(DE.Clus,DE.If))

draw.pairwise.venn(A1, A2, CA, 
                   category = c("Cluster 1\nvs\nCluster 2","Infected\nvs\nnonInfected"),
                   cat.pos = c(0, 0),
                   cat.dist = c(.08,.08),
                   fill=c("blue","yellow"))

write(DE.Clus,"clusterComp_DESeq.txt")
write(DE.If,"infectStatusComp_DEeq.txt")

