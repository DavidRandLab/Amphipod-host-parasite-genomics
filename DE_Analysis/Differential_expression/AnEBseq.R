BiocManager::install("EBSeq")
library(EBSeq)

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
Ogril.cnts<-as.matrix(Ogril.cnts)

pheno$Test<- sub(".*uninf","Uf",pheno$Test)
pheno$Test<- sub(".*inf","If",pheno$Test)
colnames(pheno)<- c("Cluster","InfectionStatus")
pheno$Condition<- paste(pheno$InfectionStatus,pheno$Cluster,sep="_")


#row.names(Ogril.cnts.m)<- row.names(Ogril.cnts)
#colnames(Ogril.cnts.m)<- colnames(Ogril.cnts)


###### Multi #########
# # Normalize
# Sizes=MedianNorm(Ogril.cnts)
# Ogril.cnts<- as.matrix(Ogril.cnts)
# 
# # setting up the levels for our count table 
# # the way we are going to subset it 
# Cond<- as.factor(paste(sub(".*_.*\\.","",row.names(pheno)),pheno$V1,sep="_"))
# 
# # this gets the patterns given the conditions
# # ie whats different than what 
# PosParti<- GetPatterns(Cond)
# # vizual plot to aid with determining patterns of interest 
# PlotPattern(PosParti)
# 
# MultiSize=MedianNorm(Ogril.cnts)
# MultiOut=EBMultiTest(Ogril.cnts,NgVector=NULL,Conditions=Cond,
#                      AllParti=PosParti[4:6,], sizeFactors=MultiSize,maxround=5)
# 
# 
# # This gets the probability that a gene belongs to a given pattern 
# # of DE 
# MultiPP=GetMultiPP(MultiOut)
# names(MultiPP)
# 
# # This gets the FC values 
# MultiFC=GetMultiFC(MultiOut)
# str(MultiFC)
# 
# QQP(MultiOut)
# MultiPP$MAP[which(MultiPP$MAP=="Pattern5")]
# 


#### Pairwise #####
# DE between clust 1 vs 2
Sizes=MedianNorm(Ogril.cnts)
Cond<- as.factor(pheno$Cluster)
EBOut=EBTest(Data=Ogril.cnts,Conditions=Cond,sizeFactors=Sizes, maxround=5)
ResF=GetDEResults(EBOut, FDR=0.05)
clusterComp<- ResF$DEfound

length(clusterComp)

# DE between Infected vs noninfected
Sizes=MedianNorm(Ogril.cnts)
Cond<- as.factor(pheno$InfectionStatus)
EBOut=EBTest(Data=Ogril.cnts,Conditions=Cond,sizeFactors=Sizes, maxround=5)
ResF=GetDEResults(EBOut, FDR=0.05,Threshold_FC = 1)
infectStatusComp<- ResF$DEfound
length(infectStatusComp)

# plot a venn diagram to determine overlap

A1= length(clusterComp)
A2= length(infectStatusComp)
CA= length(intersect(clusterComp,infectStatusComp))

draw.pairwise.venn(A1, A2, CA, 
                   category = c("Cluster 1\nvs\nCluster 2","Infected\nvs\nnonInfected"),
                   cat.pos = c(0, 0),
                   cat.dist = c(.08,.08),
                   fill=c("blue","yellow"))

write(clusterComp,"clusterComp_EBseq.txt")
write(infectStatusComp,"infectStatusComp_EBseq.txt")
