# Creating Consolidated Lists 

ConsolidateDE_cluster<- Reduce(intersect,list(DECluster_edgeR,DE.Clus,clusterComp))
ConsolidatedDE_infectionStatus<- Reduce(intersect,list(DEInfectionStatus_edgeR,DE.If,infectStatusComp))

A1=length(infectStatusComp) # EBseq
A2=length(DE.If)# DESeq
A3=length(DEInfectionStatus_edgeR) #edgeR
n.12<- length(intersect(infectStatusComp,DE.If))
n.23<- length(intersect(DE.If,DEInfectionStatus_edgeR))
n.13<- length(intersect(infectStatusComp,DEInfectionStatus_edgeR))
n.123<- length(ConsolidatedDE_infectionStatus)


draw.triple.venn(A1, A2, A3,n12 =n.12 ,n23 =n.23 ,n13 =n.13 ,n123 = n.123,
                   category = c("EBSeq","DESeq2","edgeR"),
                   cat.pos = c(0, 0,180),
                   cat.dist = c(.08,.08,0.08),
                   fill=c("blue","yellow","red"))

###### heatmap ######

ncnts[which(row.names(ncnts) %in% ConsolidatedDE_infectionStatus),]
my_palette <- colorRampPalette(c("red", "white", "blue"))(n = 299)
heatmap.2(ncnts[which(row.names(ncnts) %in% ConsolidatedDE_infectionStatus),order(coldata$Condition)],col=my_palette,
          density.info="none",  # turns off density plot inside color legend
          trace="none",
          key = F,
          Colv = NA)

##### Cluster ######
ncnts[which(row.names(ncnts) %in% ConsolidateDE_cluster),]
my_palette <- colorRampPalette(c("red", "white", "blue"))(n = 299)

heatmap.2(ncnts[which(row.names(ncnts) %in% ConsolidateDE_cluster),order(coldata$Cluster)],col=my_palette,
          density.info="none",  # turns off density plot inside color legend
          trace="none",
          key = F,
          Colv = NA)




##### PCA Plots
# Variance stabilizing 

vsd<- vst(dds,blind = T)
pcData<- plotPCA(vsd, intgroup=c("InfectionStatus", "Cluster"),returnData=T)
percentVar <- round(100 * attr(pcData, "percentVar"))
ggplot(pcData, aes(PC1, PC2, color=InfectionStatus, shape=Cluster)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()+
  scale_color_manual(values=c("navyblue","darkgreen"))+
  theme_bw(base_size = 15)


##### GO #####
# Do GO analysis on high confidence hits using JCBN blast scores
# First separate consolidated DE into up and down regulated 

consolidated.table<- results_fitQL$table[ConsolidatedDE_infectionStatus,]

upGenesindx<- which(consolidated.table$logFC >=1)
downGenesindx<- which(consolidated.table$logFC <= -1)

upGenes<- row.names(consolidated.table[upGenesindx,])
downGenes<- row.names(consolidated.table[downGenesindx,])
Universe<- row.names(cnts)

library(drosophila2.db)
blast<- read.csv("https://www.dropbox.com/s/e7dd1oh6i3o6pyk/Ogril1_DMEL_BLAST.txt?dl=1",sep="\t")
blast_FB<- read.csv("https://www.dropbox.com/s/dcx5dlzumd8opaj/DMEL-FBgn-uniprot.txt?dl=1",sep="\t")

blast_FB$Entry
blast$GENEID<- sub(".t.*","",blast$GENEID)

# Universe Set Up 
bindx<- which(blast$GENEID %in% Universe)
blastU<- blast[bindx,1:2]

# DEup 
bindx<- which(blast$GENEID %in% upGenes)
blastup<- blast[bindx,1:2]

# DEdown
bindx<- which(blast$GENEID %in% downGenes)
blastdown<- blast[bindx,1:2]

### Some non overlap, need to check 

# Convert JCBN files into a usable form 
blastU$FB=NA
for (indx in blast_FB[,"Entry"]){
  if (as.character(indx) %in% as.character(blastU$Entry)){
    print("Tr")
  blastU$FB<- blast_FB[indx,"Cross.reference_FlyBase"]
  }
}


