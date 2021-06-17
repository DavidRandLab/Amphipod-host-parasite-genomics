reactome.amnat <- function(GeneGroup){
  require(drosophila2.db)
  require(ReactomePA)
  de.genes<- na.omit(geneinfo[match(GeneGroup,geneinfo$GENEID),'PrimaryFBgn'])
  de.genes<-AnnotationDbi::mapIds(drosophila2.db, de.genes, 'ENTREZID', 'FLYBASE',multiVals='first')
  de.genes<-unlist(de.genes,use.names = F)
  de.genes<-na.omit(unique(de.genes))
  
  assayed.genes<- na.omit(geneinfo[match(row.names(Ogril.cnts),geneinfo$GENEID),'PrimaryFBgn'])
  assayed.genes<-AnnotationDbi::mapIds(drosophila2.db, assayed.genes, 'ENTREZID', 'FLYBASE',multiVals='first')
  assayed.genes<-unlist(assayed.genes,use.names = F)
  assayed.genes<-na.omit(unique(assayed.genes))
  
  
  x <- enrichPathway(gene=de.genes,pvalueCutoff=0.05, organism = 'fly', universe = assayed.genes, readable=T)
  return(x)
}