GOamnat<- function(GeneGroup){
  require(drosophila2.db)
  de.genes<- na.omit(geneinfo[match(GeneGroup,geneinfo$GENEID),'PrimaryFBgn'])
  de.genes<-AnnotationDbi::mapIds(drosophila2.db, de.genes, 'SYMBOL', 'FLYBASE',multiVals='first')
  de.genes<-unlist(de.genes,use.names = F)
  de.genes<-na.omit(unique(de.genes))
  
  assayed.genes<- na.omit(geneinfo[match(row.names(Ogril.cnts),geneinfo$GENEID),'PrimaryFBgn'])
  assayed.genes<-AnnotationDbi::mapIds(drosophila2.db, assayed.genes, 'SYMBOL', 'FLYBASE',multiVals='first')
  assayed.genes<-unlist(assayed.genes,use.names = F)
  assayed.genes<-na.omit(unique(assayed.genes))
  
  gene.vector=as.integer(assayed.genes %in% de.genes)
  names(gene.vector)=assayed.genes
  
  
  pwf=nullp(gene.vector,'dm3','geneSymbol',plot.fit = F)
  GO.wall=goseq(pwf,'dm3','geneSymbol',test.cats=c("GO:BP"))
  GO.wall= GO.wall[which(GO.wall$numDEInCat > 3),] # changed to enforce more than 3 genes constituting a category 
  
  GO.wall$padj = p.adjust(GO.wall$over_represented_pvalue, method="BH")
  enriched.GO=GO.wall$category[p.adjust(GO.wall$over_represented_pvalue, method="BH") <0.05]
  GO.sum<- GO.wall[which(GO.wall$category %in% enriched.GO),]
  
  
  out.genes<- c()
  go_id<- enriched.GO
  for (go.term in go_id){
    allegs<-get(go.term, org.Dm.egGO2ALLEGS)
    genes = unlist(mget(allegs,org.Dm.egSYMBOL),use.names = FALSE)
    
    term.cg<- de.genes[which(de.genes %in% genes)]
    out.genes <- c(out.genes, toString(term.cg))
  }
  GO.sum$Genes <- out.genes
  return(GO.sum)
  
  
}

