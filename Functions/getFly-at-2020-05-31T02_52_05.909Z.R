
getFLY <- function(geneid){
  require(drosophila2.db)
  genes <- geneinfo[match(geneid,geneinfo$GENEID),'PrimaryFBgn']
    genes <-AnnotationDbi::select(drosophila2.db, genes, 'SYMBOL', 'FLYBASE')
    return (genes$SYMBOL)
}
