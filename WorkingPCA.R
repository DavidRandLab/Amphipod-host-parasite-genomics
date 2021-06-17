pca.expnd <- function (object, intgroup = "condition", ntop = 500, returnData = FALSE, PC = c(1,2)) 
{
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, 
                                               drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = ":"))
  }
  else {
    colData(object)[[intgroup]]
  }
  d <- data.frame(PC1 = pca$x[, PC[1]], PC2 = pca$x[, PC[2]], group = group, 
                  intgroup.df, name = colnames(object))
  if (returnData) {
    attr(d, "percentVar") <- percentVar[c(PC[1],PC[2])]
    return(d)
  }
  ggplot(data = d, aes_string(x = paste("PC",PC[1]), y = paste("PC",PC[2]), color = "group")) + 
    geom_point(size = 3) + 
    xlab(paste0("PC", PC[1],": ", round(percentVar[PC[1]] *100), "% variance")) +
    ylab(paste0("PC",PC[2],": ", round(percentVar[PC[2]] *  100), "% variance")) + 
    coord_fixed()
}

pcData<- pca.expnd(vsd,  intgroup = c("colnames.Ogril.cnts.", "infection.status"),returnData = T, PC= c(1,2))


percentVar <- round(100 * attr(pcData, "percentVar"))

ggplot(pcData, aes(PC1, PC2, color = infection.status, label = colnames.Ogril.cnts.)) +
  geom_point(size=6,alpha = .8) +
  geom_text(size=2, nudge_y = 6)+
  xlab(paste0("PC: ",percentVar[1],"% variance")) +
  ylab(paste0("PC: ",percentVar[2],"% variance")) +
  theme_pubr()

