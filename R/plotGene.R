library(reshape2)
plotGene <- function(dat,
                     gene) {
  subsample <- dat[,genes]
  if (length(gene) > 1) {
    subsample <- melt(subsample)
  }
  plot <- ggplot(subsample, aes(x = subtype, y = subsample[,gene], fill = subtype)) + 
    geom_violin() + 
    geom_boxplot(width = 0.25, fill = "white") + 
    theme_classic() + 
    theme(axis.title.x = element_blank())
  if(length(gene) > 1) {
    plot <- plot + facet_wrap(~subsample[,gene]) + 
      ylab("log2 mRNA Level")
  } else {
    plot + ylab(subsample[,gene])
  }
  
  return(plot)
}