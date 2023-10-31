Mbarplot = function(height, x="Count", colorBy='pvalue', showCategory=5, font.size=12, title="", ...){ 
  object <- height 
  colorBy <- match.arg(colorBy, c("pvalue", "p.adjust", "qvalue")) 
  if (x == "geneRatio" || x == "GeneRatio") { 
    x <- "GeneRatio" 
  } 
  else if (x == "count" || x == "Count") { 
    x <- "Count" 
  } 
  Description <- Count <- NULL # to satisfy codetools 
  df <- fortify_enrichResult(object, showCategory=showCategory, by=x, ...) 
 
  p <- ggplot(df, aes_string(x = "Description", y = x)) 
  p <- p + geom_bar(stat = "identity") + coord_flip() + theme_dose(font.size) 
 
  if("pvalue" %in% colnames(p$data)) { 
    pvalue <- NULL # to satisfy codetools 
    p <- p + aes_string(fill=colorBy) + 
      scale_fill_continuous(low="red", high="blue") 
  } else { 
    p <- p+aes(fill=Description) + 
      theme(legend.position="none") 
  } 
  p <- p + ggtitle(title) + xlab("") + ylab("") 
  return(p) 
} 
 
