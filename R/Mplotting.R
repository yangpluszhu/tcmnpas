Mplotting = function(clProf.reshape.df,x = ~Cluster,type = "dot",colorBy = "p.adjust",by = "geneRatio",title="",font.size=12){ 
  Description <- Percentage <- Count <- Cluster <- GeneRatio <- p.adjust <- pvalue <- NULL # to satisfy codetools 
  if (type == "bar") { 
    if (by == "percentage") { 
      p <- ggplot(clProf.reshape.df, 
                  aes(x=Description, y = Percentage, fill=Cluster)) 
    } else if (by == "count") { 
      p <- ggplot(clProf.reshape.df, 
                  aes(x=Description, y = Count, fill=Cluster)) 
    } else { 
 
    } 
    p <- p + 
      geom_bar() + 
      coord_flip() 
  } 
  if (type == "dot") { 
    if (by == "rowPercentage") { 
      p <- ggplot(clProf.reshape.df, 
                  aes_(x = x, y = ~Description, size = ~Percentage)) 
    } else if (by == "count") { 
      p <- ggplot(clProf.reshape.df, 
                  aes_(x = x, y = ~Description, size = ~Count)) 
    } else if (by == "geneRatio") { 
      p <- ggplot(clProf.reshape.df, 
                  aes_(x = x, y = ~Description, size = ~GeneRatio)) 
    } else { 
      ## nothing here 
    } 
    if (any(colnames(clProf.reshape.df) == colorBy)) { 
      p <- p + 
        geom_point() + 
        aes_string(color=colorBy) + 
        scale_colour_gradient(low="red", high="blue") 
    } else { 
      p <- p + geom_point(colour="steelblue") 
    } 
  } 
  p <- p + xlab("") + ylab("") + ggtitle(title) + 
    theme_dose(font.size) 
  return(p) 
} 
 
