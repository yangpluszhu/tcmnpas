fortify_enrichResult = function(model, data, showCategory=5, by = "Count", order=FALSE, drop=FALSE, ...){
  res <- as.data.frame(model)
  if (drop) {
    res <- res[res$Count != 0, ]
  }
  ####
  gsize <- as.numeric(sub("/\\d+$", "", as.character(res$GeneRatio)))
  gcsize <- as.numeric(sub("^\\d+/", "", as.character(res$GeneRatio)))
  res$GeneRatio = gsize/gcsize
  ####
  #res$GeneRatio <- parse_ratio(res$GeneRatio)

  if (order) {
    if (by == "Count") {
      idx <- order(res$Count, decreasing=TRUE)
    } else {
      idx <- order(res$GeneRatio, decreasing=TRUE)
    }
    res <- res[idx,]
  }

  if ( is.numeric(showCategory) ) {
    if ( showCategory <= nrow(res) ) {
      res <- res[1:showCategory,]
    }
  } else { ## selected categories
    res <- res[res$ID %in% showCategory,]
  }

  res$Description <- factor(res$Description,
                            levels=rev(res$Description))

  return(res)
}

