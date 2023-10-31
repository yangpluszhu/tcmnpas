ToDataFrame = function(x, itemSep = ",", setStart = "{", setEnd = "}", linebreak = NULL,...) { 
  require(arules) 
  n_of_itemsets <- length(x) 
  if(n_of_itemsets == 0) return(invisible(NULL)) 
  ## Nothing to inspect here ... 
 
  l <- labels(x, itemSep, setStart, setEnd) 
  if(is.null(linebreak)) 
    linebreak <- any(nchar(l) > options("width")$width*2/3) 
 
  if(!linebreak) { 
    out <- data.frame(items = l) 
    if(nrow(quality(x)) > 0) out <- cbind(out, quality(x)) 
    #print(out, right = FALSE) 
 
  } else { 
 
 
    ## number of rows + fix empty itemsets 
    items <- unlist(lapply(as(items(x), "list"), 
                           FUN = function(x) if(length(x) == 0) "" else x)) 
    n_of_items <- size(items(x)) 
    n_of_items[n_of_items == 0] <- 1 
 
    ## calculate begin and end positions 
    entry_beg_pos <- cumsum(c(1, n_of_items[-n_of_itemsets])) 
    entry_end_pos <- entry_beg_pos+n_of_items-1 
 
    ## prepare output 
    n_of_rows <- sum(n_of_items) 
    quality <- quality(x) 
    ## Output. 
    out <- matrix("", nrow = n_of_rows+1, ncol = 2 + NCOL(quality)) 
 
    ## Column 1: itemset nr. 
    tmp <- rep.int("", n_of_rows + 1) 
    tmp[entry_beg_pos+1] <- c(1:n_of_itemsets) 
    out[,1] <- format(tmp) 
 
    ## Column 2: items in the item sets, one per line. 
    pre <- rep.int(" ", n_of_rows) 
    pre[entry_beg_pos] <- rep.int(setStart, length(entry_beg_pos)) 
    post <- rep.int(itemSep, n_of_rows) 
    post[entry_end_pos] <- rep.int(setEnd, length(entry_end_pos)) 
    out[, 2] <- format(c("items", 
                         paste(pre, unlist(items), post, sep = ""))) 
 
 
    ## Remaining columns: quality measures. 
    for(i in seq(length = NCOL(quality))) { 
      tmp <- rep.int("", n_of_rows + 1) 
      tmp[1] <- names(quality)[i] 
      tmp[entry_end_pos + 1] <- format(quality[[i]]) 
      out[, i + 2] <- format(tmp, justify = "right") 
    } 
 
    ## Output. 
    cat(t(out), sep = c(rep.int(" ", NCOL(out) - 1), "\n")) 
  } 
  return(out) 
} 
 
