##' Make the matrix Symmetric
##'
##' @title asSymmetric
##' @param x a matrix
##' @rule Symmetric rule:'upper'
##' @return a matrix object
##' @export
##' @author Yang Ming
asSymmetric = function(x,rule='upper'){ 
  require(sna) 
  y=symmetrize(x,rule=rule) 
  detach('package:sna',unload=TRUE) 
  rownames(y)=rownames(x) 
  colnames(y)=colnames(x) 
  return(y) 
} 
 
