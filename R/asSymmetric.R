asSymmetric = function(x,rule='upper'){ 
  require(sna) 
  y=symmetrize(x,rule=rule) 
  detach('package:sna',unload=TRUE) 
  rownames(y)=rownames(x) 
  colnames(y)=colnames(x) 
  return(y) 
} 
 
