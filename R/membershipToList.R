membershipToList = function(Members){ 
  require(igraph) 
  clusterName=unique(Members) 
  Nnodes=length(Members) 
  LM=list() 
  length(LM)=length(clusterName) 
  names(LM)=clusterName 
  for (k in clusterName){ 
    LM[[k]]=names(Members)[which(Members==k)] 
  } 
  return(LM) 
} 
 
