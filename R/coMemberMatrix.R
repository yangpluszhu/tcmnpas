coMemberMatrix = function(BF_result){
  require(data.table)
  require(igraph)
  BF_result2=BF_result$BF_result
  NBF=nrow(BF_result2)
  CoMatrix=matrix(0,nrow=NBF,ncol=NBF)
  colnames(CoMatrix)=BF_result2$id
  rownames(CoMatrix)=BF_result2$id
  for (i in 1:NBF){
    for (j in 1:NBF){
      tempi=BF_result2$id[i]
      tempj=BF_result2$id[j]
      makeupi=BF_result2$makeup[BF_result2$id==tempi]
      makeupi2=unlist(str_split(makeupi,','))
      makeupj=BF_result2$makeup[BF_result2$id==tempj]
      makeupj2=unlist(str_split(makeupj,','))
      CoMatrix[tempi,tempj]=length(intersect(makeupi2,makeupj2))
    }
  }
  CoMatrix=as.data.frame(CoMatrix)
  return(CoMatrix)
}

