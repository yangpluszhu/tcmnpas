checkNodeType = function(id){
  library(stringr)
  idTrans=as.numeric(id)
  idTrans2=str_detect(id,'c')
  geneType=rep('gene',length(id))
  chemType=rep('chem',length(id))
  herbType=rep('herb',length(id))
  nodeType=ifelse(is.na(idTrans),ifelse(idTrans2,chemType,herbType),geneType)
  return(nodeType)
}

