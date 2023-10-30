getDiseaseGene = function(keyWords,DGdatabase=c('OMIM','TTD','DC','GeneCard'),GeneCardScore=10){
  require(data.table)
  library(RCurl)
  library(stringr)
  library(XML)
  library(readr)
  require(clusterProfiler)
  GeneCardScore=as.numeric(GeneCardScore)
  if ('OMIM'%in%DGdatabase){
    result1=getOMIMdiseaseGene(keyWords)
	resultFromDC=getDCdiseaseGene(keywords=keyWords,Type='OMIM')
	result1$geneNames=unique(c(result1$geneNames,resultFromDC$geneNames))
	result1$geneID=unique(c(result1$geneID,resultFromDC$geneID))
  }else{
    result1=list(geneNames=character(),geneID=character(),result=data.frame())
  }
  ##
  if ('TTD'%in%DGdatabase){
    result2=getTTDdiseaseGene(keyWords)
  }else{
    result2=list(geneNames=character(),geneID=character(),result=data.frame())
  }
  ##
  if ('DC'%in%DGdatabase){
    result3=getDCdiseaseGene(keyWords)
  }else{
    result3=list(geneNames=character(),geneID=character(),result=data.frame())
  }
  #
  if ('GeneCard'%in%DGdatabase){
    result4=getGeneCardGene(keywords=keyWords,score=GeneCardScore)
  }else{
    result4=list(geneNames=character(),geneID=character(),result=data.frame())
  }
  ##
  result=list(geneNames=unique(c(result1$geneNames,result2$geneNames,result3$geneNames,result4$geneNames)),geneID=unique(c(result1$geneID,result2$geneID,result3$geneID,result4$geneID)),result=list(OMIM=result1$result,TTD=result2$result,DC=result3$result,GeneCard=result4$result))
  return(result)
}

