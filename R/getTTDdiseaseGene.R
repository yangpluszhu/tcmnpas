getTTDdiseaseGene = function(keywords){
  require(data.table)
  require(stringr)
  require(clusterProfiler)
  if (!exists('TTDdiseasedatabase')){
    load('db/TTDdiseasedatabase.RData')
  }
  keywords=tolower(keywords)
  result=TTDdiseasedatabase[TTDdiseasedatabase$disease%like%keywords,]
  if (nrow(result)!=0){
    geneID=unique(result$geneID)
    #genesTran=bitr(geneID, fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Hs.eg.db")
    #genes=genesTran$SYMBOL
	genes=unique(result$geneName)
  }else{
    geneID=''
    genes=''
  }
  return(list(geneNames=genes,geneID=geneID,result=result))
}

