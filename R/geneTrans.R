geneTrans = function(geneTerms,type){
  #geneTerms:geneID or geneName
  #type:'ID'-->ID to Name;'NAME'-->Name to ID
  require(clusterProfiler)
  require(data.table)
  require(plyr)
  geneTerms=as.character(geneTerms)
  if (type=='ID'){
    eg = bitr(geneTerms, fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Hs.eg.db")
    eg=as.data.table(eg)
    setkey(eg,ENTREZID)
    SYMBOL=eg[.(geneTerms),SYMBOL,mult='first']
    return(SYMBOL)
  }else if (type=='NAME'){
    eg = bitr(geneTerms, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
    eg=as.data.table(eg)
    setkey(eg,SYMBOL)
    ENTREZID=eg[.(geneTerms),ENTREZID,mult='first']
    return(ENTREZID)
  }
}

