geneTrans = function(geneTerms,type="ID"){ 
  ##geneTerms:geneID or geneName 
  ##type:'ID'-->ID to Name;'NAME'-->Name to ID 
  require(clusterProfiler) 
  require(data.table) 
  require(plyr) 
  geneTerms=as.character(geneTerms) 
  type_ID="ID" 
  type_Name="NAME" 
  if (type == type_ID){ 
    eg = bitr(geneTerms, fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Hs.eg.db") 
    #eg=as.data.table(eg) 
    #setkey(eg,ENTREZID) 
    #SYMBOLname=eg[.(geneTerms),SYMBOL,mult='first'] 
	SYMBOLname=eg$SYMBOL[eg$ENTREZID%in%geneTerms][1] 
    result=SYMBOLname 
   }else if (type == type_Name){ 
    eg = bitr(geneTerms, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db") 
    #eg=as.data.table(eg) 
    #setkey(eg,SYMBOL) 
    #ENTREZID0=eg[.(geneTerms),ENTREZID,mult='first'] 
	ENTREZID0=eg$ENTREZID[eg$SYMBOL%in%geneTerms][1] 
	result=ENTREZID0 
   } 
  return(result) 
} 
 
