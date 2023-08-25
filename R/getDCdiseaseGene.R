getDCdiseaseGene = function(keywords,Type=c('DEG','OMIM','GWAS')){ 
  require(data.table) 
  require(stringr) 
  require(clusterProfiler) 
  ##get "http://disease-connect.org/" disease Genes 
#DCdatabase=fread('E:/DATABASE/Disease-Gene/Disease-Gene.csv',colClasses = 'character') 
#DCdatabase$Diseas=tolower(DCdatabase$Diseas) 
if (!exists('DCdatabase')){ 
  load('db/DCdatabase.RData') 
} 
keywords=tolower(keywords) 
result=DCdatabase[(DCdatabase$Diseas%like%keywords)&(DCdatabase$Type%in%Type),] 
if (nrow(result)!=0){ 
  genes=unique(result$GeneName) 
  #genesTran=bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db") 
  #geneID=genesTran$ENTREZID 
  geneID=unique(result$geneID) 
}else{ 
  genes='' 
  geneID='' 
} 
return(list(geneNames=genes,geneID=geneID,result=result)) 
} 
 
