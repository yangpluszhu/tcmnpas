getGeneCardGene = function(keywords,score=10){ 
  require(data.table) 
  require(stringr) 
  require(RCurl) 
  require(clusterProfiler) 
  require(readr) 
  #http://www.genecards.org/Search/Export?queryString=%27fat%20liver%27&searchType=Keywords 
  url1='http://www.genecards.org/Search/Export?queryString=%27' 
  url2='%27&searchType=Keywords' 
  keywords=str_replace_all(keywords,' ','%20') 
  url=paste(url1,keywords,url2,sep='') 
  ## 
  myHttpheader=c('User-Agent'='Mozilla/5.0(Windows;U;Windows NT 6.1;zh-CN;rv:1.9.1.6)','Accept'='text/html,application/xhyml+xml,application/xml;q=0.9,*/*;q=0.8','Accept-Language'='en-us','Connection'='keep-alive','Accept-Charset'='GB2312,utf-8;q=0.7,*;q=0.7') 
  d2 =debugGatherer() 
  cHandle2<- getCurlHandle(httpheader=myHttpheader,followlocation=1, 
                           debugfunction=d2$update,verbose=TRUE) 
  temp2<- getURLContent(url,curl=cHandle2,.encoding="UTF-8") 
  result=read_csv(temp2) 
  result=result[result$`Relevance score`>score,] 
  if (nrow(result)!=0){ 
    geneName0=unique(result$`Gene Symbol`) 
    geneID=geneTrans(geneName0,type = 'NAME') 
    geneID=unique(geneID[!is.na(geneID)]) 
    geneName=geneTrans(geneID,type = 'ID') 
    #genesTran=bitr(geneID, fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Hs.eg.db") 
    #genes=genesTran$SYMBOL 
  }else{ 
    geneID='' 
    geneName='' 
  } 
  return(list(geneNames=geneName,geneID=geneID,result=result)) 
} 
 
