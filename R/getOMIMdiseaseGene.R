getOMIMdiseaseGene = function(keyWords){ 
  require(RCurl) 
  require(stringr) 
  require(XML) 
  require(readr) 
  require(clusterProfiler) 
  keyWords=str_replace_all(keyWords,' ','+') 
  url=paste('http://www.omim.org/search/?index=geneMap&search=','"',keyWords,'"','&start=1&limit=20000&format=tsv',sep='') 
  myHttpheader=c('User-Agent'='Mozilla/5.0(Windows;U;Windows NT 6.1;zh-CN;rv:1.9.1.6)','Accept'='text/html,application/xhyml+xml,application/xml;q=0.9,*/*;q=0.8','Accept-Language'='en-us','Connection'='keep-alive','Accept-Charset'='GB2312,utf-8;q=0.7,*;q=0.7') 
  d2 =debugGatherer() 
  cHandle2<- getCurlHandle(httpheader=myHttpheader,followlocation=1, 
                           debugfunction=d2$update,verbose=TRUE) 
  temp2<- getURL(url,curl=cHandle2,.encoding="UTF-8") 
  content=htmlTreeParse(temp2, asText = TRUE) 
  a=content$children$html[['body']] 
  b=a[['p']][['text']] 
  b=as.character(b) 
  write_lines(b,'temp.txt') 
  bb=read_lines('temp.txt') 
  blank=which(bb=='') 
  Copyright=which(str_detect(bb,'Copyright')) 
  lo1=blank[which(blank>=Copyright)[1]] 
  lo2=blank[which(blank>=Copyright)[2]] 
  result=character() 
  for (i in (lo1+1):(lo2-1)){ 
    tempD=unlist(str_split_fixed(bb[i],'\t',n=11)) 
    result=rbind(result,tempD) 
  } 
  colnames(result)=unlist(result[1,]) 
  result=result[-1,,drop=F] 
  #result=apply(result,2,as.character) 
  result=as.data.frame(result,stringsAsFactors = F) 
  if (nrow(result)==0){ 
    genes='' 
    geneID='' 
  }else{ 
    genes=paste(result$`Gene/Locus`,collapse = ',') 
    genes=unique(str_trim(unlist(str_split(genes,',')),side='both')) 
    #genesTran=bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db") 
    #geneID=genesTran$ENTREZID 
    if (class(try({geneID=geneTrans(genes,type = 'ID')},silent = T))=='try-error')geneID=NA 
  } 
 
  return(list(geneNames=genes,geneID=geneID,result=result)) 
} 
 
