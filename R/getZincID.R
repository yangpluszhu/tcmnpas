getZincID = function(name,cookiefileS=NULL){###chemical_name-->get ZincID 
  require(XML) 
  require(stringr) 
  require(RCurl) 
  name=tolower(name) 
  name=str_replace_all(name,' ','%20') 
  a=paste('http://zinc.docking.org/synonym/',name,sep='') 
  myheader=c('User-Agent'='Mozilla/5.0 (Windows NT 6.1; WOW64)','Accept'='text/html,application/xhyml+xml,application/xml;q=0.9,*/*;q=0.8','Accept-Language'='zh-CN,zh;q=0.8,en;q=0.6','Connection'='keep-alive','Accept-Charset'='GB2312,utf-8;q=0.7,*;q=0.7') 
  d2 =debugGatherer() 
  cHandle2<- getCurlHandle(httpheader=myheader,followlocation=1, 
                           debugfunction=d2$update,verbose=TRUE, 
                           cookiefile=cookiefileS) 
  a<- getURL(a,curl=cHandle2,.encoding="UTF-8") 
  #a=getURL(a,httpheader=myheader,followlocation=T) 
  b=htmlParse(a) 
  c=getNodeSet(b,"//input[@type='checkbox']") 
  if (class(try({d=xmlToList(c[[1]])['value']},silent =T))!='try-error'){ 
    #d=xmlToList(c[[1]])['value'] 
    chemical_Zid=paste('ZINC',d,sep='') 
  }else{ 
    chemical_Zid=NA 
  } 
  return(ZincID=chemical_Zid) 
} 
 
