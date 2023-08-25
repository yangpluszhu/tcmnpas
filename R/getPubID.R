getPubID = function(name,cookiefileS=NULL){###chemical_name-->get PubmedChemID 
  require(RCurl) 
  require(stringr) 
  require(rvest) 
  require(httr) 
  require(XML) 
  name=tolower(name) 
  #name='protocatechualdehyde' 
  if (str_detect(name,'\\s')){ 
    name2=str_replace_all(name,' ','+') 
    #name3=paste('%22',name2,'%22',sep='') 
    name3=name2 
  }else{ 
    name3=name 
  } 
  url=paste('http://www.ncbi.nlm.nih.gov/pccompound/?term=',name3,sep='') 
  myheader=c('User-Agent'='Mozilla/5.0 (Windows NT 6.1; WOW64)','Accept'='text/html,application/xhyml+xml,application/xml;q=0.9,*/*;q=0.8','Accept-Language'='zh-CN,zh;q=0.8,en;q=0.6','Connection'='keep-alive','Accept-Charset'='GB2312,utf-8;q=0.7,*;q=0.7') 
  d2 =debugGatherer() 
  cHandle2<- getCurlHandle(httpheader=myheader,followlocation=1, 
                           debugfunction=d2$update,verbose=TRUE, 
                           cookiefile=cookiefileS) 
  au<- getURL(url,curl=cHandle2,.encoding="UTF-8") 
  #au=getURL(url,httpheader=myheader,followlocation=T) 
  if (str_detect(au,'title=\"CID:')){ 
    #txt=str_sub(au,str_locate(au,'title=\"CID:')[1,2]+1,str_locate(au,'title=\"CID:')[1,2]+15) 
    txt=str_extract(au,'CID[0-9]+') 
    cid=str_extract(txt,pattern='[0-9]+') 
    sure='ok' 
  }else if (str_detect(au,'No items found')){ 
    cid=NA 
    sure=NA 
  }else if (str_detect(au,'Pccompound_ResultsPanel')){ 
    #   txt=str_sub(au,str_locate(au,'CID:')[1,2]+1,str_locate(au,'CID:')[1,2]+25) 
    #   cid=str_extract(txt,pattern='[0-9]+') 
    #   sure='multi' 
    #aurl=getURL(url,httpheader=myheader,followlocation=T) 
    bparse=htmlParse(au) 
    bparse2=xpathSApply(bparse,"//p[@class='title']",fun=xmlValue) 
    bparse0=xpathSApply(bparse,"//p[@class='title']") 
    attrtemp=lapply(bparse0,xmlToList) 
    href=sapply(attrtemp,tempcatch) 
    tempsplit=lapply(bparse2,str_split,pattern=';') 
    tempsplit2=lapply(tempsplit,unlist) 
    exact_sum=unlist(lapply(tempsplit2,function(x)sum(tolower(x)%in%name))) 
    if (sum(exact_sum)>0){ 
      cid=href[which(exact_sum==1)[1]] 
      sure='ok' 
    }else{ 
      cid=href[1] 
      sure='multi' 
    } 
  } 
  return(c(cid=cid,sure=sure))##cid-->PubmedChemID;sure-->ok:only one match;multi-->multi match;NA-->not available 
} 
 
