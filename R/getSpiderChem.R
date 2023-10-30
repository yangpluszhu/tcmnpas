getSpiderChem = function(name,cookiefileS=NULL){####chemical_name-->get SpederID&Smiles
  library(XML)
  library(stringr)
  library(RCurl)
  name=tolower(name)
  #name='Methyl palmitate'
  name=str_replace_all(name,' ','%20')
  if (str_detect(name,"'")){
    name2=str_replace(name,"'","")
    a0=paste('http://rdf.chemspider.com/search/',"%22",name2,"%22",sep='')
  }else{
    a0=paste('http://rdf.chemspider.com/search/',name,sep='')
  }
  myheader=c('User-Agent'='Mozilla/5.0 (Windows NT 6.1; WOW64)','Accept'='text/html,application/xhyml+xml,application/xml;q=0.9,*/*;q=0.8','Accept-Language'='zh-CN,zh;q=0.8,en;q=0.6','Connection'='keep-alive','Accept-Charset'='GB2312,utf-8;q=0.7,*;q=0.7')
  d2 =debugGatherer()
  cHandle2<- getCurlHandle(httpheader=myheader,followlocation=1,
                           debugfunction=d2$update,verbose=TRUE,
                           cookiefile=cookiefileS)
  a<- getURL(a0,curl=cHandle2,.encoding="UTF-8")
  #a=getURL(a,httpheader=myheader,followlocation=T)
  if (str_detect(a,'(no matches)|(Error)|(Bad Request)')){
    chemspider_id=NA
    smiles=NA
    sure=NA
  }else{
    if(str_detect(a,'ChemSpider ID:')){
      b1=htmlParse(a)
      c1=xpathSApply(b1,"//p/span",fun=xmlValue)
      chemspider_id=c1[[1]]
      smiles=c1[[4]]
      sure='multi'
    }else{
      tempa=str_extract(a,'Chemical-Structure\\.[0-9]+\\.rdf')
      chemspider_id=str_extract(tempa,'[0-9]+')
      b1=htmlParse(a)
      c1=xpathSApply(b1,"//smiles",fun=xmlValue)
      smiles=c1
      sure='ok'
    }
  }
  return(c(SpiderID=chemspider_id,smiles=smiles,sure=sure))
}

