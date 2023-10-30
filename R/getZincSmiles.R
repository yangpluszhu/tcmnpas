getZincSmiles = function(zinc_id){
  library(XML)
  library(stringr)
  library(RCurl)
  myheader=c('User-Agent'='Mozilla/5.0 (Windows NT 6.1; WOW64)','Accept'='text/html,application/xhyml+xml,application/xml;q=0.9,*/*;q=0.8','Accept-Language'='zh-CN,zh;q=0.8,en;q=0.6','Connection'='keep-alive','Accept-Charset'='GB2312,utf-8;q=0.7,*;q=0.7')
  name=zinc_id
  a=paste('http://zinc.docking.org/substance/',name,sep='')
  a=htmlParse(a)
  b=getNodeSet(a,"//input[@id='item-smiles']")
  smiles=xmlToList(b[[1]])['value']
  return(smiles)
}

