getCTDId_and_pubID = function(name,cookiefileS=NULL){####chemical_name--->CTDMeshID&pubID
  library(RCurl)
  library(stringr)
  library(rvest)
  library(httr)
  library(data.table)
  #name='DDFG'
  name=tolower(name)
  name2=str_replace_all(name,' ','%20')
  name3=paste('%22',name2,'%22',sep='')
  myheader=c('User-Agent'='Mozilla/5.0 (Windows NT 6.1; WOW64)','Accept'='text/html,application/xhyml+xml,application/xml;q=0.9,*/*;q=0.8','Accept-Language'='zh-CN,zh;q=0.8,en;q=0.6','Connection'='keep-alive','Accept-Charset'='GB2312,utf-8;q=0.7,*;q=0.7')
  chemCSV_url=paste('http://ctdbase.org/basicQuery.go?bq=',name3,'&6578706f7274=1&d-4029212-e=1&bqCat=chem',sep='')
  if (c(class(try({chemTable=suppressWarnings(fread(chemCSV_url))},silent=T))=='try-error')[1]){
    MeshID=NA
    pubID=NA
    sure=NA
  }else{
    #chemTable=suppressWarnings(fread(chemCSV_url))
    if (sum(str_detect(chemTable,'<meta'))>0){
      d2 =debugGatherer()
      cHandle2<- getCurlHandle(httpheader=myheader,followlocation=1,
                               debugfunction=d2$update,verbose=TRUE,
                               cookiefile=cookiefileS)
      au=getURL(chemCSV_url,curl=cHandle2,.encoding="UTF-8")
      MeshID0=str_extract(au,'term=[A-Z0-9]+')
      MeshID=str_replace(MeshID0,'term=','')
    }else{
      chemTable[,Chemical:=tolower(Chemical)]
      setkey(chemTable,Chemical)
      if (is.na(as.character(chemTable[name,'Accession ID',with=F]))){
        MeshID=as.character(chemTable[1,'Accession ID',with=FALSE])
      }else{
        MeshID=as.character(chemTable[name,'Accession ID',with=F])
      }
    }
    pubSurl=paste('http://www.ncbi.nlm.nih.gov/pcsubstance/?term=',MeshID,sep='')
    d2 =debugGatherer()
    cHandle2<- getCurlHandle(httpheader=myheader,followlocation=1,
                             debugfunction=d2$update,verbose=TRUE,
                             cookiefile=cookiefileS)
    au<- getURL(pubSurl,curl=cHandle2,.encoding="UTF-8")
    #au=getURL(pubSurl,httpheader=myheader,followlocation=T)
    if (str_detect(au,'No items found')){
      pubID=NA
      sure=NA
    }else{
      sid0=str_extract(au,'SID:[0-9]+')
      sid=str_extract(sid0,'[0-9]+')
      url_sid=paste('pubchem.ncbi.nlm.nih.gov/rest/rdf/substance/SID',sid,'.html',sep='')
      a=getURL(url_sid,curl=cHandle2,.encoding="UTF-8")
      #       lo=str_locate(au,'data-pubchem-title=')
      #       lo2=str_locate(au,'data-pubchem-record')
      #       temp=str_sub(au,lo[1,2]+1,lo2[1,1]-3)
      #       tempName=str_replace_all(temp,'\"','')
      #       result=getPubID(tempName)
      tempcid=str_extract(a,'CID[0-9]+')
      pubID=str_extract(tempcid,'[0-9]+')
      sure='ok'
    }

  }
  Re=c(MeshID=MeshID,pubID=pubID,sure=sure)
  names(Re)=c('MeshID','pubID','sure')
  return(Re)
}

