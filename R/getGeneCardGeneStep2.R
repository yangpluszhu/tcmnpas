getGeneCardGeneStep2 = function(MCID){
  require(data.table)
  require(stringr)
  require(RCurl)
  require(clusterProfiler)
  require(readr)
  require(rvest)
  require(stringdist)
  #library(reticulate)
  #library(jsonlite)
  #library(xml2)
  library("httr")
  options(stringsAsFactors = F)
  Symbol_detect=function(x){
    require(stringr)
    y=numeric()
    temp=unlist(x)
    y=sum(str_detect(unlist(temp),'Symbol'),na.rm =T)
    return(y)
  }
  #setwd('D:/PA2.1Plus/')
  MCID2=unlist(str_split(MCID,',|，|;'))
  MCID2=str_trim(MCID2,'both')
  MCID2=unique(MCID2)
  url1='https://www.malacards.org/card/'
  urlEND='?limit[RelatedGenes]=0#RelatedGenes-table'
  UA <- "Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:93.0) Gecko/20100101 Firefox/93.0"###update 20230212
  searchResult=data.frame()
  for (i in MCID2){
    url=paste0(url1,i,urlEND)
    a=GET(url,add_headers(`Connection` = "keep-alive", `User-Agent` = UA),config(sslversion=6,ssl_verifypeer=1))###update 20230212
    b=httr::content(a,'text')
    bb=read_html(b)
    tableNOde=html_node(x=bb,xpath='//*[@id="RelatedGenes-table"]')
    if (!is.na(tableNOde)){
      tablesBB=html_table(tableNOde)
      #GeneID=geneTrans(geneTerms=tablesBB$Symbol,type='NAME')
      if (class(try({GeneID=geneTrans(geneTerms=tablesBB$Symbol,type='NAME')},silent = T))!="try-error"){
        tablesBB$GeneID=GeneID
        tablesBB$MCID=i
        #tablesBB=tablesBB[,c('MCID','Symbol','GeneID','Description','Category','Score','Evidence','PubmedIds')]###Before 20230212
        tablesBB=tablesBB[,c('MCID','Symbol','GeneID','Description','Category','Score','Evidence','PubMed IDs')]###update 20230212
      }else{
        tablesBB$GeneID=NA
        tablesBB$MCID=i
        #tablesBB=tablesBB[,c('MCID','Symbol','GeneID','Description','Category','Score','Evidence','PubmedIds')]###Before 20230212
        tablesBB=tablesBB[,c('MCID','Symbol','GeneID','Description','Category','Score','Evidence','PubMed IDs')]###update 20230212
      }
    }else{
      tablesBB=data.frame()
    }
    searchResult=rbind(searchResult,tablesBB)
  }
  if (nrow(searchResult)>0){
    CombineResult=searchResult[,c('Symbol','GeneID','Description','Category','Score','Evidence')]
    CombineResult=as.data.table(CombineResult)
    CombineResult=unique(CombineResult,by='GeneID')
  }else{
    CombineResult=searchResult
  }
  return(list(searchResult=searchResult,CombineResult=CombineResult))
}

