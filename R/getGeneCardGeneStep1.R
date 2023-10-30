getGeneCardGeneStep1 = function(keywords){
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
  #setwd('D:/PA2.1Plus/')
  url1='https://www.malacards.org/search/results?query='
  #url2='%27&searchType=Keywords'
  #keywords='fatty liver'
  keywords2=str_replace_all(keywords,' ','+')
  url=paste(url1,keywords2,sep='')
  UA <- "Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:93.0) Gecko/20100101 Firefox/93.0"###update 20230212
  a=GET(url,add_headers(`Connection` = "keep-alive", `User-Agent` = UA),config(sslversion=6,ssl_verifypeer=1))###update 20230212
  b=httr::content(a,'text')
  bb=read_html(b)
  tablesBB=html_table(bb)
  #tablesBB=tablesBB[[1]]###Before 20230212
  tablesBB=tablesBB[[2]]###update 20230212
  tablesBB=tablesBB[!is.na(tablesBB$MIFTS),]
  #TempsimScore=stringsim(keywords,tablesBB$Name,method='soundex')
  #tablesBB$simScore=TempsimScore
  tablesBB=as.data.table(tablesBB)
  tablesBB=tablesBB[order(-Score)]
  tablesBB=tablesBB[,-2]
  #write.csv(tablesBB,'Report/GeneCardSearchStep1.csv',row.names = F)
  return(tablesBB)
}

