drug_Onlie_search = function(drug){
  url0="http://baike.baidu.com/search/word?word="
  url=paste(url0,drug,sep='')
  Order=paste('start',url,sep=' ')
  shell(Order,wait = F)
}

