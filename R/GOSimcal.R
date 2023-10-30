GOSimcal = function(genelist1,genelist2,method=c('GO','DO')){
  require(GOSemSim)
  require(DOSE)
  require(stringr)
  #####
  if (str_detect(genelist1,',| |;')){
    genelist1=as.character(unlist(str_split(genelist1,',| |;')))
    genelist1=str_replace_all(genelist1,'"','')
    genelist1=str_replace_all(genelist1,"'",'')
    genelist1=genelist1[genelist1!='']
    genelist1=unique(str_trim(genelist1,side='both'))
  }
  genelist1=as.character(unique(genelist1))
  #####
  if (str_detect(genelist2,',| |;')){
    genelist2=as.character(unlist(str_split(genelist2,',| |;')))
    genelist2=str_replace_all(genelist2,'"','')
    genelist2=str_replace_all(genelist2,"'",'')
    genelist2=genelist2[genelist2!='']
    genelist2=unique(str_trim(genelist2,side='both'))
  }
  genelist2=as.character(unique(genelist2))
  #####
  if ('GO'%in%method){
    resultGO=GOsimilarity(genelist1=genelist1,genelist2=genelist2)
  }else{
    resultGO=NULL
  }
  ######
  if ('DO'%in%method){
    resultDO=GODOsimilarity(genelist1=genelist1,genelist2=genelist2)
  }else{
    resultDO=NULL
  }
  #####
  return(list(GO=resultGO,DO=resultDO))
}

