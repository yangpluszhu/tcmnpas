PathWayScore = function(geneID1,geneID1Name='Name1',geneID2,geneID2Name='Name2',method=c('WNS','keggSim')){ 
  require(data.table) 
  require(plyr) 
  require(stringr) 
  ################################ 
  if (str_detect(geneID1,',|;| ')){ 
    geneID1=unlist(str_split(geneID1,',|;| ')) 
    geneID1=str_replace_all(geneID1,'"','') 
    geneID1=str_replace_all(geneID1,"'",'') 
    geneID1=geneID1[geneID1!=''] 
    geneID1=str_trim(geneID1,side='both') 
  } 
  geneID1=unique(geneID1) 
  Data1=data.table(id=geneID1Name,geneID=geneID1) 
  ############################## 
  if (str_detect(geneID2,',|;| ')){ 
    geneID2=unlist(str_split(geneID2,',|;| ')) 
    geneID2=str_replace_all(geneID2,'"','') 
    geneID2=str_replace_all(geneID2,"'",'') 
    geneID2=geneID2[geneID2!=''] 
    geneID2=str_trim(geneID2,side='both') 
  } 
  geneID2=unique(geneID2) 
  Data2=data.table(id=geneID2Name,geneID=geneID2) 
  ############################## 
  Data=rbind(Data1,Data2) 
  ############################# 
  if (length(method)==2){ 
    result1=WNSscore(Data) 
    result2=keggSimScore(Data) 
    result=data.frame(geneID1Name=geneID1Name,geneID2Name=geneID2Name,WNSscore=result1[geneID1Name,geneID2Name],keggSimscore=result2[geneID1Name,geneID2Name]) 
  }else if ('WNS'%in%method){ 
    result1=WNSscore(Data) 
    result=data.frame(geneID1Name=geneID1Name,geneID2Name=geneID2Name,WNSscore=result1[geneID1Name,geneID2Name]) 
  }else if ('keggSim'%in%method){ 
    result2=keggSimScore(Data) 
    result=data.frame(geneID1Name=geneID1Name,geneID2Name=geneID2Name,keggSimscore=result2[geneID1Name,geneID2Name]) 
  } 
  return(result) 
} 
 
