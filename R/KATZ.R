KATZ = function(SeedGeneID,ScoreGeneID,randomtimes=1000,testM=1){
  #SeedGeneID-->疾病靶标
  #ScoreGeneID-->药物靶标
  #testM:1-->wilcox;testM:2-->perm
  require(stringr)
  if (!exists('ppiNetMatrix')){
  load('db/PPIpath.db')
	#load('db/PPIpath.db')
  }
  ###
  if (sum(str_detect(SeedGeneID,',| |;'))>0){
    SeedGeneID=as.character(unlist(str_split(SeedGeneID,',| |;')))
    SeedGeneID=str_replace_all(SeedGeneID,'"','')
    SeedGeneID=str_replace_all(SeedGeneID,"'",'')
    SeedGeneID=SeedGeneID[SeedGeneID!='']
    SeedGeneID=unique(str_trim(SeedGeneID,side='both'))
  }
  SeedGeneID=SeedGeneID[SeedGeneID%in%colnames(ppiNetMatrix)]
  ##
  if (sum(str_detect(ScoreGeneID,',| |;'))>0){
  ScoreGeneID=as.character(unlist(str_split(ScoreGeneID,',| |;')))
  ScoreGeneID=str_replace_all(ScoreGeneID,'"','')
  ScoreGeneID=str_replace_all(ScoreGeneID,"'",'')
  ScoreGeneID=ScoreGeneID[ScoreGeneID!='']
  ScoreGeneID=unique(str_trim(ScoreGeneID,side='both'))
  }
  ScoreGeneID=ScoreGeneID[ScoreGeneID%in%colnames(ppiNetMatrix)]
  ##
  ScoreGeneData=data.frame(SeedGeneIDPPI=toString(SeedGeneID),NumSeedGeneID=length(SeedGeneID),ScoreGeneIDPPI=toString(ScoreGeneID),NumScoreGeneID=length(ScoreGeneID),overlapScore=0,Path1Score=0,Path2Score=0,Path3Score=0,TotalScore=0,RandomScoreMedian=0,Pvalue=1)
  if (length(SeedGeneID)!=0&length(ScoreGeneData)!=0){
  tempGene=SeedGeneID
  overlapScore=length(intersect(tempGene,ScoreGeneID))
  N1Score=sum(ppiNetMatrix[tempGene,ScoreGeneID])*0.001
  N2Score=sum(ppiNetMatrix2[tempGene,ScoreGeneID])*0.000001
  N3Score=sum(ppiNetMatrix3[tempGene,ScoreGeneID])*0.000000001
  ScoreGeneData[1,'overlapScore']=overlapScore
  ScoreGeneData[1,'Path1Score']=N1Score
  ScoreGeneData[1,'Path2Score']=N2Score
  ScoreGeneData[1,'Path3Score']=N3Score
  ScoreGeneData[1,'TotalScore']=overlapScore+N1Score+N2Score+N3Score
  }
  #####################
  #randomScore=numeric()
  randomScore2=numeric()
  seedNum=1234567
  for (i in 1:randomtimes){
    #randomGeneID=sample(colnames(ppiNetMatrix),length(SeedGeneID))
    #overlapScore=length(intersect(randomGeneID,ScoreGeneID))
    #N1Score=sum(ppiNetMatrix[randomGeneID,ScoreGeneID])*0.001
    #N2Score=sum(ppiNetMatrix2[randomGeneID,ScoreGeneID])*0.000001
    #N3Score=sum(ppiNetMatrix3[randomGeneID,ScoreGeneID])*0.000000001
    #randomScore[i]=overlapScore+N1Score+N2Score+N3Score
    #####
    set.seed(seedNum+i)
    randomGeneID2=sample(colnames(ppiNetMatrix),length(ScoreGeneID))
    overlapScore2=length(intersect(randomGeneID2,SeedGeneID))
    N1Score2=sum(ppiNetMatrix[randomGeneID2,SeedGeneID])*0.001
    N2Score2=sum(ppiNetMatrix2[randomGeneID2,SeedGeneID])*0.000001
    N3Score2=sum(ppiNetMatrix3[randomGeneID2,SeedGeneID])*0.000000001
    randomScore2[i]=overlapScore2+N1Score2+N2Score2+N3Score2
  }
  #Pvalue1=sum(randomScore>ScoreGeneData[1,'TotalScore'])/randomtimes
  if (testM==1){
  #Pvalue1=wilcox.test(randomScore,mu=ScoreGeneData[1,'TotalScore'],alternative='less')$p.value
  Pvalue2=wilcox.test(randomScore2,mu=ScoreGeneData[1,'TotalScore'],alternative='less')$p.value
  #ScoreGeneData[1,'Pvalue']=max(Pvalue1,Pvalue2)
  }else{
  Pvalue2=sum(randomScore2>ScoreGeneData[1,'TotalScore'])/randomtimes
  }
  ScoreGeneData[1,'Pvalue']=Pvalue2
  #Pvalue2=sum(randomScore2>ScoreGeneData[1,'TotalScore'])/randomtimes
  ScoreGeneData[1,'RandomScoreMedian']=paste(mean(randomScore2),' [',quantile(randomScore2,0.025),',',quantile(randomScore2,0.975),']',sep='')
  rm('ppiNetMatrix')
  rm('ppiNetMatrix2')
  rm('ppiNetMatrix3')
  gc()
  return(ScoreGeneData)
}

