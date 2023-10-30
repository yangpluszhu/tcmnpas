DisAccess = function(SeedGeneID,ScoreGeneID,randomtimes=1000,testM=1){
  require(data.table)
  require(plyr)
  require(stringr)
  require(igraph)
  #testM:1-->wilcox;testM:2-->perm
  #require("org.Hs.eg.db")
  require(pracma)
  if (!exists('ppiBinaryNet')){
   load('db/ppiNetData.db')
	 #load('db/ppiNetData.db')
  }
  ##
  if (sum(str_detect(SeedGeneID,',| |;'))>0){
    SeedGeneID=as.character(unlist(str_split(SeedGeneID,',| |;')))
    SeedGeneID=str_replace_all(SeedGeneID,'"','')
    SeedGeneID=str_replace_all(SeedGeneID,"'",'')
    SeedGeneID=SeedGeneID[SeedGeneID!='']
    SeedGeneID=unique(str_trim(SeedGeneID,side='both'))
  }
  SeedGeneID=as.character(unique(SeedGeneID))
  SeedGeneID=SeedGeneID[SeedGeneID%in%V(ppiBinaryNet)$name]
  ##
  if (sum(str_detect(ScoreGeneID,',| |;'))>0){
    ScoreGeneID=as.character(unlist(str_split(ScoreGeneID,',| |;')))
    ScoreGeneID=str_replace_all(ScoreGeneID,'"','')
    ScoreGeneID=str_replace_all(ScoreGeneID,"'",'')
    ScoreGeneID=ScoreGeneID[ScoreGeneID!='']
    ScoreGeneID=unique(str_trim(ScoreGeneID,side='both'))
  }
  ScoreGeneID=as.character(unique(ScoreGeneID))
  ScoreGeneID=ScoreGeneID[ScoreGeneID%in%V(ppiBinaryNet)$name]
  ##
  ppiDiameter=diameter(ppiBinaryNet)
  ##############
  Yan_inPPI_geneID=SeedGeneID
  BXXXT_inPPI_geneID=ScoreGeneID
  Dis_Yan=distances(ppiBinaryNet,BXXXT_inPPI_geneID,Yan_inPPI_geneID)
  dtimes=randomtimes
  #BXXXT_DRandomYan=matrix(0,nrow=dtimes,ncol=length(BXXXT_inPPI_geneID))
  ##
  BXXXT_DRandomYan=numeric()
  seedNum=12345678
  for (i in 1:dtimes){
    set.seed(seedNum+i)
    RandomNodes=sample(V(ppiBinaryNet)$name,length(BXXXT_inPPI_geneID))
    ##
    RdisY=distances(ppiBinaryNet,RandomNodes,Yan_inPPI_geneID)
    RdisY[RdisY==Inf]<-ppiDiameter
    #RdisY=rowMeans(RdisY)
    RdisY=mean(RdisY)
    #BXXXT_DRandomYan[i,]=RdisY
	BXXXT_DRandomYan[i]=RdisY
    ##
  }
  #########
  rm('ppiBinaryNet')
  gc()
  BXXXT_DRandomYan2=BXXXT_DRandomYan
  #BXXXT_DRandomYan2[which(BXXXT_DRandomYan2==Inf)]=NA
  #Dis_YanScore=rowMeans(Dis_Yan)
  Dis_YanScore=mean(Dis_Yan)
  #BXXXT_DRandomYan2Score=rowMeans(BXXXT_DRandomYan2)
  BXXXT_DRandomYan2Score=BXXXT_DRandomYan2
  if (testM==1){
  Pvalue=wilcox.test(Dis_YanScore,BXXXT_DRandomYan2Score)$p.value
  #Pvalue=sum(BXXXT_DRandomYan2Score<median(Dis_YanScore))/randomtimes
  }else{
  Pvalue=sum(BXXXT_DRandomYan2Score<median(Dis_YanScore))/randomtimes
  }
  
  return(list(SeedGeneID=SeedGeneID,NumSeedGeneID=length(SeedGeneID),ScoreGeneID=ScoreGeneID,NumScoreGeneID=length(ScoreGeneID),MedianDistance=paste(Dis_YanScore,' [',quantile(Dis_Yan,0.025),',',quantile(Dis_Yan,0.975),']',sep=''),RandomScoreMedian=paste(mean(BXXXT_DRandomYan2Score),' [',quantile(BXXXT_DRandomYan2Score,0.025),',',quantile(BXXXT_DRandomYan2Score,0.975),']',sep=''),Pvalue=Pvalue))
}

