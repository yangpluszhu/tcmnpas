RWanalysis = function(SeedGeneID,ScoreGeneID,scoreStand=1000,randomtimes=1000,EdgeWeight=FALSE,gamma=0.7,TopNCal=FALSE,TopN=1:10,Plot=FALSE,testM=1,ppiBinaryNet=ppiBinaryNet){ 
  #ppiBinaryNet 
  #testM:1-->wilcox;testM:2-->perm 
  require(data.table) 
  require(plyr) 
  require(stringr) 
  require(igraph) 
  require("org.Hs.eg.db") 
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
  if (!TopNCal){ 
    RWR_Seed=randomWalk(ppiBinaryNet,SeedGeneID,EdgeWeight=EdgeWeight,gamma=gamma)*scoreStand 
    TopNN=1:length(ScoreGeneID) 
    RWR_Score=numeric() 
    for (i in TopNN){ 
      tempN=TopNN[i] 
      RWR_Score[i]=sum(RWR_Seed[ScoreGeneID[1:tempN]],na.rm=T)/tempN 
      } 
    ####### 
    TopNNauc=trapz(TopNN,RWR_Score) 
    GeneScore=RWR_Seed[ScoreGeneID] 
    TotalScore=sum(RWR_Seed[ScoreGeneID])*length(unique(SeedGeneID)) 
    if (Plot){ 
      RWRscore=data.frame(TopN=TopNN,geneScore=RWR_Score) 
      #RWRscoreMelt=melt(RWRscore,id='TopN') 
      colnames(RWRscore)=c('TopN','Score') 
      ggplot(RWRscore)+aes(x=TopN,y=Score)+geom_point()+geom_line()+xlab('TopN Genes')+ylab('RWR Score') 
      ggsave('RWRscore.png') 
    } 
 
  }else{ 
    RWR_Seed=randomWalk(ppiBinaryNet,SeedGeneID,EdgeWeight=EdgeWeight,gamma=gamma)*scoreStand 
    #TopNN=TopN 
    RWR_Score=numeric() 
    for (i in 1:length(TopN)){ 
      tempN=TopN[i] 
      RWR_Score[i]=sum(RWR_Seed[ScoreGeneID[1:tempN]],na.rm=T)/tempN 
    } 
    ####### 
    TopNNauc=trapz(TopN,RWR_Score) 
    GeneScore=RWR_Seed[ScoreGeneID] 
    TotalScore=sum(RWR_Seed[ScoreGeneID])*length(unique(SeedGeneID)) 
    if (Plot){ 
      RWRscore=data.frame(TopN=TopN,geneScore=RWR_Score) 
      #RWRscoreMelt=melt(RWRscore,id='TopN') 
      colnames(RWRscore)=c('TopN','Score') 
      ggplot(RWRscore)+aes(x=TopN,y=Score)+geom_point()+geom_line()+xlab('TopN Genes')+ylab('RWR Score') 
      ggsave('RWRscore.png') 
    } 
  } 
  ################### 
  tempRandomScore=numeric() 
  #tempRandomScore2=numeric() 
  #RWR_ScoregeneID=randomWalk(ppiBinaryNet,ScoreGeneID,EdgeWeight=EdgeWeight,gamma=gamma)*scoreStand 
  #TotalScore_ScoreID=sum(RWR_ScoregeneID[SeedGeneID])*length(unique(ScoreGeneID)) 
  seedNum=123456 
  for (j in 1:randomtimes){ 
    set.seed(seedNum+j) 
    tempRandomSeed1=sample(V(ppiBinaryNet)$name,length(ScoreGeneID)) 
    #tempRandomSeed2=sample(V(ppiBinaryNet)$name,length(SeedGeneID)) 
    tempRandomScore[j]=sum(RWR_Seed[tempRandomSeed1])*length(unique(SeedGeneID)) 
    #tempRandomScore2[i]=sum(RWR_ScoregeneID[tempRandomSeed2])*length(unique(ScoreGeneID)) 
  } 
  #tempRandomScore2[is.na(tempRandomScore2)]=0 
  tempRandomScore[is.na(tempRandomScore)]=0 
  if (testM==1){ 
  Pvalue1=wilcox.test(tempRandomScore,mu=TotalScore,alternative='less')$p.value 
  }else{ 
  Pvalue1=sum(tempRandomScore>TotalScore)/randomtimes   
  } 
  #Pvalue1=sum(tempRandomScore>TotalScore)/randomtimes   
  #Pvalue2=sum(tempRandomScore2>TotalScore_ScoreID)/randomtimes 
  #Pvalue2=wilcox.test(tempRandomScore2,mu=TotalScore,alternative='less')$p.value 
  #Pvalue=max(Pvalue1,Pvalue2) 
  Pvalue=Pvalue1 
  rm('ppiBinaryNet') 
  gc()   
  return(list(SeedGeneID=SeedGeneID,NumSeedGeneIDPPI=length(SeedGeneID),ScoreGeneID=ScoreGeneID,NumScoreGeneIDPPI=length(ScoreGeneID),GeneScore=GeneScore,TopNauc=TopNNauc,TotalScore=TotalScore,Pvalue=Pvalue,RandomScoreMedian=paste(mean(tempRandomScore),' [',quantile(tempRandomScore,0.025),',',quantile(tempRandomScore,0.975),']',sep=''))) 
} 
 
