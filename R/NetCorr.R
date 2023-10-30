NetCorr = function(SeedGeneID,ScoreGeneID,Methodnet,randomtimes=1000,gamma=0.7,testM=1){
  if (Methodnet=='RW'){
    resultRW=RWanalysis(SeedGeneID=SeedGeneID,ScoreGeneID=ScoreGeneID,randomtimes=randomtimes,gamma=gamma,testM=testM)
    resultRWTable=data.frame(SeedGeneIDPPI=toString(resultRW$SeedGeneID),NumSeedGeneID=resultRW$NumSeedGeneIDPPI,ScoreGeneIDPPI=toString(resultRW$ScoreGeneID),NumScoreGeneID=resultRW$NumScoreGeneIDPPI,TotalScore=resultRW$TotalScore,RandomScoreMean=resultRW$RandomScoreMedian,Pvalue=resultRW$Pvalue)
  }else{
    resultRW=NULL
    resultRWTable=NULL
  }
  #####
  if (Methodnet=='KATZ'){
    resultKATZ=KATZ(SeedGeneID,ScoreGeneID,randomtimes=randomtimes,testM=testM)
  }else{
    resultKATZ=NULL
  }
  #####
  if (Methodnet=='DIS'){
    resultDis=DisAccess(SeedGeneID,ScoreGeneID,randomtimes=randomtimes,testM=testM)
    resultDisTable=data.frame(SeedGeneIDPPI=toString(resultDis$SeedGeneID),NumSeedGeneID=resultDis$NumSeedGeneID,ScoreGeneIDPPI=toString(resultDis$ScoreGeneID),NumScoreGeneID=resultDis$NumScoreGeneID,MeanDistance=resultDis$MedianDistance,RandomScoreMean=resultDis$RandomScoreMedian,Pvalue=resultDis$Pvalue)
  }else{
    resultDis=NULL
    resultDisTable=NULL
  }
  #####
return(list(resultRW=resultRW,resultRWTable=resultRWTable,resultKATZ=resultKATZ,resultDis=resultDis,resultDisTable=resultDisTable))
}

