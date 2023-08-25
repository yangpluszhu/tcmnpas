ChemTargetNetScore = function(ChemData,diseaseGeneID,Methodnet='KATZ',IFcombine=F){ 
  require(data.table) 
  require(stringr) 
  require(RSQLite) 
  require(plyr) 
  require(ggplot2) 
  require(pracma) 
  options(stringsAsFactors = F) 
  ChemData$geneID=as.character(ChemData$geneID) 
  if (IFcombine){ 
    queryGeneID=unique(ChemData$geneID) 
    NetCorrResult=NetCorr(SeedGeneID=diseaseGeneID,ScoreGeneID=queryGeneID,Methodnet=Methodnet,randomtimes=1,gamma=0.7,testM=2) 
    ScoreGeneData=data.frame(SeedGeneIDPPI=0,NumSeedGeneID=0,ScoreGeneIDPPI=0,NumScoreGeneID=0,overlapScore=0,Path1Score=0,Path2Score=0,Path3Score=0,TotalScore=0,AdjTotalScore=0) 
    if (Methodnet=='KATZ'){ 
      KATZresult=NetCorrResult$resultKATZ 
      #KATZresult=data.frame(SeedGeneIDPPI=toString(SeedGeneID),NumSeedGeneID=length(SeedGeneID),ScoreGeneIDPPI=toString(ScoreGeneID),NumScoreGeneID=length(ScoreGeneID),overlapScore=0,Path1Score=0,Path2Score=0,Path3Score=0,TotalScore=0,RandomScoreMedian=0,Pvalue=1) 
      ScoreGeneData$SeedGeneIDPPI=KATZresult$SeedGeneIDPPI 
      ScoreGeneData$NumSeedGeneID=KATZresult$NumSeedGeneID 
      ScoreGeneData$ScoreGeneIDPPI=KATZresult$ScoreGeneIDPPI 
      ScoreGeneData$NumScoreGeneID=KATZresult$NumScoreGeneID 
      ScoreGeneData$overlapScore=KATZresult$overlapScore 
      ScoreGeneData$Path1Score=KATZresult$Path1Score 
      ScoreGeneData$Path2Score=KATZresult$Path2Score 
      ScoreGeneData$Path3Score=KATZresult$Path3Score 
      ScoreGeneData$TotalScore=KATZresult$TotalScore 
      ##score/n0*exp(score/n1-1)*1000 
      #ScoreGeneData$AdjTotalScore=KATZresult$TotalScore/KATZresult$NumSeedGeneID*exp(KATZresult$TotalScore/KATZresult$NumScoreGeneID-1)*1000 
      #score/n0*log10(1+n0/n1/2) 
      ScoreGeneData$AdjTotalScore=KATZresult$TotalScore/KATZresult$NumSeedGeneID*log10(1+KATZresult$NumSeedGeneID/KATZresult$NumScoreGeneID/2)*1000 
    }else if (Methodnet=='RW'){ 
      RWresult=NetCorrResult$resultRWTable 
      #resultRWTable:data.frame(SeedGeneIDPPI=toString(resultRW$SeedGeneID),NumSeedGeneID=resultRW$NumSeedGeneIDPPI,ScoreGeneIDPPI=toString(resultRW$ScoreGeneID),NumScoreGeneID=resultRW$NumScoreGeneIDPPI,TotalScore=resultRW$TotalScore,RandomScoreMean=resultRW$RandomScoreMedian,Pvalue=resultRW$Pvalue) 
      ScoreGeneData$SeedGeneIDPPI=RWresult$SeedGeneIDPPI 
      ScoreGeneData$NumSeedGeneID=RWresult$NumSeedGeneID 
      ScoreGeneData$ScoreGeneIDPPI=RWresult$ScoreGeneIDPPI 
      ScoreGeneData$NumScoreGeneID=RWresult$NumScoreGeneID 
      ScoreGeneData$TotalScore=RWresult$TotalScore 
      ##score/n0*exp(score/n1-1)*1000 
      #ScoreGeneData$AdjTotalScore=RWresult$TotalScore/RWresult$NumSeedGeneID*exp(RWresult$TotalScore/RWresult$NumScoreGeneID-1) 
      #score/n0*log10(1+n0/n1/2) 
      ScoreGeneData$AdjTotalScore=RWresult$TotalScore/RWresult$NumSeedGeneID*log10(1+RWresult$NumSeedGeneID/RWresult$NumScoreGeneID/2) 
    } 
  }else{ 
    queryGeneIDdata=ChemData 
    queryCID=unique(queryGeneIDdata$cid) 
    ScoreGeneData=data.frame() 
    for (i in 1:length(queryCID)){ 
      tempqueryGeneID=queryGeneIDdata$geneID[queryGeneIDdata$cid==queryCID[i]] 
      tempNetCorrResult=NetCorr(SeedGeneID=diseaseGeneID,ScoreGeneID=tempqueryGeneID,Methodnet=Methodnet,randomtimes=1,gamma=0.7,testM=1) 
      tempScoreGeneData=data.frame(cid=0,SeedGeneIDPPI=0,NumSeedGeneID=0,ScoreGeneIDPPI=0,NumScoreGeneID=0,overlapScore=0,Path1Score=0,Path2Score=0,Path3Score=0,TotalScore=0,AdjTotalScore=0) 
      if (Methodnet=='KATZ'){ 
        tempKATZresult=tempNetCorrResult$resultKATZ 
        tempScoreGeneData$cid=queryCID[i] 
        #tempScoreGeneData$chemical_name=queryGeneIDdata$chemical_name[queryGeneIDdata$cid==queryCID[i]][1] 
        #tempScoreGeneData$inchikey=queryGeneIDdata$inchikey[queryGeneIDdata$cid==queryCID[i]][1] 
        tempScoreGeneData$SeedGeneIDPPI=tempKATZresult$SeedGeneIDPPI 
        tempScoreGeneData$NumSeedGeneID=tempKATZresult$NumSeedGeneID 
        tempScoreGeneData$ScoreGeneIDPPI=tempKATZresult$ScoreGeneIDPPI 
        tempScoreGeneData$NumScoreGeneID=tempKATZresult$NumScoreGeneID 
        tempScoreGeneData$overlapScore=tempKATZresult$overlapScore 
        tempScoreGeneData$Path1Score=tempKATZresult$Path1Score 
        tempScoreGeneData$Path2Score=tempKATZresult$Path2Score 
        tempScoreGeneData$Path3Score=tempKATZresult$Path3Score 
        tempScoreGeneData$TotalScore=tempKATZresult$TotalScore 
        ##score/n0*exp(score/n1-1)*1000 
        ###score/n0*log10(1+n0/n1/2) 
        #tempScoreGeneData$AdjTotalScore=tempKATZresult$TotalScore/tempKATZresult$NumSeedGeneID*exp(tempKATZresult$TotalScore/tempKATZresult$NumScoreGeneID-1)*1000 
        tempScoreGeneData$AdjTotalScore=tempKATZresult$TotalScore/tempKATZresult$NumSeedGeneID*log10(1+tempKATZresult$NumSeedGeneID/tempKATZresult$NumScoreGeneID/2)*1000 
        ScoreGeneData=rbind(ScoreGeneData,tempScoreGeneData) 
      }else if (Methodnet=='RW'){ 
        tempKATZresult=tempNetCorrResult$resultRWTable 
        tempScoreGeneData$cid=queryCID[i] 
        #tempScoreGeneData$chemical_name=queryGeneIDdata$chemical_name[queryGeneIDdata$cid==queryCID[i]][1] 
        #tempScoreGeneData$inchikey=queryGeneIDdata$inchikey[queryGeneIDdata$cid==queryCID[i]][1] 
        tempScoreGeneData$SeedGeneIDPPI=tempKATZresult$SeedGeneIDPPI 
        tempScoreGeneData$NumSeedGeneID=tempKATZresult$NumSeedGeneID 
        tempScoreGeneData$ScoreGeneIDPPI=tempKATZresult$ScoreGeneIDPPI 
        tempScoreGeneData$NumScoreGeneID=tempKATZresult$NumScoreGeneID 
        ##score/n0*exp(score/n1-1)*1000 
        ###score/n0*log10(1+n0/n1/2) 
        #tempScoreGeneData$AdjTotalScore=tempKATZresult$TotalScore/tempKATZresult$NumSeedGeneID*exp(tempKATZresult$TotalScore/tempKATZresult$NumScoreGeneID-1)*1000 
        tempScoreGeneData$AdjTotalScore=tempKATZresult$TotalScore/tempKATZresult$NumSeedGeneID*log10(1+tempKATZresult$NumSeedGeneID/tempKATZresult$NumScoreGeneID/2) 
        ScoreGeneData=rbind(ScoreGeneData,tempScoreGeneData) 
      } 
    } 
  } 
  return(ScoreGeneData) 
} 
 
