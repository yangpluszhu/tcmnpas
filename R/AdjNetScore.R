AdjNetScore = function(TotalScore,NumSeedGene,NumScoreGene,multiple=1000){
#ScoreGeneData$AdjTotalScore=KATZresult$TotalScore/KATZresult$NumSeedGeneID*log10(1+KATZresult$NumSeedGeneID/KATZresult$NumScoreGeneID/2)*1000
Score=TotalScore/NumSeedGene*log10(1+NumSeedGene/NumScoreGene/2)*multiple
return(Score)
}

