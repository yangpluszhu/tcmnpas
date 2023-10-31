##' Adjusted Network based relevance score
##' Adjusted based on the number of genes: Score=TotalScore/NumSeedGene*log10(1+NumSeedGene/NumScoreGene/2)*multiple
##'
##' @title AdjNetScore
##' @param TotalScore Total Score
##' @param NumSeedGene Number of Seed Gene
##' @param NumScoreGene Seed Gene Score
##' @param multiple Adjusted multiple constant
##' @return a numeric object
##' @export
##' @author Yang Ming
AdjNetScore = function(TotalScore,NumSeedGene,NumScoreGene,multiple=1000){ 
Score=TotalScore/NumSeedGene*log10(1+NumSeedGene/NumScoreGene/2)*multiple 
return(Score) 
} 
 
