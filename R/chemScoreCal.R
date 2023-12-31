##' Compound Score 
##' @title chemScoreCal 
##' @param chem_target data.frame-->colnames(chem_target)=c('id','geneID') 
##' @param geneScore gene Score 
##' @param selectGeneId Gene Id 
##' @param method 1: Use the cutoffValue; 2: SelectGeneId with a minimum of 100% coverage. 
##' @param cutoffValue cutoff Value 
##' @return a list Object 
##' @export 
##' @author Yang Ming 
chemScoreCal = function(chem_target,geneScore,selectGeneId,method=1,cutoffValue=0){ 
  require(data.table) 
  require(plyr) 
  require(stringr) 
  options(stringsAsFactors = F) 
  chem_target=as.data.frame(chem_target) 
  chem_target$geneID=as.character(chem_target$geneID) 
  chemID=unique(chem_target$id) 
  geneID=unique(chem_target$geneID) 
  chemScore=numeric() 
  length(chemScore)=length(chemID) 
  names(chemScore)=chemID 
  IdentitySelectGene=rep(0,length=length(geneID)) 
  names(IdentitySelectGene)=geneID 
  IdentitySelectGene[selectGeneId]=1 
  for (i in chemID){ 
    tempD=chem_target[chem_target$id==i,] 
    tempGene=as.character(tempD$geneID) 
    if (sum(tempGene%in%selectGeneId)>0){ 
      tempGene=tempGene[tempGene%in%selectGeneId] 
      chemScore[i]=mean(geneScore[tempGene]*IdentitySelectGene[tempGene]) 
    }else{ 
      chemScore[i]=0 
    } 
  } 
  chemScore=sort(chemScore,decreasing = T) 
  if (method==1){ 
    selectChemId=names(chemScore)[chemScore>cutoffValue] 
    return(list(chemScore=chemScore,selectChemId=selectChemId)) 
  }else{ 
    n=1 
    percent=0 
    while (percent<cutoffValue){ 
      percent=sum(unique(chem_target$geneID[chem_target$id%in%names(chemScore)[1:n]])%in%selectGeneId)/length(selectGeneId) 
      n=n+1 
    } 
    selectChemId=names(chemScore)[1:n] 
    return(list(chemScore=chemScore,selectChemId=selectChemId)) 
  } 
} 
 
