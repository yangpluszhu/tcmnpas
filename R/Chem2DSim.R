##' Calculate Compound 2D Similarity
##'
##' @title Chem2DSim
##' @param SMI molecular file: a dataframe with smiles
##' @return a list object
##' @export
##' @author Yang Ming
Chem2DSim = function(SMI){ 
  require(data.table) 
  require(plyr) 
  require(stringr) 
  require(igraph) 
  require(ggplot2) 
  require(ChemmineR) 
  require(PaDEL) 
  require(fingerprint) 
  if (!dir.exists('MOL')){ 
    dir.create('MOL') 
  } 
  SMI=as.data.table(SMI) 
  SMI2=SMI[smiles!='',] 
  smiles=SMI2$smiles 
  names(smiles)=SMI2$id 
  SMISdfSets=ChemmineR::smiles2sdf(smiles) 
  SMISdfSets=SMISdfSets[validSDF(SMISdfSets)] 
  #sum(validSDF(SMISdfSets))==length(SMISdfSets) 
  write.SDF(SMISdfSets,'MOL/SMISdfSets.sdf',cid=T) 
  FP=CalDescriptors('MOL/SMISdfSets.sdf') 
  FPv=setdiff(colnames(FP),'Name') 
  FPlist=list() 
  CID=FP$Name 
  for (i in 1:nrow(FP)){ 
    temp=FP[i,FPv,with=F] 
    tempObject=new("fingerprint", nbit=length(FPv), bits=which(temp!=0)) 
    FPlist=c(FPlist,tempObject) 
  } 
  length(FPlist)==length(CID) 
  ChemSimMatrix=fp.sim.matrix(FPlist,method='tanimoto') 
  rownames(ChemSimMatrix)=CID 
  colnames(ChemSimMatrix)=CID 
  return(list(ChemSimMatrix=ChemSimMatrix,FP=FP)) 
} 
 
