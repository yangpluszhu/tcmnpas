Chem2DSim = function(SMI){
  #SMI:data.frame-->colnames(SMI)=c('id','smiles')
  #calMolecularDes--->set padel_types.xml!
  library(data.table)
  library(plyr)
  library(stringr)
  library(igraph)
  library(ggplot2)
  library(ChemmineR)
  library(PaDEL)
  library(fingerprint)
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

