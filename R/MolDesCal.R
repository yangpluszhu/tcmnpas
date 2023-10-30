MolDesCal = function(SMI){
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
  return(FP)
}

