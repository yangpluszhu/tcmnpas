Chem3DSim = function(SMI){
  #SMI:data.frame-->colnames(SMI)=c('id','smiles')
  library(data.table)
  library(plyr)
  library(stringr)
  library(clipr)
  library(igraph)
  library(ggplot2)
  library(ChemmineR)
  if (!dir.exists('MOL'))dir.create('MOL')
  SMI=as.data.table(SMI)
  SMI2=SMI[smiles!='',]
  smiles=SMI2$smiles
  names(smiles)=SMI2$id
  SMISdfSets=ChemmineR::smiles2sdf(smiles)
  SMISdfSets=SMISdfSets[validSDF(SMISdfSets)]
  #sum(validSDF(SMISdfSets))==length(SMISdfSets)
  datablock(SMISdfSets)<-data.frame(id=SMISdfSets@ID)
  write.SDF(SMISdfSets,'MOL/SMISdfSets.sdf',cid=T)
  ###
  setwd(paste(getwd(),'/MOL',sep=''))
  TransMol('SMISdfSets.sdf')
  DrugId=unique(SMISdfSets@ID)
  ChemShapeSimMatrix=matrix(0,nrow=length(DrugId),ncol=length(DrugId))
  colnames(ChemShapeSimMatrix)=DrugId
  rownames(ChemShapeSimMatrix)=DrugId
  ChemFeatureSimMatrix=matrix(0,nrow=length(DrugId),ncol=length(DrugId))
  colnames(ChemFeatureSimMatrix)=DrugId
  rownames(ChemFeatureSimMatrix)=DrugId
  ChemHybridSimMatrix=matrix(0,nrow=length(DrugId),ncol=length(DrugId))
  colnames(ChemHybridSimMatrix)=DrugId
  rownames(ChemHybridSimMatrix)=DrugId
  for (i in 1:length(DrugId)){
    cid1=DrugId[i]
    molFilename=paste(cid1,'mol2',sep='.')
    print(i)
    j=i+1
    while (j<=length(DrugId)){
      cid2=DrugId[j]
      molFilename2=paste(cid2,'mol2',sep='.')
      dos=paste('Cynthia -q',molFilename,'-t',molFilename2)
      shell(dos,intern=T)
      if (class(try(fread('Result.list'),silent = T))=='try-error'){
        ChemHybridSimMatrix[cid1,cid2]=NA
        ChemShapeSimMatrix[cid1,cid2]=NA
        ChemFeatureSimMatrix[cid1,cid2]=NA
      }else{
        tempResult=fread('Result.list')
        if (nrow(tempResult)==0){
          ChemHybridSimMatrix[cid1,cid2]=NA
          ChemShapeSimMatrix[cid1,cid2]=NA
          ChemFeatureSimMatrix[cid1,cid2]=NA
        }else{
          ChemHybridSimMatrix[cid1,cid2]=tempResult$HybridScore
          ChemShapeSimMatrix[cid1,cid2]=tempResult$ShapeScore
          ChemFeatureSimMatrix[cid1,cid2]=tempResult$FeatureScore
        }
      }
      j=j+1
    }
  }
  setwd(str_replace(getwd(),'/MOL',''))
  ChemHybridSimMatrix=asSymmetric(ChemHybridSimMatrix)
  ChemShapeSimMatrix=asSymmetric(ChemShapeSimMatrix)
  ChemFeatureSimMatrix=asSymmetric(ChemFeatureSimMatrix)
 return(list(ChemHybridSimMatrix=ChemHybridSimMatrix,ChemShapeSimMatrix=ChemShapeSimMatrix,ChemFeatureSimMatrix=ChemFeatureSimMatrix))
}

