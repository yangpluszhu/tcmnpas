SmiToSDF = function(SMI,type='C'){ 
  #SMI='SMI.csv'-->colnames=c('id','smiles') or SMItext-->Named character vector 
  #type='F'-->file for SMI or 'C'-->character for SMItext 
  #outPut:'Molelular.sdf' 
  require(data.table) 
  require(plyr) 
  require(stringr) 
  require(clipr) 
  require(igraph) 
  require(readr) 
  require(ggplot2) 
  require(ChemmineR) 
  if (type=='F'){ 
    SMI1=fread(SMI) 
    SMI2=SMI1[smiles!='',] 
    smiles=SMI2$smiles 
    names(smiles)=SMI2$id 
    SMISdfSets=ChemmineR::smiles2sdf(smiles) 
    SMISdfSets=SMISdfSets[validSDF(SMISdfSets)] 
    #sum(validSDF(SMISdfSets))==length(SMISdfSets) 
    datablock(SMISdfSets)<-data.frame(id=SMISdfSets@ID) 
    write.SDF(SMISdfSets,'Molelular.sdf',cid=T) 
  }else if (type=='C'){ 
    SMISdfSets=ChemmineR::smiles2sdf(SMI) 
    SMISdfSets=SMISdfSets[validSDF(SMISdfSets)] 
    #sum(validSDF(SMISdfSets))==length(SMISdfSets) 
    datablock(SMISdfSets)<-data.frame(id=SMISdfSets@ID) 
    write.SDF(SMISdfSets,'Molelular.sdf',cid=T) 
  }else{ 
    print('ERROR:type not found!') 
  } 
} 
 
