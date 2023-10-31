keggSimScore = function(chem_target,keggPath){ 
  #chem_target:data.frame-->colnames(chem_target)=c('id','geneID') 
  require(data.table) 
  require(plyr) 
  require(stringr) 
  require(clipr) 
  require(igraph) 
  require(ggplot2) 
  require(ChemmineR) 
  require(pracma) 
  #chem_target$geneID=as.character(chem_target$geneID) 
  if (!exists('keggPath')){ 
    #keggPath=fread('db/KEGG_data.csv') 
    keggPath=fread('db/KEGG_data.db',encoding='UTF-8') 
  } 
  #keggPath=fread('KEGG_data160517.csv') 
  keggPath$geneID=as.character(keggPath$geneID) 
  DrugId=unique(chem_target$id) 
  chem_target$geneID=as.character(chem_target$geneID) 
  keggSimMatrix=matrix(0,nrow=length(DrugId),ncol=length(DrugId)) 
  colnames(keggSimMatrix)=DrugId 
  rownames(keggSimMatrix)=DrugId 
  for (i in 1:length(DrugId)){ 
    print(i) 
    j=i+1 
    while (j<=length(DrugId)){ 
      cid1=DrugId[i] 
      cid1Genetemp=chem_target[id==cid1,geneID] 
      keggid1=keggPath[geneID%in%cid1Genetemp,unique(PathwayID)] 
      cid2=DrugId[j] 
      cid2Genetemp=chem_target[id==cid2,geneID] 
      keggid2=keggPath[geneID%in%cid2Genetemp,unique(PathwayID)] 
      if (length(keggid1)==0|length(keggid2)==0){ 
        keggSimMatrix[cid1,cid2]=0 
      }else{ 
        keggSimMatrix[cid1,cid2]=length(intersect(keggid1,keggid2))/length(union(keggid1,keggid2)) 
      } 
      j=j+1 
    } 
  } 
  keggSimMatrix=asSymmetric(keggSimMatrix) 
  return(keggSimMatrix) 
} 
 
