makeHerb_ChemNet = function(NetworkData,NetworkType='HerbChem',INppi=T){
  library(igraph)
  library(rio)
  library(data.table)
  #library(visNetwork)
  ###NetworkData-->from TCMNPAS outPut[colnames:"herb","cid","chemical_name","geneID","score","inchikey", "smiles","database","QED_DES","Lipinski_Vios","Veber_Vios"]
  ###NetworkType:'HerbChem','ChemGene',or 'HerbChemGene'
  ###list(outputNetwork=outputNetwork,ChemCid_ChemName=ChemCid_ChemName):outputNetwork-->igraphObj;ChemCid_ChemName--dataFrame[colname->('cid','chemical_name')]
  options(stringsAsFactors = F)
  #dd=import('herb-chem-target.csv')
  dd=NetworkData
  dd$geneID=as.character(dd$geneID)
  ChemCid_ChemName=dd[,c('cid','chemical_name')]
  ChemCid_ChemName=as.data.table(ChemCid_ChemName)
  ChemCid_ChemName=unique(ChemCid_ChemName,by='cid')
  ChemCid_ChemName=as.data.frame(ChemCid_ChemName)
  ###
  Herb_Chem=dd[,c('herb','cid')]
  Herb_Chem=unique(Herb_Chem)
  Nodes_Herb_Chem=c(Herb_Chem$herb,Herb_Chem$cid)
  Nodes_Herb_Chem_Group=c(rep('Herb',length(Herb_Chem$herb)),rep('Chemical',length(Herb_Chem$cid)))
  Nodes_Herb_ChemData=data.frame(name=Nodes_Herb_Chem,group=Nodes_Herb_Chem_Group)
  Nodes_Herb_ChemData=unique(Nodes_Herb_ChemData)
  Herb_ChemNet=graph_from_data_frame(d=Herb_Chem, directed = F, vertices = Nodes_Herb_ChemData)#####Herb_ChemNet
  ####################
  Chem_Gene=dd[,c('cid','geneID')]
  #Chem_Gene=as.data.table(Chem_Gene)
  Chem_Gene=unique(Chem_Gene)
  #Chem_Gene$geneID=as.character(Chem_Gene$geneID)
  Nodes_Chem_Gene=c(Chem_Gene$cid,Chem_Gene$geneID)
  Nodes_Chem_Gene_Group=c(rep('Chemical',length(Chem_Gene$cid)),rep('Gene',length(Chem_Gene$geneID)))
  Nodes_Chem_GeneData=data.frame(name=Nodes_Chem_Gene,group=Nodes_Chem_Gene_Group)
  Nodes_Chem_GeneData=unique(Nodes_Chem_GeneData)
  Chem_GeneNet=graph_from_data_frame(d=Chem_Gene, directed = F, vertices = Nodes_Chem_GeneData)####Chem_GeneNet
  #####################
  #Herb_Chem_GeneNet=combineIgraphOB(Herb_ChemNet,Chem_GeneNet)
  Herb_Chem_GeneNet=Herb_ChemNet+Chem_GeneNet
  g1_group=V(Herb_Chem_GeneNet)$group_1
  g2_group=V(Herb_Chem_GeneNet)$group_2
  g_group=g1_group
  g_group[is.na(g1_group)]<-g2_group[is.na(g1_group)]
  V(Herb_Chem_GeneNet)$group<-g_group######Herb_Chem_GeneNet
  Herb_Chem_GeneNet=delete_vertex_attr(Herb_Chem_GeneNet,'group_1')
  Herb_Chem_GeneNet=delete_vertex_attr(Herb_Chem_GeneNet,'group_2')
  if (INppi){
    Chem_GeneNetCore=ExpandNet(geneID=unique(Chem_Gene$geneID),method=0)
    Chem_GeneNetCombinePPI=Chem_GeneNet+Chem_GeneNetCore
    Herb_Chem_GeneNetCombinePPI=Herb_Chem_GeneNet+Chem_GeneNetCore
    Herb_Chem_GeneNet=Herb_Chem_GeneNetCombinePPI###
    Chem_GeneNet=Chem_GeneNetCombinePPI#####
  }
  if (NetworkType=='HerbChem'){
    outputNetwork=Herb_ChemNet
  }else if (NetworkType=='ChemGene'){
    outputNetwork=Chem_GeneNet
  }else{
    outputNetwork=Herb_Chem_GeneNet
  }
  return(list(outputNetwork=outputNetwork,ChemCid_ChemName=ChemCid_ChemName))
}

