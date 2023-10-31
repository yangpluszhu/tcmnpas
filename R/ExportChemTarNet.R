ExportChemTarNet = function(chemTar,geneListName='geneList',diseaseGeneID=NULL,diseaseName='Disease',expandMethod=c('Core','EP','EN'),cutoffValue=2,IncludeHerb=T,IncludeChem=T,IncludePPI=T,Savefilename='EpChemTarNet'){ 
  ##chemTar:colnames:herb,cid,chemical_name,geneID 
  require(data.table) 
  require(plyr) 
  require(igraph) 
  if (!dir.exists('Report')) dir.create('Report') 
  chemTar=as.data.table(chemTar) 
  if (!'geneName'%in%colnames(chemTar))chemTar=transform(chemTar,geneName=geneTrans(geneID,type='ID')) 
  if (!'chemical_name'%in%colnames(chemTar))chemTar=transform(chemTar,chemical_name=cid) 
  if (!'herb'%in%colnames(chemTar))chemTar=transform(chemTar,herb='drug') 
  chemTar=chemTar[!is.na(geneName),] 
  chemTar=unique(chemTar,by=c('herb','cid','geneID')) 
  #################chemHerbNet 
  chemHerb=chemTar[,.N,by=c('herb','cid')] 
  chemHerb=rename(chemHerb,c('N'='weight')) 
  chemHerbNodeName1=unique(data.table(name=unique(chemHerb$herb),symbol=unique(chemHerb$herb),class='herb')) 
  chemHerbNodeName2=unique(chemTar[,.(name=cid,symbol=chemical_name,class='chemical')],by='name') 
  chemHerbNodeName=rbind(chemHerbNodeName1,chemHerbNodeName2) 
  chemHerbNet=graph_from_data_frame(chemHerb,directed = F,vertices =chemHerbNodeName ) 
  ###############chemGeneNet 
  chemGene=chemTar[,.N,by=c('cid','geneID')] 
  chemGene=rename(chemGene,c('N'='weight')) 
  chemGeneNodeName1=unique(chemTar[,.(name=cid,symbol=chemical_name,class='chemical')],by='name') 
  chemGeneNodeName2=unique(chemTar[,.(name=geneID,symbol=geneName,class='gene')],by='name') 
  chemGeneNodeName=rbind(chemGeneNodeName1,chemGeneNodeName2) 
  chemGeneNet=graph_from_data_frame(chemGene,directed = F,vertices =chemGeneNodeName ) 
  ###############herbGeneNet### 
  herbGene=chemTar[,.N,by=c('herb','geneID')] 
  herbGene=rename(herbGene,c('N'='weight')) 
  herbGeneNodeName1=unique(data.table(name=unique(herbGene$herb),symbol=unique(herbGene$herb),class='herb')) 
  herbGeneNodeName2=unique(chemTar[,.(name=geneID,symbol=geneName,class='gene')],by='name') 
  herbGeneNodeName=rbind(herbGeneNodeName1,herbGeneNodeName2) 
  herbGeneNet=graph_from_data_frame(herbGene,directed = F,vertices =herbGeneNodeName ) 
  #############PPI 
  #if ('Core'%in%expandMethod)method=0 
  #if ('EP'%in%expandMethod)method=1 
  #if ('EN'%in%expandMethod)method=2 
  if (is.null(Savefilename))Savefilename='EpChemTarNet' 
  ExNet=ExportNetwithin(geneID=unique(chemTar$geneID),geneIDName=geneListName,diseaseGeneID=diseaseGeneID,diseaseName=diseaseName,expandMethod=expandMethod,cutoffValue=cutoffValue) 
  E(ExNet)$weight<-1 
  if (IncludeHerb&IncludeChem&IncludePPI){ 
    NetALL=combineIgraphOB(chemHerbNet,chemGeneNet) 
    NetALL=combineIgraphOB(NetALL,ExNet) 
  } 
  if (!IncludeHerb&IncludeChem&IncludePPI)NetALL=combineIgraphOB(chemGeneNet,ExNet) 
  if (IncludeHerb&!IncludeChem&IncludePPI)NetALL=combineIgraphOB(herbGeneNet,ExNet) 
  if (IncludeHerb&IncludeChem&!IncludePPI)NetALL=combineIgraphOB(chemHerbNet,chemGeneNet) 
  if (!IncludeHerb&!IncludeChem&IncludePPI)NetALL=ExNet 
  if (!IncludeHerb&IncludeChem&!IncludePPI)NetALL=chemGeneNet 
  if (IncludeHerb&!IncludeChem&!IncludePPI)NetALL=herbGeneNet 
  if (!IncludeHerb&!IncludeChem&!IncludePPI)NetALL=make_empty_graph() 
  file=paste('Report/',Savefilename,'.gml',sep='') 
  write_graph(NetALL, file=file, format='gml') 
  return(NetALL) 
} 
 
