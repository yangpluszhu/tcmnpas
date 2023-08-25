ExportNetwithin = function(geneID,geneIDName='geneList',diseaseGeneID,diseaseName='Disease',expandMethod=c('Core','EP','EN'),cutoffValue=2){ 
  #return(list(overlapGene=EXoverlapNodes,overlapGeneCount=length(EXoverlapNodes),ExGene=V(ExNet)$name,Exdisease=V(ExDiseaseNet)$name,OverlapResult=OverlapResult,NetStat=NetStat)) 
  require(data.table) 
  require(plyr) 
  require(stringr) 
  require(igraph) 
  #require(NetPathMiner) 
  if (sum(str_detect(geneID,',|;'))){ 
    geneID=as.character(unlist(str_split(geneID,',|;'))) 
    geneID=str_replace_all(geneID,'"','') 
    geneID=str_replace_all(geneID,"'",'') 
    geneID=geneID[geneID!=''] 
    geneID=unique(str_trim(geneID,side='both')) 
  } 
  ####### 
  if (is.null(geneIDName))geneIDName='geneList' 
  if (is.null(diseaseName))diseaseName='Disease' 
  if ('Core'%in%expandMethod)method=0 
  if ('EP'%in%expandMethod)method=1 
  if ('EN'%in%expandMethod)method=2 
  #### 
  if (!is.null(diseaseGeneID)){ 
    if (sum(str_detect(diseaseGeneID,',|;'))){ 
      diseaseGeneID=as.character(unlist(str_split(diseaseGeneID,',|;'))) 
      diseaseGeneID=str_replace_all(diseaseGeneID,'"','') 
      diseaseGeneID=str_replace_all(diseaseGeneID,"'",'') 
      diseaseGeneID=diseaseGeneID[diseaseGeneID!=''] 
      diseaseGeneID=unique(str_trim(diseaseGeneID,side='both')) 
    } 
    #### 
    CoreOverlap=intersect(geneID,diseaseGeneID) 
    ExNet=ExpandNet(geneID=geneID,method=method,cutoffValue=cutoffValue) 
    ExDiseaseNet=ExpandNet(geneID=diseaseGeneID,method=method,cutoffValue=cutoffValue) 
    NetAll=ExNet+ExDiseaseNet 
    overlapNodes=intersect(V(ExNet)$name,V(ExDiseaseNet)$name) 
    EXoverlapNodes=setdiff(intersect(V(ExNet)$name,V(ExDiseaseNet)$name),CoreOverlap) 
    NonlapDiseaseNodes=setdiff(V(ExDiseaseNet)$name,c(EXoverlapNodes,diseaseGeneID)) 
    NongeneIDNods=setdiff(V(ExNet)$name,c(EXoverlapNodes,geneID)) 
    netGeneID=V(NetAll)$name 
    netGeneName=geneTrans(netGeneID,type='ID') 
    V(NetAll)$symbol=netGeneName 
    V(NetAll)$type<-'NA' 
    V(NetAll)$type[V(NetAll)$name%in%geneID]<-geneIDName 
    V(NetAll)$type[V(NetAll)$name%in%diseaseGeneID]<-diseaseName 
    V(NetAll)$type[V(NetAll)$name%in%CoreOverlap]<-'CoreOverlap' 
    V(NetAll)$type[V(NetAll)$name%in%EXoverlapNodes]<-'ExpandOverlap' 
    V(NetAll)$type[V(NetAll)$name%in%NonlapDiseaseNodes]<-paste('Expand',diseaseName,sep='_') 
    V(NetAll)$type[V(NetAll)$name%in%NongeneIDNods]<-paste('Expand',geneIDName,sep='_') 
    #V(NetAll)$degree=igraph::degree(NetAll) 
    return(NetAll) 
  }else{ 
    ################### 
    ExNet=ExpandNet(geneID=geneID,method=method,cutoffValue=cutoffValue) 
    NetAll=ExNet 
    netGeneID=V(NetAll)$name 
    netGeneName=geneTrans(netGeneID,type='ID') 
    V(NetAll)$symbol=netGeneName 
    #V(NetAll)$degree=igraph::degree(NetAll) 
    return(NetAll) 
  } 
} 
 
