NetStatCal = function(NetList){
  #NetList:list-->names(list)=NetName(igraphObject)
  library(data.table)
  library(plyr)
  library(stringr)
  library(igraph)
  library("org.Hs.eg.db")
  metric=c('Mean_degree','Density','Mean_betweenness','Count_nodes','Cluster_coefficient','Diameter','Shortest_path')
  #NetList=list(PPI=ppiBinaryNet,BXXXT=BXXXT_Net,Colitis=Yan_Net,Diabet=Diabet_Net,Cancer=Cancer_Net)
  statistics_matrix=matrix(0,nrow=length(metric),ncol=length(NetList))
  colnames(statistics_matrix)=names(NetList)
  rownames(statistics_matrix)=metric
  ###
  statistics_matrix['Mean_degree',]=sapply(NetList,function(x)paste(round(mean(igraph::degree(x)),2),round(sd(igraph::degree(x)),2),sep= '±'))
  statistics_matrix['Mean_betweenness',]=sapply(NetList,function(x)paste(round(mean(igraph::betweenness(x,directed=F)),2),round(sd(igraph::betweenness(x,directed=F)),2),sep= '±'))
  statistics_matrix['Density',]=sapply(NetList,igraph::edge_density)
  statistics_matrix['Count_nodes',]=sapply(NetList,igraph::vcount)
  statistics_matrix['Cluster_coefficient',]=sapply(NetList,igraph::transitivity,type='global')
  statistics_matrix['Diameter',]=sapply(NetList,igraph::diameter,directed=F)
  statistics_matrix['Shortest_path',]=sapply(NetList,igraph::mean_distance,directed=F)
  statistics_matrix=as.data.frame(statistics_matrix)
  rownames(statistics_matrix)=metric
  return(statistics_matrix)
}

