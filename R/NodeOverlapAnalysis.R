NodeOverlapAnalysis = function(NetList,Plot=TRUE){
  #NetList:list-->names(list)=NetName(igraphObject)
  library(data.table)
  library(plyr)
  library(stringr)
  library(igraph)
  library(corrplot)
  library("org.Hs.eg.db")
  #NetList=list(BXXXT=BXXXT_kpath_Net,Colitis=Yan_kpath_Net,Diabet=Diabet_kpath_Net,Cancer=Cancer_kpath_Net)
  overlap_DISmatrix=matrix(0,nrow=length(NetList),ncol=length(NetList))
  colnames(overlap_DISmatrix)=names(NetList)
  rownames(overlap_DISmatrix)=names(NetList)
  for (i in names(NetList)){
    for (j in names(NetList)){
      tempNet_i=NetList[[i]]
      tempNet_j=NetList[[j]]
      overlap_DISmatrix[i,j]=length(intersect(V(tempNet_i)$name,V(tempNet_j)$name))/min(c(vcount(tempNet_i),vcount(tempNet_j)))
    }
  }
  diag(overlap_DISmatrix)=0
  colnames(overlap_DISmatrix)=names(NetList)
  rownames(overlap_DISmatrix)=names(NetList)
  if (Plot){
    png('NetNodeOverlap_matrix.png',width=1280,height=1024)
    corrplot(overlap_DISmatrix,method='pie',diag=F,cl.lim=c(0,1),tl.cex=3,cl.cex=3)
    dev.off()
  }
  return(overlap_DISmatrix)
}

