ExpandNet = function(geneID,method=1,cutoffValue=2,ppiBinaryNet=ppiBinaryNet){ 
  require(data.table) 
  require(plyr) 
  require(stringr) 
  require(igraph) 
  require("org.Hs.eg.db") 
  if (!exists('ppiBinaryNet')){ 
    load('db/ppiNetData.db') 
  } 
  geneID=as.character(geneID) 
  geneID=geneID[geneID%in%V(ppiBinaryNet)$name] 
  if (method==1){ 
    BXXXT_Dis=distances(ppiBinaryNet,geneID,geneID) 
    BXXXT_Dis[BXXXT_Dis>cutoffValue]<-0 
    Index_dis=which(BXXXT_Dis!=0,arr.ind=T) 
    kpath_BXXXT_Dis=list() 
    length(kpath_BXXXT_Dis)=nrow(Index_dis) 
    for (i in 1:nrow(Index_dis)){ 
      templocation1=Index_dis[i,1] 
      templocation2=Index_dis[i,2] 
      tempPath=shortest_paths(ppiBinaryNet,rownames(BXXXT_Dis)[templocation1],colnames(BXXXT_Dis)[templocation2],output='vpath') 
      kpath_BXXXT_Dis[[i]]=as_ids(tempPath$vpath[[1]]) 
    } 
    BXXXT_kpath_nodes=unique(unlist(kpath_BXXXT_Dis)) 
    BXXXT_kpath_Net=induced_subgraph(ppiBinaryNet,vids=BXXXT_kpath_nodes) 
    BXXXT_kpath_Net=igraph::simplify(BXXXT_kpath_Net) 
    return(BXXXT_kpath_Net) 
  }else if (method==2){ 
    BXXXT_Neighbor=neighborhood(ppiBinaryNet,1,nodes=geneID[1]) 
    for (i in 2:length(geneID)){ 
      tempID=geneID[i] 
      temp=neighborhood(ppiBinaryNet,1,nodes=tempID) 
      if (i==2){ 
        temp2=c(BXXXT_Neighbor[[1]],temp[[1]]) 
      }else{ 
        temp2=c(BXXXT_Neighbor,temp[[1]]) 
      } 
      BXXXT_Neighbor=unique(temp2) 
    } 
    #### 
    BXXXT_Net=induced_subgraph(ppiBinaryNet,vids=BXXXT_Neighbor) 
    return(BXXXT_Net) 
  }else{ 
    BXXXT_CoreNet=induced_subgraph(ppiBinaryNet,vids=geneID) 
    return(BXXXT_CoreNet) 
  } 
} 
 
