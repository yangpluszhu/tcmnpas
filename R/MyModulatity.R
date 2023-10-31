MyModulatity = function(g,LisMembership){ 
  ##g-->igraph object,Must have node name;LisMembership-->list,names(list):name of cluster 
  require(igraph) 
  require(clValid) 
  Ncluster=length(LisMembership) 
  mem=unique(unlist(LisMembership)) 
  if (length(mem)==vcount(g)&(sum(mem%in%V(g)$name)==length(mem))){ 
    membership_Mat=annotationListToMatrix(LisMembership,V(g)$name) 
    Co_membership=membership_Mat%*%t(membership_Mat) 
    adj=as_adj(g) 
    D=degree(g) 
    D1=replicate(length(D),as.matrix(D),simplify ='matrix' ) 
    D2=replicate(length(D),as.matrix(D),simplify ='matrix' ) 
    D2=t(D2) 
    DD=D1*D2 
    EC=ecount(g) 
    MQ=(adj-DD/(2*EC))*Co_membership 
    Q=sum(MQ)/(2*EC) 
    return(Q) 
  }else{ 
    g_sub=induced_subgraph(g,mem) 
    membership_Mat=annotationListToMatrix(LisMembership,V(g_sub)$name) 
    Co_membership=membership_Mat%*%t(membership_Mat) 
    adj=as_adj(g_sub) 
    D=degree(g_sub) 
    D1=replicate(length(D),as.matrix(D),simplify ='matrix' ) 
    D2=replicate(length(D),as.matrix(D),simplify ='matrix' ) 
    D2=t(D2) 
    DD=D1*D2 
    EC=ecount(g_sub) 
    MQ=(adj-DD/(2*EC))*Co_membership 
    Q=sum(MQ)/(2*EC) 
    return(Q) 
  } 
} 
 
