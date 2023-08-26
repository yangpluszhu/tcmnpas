##' Combine igraph Object
##'
##' @title combineIgraphOB
##' @param g1 a igraph Object
##' @param g2 a igraph Object
##' @return a igraph object
##' @export
##' @author Yang Ming
combineIgraphOB = function(g1,g2){ 
  require(igraph) 
  Node1_attr=vertex_attr_names(g1) 
  Node2_attr=vertex_attr_names(g2) 
  common_attr=intersect(Node1_attr,Node2_attr) 
  common_attr=setdiff(common_attr,'name') 
  Node1_attr_1=paste(common_attr,'1',sep='_') 
  Node1_attr_2=paste(common_attr,'2',sep='_') 
  ########## 
  E1_attr=edge_attr_names(g1) 
  E2_attr=edge_attr_names(g2) 
  Ecommon_attr=intersect(E1_attr,E2_attr) 
  E1_attr_1=paste(Ecommon_attr,'1',sep='_') 
  E1_attr_2=paste(Ecommon_attr,'2',sep='_') 
  ######### 
  g=g1+g2 
  ########## 
 
  if (length(common_attr)!=0){ 
    for (i in 1:length(common_attr)){ 
      tempText1=paste('V(g)$',Node1_attr_2[i],'[is.na(V(g)$',Node1_attr_2[i],')]','<-','V(g)$',Node1_attr_1[i],'[is.na(V(g)$',Node1_attr_2[i],')]',sep='') 
      eval(parse(text=tempText1)) 
      order1=paste('V(g)$',common_attr[i],'=','V(g)$',Node1_attr_2[i],sep='') 
      eval(parse(text=order1)) 
      ###################### 
      TTorder1=paste('g=delete_vertex_attr(g,','"',Node1_attr_1[i],'"',')',sep='') 
      TTorder2=paste('g=delete_vertex_attr(g,','"',Node1_attr_2[i],'"',')',sep='') 
      eval(parse(text=TTorder1)) 
      eval(parse(text=TTorder2)) 
    } 
    ### 
  } 
  ######################## 
  if (length(Ecommon_attr)!=0){ 
    for (i in 1:length(Ecommon_attr)){ 
      tempText2=paste('E(g)$',E1_attr_2[i],'[is.na(E(g)$',E1_attr_2[i],')]','<-','E(g)$',E1_attr_1[i],'[is.na(E(g)$',E1_attr_2[i],')]',sep='') 
      eval(parse(text=tempText2)) 
      ### 
      order2=paste('E(g)$',Ecommon_attr[i],'=','E(g)$',E1_attr_2[i],sep='') 
      eval(parse(text=order2)) 
      ### 
      ETTorder1=paste('g=delete_edge_attr(g,','"',E1_attr_1[i],'"',')',sep='') 
      ETTorder2=paste('g=delete_edge_attr(g,','"',E1_attr_2[i],'"',')',sep='') 
      eval(parse(text=ETTorder1)) 
      eval(parse(text=ETTorder2)) 
    } 
    ####### 
  } 
  return(g) 
} 
 
