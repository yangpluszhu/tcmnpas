toVnaNet = function(net,filename='net.vna',NodeAtr=F,NodeData){ 
  ##net:igraphObject,V(net)$name,E(net)$weight 
  ##filename:*.vna 
  ##NodeAtr:whether to import Node Attributes DataFrame 
  ##NodeData:Node Attributes DataFrame 
  require(igraph) 
  require(readr) 
  require(data.table) 
  a=as_data_frame(net,  what  =  c("both")) 
  if (NodeAtr){ 
    Node=NodeData 
  }else{ 
    if (ncol(a$vertices)==0){ 
      a$vertices=data.frame(Node=as.character(1:nrow(a$vertices))) 
    } 
    Node=a$vertices 
  } 
  ##################### 
  if (ncol(a$edges)==2){ 
    a$edges$weight=1 
  } 
  Edge=as.data.table(a$edges) 
  Edge=unique(Edge) 
  Edge=as.data.frame(Edge) 
  Edge$from=as.character(Edge$from) 
  Edge$to=as.character(Edge$to) 
  write_lines('*Node data',filename) 
  write_csv(Node,filename,append=T,col_names=T) 
  write_lines('*Tie data',filename,append=T) 
  write_csv(Edge,filename,append=T,col_names=T) 
} 
 
