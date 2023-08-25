randomWalk = function(igraphM,queryGenes,EdgeWeight = FALSE, gamma = 0.7){ 
  require(PAGI) 
  require(igraph) 
  ##igraphM-->igraph should has vertex name! 
  ##queryGenes-->seedGeneNames_vector,vertex name 
  ##return-score--->Named Vector 
  #g=make_ring(10) 
  #V(g)$name <- letters[1:10] 
  vg=numeric(length(V(igraphM))) 
  names(vg)=V(igraphM)$name 
  index=queryGenes 
  vg[index]=1 
  score=RandomWalk2igraph(igraphM,vg,EdgeWeight,gamma) 
  names(score)=names(vg) 
  return(score) 
} 
 
