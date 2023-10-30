GODOsimilarity = function(genelist1,genelist2){
  require(DOSE)
  DOGOsimilarity=DOSE::clusterSim(genelist1,genelist2)
  return(DOGOsimilarity)
}

