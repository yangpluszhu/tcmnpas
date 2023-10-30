GOsimilarity = function(genelist1,genelist2){
  require(GOSemSim)
  d <- godata('org.Hs.eg.db', ont="BP", computeIC=FALSE)
  resultBP=GOSemSim::clusterSim(genelist1,genelist2,semData=d, measure="Wang")
  d <- godata('org.Hs.eg.db', ont="CC", computeIC=FALSE)
  resultCC=GOSemSim::clusterSim(genelist1,genelist2,semData=d, measure="Wang")
  d <- godata('org.Hs.eg.db', ont="MF", computeIC=FALSE)
  resultMF=GOSemSim::clusterSim(genelist1,genelist2,semData=d, measure="Wang")
  return(data.frame(BP=resultBP,MF=resultMF,CC=resultCC,MeanScore=mean(c(resultBP,resultMF,resultCC),na.rm=T)))
}

