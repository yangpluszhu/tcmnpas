matchClassic = function(TopN=5){
  require(data.table)
  require(stringr)
  prescriptionC=prescription
  setkey(prescriptionC,yaoId)
  setkey(mine_data,idm)
  prescriptionC$class=mine_data[prescriptionC$yaoId,class]
  prescriptionC=prescriptionC[prescriptionC$class%like%'h',]
  prescriptionC$mine_Cname=mine_data[prescriptionC$yaoId,mine_Cname]
  prescriptionCs=prescriptionC[,.(id,mine_Cname)]
  prescriptionCsHerbCount=unique(prescriptionCs,by=c('id','mine_Cname'))
  prescriptionCsHerbCount2=prescriptionCsHerbCount[,.(HerbCount=length(mine_Cname)),by='id']
  prescriptionCs=prescriptionCs[prescriptionCs$mine_Cname%in%classicYao,]
  prescriptionCs=unique(prescriptionCs,by=c('id','mine_Cname'))
  Matched=prescriptionCs[,.(Matched=MClassic(mine_Cname,classicF_Matrix,TopN)$Sim,MaxSim=MClassic(mine_Cname,classicF_Matrix,TopN)$MaxSim,MeanSim=MClassic(mine_Cname,classicF_Matrix,TopN)$MeanSim),by='id']
  Matched=Matched[order(-MaxSim,-MeanSim)]
  Matched2=merge(Matched,prescriptionCsHerbCount2,all.x = T, all.y = F,by='id')
  Matched2=Matched2[,.(id,HerbCount,Matched,MaxSim,MeanSim)]
  write.csv(Matched2,'Report/Tab15经典方匹配结果.csv',row.names=F)
  return(Matched2)
}

