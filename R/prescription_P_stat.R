prescription_P_stat = function(){
  require(data.table)
  prescription$yaoId=as.character(prescription$yaoId)
  mine_data$idm=as.character(mine_data$idm)
  herb_id=unique(mine_data[class%in%c('h','hf'),idm])
  z_id=unique(mine_data[class%in%c('z'),idm])
  x_id=unique(mine_data[class%in%c('x'),idm])
  Result_P_stat=unique(prescription,by=c('id','drug'))
  P_stat=Result_P_stat[,.(herb_count=sum(yaoId%in%herb_id),z_count=sum(yaoId%in%z_id),drug_count=sum(yaoId%in%x_id)),by='id']
  P_stat=P_stat[order(-herb_count,-z_count,-drug_count)]
  write.csv(P_stat,'Report/Tab1各处方草药、成药、西药的使用个数统计列表.csv')
  #P_stat统计各处方草药、成药、西药的使用个数
  return(P_stat)
}

