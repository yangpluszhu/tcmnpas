overlap_analysis = function(){
  ##overlap_analysis=function(prescription,mine_data)
  require(data.table)
  require(plyr)
  require(stringr)
  require(foreach)
  require(igraph)
  require(doParallel)
  require(iterators)
  require(Matrix)
  prescription_transData=prescription_transDBF(prescriptionDBF_data=prescription,mine_data=mine_data)
  overlap_result=overlapDBF(prescription_transDBF_data=prescription_transData,mine_data=mine_data)
  overlap_result2=as.data.table(overlap_result)
  overlap_result2=overlap_result2[fufang_overlap_sim>0|herb_fufang_overlap_sim>0][order(-herb_fufang_overlap_sim)]
  overlap_stat_all=merge(overlap_result2,prescription,by='id',all.x=T,all.y=F)
  overlap_stat_all$fufang_overlap_sim=as.numeric(overlap_stat_all$fufang_overlap_sim)
  overlap_stat_all$herb_fufang_overlap_sim=as.numeric(overlap_stat_all$herb_fufang_overlap_sim)
  overlap_stat_report=list()
  overlap_stat_report[['doctor']]=overlap_stat_all[,.(fufang_overlap_sim=mean(fufang_overlap_sim),herb_fufang_overlap_sim=mean(herb_fufang_overlap_sim)),by='doctor']
  overlap_stat_report[['diagnosis']]=overlap_stat_all[,.(fufang_overlap_sim=mean(fufang_overlap_sim),herb_fufang_overlap_sim=mean(herb_fufang_overlap_sim)),by='diagnosis']
  overlap_stat_report[['TCMdiagnosis']]=overlap_stat_all[,.(fufang_overlap_sim=mean(fufang_overlap_sim),herb_fufang_overlap_sim=mean(herb_fufang_overlap_sim)),by='TCMdiagnosis']
  overlap_stat_report[['patient']]=overlap_stat_all[,.(fufang_overlap_sim=mean(fufang_overlap_sim),herb_fufang_overlap_sim=mean(herb_fufang_overlap_sim)),by='patient']
  #overlap_result3=overlap_result2[,!c('fufang_overlap_sim','herb_fufang_overlap_sim'),with=FALSE]
  overlap_result3=overlap_result2
  if (nrow(overlap_result3)==0){
    overlap_result4=as.data.frame(matrix(NA,nrow=1,ncol=ncol(overlap_result3)))
    colnames(overlap_result4)=colnames(overlap_result3)
  }else{
    overlap_result4=overlap_result3
  }
  write.csv(overlap_result4,'Report/Tab3重复用药审查结果.csv',row.names = F)
  write.csv(overlap_stat_report$doctor,'Report/Tab3_1基于开方医生的重复用药情况.csv',row.names = F)
  write.csv(overlap_stat_report$diagnosis,'Report/Tab3_2基于诊断疾病的重复用药情况.csv',row.names = F)
  write.csv(overlap_stat_report$TCMdiagnosis,'Report/Tab3_3基于中医诊断的重复用药情况.csv',row.names = F)
  write.csv(overlap_stat_report$patient,'Report/Tab3_4基于患者的重复用药情况.csv',row.names = F)
}

