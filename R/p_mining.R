p_mining = function(tS0.9=0.03,min_yao=8,threshold_leastYao=10,onlyHerb=TRUE){
  ###p_mining=function(prescription,mine_data,TopDisease=5,TopDrug=10,tS0.9=0.03,min_yao=8,threshold_leastYao=10,onlyHerb=TRUE,Drug_pair_analysis=FALSE)
  require(data.table)
  require(plyr)
  require(stringr)
  require(foreach)
  require(igraph)
  require(doParallel)
  require(iterators)
  require(Matrix)
  prescription_basicData=basic_transDBF(prescriptionDBF_data=prescription,mine_data=mine_data)
  P_analysis_result=P_analysis(prescription_basic=prescription_basicData$prescription_basic,yao=prescription_basicData$yao,mine_data=mine_data,tS0.9=tS0.9,min_yao=min_yao,threshold_leastYao=threshold_leastYao,onlyHerb=onlyHerb)
  #write.csv(P_analysis_result$noherb_fre,'Report/Tab4_1西药及中成药使用频次统计.csv',row.names = F)
  #write.csv(P_analysis_result$herb_fre,'Report/Tab4_2草药使用频次统计.csv',row.names = F)
  #write.csv(P_analysis_result$disease_fre,'Report/Tab4_3疾病频次统计.csv',row.names = F)
  #write.csv(P_analysis_result$P_stat,'Report/Tab4_4处方中草药_中成药_西药的使用个数统计.csv',row.names = F)
  write.csv(P_analysis_result$BF,'Report/Tab4_2核心处方分析结果.csv',row.names = F)
  write.csv(P_analysis_result$Degree_yao,'Report/Tab4_1核心药物分析结果.csv',row.names = F)

}

