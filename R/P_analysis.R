P_analysis = function(prescription_basic,yao,mine_data,tS0.9=0.03,min_yao=8,threshold_leastYao=10,onlyHerb=FALSE){
  ##P_analysis(prescription_basic,yao,mine_data,TopDisease=5,TopDrug=10,tS0.9=0.03,min_yao=8,threshold_leastYao=10,onlyHerb=FALSE,Drug_pair_analysis=FALSE,local=TRUE)
  ##prescription_basic-->"id","time", "patient","diagnosis", "doctor","","drug"[binary]
  ##mine_data--->"num","drug_id","mine_id","mine_Cname", "drug_commonName", "mine_name" ,"class","fufang", "jixing","guige""生产企业"         "规格", "jixing2" , "innComponentName","type"合并了<--herb_data
  ##yao-->drug向量包含西药及成药，由basic_trans输出得到
  ##TopDisease-->统计前几位的疾病
  ##TopDrug--->统计关联度最强的前几位的药物
  ##min_yao-->基本方发现时需要指定的最小含药数
  ##threshold_leastYao-->自适应优化二值化阈值时采用的预期含药数
  ##herb_name-->由herb_data得到unique(herb_data$mine_name+herb_data$innComponentName)
  ##z_name-->表示成药向量
  ##x_name-->表示西药向量
  require(plyr)
  require(igraph)
  require(data.table)
  require(Matrix)
  prescription_basic=as.data.frame(prescription_basic)
  herb_name=unique(mine_data[class=='h',mine_Cname])
  ##########################################
  if (onlyHerb){
    noherb=yao[!yao%in%herb_name]
    Pyao=yao[yao%in%herb_name]
    P=subset(prescription_basic,select=Pyao)
    P=as.matrix(P)
    P=Matrix(P,sparse=T)
    rownames(P)=prescription_basic$id
    # if (Drug_pair_analysis){
    #   Drug_pair_result=drug_pair(P)
    # }else{
    #   Drug_pair_result=na_table
    # }
    #####
    BF_result=BF_analysis(P,tS0.9=tS0.9,min_yao=min_yao,threshold_leastYao=threshold_leastYao)
    BF=BF_result$BF_result####BF基本方信息
    }else{
    P=subset(prescription_basic,select=yao)
    P=as.matrix(P)
    P=Matrix(P,sparse=T)
    rownames(P)=prescription_basic$id
    # if (Drug_pair_analysis){
    #   Drug_pair_result=drug_pair(P)
    # }else{
    #   Drug_pair_result=na_table
    # }
    #####
    BF_result=BF_analysis(P,tS0.9=tS0.9,min_yao=min_yao,threshold_leastYao=threshold_leastYao)
    BF=BF_result$BF_result####BF基本方信息

  }
  ###########################################

  #############################################
  #############################################
  #return(list(BF=BF,herb_fre=herb_fre,opt_a=BF_result$opt_a,Degree_yao=BF_result$Degree_yao,Drug_pair_result=Drug_pair_result,noherb_fre=noherb_fre,disease_fre=disease_fre,P_stat=P_stat,disease_yao=disease_yao))
  return(list(BF=BF,opt_a=BF_result$opt_a,Degree_yao=BF_result$Degree_yao))
}

