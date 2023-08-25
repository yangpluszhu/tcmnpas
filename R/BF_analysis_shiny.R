BF_analysis_shiny = function(Pdata,tS0.9=0.03,min_yao=8,threshold_leastYao=10,exact=FALSE,combineSim=FALSE,combine_method=2,combine_threshold=0.7,csvInput=T,calRs=F){ 
  ##output:list$BF_result:makeup+statstatistic;list$yao_fre 
  require(igraph) 
  require(plyr) 
  require(data.table) 
  require(stringr) 
  options(stringsAsFactors = F) 
  BF_analysisOB=BF_analysis(Pdata=Pdata,tS0.9=tS0.9,min_yao=min_yao,threshold_leastYao=threshold_leastYao,exact=exact,combineSim=combineSim,combine_method=combine_method,combine_threshold=combine_threshold,csvInput=csvInput) 
  ##return(list(BF_result=BF_result,yao_fre=yao_fre,opt_a=k_opt,Degree_yao=Degree_yao2,network=g2,dataP0=tempP0)) 
  if (nrow(BF_analysisOB$BF_result)>1){ 
    coMemberMatrixResult=coMemberMatrix(BF_analysisOB) 
  }else{ 
    coMemberMatrixResult=NULL 
  } 
  ##### 
  if (calRs){ 
    BF_evalResult=BF_eval(BF_result=BF_analysisOB,yao=BF_analysisOB$dataP0) 
  }else{ 
    BF_evalResult=NULL 
  } 
  return(list(BF_result=BF_analysisOB$BF_result,yao_fre=BF_analysisOB$yao_fre,opt_a=BF_analysisOB$opt_a,Degree_yao=BF_analysisOB$Degree_yao,BFnetwork=BF_analysisOB$BFnetwork,dataP0=BF_analysisOB$dataP0,BF_evalResult=BF_evalResult,coMemberMatrixResult=coMemberMatrixResult,yao_yaoNetwork=BF_analysisOB$yao_yaoNetwork)) 
} 
 
