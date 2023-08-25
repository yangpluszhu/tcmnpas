P_analysis = function(prescription_basic,yao,mine_data,tS0.9=0.03,min_yao=8,threshold_leastYao=10,onlyHerb=FALSE){ 
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
    BF=BF_result$BF_result####BF 
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
    BF=BF_result$BF_result####BF 
 
  }  #return(list(BF=BF,herb_fre=herb_fre,opt_a=BF_result$opt_a,Degree_yao=BF_result$Degree_yao,Drug_pair_result=Drug_pair_result,noherb_fre=noherb_fre,disease_fre=disease_fre,P_stat=P_stat,disease_yao=disease_yao)) 
  return(list(BF=BF,opt_a=BF_result$opt_a,Degree_yao=BF_result$Degree_yao)) 
} 
 
