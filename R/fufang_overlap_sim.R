fufang_overlap_sim = function(drug_trans,drug,fufang_name){
  require(data.table)
  require(stringr)
  require(lubridate)
  temp_fufang=unique(drug[drug%in%fufang_name])
  fufang_overlap_sim_result=0
  Nfufang=length(temp_fufang)
  if (Nfufang<2){
    fufang_overlap_sim_result=0
  }else{
    #temp_result=matrix(0,nrow=Nfufang,ncol=Nfufang)
    #name_result=matrix(0,nrow=Nfufang,ncol=Nfufang)
    temp_result_score=matrix(0,nrow=Nfufang,ncol=Nfufang)
    for (i in 1:(Nfufang-1)){
      f1=temp_fufang[i]
      j=i+1
      while (j<=Nfufang){
        f2=temp_fufang[j]
        sim_result=sim(drug_trans[drug%in%f1],drug_trans[drug%in%f2])
        #temp_result[i,j]=sim_result$result
        #name_result[i,j]=paste(f1,f2,sep='_')
        temp_result_score[i,j]=sim_result$sim_score
        j=j+1
      }
    }
    #temp_result_sub=temp_result[upper.tri(temp_result,diag = FALSE)]
    #name_result_sub=name_result[upper.tri(name_result,diag = FALSE)]
    temp_result_score_sub=temp_result_score[upper.tri(temp_result_score,diag = FALSE)]
    fufang_overlap_sim_result=max(temp_result_score_sub,na.rm=T)
  }
  return(fufang_overlap_sim_result)
}

