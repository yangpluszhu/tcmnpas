sim = function(x,y){ 
  ##x,y--->two vectors not continurous! 
  ##output-->list(result=result,sim_score=sim_score) 
  common_item=intersect(x,y) 
  sim_score=length(common_item)/min(length(x),length(y)) 
  sim_score=round(sim_score,digits = 3) 
  if (length(common_item)==0){ 
    result_txt='non-overlapping' 
    sim_score=0 
    result=result_txt 
  }else{ 
    result_txt=paste(common_item,collapse = ',') 
    result=paste('S=',sim_score,'->',result_txt,sep='') 
  } 
  return(list(result=result,sim_score=sim_score)) 
} 
 
