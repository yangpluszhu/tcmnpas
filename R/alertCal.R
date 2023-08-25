alertCal = function(sdfFile){ 
  require(ChemmineR) 
  require(ChemmineOB) 
  unwant=read.csv('db/unwanted_structure.csv',head=T,stringsAsFactors=F) 
  sdf=read.SDFset(sdfFile) 
  valid <- validSDF(sdf) 
  #sum(valid) 
  # datablocktag(sdf[which(!valid)],'id') 
  sdf <- sdf[valid] 
  cid(sdf)=datablocktag(sdf,'id') 
  alert=matrix(0,nrow=length(sdf),ncol=nrow(unwant)) 
  rownames(alert)=datablocktag(sdf,'id') 
  colnames(alert)=unwant$Name 
  for (i in cid(sdf)){ 
    for (j in 1:nrow(unwant)){ 
      alert[i,j]=smartsSearchOB(sdf[i],unwant$smarts[j]) 
    } 
  } 
  alert_value=rowSums(alert) 
  alert=as.data.frame(alert) 
  alert$id=datablocktag(sdf,'id') 
  alert$ALERTS=alert_value 
  return(alert) 
} 
 
