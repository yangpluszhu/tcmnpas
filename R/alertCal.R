##' Calculate ALERTS score in QED
##'
##' @title alertCal
##' @param sdfFile molecular file *.sdf
##' @param unwanted_structureFile a csv file with colnames:Name, smarts
##' @return a dataframe object
##' @export
##' @author Yang Ming
alertCal = function(sdfFile,unwanted_structureFile='db/unwanted_structure.csv'){ 
  require(ChemmineR) 
  require(ChemmineOB) 
  unwant=read.csv(unwanted_structureFile,head=T,stringsAsFactors=F) 
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
 
