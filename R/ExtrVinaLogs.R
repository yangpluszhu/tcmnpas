ExtrVinaLogs = function(LogFile='result.log'){ 
  require(stringr) 
  options(stringsAsFactors = F) 
  ExtractAffvalue=function(x){ 
    x=unlist(str_split(x,' ')) 
    x=x[x!=''] 
    x=as.numeric(x) 
    return(x) 
  } 
  #### 
  a=readLines(LogFile,skipNul=T) 
  a=a[(which(str_detect(a,'mode'))+3):(length(a)-1)] 
  Affvalue=numeric() 
  for (i in 1:length(a)){ 
    temp=ExtractAffvalue(a[i]) 
    Affvalue=rbind(Affvalue,temp) 
  } 
  colnames(Affvalue)=c('mode','affinity','dist1','dist2') 
  Affvalue=as.data.frame(Affvalue) 
  rownames(Affvalue)<-NULL 
  #### 
  return(Affvalue) 
} 
 
