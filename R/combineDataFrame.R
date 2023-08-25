combineDataFrame = function(D1,D2){ 
  name1=colnames(D1) 
  name2=colnames(D2) 
  name=unique(c(name1,name2)) 
  D1out=setdiff(name,name1) 
  if (length(D1out)>0){ 
    D1outData=matrix(0,nrow=nrow(D1),ncol=length(D1out)) 
    D1outData=as.data.frame(D1outData) 
    colnames(D1outData)=D1out 
    D1C=cbind(D1,D1outData) 
  }else{ 
    D1C=D1 
  } 
  D1C=D1C[,name] 
  ###### 
  D2out=setdiff(name,name2) 
  if (length(D1out)>0){ 
    D2outData=matrix(0,nrow=nrow(D2),ncol=length(D2out)) 
    D2outData=as.data.frame(D2outData) 
    colnames(D2outData)=D2out 
    D2C=cbind(D2,D2outData) 
  }else{ 
    D2C=D2 
  } 
  D2C=D2C[,name] 
  result=rbind(D1C,D2C) 
} 
 
