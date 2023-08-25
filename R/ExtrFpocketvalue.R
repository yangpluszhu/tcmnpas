ExtrFpocketvalue = function(pqrFile='966c_out/966c_pockets.pqr',offsetvalue=5){ 
  require(stringr) 
  require(data.table) 
  options(stringsAsFactors = F) 
  a=readLines(pqrFile) 
  a=a[which(str_detect(a,'ATOM'))] 
  ExtrFpocket=function(x){ 
    x=unlist(str_split(x,' ')) 
    x=x[x!=''] 
    return(x) 
  } 
  Fpocketvalue=character() 
  for (i in 1:length(a)){ 
    temp=ExtrFpocket(a[i]) 
    Fpocketvalue=rbind(Fpocketvalue,temp) 
  } 
  colnames(Fpocketvalue)=c('atom','n1','n2','n3','pockID','x','y','z','d','r') 
  Fpocketvalue=as.data.frame(Fpocketvalue) 
  Fpocketvalue$x=as.numeric(Fpocketvalue$x) 
  Fpocketvalue$y=as.numeric(Fpocketvalue$y) 
  Fpocketvalue$z=as.numeric(Fpocketvalue$z) 
  Fpocketvalue$pockID=as.numeric(Fpocketvalue$pockID) 
  rownames(Fpocketvalue)<-NULL 
  Fpocketvalue=as.data.table(Fpocketvalue) 
  Fpocketvalue$pockID=as.character(Fpocketvalue$pockID) 
  MaxPocketIDtable=table(Fpocketvalue$pockID) 
  MaxPocketID=names(MaxPocketIDtable)[MaxPocketIDtable==max(MaxPocketIDtable)] 
  FpocketvalueStat=Fpocketvalue[,.(x=mean(x),y=mean(y),z=mean(z),sizex=offsetvalue+max(x)-min(x),sizey=offsetvalue+max(y)-min(y),sizez=offsetvalue+max(z)-min(z)),by='pockID'] 
  FpocketvalueStat=as.data.frame(FpocketvalueStat) 
  firstPid=FpocketvalueStat[FpocketvalueStat$pockID==MaxPocketID,] 
  return(list(firstPid=firstPid,FpocketvalueStat=FpocketvalueStat)) 
} 
 
