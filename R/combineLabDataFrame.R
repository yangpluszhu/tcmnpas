combineLabDataFrame = function(labData1,labData2){ 
  Name1=unlist(labData1[1,]) 
  Name2=unlist(labData2[1,]) 
  NameAll=unique(c(Name1,Name2)) 
  #First6Names=Name2[1:6] 
  #NameAllDel6name=setdiff(NameAll,First6Names) 
  NameAllD=data.frame(ID=paste('V',1:length(NameAll),sep=''),Name=NameAll,stringsAsFactors=F) 
  Name1Out=setdiff(NameAll,Name1) 
  if (length(Name1Out)>0){ 
    Name1OutD=as.data.frame(matrix('*',nrow=nrow(labData1),ncol=length(Name1Out)),stringsAsFactors=F) 
    Name1OutD[1,]=Name1Out 
    Name1D=cbind(labData1,Name1OutD) 
  }else{ 
    Name1D=labData1 
  } 
  Name1List=unlist(Name1D[1,]) 
  colName1=NameAllD$ID[match(Name1List,NameAllD$Name)] 
  colnames(Name1D)=colName1 
  ### 
  Name2Out=setdiff(NameAll,Name2) 
  if (length(Name2Out)>0){ 
    Name2OutD=as.data.frame(matrix('*',nrow=nrow(labData2),ncol=length(Name2Out)),stringsAsFactors=F) 
    Name2OutD[1,]=Name2Out 
    Name2D=cbind(labData2,Name2OutD) 
  }else{ 
    Name2D=labData2 
  } 
  Name2List=unlist(Name2D[1,]) 
  colName2=NameAllD$ID[match(Name2List,NameAllD$Name)] 
  colnames(Name2D)=colName2 
  Name2D=Name2D[,colName1] 
  Name2D=Name2D[-1,] 
  ### 
  dataAll=rbind(Name1D,Name2D) 
  #dataAll=dataAll[,c(First6Names,NameAllDel6name)] 
  return(dataAll) 
} 
 
