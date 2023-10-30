combineDataBase = function(path=NULL,filename1,filename2,combinefileName,removeDuplicate='TRUE'){
  require(data.table)
  require(stringr)
  if (is.null(path)){
    data1=fread(filename1)
    data2=fread(filename2)
    dataAll=rbind(data1,data2)
    if (removeDuplicate){dataAll=unique(dataAll)}
    outPut=paste('Report/',combinefileName,sep='')
    write.csv(dataAll,outPut,row.names=F)
  }else{
    filesCSV=list.files(path,pattern='.csv')
    dataAll=data.table()
    for (i in filesCSV){
      tempName=paste(path,i,sep='/')
      tempData=fread(tempName)
      dataAll=rbind(dataAll,tempData)
    }
    if (removeDuplicate){dataAll=unique(dataAll)}
    outPut=paste('Report/',combinefileName,sep='')
    write.csv(dataAll,outPut,row.names=F)
  }
}

