DotPlot = function(Data,cluster,item,colorBy,sizeby,showCategory=5,ascending=T,xlabSize=8,ylabSize=8,plotWith=20,plotHeight=18,filename=NULL){ 
  require(data.table) 
  require(ggplot2) 
  if (is.null(filename))filename='DotPlot' 
  Data=as.data.table(Data) 
  Data=Data[,c(cluster,item,colorBy,sizeby),with=F] 
  if (ascending){ 
    DataSub=Data[order(get(cluster),get(colorBy))] 
  }else{ 
    DataSub=Data[order(get(cluster),-get(colorBy))] 
  } 
  DataSub=DataSub[,head(.SD, showCategory),.SDcols=setdiff(colnames(DataSub),cluster),by=cluster] 
  #DataSub[,c('sizebyVar'):=list((get(sizeby)-mean(get(sizeby)))/mean(get(sizeby))+1.5)] 
  #p=ggplot(DataSub)+aes(x=get(cluster),y=get(item),colour=get(colorBy))+geom_point(aes(size=get(sizeby)))+labs(size=sizeby,colour=colorBy)+scale_size_continuous(range(1,4))+scale_colour_gradient(low="green", high="red") 
  p=ggplot(DataSub)+aes(x=get(cluster),y=get(item),colour=get(colorBy))+geom_point(aes(size=get(sizeby)))+labs(size=sizeby,colour=colorBy)+scale_colour_gradient(low="red", high="green") 
  p=p+xlab('')+ylab('') 
  p=p+theme(axis.text.y=element_text(size=ylabSize,colour='black'),axis.text.x=element_text(size=xlabSize,colour='black')) 
  print(p) 
  if (!dir.exists('Report'))dir.create('Report') 
  ggsave(paste('Report/',filename,'.png',sep=''),plot=p,width=plotWith,height=plotHeight,units='cm',dpi=600) 
} 
 
