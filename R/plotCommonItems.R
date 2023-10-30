plotCommonItems = function(DataFrame,Ca,Cb,TopNIDCurve=30,pointSize=2,legendTextSize=8,axisTitleSize=8,axisTextSize=6,plotName='Commonplot.png',plotWidth=10,plotHeight=8,xLab='Top N enriched-Item',yLab='Percentage of common enriched-Items'){
library(ggplot2)
library(pracma)
library(data.table)
options(stringsAsFactors = F)
tempD=DataFrame
tempD$Cluster=as.character(tempD$Cluster)
tempD$ID=as.character(tempD$ID)
item=unique(tempD$Cluster)
tempD=as.data.table(tempD)
tempD=tempD[order(Cluster,qvalue)]
herbItem=tempD$ID[tempD$Cluster==Ca]
diseaseItem=tempD$ID[tempD$Cluster==Cb]
CommonItem=intersect(herbItem,diseaseItem)
if (length(diseaseItem)>0){
  TopNIDCurve=min(TopNIDCurve,length(herbItem))
  common_CommonItem=numeric()
  for (k in 1:TopNIDCurve){
    tempC=herbItem[1:k]
    common_CommonItem[k]=length(intersect(tempC,diseaseItem))/k
  }
  common_data=data.frame(ID=1:TopNIDCurve,Disease=common_CommonItem)
  common_data=as.data.table(common_data)
  common_data=melt(common_data,id='ID',variable='Disease',value='Percentage')
  common_data=as.data.frame(common_data)
  Commonauc=trapz(1:TopNIDCurve,common_CommonItem)
  PlotCommon=ggplot(common_data)+aes(x=ID,y=Percentage,group=Disease,colour=Disease,shape=Disease,fill=Disease)+geom_line()+geom_point(size=pointSize)+xlab(xLab)+ylab(yLab)+theme(legend.text=element_text(size=legendTextSize),legend.title=element_text(size=legendTextSize))+theme(axis.title.x=element_text(size=axisTitleSize),axis.title.y=element_text(size=axisTitleSize))+theme(axis.text.x=element_text(size=axisTextSize),axis.text.y=element_text(size=axisTextSize))+theme(legend.position = 'none')
ggsave(filename=plotName,plot=PlotCommon,width=plotWidth,height=plotHeight,dpi=600,units='in')
  }else{
  PlotCommon=NULL
  Commonauc=NULL
}
return(list(PlotObject=PlotCommon,AUC=Commonauc))
}

