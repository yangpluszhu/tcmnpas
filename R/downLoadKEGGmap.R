downLoadKEGGmap = function(geneID=NULL,keggid){
  require(pathview)
  require(stringr)
  #############
  pathWayIDtoPlot=unlist(str_split(keggid,',|;| '))
  pathWayIDtoPlot=str_replace_all(pathWayIDtoPlot,'"','')
  pathWayIDtoPlot=str_replace_all(pathWayIDtoPlot,"'",'')
  pathWayIDtoPlot=pathWayIDtoPlot[pathWayIDtoPlot!='']
  pathWayIDtoPlot=str_trim(pathWayIDtoPlot,side='both')
  ############
  if (str_detect(geneID,',|;| ')){
    geneID=unlist(str_split(geneID,',|;| '))
    geneID=str_replace_all(geneID,'"','')
    geneID=str_replace_all(geneID,"'",'')
    geneID=geneID[geneID!='']
    geneID=str_trim(geneID,side='both')
  }
  geneID=unique(geneID)
  ############
  NumPath=length(unique(pathWayIDtoPlot))
  pathview(gene.data=geneID,pathway.id=pathWayIDtoPlot,species="hsa",out.suffix='PA')
  for (i in 1:NumPath){
    tempSheetName=pathWayIDtoPlot[i]
    #sheetTemp=createSheet(wb, sheetName =tempSheetName)
    picfile=paste(tempSheetName,'.PA.png',sep='')
    #addPicture(picfile,sheetTemp,scale=5)
    file.copy(picfile,to='Report/',overwrite=T)
    file.remove(picfile)
    file.remove(paste(tempSheetName,'.png',sep=''))
    file.remove(paste(tempSheetName,'.xml',sep=''))
  }
}

