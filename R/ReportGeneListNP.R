ReportGeneListNP = function(GeneListFile=NULL,GeneList=NULL,GeneListName='CUSTOM',diseaseGeneID=NULL,GOMFenrich=T,GOBPenrich=T,GOCCenrich=T,KEGGenrich=T,PAenrich=T,DOenrich=T,MESHenrich=T,meshcategory='C',qvalueCutoff=0.05,TopNIDCurve=10,TOPNkeggPlot,showCategory,PLOTfontsize,PLOTwidth,PLOTheight,Savefilename=NULL){ 
  require(r2excel) 
  require(xlsx) 
  require(data.table) 
  require(pathview) 
  require(stringr) 
  if (!is.null(GeneListFile)){ 
    GeneList=readLines(GeneListFile) 
  } 
  ###### 
  if (sum(str_detect(GeneList,',|;| '))){ 
    GeneList=unlist(str_split(GeneList,',|;| ')) 
    GeneList=str_replace_all(GeneList,'"','') 
    GeneList=str_replace_all(GeneList,"'",'') 
    GeneList=GeneList[GeneList!=''] 
    GeneList=str_trim(GeneList,side='both') 
  } 
  GeneList=unique(GeneList) 
  ###### 
  ###### 
  TCMNPobj=GeneListNP(GeneList=GeneList,GeneListName=GeneListName,diseaseGeneID=diseaseGeneID,GOMFenrich=GOMFenrich,GOBPenrich=GOBPenrich,GOCCenrich=GOCCenrich,KEGGenrich=KEGGenrich,PAenrich=PAenrich,DOenrich=DOenrich,MESHenrich=MESHenrich,meshcategory=meshcategory,qvalueCutoff=qvalueCutoff,TopNIDCurve=TopNIDCurve) 
  wb <- createWorkbook(type="xlsx") 
  ################## 
  if (GOMFenrich&(!is.null(TCMNPobj$enrich))){ 
    sheetGOMF=createSheet(wb, sheetName = "基因本体分子功能富集") 
    xlsx.addTable(wb, sheetGOMF,TCMNPobj$enrich$GOmf$Enrichs) 
    sheetGOMFPIC=createSheet(wb, sheetName = "基因本体分子功能富集PIC") 
    plotFunction<-function(){ 
      p=Mdotplot(TCMNPobj$enrich$GOmf$OBject,showCategory=showCategory,font.size=PLOTfontsize) 
      print(p) 
    } 
    xlsx.addPlot(wb, sheetGOMFPIC, plotFunction=plotFunction,width=PLOTwidth,height=PLOTheight) 
  } 
  ### 
  if (GOBPenrich&(!is.null(TCMNPobj$enrich))){ 
    sheetGOBP=createSheet(wb, sheetName = "基因本体生物过程富集") 
    xlsx.addTable(wb, sheetGOBP,TCMNPobj$enrich$GObp$Enrichs) 
    sheetGOBPPIC=createSheet(wb, sheetName = "基因本体生物过程富集PIC") 
    plotFunction<-function(){ 
      p=Mdotplot(TCMNPobj$enrich$GObp$OBject,showCategory=showCategory,font.size=PLOTfontsize) 
      print(p) 
    } 
    xlsx.addPlot(wb, sheetGOBPPIC, plotFunction=plotFunction,width=PLOTwidth,height=PLOTheight) 
  } 
  ### 
  if (GOCCenrich&(!is.null(TCMNPobj$enrich))){ 
    sheetGOCC=createSheet(wb, sheetName = "基因本体细胞组件富集") 
    xlsx.addTable(wb, sheetGOCC,TCMNPobj$enrich$GOcc$Enrichs) 
    sheetGOCCPIC=createSheet(wb, sheetName = "基因本体细胞组件富集PIC") 
    plotFunction=function(){ 
      p=Mdotplot(TCMNPobj$enrich$GOcc$OBject,showCategory=showCategory,font.size=PLOTfontsize) 
      print(p) 
    } 
    xlsx.addPlot(wb, sheetGOCCPIC, plotFunction=plotFunction,width=PLOTwidth,height=PLOTheight) 
  } 
  #### 
  if ((GOMFenrich|GOBPenrich|GOCCenrich)&(!is.null(diseaseGeneID))){ 
    sheetGODisease=createSheet(wb, sheetName = "复方与疾病基因共关联GO") 
    txt1=paste('复方与疾病基因共关联GO相似度：',TCMNPobj$Disease$CommonDiseaseGOSim,sep='') 
    txt2=paste('复方与疾病基因共关联前',TopNIDCurve,'个GO本体曲线下面积AUC：',TCMNPobj$Disease$CommonGOauc,sep='') 
    xlsx.addParagraph(wb, sheetGODisease, txt1, fontSize=12, isItalic=TRUE, 
                      fontColor="darkred", backGroundColor="gray",startRow=1, 
                      colSpan=10, rowSpan=1) 
    xlsx.addParagraph(wb, sheetGODisease, txt2, fontSize=12, isItalic=TRUE, 
                      fontColor="darkred", backGroundColor="gray",startRow=2, 
                      colSpan=10, rowSpan=1) 
    plotFunction<-function(){ 
      print(TCMNPobj$Disease$PlotCommonGO) 
    } 
    xlsx.addPlot(wb, sheetGODisease, plotFunction=plotFunction,startRow=5) 
  } 
  #### 
  ### 
  if (KEGGenrich&(!is.null(TCMNPobj$enrich))){ 
    sheetKEGG=createSheet(wb, sheetName = "KEGG信号通路富集") 
    xlsx.addTable(wb, sheetKEGG,TCMNPobj$enrich$KEGG$Enrichs) 
    sheetKEGGPIC=createSheet(wb, sheetName = "KEGG信号通路富集PIC") 
    plotFunction<-function(){ 
      p=Mdotplot(TCMNPobj$enrich$KEGG$OBject,showCategory=showCategory,font.size=PLOTfontsize) 
      print(p) 
    } 
    xlsx.addPlot(wb, sheetKEGGPIC, plotFunction=plotFunction,width=PLOTwidth,height=PLOTheight) 
    ## 
    if (TOPNkeggPlot!=0){ 
      TOPNkeggPlot=as.numeric(TOPNkeggPlot) 
      pathWayID=unique(TCMNPobj$enrich$KEGG$Enrichs$ID) 
      TOPNkeggPlot=min(TOPNkeggPlot,length(pathWayID)) 
      pathWayIDtoPlot=pathWayID[1:TOPNkeggPlot] 
      NumPath=length(pathWayIDtoPlot) 
      pathview(gene.data=TCMNPobj$selectGene$GeneID,pathway.id=pathWayIDtoPlot,species="hsa",out.suffix='PA') 
      for (i in 1:NumPath){ 
        tempSheetName=pathWayIDtoPlot[i] 
        sheetTemp=createSheet(wb, sheetName =tempSheetName) 
        picfile=paste(tempSheetName,'.PA.png',sep='') 
        addPicture(picfile,sheetTemp,scale=5) 
        file.remove(picfile) 
        file.remove(paste(tempSheetName,'.png',sep='')) 
        file.remove(paste(tempSheetName,'.xml',sep='')) 
      } 
    } 
  } 
  ######## 
  if (PAenrich&(!is.null(TCMNPobj$enrich))){ 
    sheetPA=createSheet(wb, sheetName = "Reactome信号通路富集") 
    xlsx.addTable(wb, sheetPA,TCMNPobj$enrich$PA$Enrichs) 
    sheetPAPIC=createSheet(wb, sheetName = "Reactome信号通路富集PIC") 
    plotFunction<-function(){ 
      p=Mdotplot(TCMNPobj$enrich$PA$OBject,showCategory=showCategory,font.size=PLOTfontsize) 
      print(p) 
    } 
    xlsx.addPlot(wb, sheetPAPIC, plotFunction=plotFunction,width=PLOTwidth,height=PLOTheight) 
    ## 
  } 
  ######## 
  ## 
  if ((KEGGenrich)&(!is.null(diseaseGeneID))){ 
    sheetKEGGDisease=createSheet(wb, sheetName = "复方与疾病基因共关联KEGG通路") 
    txt1=paste('复方与疾病基因共关联KEGG信号通路相似度：',TCMNPobj$Disease$CommonDiseaseKEGGSim,sep='') 
    txt2=paste('复方与疾病基因共关联前',TopNIDCurve,'个KEGG信号通路曲线下面积AUC：',TCMNPobj$Disease$CommonKEGGauc,sep='') 
    xlsx.addParagraph(wb, sheetKEGGDisease, txt1, fontSize=12, isItalic=TRUE, 
                      fontColor="darkred", backGroundColor="gray",startRow=1, 
                      colSpan=10, rowSpan=1) 
    xlsx.addParagraph(wb, sheetKEGGDisease, txt2, fontSize=12, isItalic=TRUE, 
                      fontColor="darkred", backGroundColor="gray",startRow=2, 
                      colSpan=10, rowSpan=1) 
    plotFunction<-function(){ 
      print(TCMNPobj$Disease$PlotCommonKEGG) 
    } 
    xlsx.addPlot(wb, sheetKEGGDisease, plotFunction=plotFunction,startRow=5) 
  } 
  ########### 
  if ((PAenrich)&(!is.null(diseaseGeneID))){ 
    sheetPADisease=createSheet(wb, sheetName = "复方与疾病基因共关联Reactome通路") 
    txt1=paste('复方与疾病基因共关联Reactome信号通路相似度：',TCMNPobj$Disease$CommonDiseasePASim,sep='') 
    txt2=paste('复方与疾病基因共关联前',TopNIDCurve,'个Reactome信号通路曲线下面积AUC：',TCMNPobj$Disease$CommonPAauc,sep='') 
    xlsx.addParagraph(wb, sheetPADisease, txt1, fontSize=12, isItalic=TRUE, 
                      fontColor="darkred", backGroundColor="gray",startRow=1, 
                      colSpan=10, rowSpan=1) 
    xlsx.addParagraph(wb, sheetPADisease, txt2, fontSize=12, isItalic=TRUE, 
                      fontColor="darkred", backGroundColor="gray",startRow=2, 
                      colSpan=10, rowSpan=1) 
    plotFunction<-function(){ 
      print(TCMNPobj$Disease$PlotCommonPA) 
    } 
    xlsx.addPlot(wb, sheetPADisease, plotFunction=plotFunction,startRow=5) 
  } 
  ########### 
  if (DOenrich&(!is.null(TCMNPobj$enrich))){ 
    sheetDO=createSheet(wb, sheetName = "疾病本体富集") 
    xlsx.addTable(wb, sheetDO,TCMNPobj$enrich$DO$Enrichs) 
    sheetDOPIC=createSheet(wb, sheetName = "疾病本体富集PIC") 
    plotFunction<-function(){ 
      p=Mdotplot(TCMNPobj$enrich$DO$OBject,showCategory=showCategory,font.size=PLOTfontsize) 
      print(p) 
    } 
    xlsx.addPlot(wb, sheetDOPIC, plotFunction=plotFunction,width=PLOTwidth,height=PLOTheight) 
  } 
  ### 
  if (MESHenrich&(!is.null(TCMNPobj$enrich))){ 
    sheetMESH=createSheet(wb, sheetName = "主题词富集") 
    xlsx.addTable(wb, sheetMESH,TCMNPobj$enrich$MESH$Enrichs) 
    sheetMESHPIC=createSheet(wb, sheetName = "主题词富集PIC") 
    plotFunction<-function(){ 
      p=Mbarplot(TCMNPobj$enrich$MESH$OBject,showCategory=showCategory,font.size=PLOTfontsize) 
      print(p) 
    } 
    xlsx.addPlot(wb, sheetMESHPIC, plotFunction=plotFunction,width=PLOTwidth,height=PLOTheight) 
  } 
  ############################# 
  if (!dir.exists('Report'))dir.create('Report') 
  if (is.null(Savefilename))Savefilename='指定靶标群的机制分析' 
  saveWorkbook(wb, paste('Report/',Savefilename,'.xlsx',sep='')) 
} 
 
