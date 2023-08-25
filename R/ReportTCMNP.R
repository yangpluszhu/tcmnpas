ReportTCMNP = function(CustomherbTargetDat=NULL,herbList,herbListName='TCM',QEDset=0.2,geneScore=800,targetDatabase=c('HIT','TCMID','STITCH','TCMSP','CUSTOM'),diseaseGeneID=NULL,geneSelectPv=0.01,GOMFenrich=T,GOBPenrich=T,GOCCenrich=T,KEGGenrich=T,PAenrich=T,DOenrich=T,MESHenrich=T,meshcategory='C',qvalueCutoff=0.05,TOPNkeggPlot=1,showCategory=5,font.size=12,width=1024,height=950,TopNIDCurve=10,Savefilename){ 
  require(r2excel) 
  require(xlsx) 
  require(data.table) 
  require(pathview) 
  if (!is.null(CustomherbTargetDat)){ 
    CustomherbTargetDat=fread(CustomherbTargetDat,encoding='UTF-8') 
  } 
  TCMNPobj=TCMNP(CustomherbTargetDat=CustomherbTargetDat,herbList=herbList,herbListName=herbListName,QEDset=QEDset,geneScore=geneScore,targetDatabase=targetDatabase,diseaseGeneID=diseaseGeneID,geneSelectPv=geneSelectPv,GOMFenrich=GOMFenrich,GOBPenrich=GOBPenrich,GOCCenrich=GOCCenrich,KEGGenrich=KEGGenrich,PAenrich=PAenrich,DOenrich=DOenrich,MESHenrich=MESHenrich,meshcategory=meshcategory,qvalueCutoff=qvalueCutoff,TopNIDCurve=TopNIDCurve) 
  wb <- createWorkbook(type="xlsx") 
  sheet1 <- createSheet(wb, sheetName = "复方靶标筛选结果") 
  sheet2 <- createSheet(wb, sheetName = "复方活性化合物筛选结果") 
  sheet3 <- createSheet(wb, sheetName = "药物-成分-靶标检索结果") 
  ################## 
  xlsx.addHeader(wb, sheet1, " 复方靶标筛选结果", level=2, underline=1) 
  xlsx.addLineBreak(sheet1, 2) 
  paragraph1=paste('输入的药物列表：',toString(herbList),'||检索到的药物列表：',toString(TCMNPobj$actHerb),'||检索的靶标数据库：',toString(targetDatabase),'||设定的基因靶标关联阈值：',geneScore,'||设定的基因靶标显著性阈值：',geneSelectPv) 
  xlsx.addParagraph(wb, sheet1, paragraph1, fontSize=12, isItalic=TRUE, 
                    fontColor="darkred", backGroundColor="gray", 
                    colSpan=10, rowSpan=3) 
  xlsx.addTable(wb, sheet1,TCMNPobj$selectGene,startRow=8) 
  ################# 
  xlsx.addHeader(wb, sheet2, " 复方活性化合物筛选结果", level=2, underline=1) 
  xlsx.addLineBreak(sheet2, 2) 
  paragraph2=paste('输入的药物列表：',toString(herbList),'||检索到的药物列表：',toString(TCMNPobj$actHerb),'||检索的靶标数据库：',toString(targetDatabase),'||设定的化合物类药性指数阈值：',QEDset,'||设定的基因靶标显著性阈值：',geneSelectPv) 
  xlsx.addParagraph(wb, sheet2, paragraph2, fontSize=12, isItalic=TRUE, 
                    fontColor="darkred", backGroundColor="gray", 
                    colSpan=10, rowSpan=3) 
  xlsx.addTable(wb, sheet2,TCMNPobj$selectChem,startRow=8) 
  ################# 
  xlsx.addHeader(wb, sheet3, " 药物-成分-靶标检索结果", level=2, underline=1) 
  xlsx.addLineBreak(sheet3, 2) 
  paragraph3=paste('输入的药物列表：',toString(herbList),'||检索到的药物列表：',toString(TCMNPobj$actHerb),'||检索的靶标数据库：',toString(targetDatabase),'||设定的基因靶标关联阈值：',geneScore,'||设定的化合物类药性指数阈值：',QEDset,'||设定的基因靶标显著性阈值：',geneSelectPv) 
  xlsx.addParagraph(wb, sheet3, paragraph3, fontSize=12, isItalic=TRUE, 
                    fontColor="darkred", backGroundColor="gray", 
                    colSpan=10, rowSpan=3) 
  xlsx.addTable(wb, sheet3,TCMNPobj$HerbTarget,startRow=8) 
  ################## 
  if (GOMFenrich&(!is.null(TCMNPobj$enrich))){ 
    sheetGOMF=createSheet(wb, sheetName = "基因本体分子功能富集") 
    xlsx.addTable(wb, sheetGOMF,TCMNPobj$enrich$GOmf$Enrichs) 
    sheetGOMFPIC=createSheet(wb, sheetName = "基因本体分子功能富集PIC") 
    plotFunction<-function(){ 
      p=Mdotplot(TCMNPobj$enrich$GOmf$OBject,showCategory=showCategory,font.size=font.size) 
      print(p) 
      } 
    xlsx.addPlot(wb, sheetGOMFPIC, plotFunction=plotFunction,width=width,height=height) 
  } 
  ### 
  if (GOBPenrich&(!is.null(TCMNPobj$enrich))){ 
    sheetGOBP=createSheet(wb, sheetName = "基因本体生物过程富集") 
    xlsx.addTable(wb, sheetGOBP,TCMNPobj$enrich$GObp$Enrichs) 
    sheetGOBPPIC=createSheet(wb, sheetName = "基因本体生物过程富集PIC") 
    plotFunction<-function(){ 
      p=Mdotplot(TCMNPobj$enrich$GObp$OBject,showCategory=showCategory,font.size=font.size) 
      print(p) 
      } 
    xlsx.addPlot(wb, sheetGOBPPIC, plotFunction=plotFunction,width=width,height=height) 
  } 
  ### 
  if (GOCCenrich&(!is.null(TCMNPobj$enrich))){ 
    sheetGOCC=createSheet(wb, sheetName = "基因本体细胞组件富集") 
    xlsx.addTable(wb, sheetGOCC,TCMNPobj$enrich$GOcc$Enrichs) 
    sheetGOCCPIC=createSheet(wb, sheetName = "基因本体细胞组件富集PIC") 
    plotFunction=function(){ 
      p=Mdotplot(TCMNPobj$enrich$GOcc$OBject,showCategory=showCategory,font.size=font.size) 
      print(p) 
      } 
    xlsx.addPlot(wb, sheetGOCCPIC, plotFunction=plotFunction,width=width,height=height) 
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
      p=Mdotplot(TCMNPobj$enrich$KEGG$OBject,showCategory=showCategory,font.size=font.size) 
      print(p) 
      } 
    xlsx.addPlot(wb, sheetKEGGPIC, plotFunction=plotFunction,width=width,height=height) 
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
      p=Mdotplot(TCMNPobj$enrich$PA$OBject,showCategory=showCategory,font.size=font.size) 
      print(p) 
    } 
    xlsx.addPlot(wb, sheetPAPIC, plotFunction=plotFunction,width=width,height=height) 
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
      p=Mdotplot(TCMNPobj$enrich$DO$OBject,showCategory=showCategory,font.size=font.size) 
      print(p) 
      } 
    xlsx.addPlot(wb, sheetDOPIC, plotFunction=plotFunction,width=width,height=height) 
  } 
  ### 
  if (MESHenrich&(!is.null(TCMNPobj$enrich))){ 
    sheetMESH=createSheet(wb, sheetName = "主题词富集") 
    xlsx.addTable(wb, sheetMESH,TCMNPobj$enrich$MESH$Enrichs) 
    sheetMESHPIC=createSheet(wb, sheetName = "主题词富集PIC") 
    plotFunction<-function(){ 
      p=Mbarplot(TCMNPobj$enrich$MESH$OBject,showCategory=showCategory,font.size=font.size) 
      print(p) 
      } 
    xlsx.addPlot(wb, sheetMESHPIC, plotFunction=plotFunction,width=width,height=height) 
  } 
  ############################# 
  if (!dir.exists('Report'))dir.create('Report') 
  saveWorkbook(wb, paste('Report/',Savefilename,'.xlsx',sep='')) 
} 
 
