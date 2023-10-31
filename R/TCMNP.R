##' Formula Mechanism 
##' 
##' @title TCMNP 
##' @param CustomherbTargetDat a dataframe with colname:herb,cid,chemical_name,geneID,score,inchikey,smiles,database,QED_DED 
##' @param herbList a character vector  
##' @param herbListName character 
##' @param QEDset QED threshold
##' @param Lipinski_Viosset Lipinski-Rules screening 
##' @param Veber_Viosset Veber-Rules screening 
##' @param geneScore gene Score threshold 
##' @param targetDatabase 'HIT','TCMID','STITCH','TCMSP', or 'CUSTOM' 
##' @param diseaseGeneID disease GeneID 
##' @param GOMFenrich True or False 
##' @param GOBPenrich True or False 
##' @param GOCCenrich True or False 
##' @param KEGGenrich True or False 
##' @param PAenrich True or False 
##' @param DOenrich True or False 
##' @param meshcategory mesh enrichment category  
##' @param qvalueCutoff enrichment adjusted P value 
##' @param TopNIDCurve The association curve of the top N items shared by diseases and the compound   
##' @return a list Object 
##' @export 
##' @author Yang Ming 
TCMNP = function(CustomherbTargetDat=NULL,herbList,herbListName='TCM',QEDset=0.2,geneScore=400,targetDatabase=c('HIT','TCMID','STITCH','TCMSP','CUSTOM'),diseaseGeneID=NULL,geneSelectPv=0.05,GOMFenrich=T,GOBPenrich=T,GOCCenrich=T,KEGGenrich=T,PAenrich=T,DOenrich=T,MESHenrich=T,meshcategory='C',qvalueCutoff=0.05,TopNIDCurve=10,Lipinski_Viosset=0,Veber_Viosset=0){
  #diseaseGeneID:data.frame('class','geneID') or character
  ##CustomherbTargetDat:colname-->c(herb,cid,chemical_name,geneID,score,inchikey,smiles,database,QED_DES,Lipinski_Vios,Veber_Vios)
  require(data.table)
  require(stringr)
  require(RSQLite)
  require(plyr)
  require(ggplot2)
  require(pracma)
  ##############################################
  if (!is.null(diseaseGeneID)){
    if (class(diseaseGeneID)%in%c('data.table','data.frame')){
      diseaseGeneID$geneID=as.character(diseaseGeneID$geneID)
      DGeneID=diseaseGeneID$geneID
      names(DGeneID)=diseaseGeneID$class
      DiseaseName=unique(names(DGeneID))
      DiseaseNum=length(DiseaseName)
    }else{
      if (is.null(names(diseaseGeneID))){
        diseaseGeneID=as.character(unlist(str_split(diseaseGeneID,',| |;')))
        diseaseGeneID=str_replace_all(diseaseGeneID,'"','')
        diseaseGeneID=str_replace_all(diseaseGeneID,"'",'')
        diseaseGeneID=diseaseGeneID[diseaseGeneID!='']
        diseaseGeneID=unique(str_trim(diseaseGeneID,side='both'))
        DGeneID=diseaseGeneID
        names(DGeneID)=rep('Disease',length(diseaseGeneID))
        DiseaseName='Disease'
        DiseaseNum=1
      }else{
        diseaseGeneID=as.character(unlist(str_split(diseaseGeneID,',| |;')))
        diseaseGeneID=str_replace_all(diseaseGeneID,'"','')
        diseaseGeneID=str_replace_all(diseaseGeneID,"'",'')
        diseaseGeneID=diseaseGeneID[diseaseGeneID!='']
        diseaseGeneID=unique(str_trim(diseaseGeneID,side='both'))
        DGeneID=diseaseGeneID
        DiseaseName=unique(names(DGeneID))
        DiseaseNum=length(DiseaseName)
      }
    }
  }
  #############################################################################################################################自定义输入靶标
  if (!is.null(CustomherbTargetDat)){
    standColname=c('herb','cid','chemical_name','geneID','score','inchikey','smiles','database','QED_DES',"Lipinski_Vios","Veber_Vios")
    colNameCustom=colnames(CustomherbTargetDat)
    if (sum(colNameCustom%in%standColname)==9){
      subD=CustomherbTargetDat
    }else{
      OutcolNameCustom=standColname[!standColname%in%colNameCustom]
      OutSideCustomData=data.frame(matrix(NA,ncol=length(OutcolNameCustom),nrow=nrow(CustomherbTargetDat)))
      colnames(OutSideCustomData)=OutcolNameCustom
      if ('score'%in%OutcolNameCustom)OutSideCustomData$score=9999
      if ('database'%in%OutcolNameCustom)OutSideCustomData$database='Custom'
      if ('QED_DES'%in%OutcolNameCustom)OutSideCustomData$QED_DES=1
	  if ('Lipinski_Vios'%in%OutcolNameCustom)OutSideCustomData$Lipinski_Vios=0
	  if ('Veber_Vios'%in%OutcolNameCustom)OutSideCustomData$Veber_Vios=0
      TempsubD=cbind(CustomherbTargetDat,OutSideCustomData)
      subD=TempsubD[,c(standColname,setdiff(colNameCustom,standColname)),with=F]
    }
    #####################################
    actHerb=unique(subD$herb)
    subD=subD[QED_DES>=QEDset&score>=geneScore&Lipinski_Vios<=Lipinski_Viosset&Veber_Vios<=Veber_Viosset,]
    subD$geneID=as.character(subD$geneID)
    tempR=geneScoreCal(unique(subD[,.(herb,id=cid,geneID=geneID)]),Pcutoff = geneSelectPv)
    selectGeneID=tempR$selectGeneId###selectGeneID
    if (length(selectGeneID)!=0){
      selectGeneName=geneTrans(selectGeneID,type='ID')
      tempR2=chemScoreCal(unique(subD[,.(id=cid,geneID=geneID)]),geneScore = tempR$geneScore,selectGeneId=selectGeneID)
      selectCID=tempR2$selectChemId#####selectCID
      subDD=unique(subD[,.(cid,chemical_name)])
      setkey(subDD,cid)
      selectCName=subDD[selectCID,chemical_name,mult="first"]
      subD2=subD[cid%in%selectCID&geneID%in%selectGeneID,]###subD2
      if (is.null(diseaseGeneID)){
        Order=paste('geneCluster=list(',herbListName,'=selectGeneID)',sep='')
      }else{
        if (DiseaseNum==1){
          Order=paste('geneCluster=list(',herbListName,'=selectGeneID,',DiseaseName,'=DGeneID[names(DGeneID)==DiseaseName])',sep='')
        }
        ##
        if (DiseaseNum==2){
          Order=paste('geneCluster=list(',herbListName,'=selectGeneID,',DiseaseName[1],'=DGeneID[names(DGeneID)==DiseaseName[1]],',DiseaseName[2],'=DGeneID[names(DGeneID)==DiseaseName[2]])',sep='')
        }
        if (DiseaseNum==3){
          Order=paste('geneCluster=list(',herbListName,'=selectGeneID,',DiseaseName[1],'=DGeneID[names(DGeneID)==DiseaseName[1]],',DiseaseName[2],'=DGeneID[names(DGeneID)==DiseaseName[2]],',DiseaseName[3],'=DGeneID[names(DGeneID)==DiseaseName[3]])',sep='')
        }
      }
      eval(parse(text=Order))
      ###################################
      if (KEGGenrich){
        if (class(try({EnrichKEGG=EnrichMentAnalysis(geneCluster,type='KEGG',qvalueCutoff=qvalueCutoff)},silent=T))=='try-error')EnrichKEGG=NULL
      }else{
        EnrichKEGG=NULL
      }
      #################################
      if (PAenrich){
        if (class(try({EnrichPA=EnrichMentAnalysis(geneCluster,type='ReactomePA',qvalueCutoff=qvalueCutoff)},silent=T))=='try-error')EnrichPA=NULL
      }else{
        EnrichPA=NULL
      }
      ############
      if (GOBPenrich){
        if (class(try({enrichGObp=EnrichMentAnalysis(geneCluster,type='GO',ontGO='BP',qvalueCutoff=qvalueCutoff)},silent=T))=='try-error')enrichGObp=NULL
      }else{
        enrichGObp=NULL
      }
      #############
      if (GOMFenrich){
        if (class(try({enrichGOmf=EnrichMentAnalysis(geneCluster,type='GO',ontGO='MF',qvalueCutoff=qvalueCutoff)},silent=T))=='try-error')enrichGOmf=NULL
      }else{
        enrichGOmf=NULL
      }
      ############
      if (GOCCenrich){
        if (class(try({enrichGOcc=EnrichMentAnalysis(geneCluster,type='GO',ontGO='CC',qvalueCutoff=qvalueCutoff)},silent=T))=='try-error')enrichGOcc=NULL
      }else{
        enrichGOcc=NULL
      }
      ############
      if (DOenrich){
        data(EG2DOLite)
        data(DOLiteTerm)
        data(DOLite2EG)
        if (class(try({enrichmentDO=EnrichMentAnalysis(geneCluster,type='DO',qvalueCutoff=qvalueCutoff)},silent=T))=='try-error')enrichmentDO=NULL
      }else{
        enrichmentDO=NULL
      }
      ############
      if (MESHenrich){
        if (class(try({enrichMESH=EnrichMentAnalysis(geneCluster,type='MESH',meshcategory=meshcategory,qvalueCutoff=qvalueCutoff)},silent=T))=='try-error')enrichMESH=NULL
      }else{
        enrichMESH=NULL
      }
      ##############################################################################
      if (!is.null(diseaseGeneID)&!is.null(EnrichKEGG)&(herbListName%in%EnrichKEGG$Enrichs$Cluster)){
        tempD=EnrichKEGG$Enrichs
        tempD$Cluster=as.character(tempD$Cluster)
        tempD$ID=as.character(tempD$ID)
        item=unique(tempD$Cluster)
        diseaseItem=item[!item%in%herbListName]
        CommonDiseaseKEGGSim=numeric()
        herbListNameItem=tempD$ID[tempD$Cluster==herbListName]
        if (length(diseaseItem)>0){
          length(CommonDiseaseKEGGSim)=length(diseaseItem)
          names(CommonDiseaseKEGGSim)=diseaseItem
          for (i in 1:length(diseaseItem)){
            diseaseItemtemp=tempD$ID[tempD$Cluster==diseaseItem[i]]
            CommonDiseaseKEGGSim[i]=length(intersect(herbListNameItem,diseaseItemtemp))/min(length(herbListNameItem),length(diseaseItemtemp))
          }
          diseaseToCurve=diseaseItem[1]
          diseasetNameItem=tempD$ID[tempD$Cluster==diseaseToCurve]
          TopNIDCurve=min(TopNIDCurve,length(herbListNameItem))
          common_KeggDisaeae=numeric()
          for (k in 1:TopNIDCurve){
            tempC=herbListNameItem[1:k]
            common_KeggDisaeae[k]=length(intersect(tempC,diseasetNameItem))/k
          }
          common_Keggdata=data.frame(ID=1:TopNIDCurve,Disease=common_KeggDisaeae)
          common_Keggdata=as.data.table(common_Keggdata)
          common_Keggdata=melt(common_Keggdata,id='ID',variable='Disease',value='Percentage')
          common_Keggdata=as.data.frame(common_Keggdata)
          CommonKEGGauc=trapz(1:TopNIDCurve,common_KeggDisaeae)
          PlotCommonKEGG=ggplot(common_Keggdata)+aes(x=ID,y=Percentage,group=Disease,colour=Disease,shape=Disease,fill=Disease)+geom_line()+geom_point(size=2)+xlab(paste('Top N enriched-Kegg pathway in ',herbListName))+ylab('Percentage of common enriched-Kegg pathway')+theme(legend.text=element_text(size=10),legend.title=element_text(size=10))+theme(axis.title.x=element_text(size=8),axis.title.y=element_text(size=8))+theme(axis.text.x=element_text(size=6),axis.text.y=element_text(size=6))
        }
        else{
          CommonDiseaseKEGGSim=NULL
          CommonKEGGauc=NULL
          PlotCommonKEGG=NULL
        }
        ######
      }else{
        CommonDiseaseKEGGSim=NULL
        CommonKEGGauc=NULL
        PlotCommonKEGG=NULL
      }
      ###############################################################################
      ###########ReactomePA#############################
      ##############################################################################
      if (!is.null(diseaseGeneID)&!is.null(EnrichPA)&(herbListName%in%EnrichPA$Enrichs$Cluster)){
        tempD=EnrichPA$Enrichs
        tempD$Cluster=as.character(tempD$Cluster)
        tempD$ID=as.character(tempD$ID)
        item=unique(tempD$Cluster)
        diseaseItem=item[!item%in%herbListName]
        CommonDiseasePASim=numeric()
        herbListNameItem=tempD$ID[tempD$Cluster==herbListName]
        if (length(diseaseItem)>0){
          length(CommonDiseasePASim)=length(diseaseItem)
          names(CommonDiseasePASim)=diseaseItem
          for (i in 1:length(diseaseItem)){
            diseaseItemtemp=tempD$ID[tempD$Cluster==diseaseItem[i]]
            CommonDiseasePASim[i]=length(intersect(herbListNameItem,diseaseItemtemp))/min(length(herbListNameItem),length(diseaseItemtemp))
          }
          diseaseToCurve=diseaseItem[1]
          diseasetNameItem=tempD$ID[tempD$Cluster==diseaseToCurve]
          TopNIDCurve=min(TopNIDCurve,length(herbListNameItem))
          common_PADisaeae=numeric()
          for (k in 1:TopNIDCurve){
            tempC=herbListNameItem[1:k]
            common_PADisaeae[k]=length(intersect(tempC,diseasetNameItem))/k
          }
          common_PAdata=data.frame(ID=1:TopNIDCurve,Disease=common_PADisaeae)
          common_PAdata=as.data.table(common_PAdata)
          common_PAdata=melt(common_PAdata,id='ID',variable='Disease',value='Percentage')
          common_PAdata=as.data.frame(common_PAdata)
          CommonPAauc=trapz(1:TopNIDCurve,common_PADisaeae)
          PlotCommonPA=ggplot(common_PAdata)+aes(x=ID,y=Percentage,group=Disease,colour=Disease,shape=Disease,fill=Disease)+geom_line()+geom_point(size=2)+xlab(paste('Top N enriched-Reactome pathway in ',herbListName))+ylab('Percentage of common enriched-Reactome pathway')+theme(legend.text=element_text(size=10),legend.title=element_text(size=10))+theme(axis.title.x=element_text(size=8),axis.title.y=element_text(size=8))+theme(axis.text.x=element_text(size=6),axis.text.y=element_text(size=6))
        }
        else{
          CommonDiseasePASim=NULL
          CommonPAauc=NULL
          PlotCommonPA=NULL
        }
        ######
      }else{
        CommonDiseasePASim=NULL
        CommonPAauc=NULL
        PlotCommonPA=NULL
      }
      ################PA END
      ################################################################################
      BlankGO=(!is.null(enrichGObp))|(!is.null(enrichGOcc))|(!is.null(enrichGOmf))
      BlankherbListName=(herbListName%in%enrichGObp$Enrichs$Cluster)|(herbListName%in%enrichGOcc$Enrichs$Cluster)|(herbListName%in%enrichGOmf$Enrichs$Cluster)
      if (!is.null(diseaseGeneID)&BlankGO&BlankherbListName){
        tempD=rbind(enrichGObp$Enrichs,enrichGOcc$Enrichs,enrichGOmf$Enrichs)
        tempD=arrange(tempD,qvalue)
        tempD$Cluster=as.character(tempD$Cluster)
        tempD$ID=as.character(tempD$ID)
        item=unique(tempD$Cluster)
        diseaseItem=item[!item%in%herbListName]
        CommonDiseaseGOSim=numeric()
        herbListNameItem=unique(tempD$ID[tempD$Cluster==herbListName])
        if (length(diseaseItem)>0){
          length(CommonDiseaseGOSim)=length(diseaseItem)
          names(CommonDiseaseGOSim)=diseaseItem
          for (i in 1:length(diseaseItem)){
            diseaseItemtemp=unique(tempD$ID[tempD$Cluster==diseaseItem[i]])
            CommonDiseaseGOSim[i]=length(intersect(herbListNameItem,diseaseItemtemp))/min(length(herbListNameItem),length(diseaseItemtemp))
          }
          diseaseToCurve=diseaseItem[1]
          diseasetNameItem=unique(tempD$ID[tempD$Cluster==diseaseToCurve])
          TopNIDCurve=min(TopNIDCurve,length(herbListNameItem))
          common_GODisaeae=numeric()
          for (k in 1:TopNIDCurve){
            tempC=herbListNameItem[1:k]
            common_GODisaeae[k]=length(intersect(tempC,diseasetNameItem))/k
          }
          common_GOdata=data.frame(ID=1:TopNIDCurve,Disease=common_GODisaeae)
          common_GOdata=as.data.table(common_GOdata)
          common_GOdata=melt(common_GOdata,id='ID',variable='Disease',value='Percentage')
          common_GOdata=as.data.frame(common_GOdata)
          CommonGOauc=trapz(1:TopNIDCurve,common_GODisaeae)
          PlotCommonGO=ggplot(common_GOdata)+aes(x=ID,y=Percentage,group=Disease,colour=Disease,shape=Disease,fill=Disease)+geom_line()+geom_point(size=2)+xlab(paste('Top N enriched-GO in ',herbListName))+ylab('Percentage of common enriched-GO')+theme(legend.text=element_text(size=10),legend.title=element_text(size=10))+theme(axis.title.x=element_text(size=8),axis.title.y=element_text(size=8))+theme(axis.text.x=element_text(size=6),axis.text.y=element_text(size=6))
        }else{
          CommonDiseaseGOSim=NULL
          CommonGOauc=NULL
          PlotCommonGO=NULL
        }
        ##############
      }else{
        CommonDiseaseGOSim=NULL
        CommonGOauc=NULL
        PlotCommonGO=NULL
      }
      
      #result=list(selectGene=data.table(GeneID=selectGeneID,GeneName=selectGeneName),selectChem=data.table(CID=selectCID,ChemName=selectCName),HerbTarget=subD2,actHerb=actHerb,enrich=list(KEGG=EnrichKEGG,PA=EnrichPA,GObp=enrichGObp,GOmf=enrichGOmf,GOcc=enrichGOcc,DO=enrichmentDO,MESH=enrichMESH),Disease=list(CommonDiseaseKEGGSim=CommonDiseaseKEGGSim,CommonKEGGauc=CommonKEGGauc,PlotCommonKEGG=PlotCommonKEGG,CommonDiseasePASim=CommonDiseasePASim,CommonPAauc=CommonPAauc,PlotCommonPA=PlotCommonPA,CommonDiseaseGOSim=CommonDiseaseGOSim,CommonGOauc=CommonGOauc,PlotCommonGO=PlotCommonGO),geneScore=tempR,chemScore=tempR2)
	  result=list(selectGene=data.table(GeneID=selectGeneID,GeneName=selectGeneName,geneScore=tempR$geneScore),selectChem=data.table(CID=selectCID,ChemName=selectCName,chemScore=tempR2$chemScore[tempR2$chemScore>0]),HerbTarget=subD2,actHerb=actHerb,enrich=list(KEGG=EnrichKEGG,PA=EnrichPA,GObp=enrichGObp,GOmf=enrichGOmf,GOcc=enrichGOcc,DO=enrichmentDO,MESH=enrichMESH),Disease=list(CommonDiseaseKEGGSim=CommonDiseaseKEGGSim,CommonKEGGauc=CommonKEGGauc,PlotCommonKEGG=PlotCommonKEGG,CommonDiseasePASim=CommonDiseasePASim,CommonPAauc=CommonPAauc,PlotCommonPA=PlotCommonPA,CommonDiseaseGOSim=CommonDiseaseGOSim,CommonGOauc=CommonGOauc,PlotCommonGO=PlotCommonGO),geneScore=tempR,chemScore=tempR2)
    }else{
      result=list(selectGene=data.table(GeneID=character(),GeneName=character()),selectChem=data.table(CID=character(),ChemName=character()),HerbTarget=data.table(),actHerb=actHerb,enrich=NULL,Disease=NULL)
    }
    ###############
    
    #####################################
  }else{###非自定义输入靶标
    #########################################################################################################
    if (!exists('AllHerbListData')){
      #ALL=fread('db/AllHerbListData.csv')
      AllHerbListData=fread('db/AllHerbListData.csv',encoding='UTF-8')
      #load('db/AllHerbListData.RData')
    }
    herbListALL=AllHerbListData$herb
    herbList2=unlist(str_split(herbList,',| |;'))
    actHerb=herbList2[herbList2%in%herbListALL]
    actHerbID=AllHerbListData[AllHerbListData$herb%in%actHerb,herbID]
    if (length(actHerbID)!=0){
      paID2=paste("'",actHerbID,"'",sep='',collapse = ',')
      paID3=paste("(",paID2,")",sep='')
      herbIDOrder=paste("herbID IN ",paID3,sep='')
      OrderherbID=paste("SELECT * FROM STITCH_HIT_TCMID_TCMSP_Tar where",herbIDOrder,sep=' ')
      tmp <- dbConnect(SQLite(), 'db/Tar.db')
      #tmp <- dbConnect(SQLite(), 'program/library/base/R/Tar.db')
      res <- dbSendQuery(tmp, OrderherbID)
      subD0 <- fetch(res, n =-1)
      dbDisconnect(tmp)
      subD=join(subD0,AllHerbListData,by='herbID')
      subD=as.data.table(subD)
      subD=subD[,!'herbID',with=F]
      subD=subD[,c('herb',setdiff(colnames(subD),'herb')),with=F]
      subD=subD[QED_DES>=QEDset&score>=geneScore&database%in%targetDatabase&Lipinski_Vios<=Lipinski_Viosset&Veber_Vios<=Veber_Viosset,]
      rm(subD0)
      gc()
      subD$geneID=as.character(subD$geneID)
      tempR=geneScoreCal(unique(subD[,.(herb,id=cid,geneID=geneID)]),Pcutoff = geneSelectPv)
      selectGeneID=tempR$selectGeneId###selectGeneID
      if (length(selectGeneID)!=0){
        selectGeneName=geneTrans(selectGeneID,type='ID')
        tempR2=chemScoreCal(unique(subD[,.(id=cid,geneID=geneID)]),geneScore = tempR$geneScore,selectGeneId=selectGeneID)
        selectCID=tempR2$selectChemId#####selectCID
        subDD=unique(subD[,.(cid,chemical_name)])
        setkey(subDD,cid)
        selectCName=subDD[selectCID,chemical_name,mult="first"]
        subD2=subD[cid%in%selectCID&geneID%in%selectGeneID,]###subD2
        if (is.null(diseaseGeneID)){
          Order=paste('geneCluster=list(',herbListName,'=selectGeneID)',sep='')
        }else{
          if (DiseaseNum==1){
            Order=paste('geneCluster=list(',herbListName,'=selectGeneID,',DiseaseName,'=DGeneID[names(DGeneID)==DiseaseName])',sep='')
          }
          ##
          if (DiseaseNum==2){
            Order=paste('geneCluster=list(',herbListName,'=selectGeneID,',DiseaseName[1],'=DGeneID[names(DGeneID)==DiseaseName[1]],',DiseaseName[2],'=DGeneID[names(DGeneID)==DiseaseName[2]])',sep='')
          }
          if (DiseaseNum==3){
            Order=paste('geneCluster=list(',herbListName,'=selectGeneID,',DiseaseName[1],'=DGeneID[names(DGeneID)==DiseaseName[1]],',DiseaseName[2],'=DGeneID[names(DGeneID)==DiseaseName[2]],',DiseaseName[3],'=DGeneID[names(DGeneID)==DiseaseName[3]])',sep='')
          }
        }
        eval(parse(text=Order))
        ###################################
        if (KEGGenrich){
          if (class(try({EnrichKEGG=EnrichMentAnalysis(geneCluster,type='KEGG',qvalueCutoff=qvalueCutoff)},silent=T))=='try-error')EnrichKEGG=NULL
        }else{
          EnrichKEGG=NULL
        }
        ##############
        if (PAenrich){
          if (class(try({EnrichPA=EnrichMentAnalysis(geneCluster,type='ReactomePA',qvalueCutoff=qvalueCutoff)},silent=T))=='try-error')EnrichPA=NULL
        }else{
          EnrichPA=NULL
        }
        ############
        if (GOBPenrich){
          if (class(try({enrichGObp=EnrichMentAnalysis(geneCluster,type='GO',ontGO='BP',qvalueCutoff=qvalueCutoff)},silent=T))=='try-error')enrichGObp=NULL
        }else{
          enrichGObp=NULL
        }
        #############
        if (GOMFenrich){
          if (class(try({enrichGOmf=EnrichMentAnalysis(geneCluster,type='GO',ontGO='MF',qvalueCutoff=qvalueCutoff)},silent=T))=='try-error')enrichGOmf=NULL
        }else{
          enrichGOmf=NULL
        }
        ############
        if (GOCCenrich){
          if (class(try({enrichGOcc=EnrichMentAnalysis(geneCluster,type='GO',ontGO='CC',qvalueCutoff=qvalueCutoff)},silent=T))=='try-error')enrichGOcc=NULL
        }else{
          enrichGOcc=NULL
        }
        ############
        if (DOenrich){
          data(EG2DOLite)
          data(DOLiteTerm)
          data(DOLite2EG)
          if (class(try({enrichmentDO=EnrichMentAnalysis(geneCluster,type='DO',qvalueCutoff=qvalueCutoff)},silent=T))=='try-error')enrichmentDO=NULL
        }else{
          enrichmentDO=NULL
        }
        ############
        if (MESHenrich){
          if (class(try({enrichMESH=EnrichMentAnalysis(geneCluster,type='MESH',meshcategory=meshcategory,qvalueCutoff=qvalueCutoff)},silent=T))=='try-error')enrichMESH=NULL
        }else{
          enrichMESH=NULL
        }
        ##############################################################################
        if (!is.null(diseaseGeneID)&!is.null(EnrichKEGG)&(herbListName%in%EnrichKEGG$Enrichs$Cluster)){
          tempD=EnrichKEGG$Enrichs
          tempD$Cluster=as.character(tempD$Cluster)
          tempD$ID=as.character(tempD$ID)
          item=unique(tempD$Cluster)
          diseaseItem=item[!item%in%herbListName]
          CommonDiseaseKEGGSim=numeric()
          herbListNameItem=tempD$ID[tempD$Cluster==herbListName]
          if (length(diseaseItem)>0){
            length(CommonDiseaseKEGGSim)=length(diseaseItem)
            names(CommonDiseaseKEGGSim)=diseaseItem
            for (i in 1:length(diseaseItem)){
              diseaseItemtemp=tempD$ID[tempD$Cluster==diseaseItem[i]]
              CommonDiseaseKEGGSim[i]=length(intersect(herbListNameItem,diseaseItemtemp))/min(length(herbListNameItem),length(diseaseItemtemp))
            }
            diseaseToCurve=diseaseItem[1]
            diseasetNameItem=tempD$ID[tempD$Cluster==diseaseToCurve]
            TopNIDCurve=min(TopNIDCurve,length(herbListNameItem))
            common_KeggDisaeae=numeric()
            for (k in 1:TopNIDCurve){
              tempC=herbListNameItem[1:k]
              common_KeggDisaeae[k]=length(intersect(tempC,diseasetNameItem))/k
            }
            common_Keggdata=data.frame(ID=1:TopNIDCurve,Disease=common_KeggDisaeae)
            common_Keggdata=as.data.table(common_Keggdata)
            common_Keggdata=melt(common_Keggdata,id='ID',variable='Disease',value='Percentage')
            common_Keggdata=as.data.frame(common_Keggdata)
            CommonKEGGauc=trapz(1:TopNIDCurve,common_KeggDisaeae)
            PlotCommonKEGG=ggplot(common_Keggdata)+aes(x=ID,y=Percentage,group=Disease,colour=Disease,shape=Disease,fill=Disease)+geom_line()+geom_point(size=2)+xlab(paste('Top N enriched-Kegg pathway in ',herbListName))+ylab('Percentage of common enriched-Kegg pathway')+theme(legend.text=element_text(size=10),legend.title=element_text(size=10))+theme(axis.title.x=element_text(size=8),axis.title.y=element_text(size=8))+theme(axis.text.x=element_text(size=6),axis.text.y=element_text(size=6))
          }
          else{
            CommonDiseaseKEGGSim=NULL
            CommonKEGGauc=NULL
            PlotCommonKEGG=NULL
          }
          ######
        }else{
          CommonDiseaseKEGGSim=NULL
          CommonKEGGauc=NULL
          PlotCommonKEGG=NULL
        }
        ###############################################################################
        ###########ReactomePA#############################
        ##############################################################################
        if (!is.null(diseaseGeneID)&!is.null(EnrichPA)&(herbListName%in%EnrichPA$Enrichs$Cluster)){
          tempD=EnrichPA$Enrichs
          tempD$Cluster=as.character(tempD$Cluster)
          tempD$ID=as.character(tempD$ID)
          item=unique(tempD$Cluster)
          diseaseItem=item[!item%in%herbListName]
          CommonDiseasePASim=numeric()
          herbListNameItem=tempD$ID[tempD$Cluster==herbListName]
          if (length(diseaseItem)>0){
            length(CommonDiseasePASim)=length(diseaseItem)
            names(CommonDiseasePASim)=diseaseItem
            for (i in 1:length(diseaseItem)){
              diseaseItemtemp=tempD$ID[tempD$Cluster==diseaseItem[i]]
              CommonDiseasePASim[i]=length(intersect(herbListNameItem,diseaseItemtemp))/min(length(herbListNameItem),length(diseaseItemtemp))
            }
            diseaseToCurve=diseaseItem[1]
            diseasetNameItem=tempD$ID[tempD$Cluster==diseaseToCurve]
            TopNIDCurve=min(TopNIDCurve,length(herbListNameItem))
            common_PADisaeae=numeric()
            for (k in 1:TopNIDCurve){
              tempC=herbListNameItem[1:k]
              common_PADisaeae[k]=length(intersect(tempC,diseasetNameItem))/k
            }
            common_PAdata=data.frame(ID=1:TopNIDCurve,Disease=common_PADisaeae)
            common_PAdata=as.data.table(common_PAdata)
            common_PAdata=melt(common_PAdata,id='ID',variable='Disease',value='Percentage')
            common_PAdata=as.data.frame(common_PAdata)
            CommonPAauc=trapz(1:TopNIDCurve,common_PADisaeae)
            PlotCommonPA=ggplot(common_PAdata)+aes(x=ID,y=Percentage,group=Disease,colour=Disease,shape=Disease,fill=Disease)+geom_line()+geom_point(size=2)+xlab(paste('Top N enriched-Reactome pathway in ',herbListName))+ylab('Percentage of common enriched-Reactome pathway')+theme(legend.text=element_text(size=10),legend.title=element_text(size=10))+theme(axis.title.x=element_text(size=8),axis.title.y=element_text(size=8))+theme(axis.text.x=element_text(size=6),axis.text.y=element_text(size=6))
          }
          else{
            CommonDiseasePASim=NULL
            CommonPAauc=NULL
            PlotCommonPA=NULL
          }
          ######
        }else{
          CommonDiseasePASim=NULL
          CommonPAauc=NULL
          PlotCommonPA=NULL
        }
        ################PA END
        ################################################################################
        ################
        BlankGO=(!is.null(enrichGObp))|(!is.null(enrichGOcc))|(!is.null(enrichGOmf))
        BlankherbListName=(herbListName%in%enrichGObp$Enrichs$Cluster)|(herbListName%in%enrichGOcc$Enrichs$Cluster)|(herbListName%in%enrichGOmf$Enrichs$Cluster)
        if (!is.null(diseaseGeneID)&BlankGO&BlankherbListName){
          tempD=rbind(enrichGObp$Enrichs,enrichGOcc$Enrichs,enrichGOmf$Enrichs)
          tempD=arrange(tempD,qvalue)
          tempD$Cluster=as.character(tempD$Cluster)
          tempD$ID=as.character(tempD$ID)
          item=unique(tempD$Cluster)
          diseaseItem=item[!item%in%herbListName]
          CommonDiseaseGOSim=numeric()
          herbListNameItem=unique(tempD$ID[tempD$Cluster==herbListName])
          if (length(diseaseItem)>0){
            length(CommonDiseaseGOSim)=length(diseaseItem)
            names(CommonDiseaseGOSim)=diseaseItem
            for (i in 1:length(diseaseItem)){
              diseaseItemtemp=unique(tempD$ID[tempD$Cluster==diseaseItem[i]])
              CommonDiseaseGOSim[i]=length(intersect(herbListNameItem,diseaseItemtemp))/min(length(herbListNameItem),length(diseaseItemtemp))
            }
            diseaseToCurve=diseaseItem[1]
            diseasetNameItem=unique(tempD$ID[tempD$Cluster==diseaseToCurve])
            TopNIDCurve=min(TopNIDCurve,length(herbListNameItem))
            common_GODisaeae=numeric()
            for (k in 1:TopNIDCurve){
              tempC=herbListNameItem[1:k]
              common_GODisaeae[k]=length(intersect(tempC,diseasetNameItem))/k
            }
            common_GOdata=data.frame(ID=1:TopNIDCurve,Disease=common_GODisaeae)
            common_GOdata=as.data.table(common_GOdata)
            common_GOdata=melt(common_GOdata,id='ID',variable='Disease',value='Percentage')
            common_GOdata=as.data.frame(common_GOdata)
            CommonGOauc=trapz(1:TopNIDCurve,common_GODisaeae)
            PlotCommonGO=ggplot(common_GOdata)+aes(x=ID,y=Percentage,group=Disease,colour=Disease,shape=Disease,fill=Disease)+geom_line()+geom_point(size=2)+xlab(paste('Top N enriched-GO in ',herbListName))+ylab('Percentage of common enriched-GO')+theme(legend.text=element_text(size=10),legend.title=element_text(size=10))+theme(axis.title.x=element_text(size=8),axis.title.y=element_text(size=8))+theme(axis.text.x=element_text(size=6),axis.text.y=element_text(size=6))
          }else{
            CommonDiseaseGOSim=NULL
            CommonGOauc=NULL
            PlotCommonGO=NULL
          }
          ##############
        }else{
          CommonDiseaseGOSim=NULL
          CommonGOauc=NULL
          PlotCommonGO=NULL
        }
        
        result=list(selectGene=data.table(GeneID=selectGeneID,GeneName=selectGeneName,geneScore=tempR$geneScore),selectChem=data.table(CID=selectCID,ChemName=selectCName,chemScore=tempR2$chemScore[tempR2$chemScore>0]),HerbTarget=subD2,actHerb=actHerb,enrich=list(KEGG=EnrichKEGG,PA=EnrichPA,GObp=enrichGObp,GOmf=enrichGOmf,GOcc=enrichGOcc,DO=enrichmentDO,MESH=enrichMESH),Disease=list(CommonDiseaseKEGGSim=CommonDiseaseKEGGSim,CommonKEGGauc=CommonKEGGauc,PlotCommonKEGG=PlotCommonKEGG,CommonDiseasePASim=CommonDiseasePASim,CommonPAauc=CommonPAauc,PlotCommonPA=PlotCommonPA,CommonDiseaseGOSim=CommonDiseaseGOSim,CommonGOauc=CommonGOauc,PlotCommonGO=PlotCommonGO),geneScore=tempR,chemScore=tempR2)
      }else{
        result=list(selectGene=data.table(GeneID=character(),GeneName=character()),selectChem=data.table(CID=character(),ChemName=character()),HerbTarget=data.table(),actHerb=actHerb,enrich=NULL,Disease=NULL)
      }
      ###############
    }else{
      result=list(selectGene=data.table(GeneID=character(),GeneName=character()),selectChem=data.table(CID=character(),ChemName=character()),HerbTarget=data.table(),actHerb=actHerb,enrich=NULL,Disease=NULL)
    }
    #######################################
  }
  return(result)
}
 
