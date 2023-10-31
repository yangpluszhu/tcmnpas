GeneListNP = function(GeneList,GeneListName='CUSTOM',diseaseGeneID=NULL,GOMFenrich=T,GOBPenrich=T,GOCCenrich=T,KEGGenrich=T,PAenrich=T,DOenrich=T,MESHenrich=T,meshcategory='C',qvalueCutoff=0.05,TopNIDCurve=10){ 
  #diseaseGeneID:data.frame('class','geneID') or character 
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
############################################################################################# 
######################################################################################################### 
    selectGeneID=GeneList 
    herbListName=GeneListName 
    selectGeneName=geneTrans(selectGeneID,type='ID') 
    ############## 
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
 
        result=list(selectGene=data.table(GeneID=selectGeneID,GeneName=selectGeneName),enrich=list(KEGG=EnrichKEGG,PA=EnrichPA,GObp=enrichGObp,GOmf=enrichGOmf,GOcc=enrichGOcc,DO=enrichmentDO,MESH=enrichMESH),Disease=list(CommonDiseaseKEGGSim=CommonDiseaseKEGGSim,CommonKEGGauc=CommonKEGGauc,PlotCommonKEGG=PlotCommonKEGG,CommonDiseasePASim=CommonDiseasePASim,CommonPAauc=CommonPAauc,PlotCommonPA=PlotCommonPA,CommonDiseaseGOSim=CommonDiseaseGOSim,CommonGOauc=CommonGOauc,PlotCommonGO=PlotCommonGO)) 
  ####################################### 
  return(result) 
} 
 
