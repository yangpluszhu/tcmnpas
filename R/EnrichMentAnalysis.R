##' Enrichment Analysis 
##' @title EnrichMentAnalysis 
##' @param geneClusterID a named list 
##' @param type Enrichment type:'KEGG' or 'GO' 'DO' or 'MESH' or 'ReactomePA' 
##' @param qvalueCutoff Adjusted P value 
##' @param filtKEGGdisease Remove KEGG disease terms 
##' @param ontGO Enrichment GO type 
##' @param ontDO Enrichment DO 
##' @param ClusterNameToMESHenrich a named list for gene cluster 
##' @param meshDatabase for MESH analysis 
##' @return a list Object 
##' @export 
##' @author Yang Ming 
EnrichMentAnalysis = function(geneClusterID,type,qvalueCutoff=0.05,filtKEGGdisease=F,ontGO='BP',ontDO='DOLite',ClusterNameToMESHenrich=names(geneClusterID)[1],meshDatabase='gene2pubmed',meshcategory='C'){ 
  #geneClusterID:list-->names(geneClusterID) 
  #type:'KEGG' or 'GO' 'DO' or 'MESH' or 'ReactomePA' 
  #ont:'BP','MF','CC' 
  #MeshEnrich: 
  ##ClusterNameToMESHenrich:names of Input List to analysis 
  #            meshDatabase:'gendoo','gene2pubmed','RBBH' 
  #            meshcategory:'C'-->disease,'D'-->Chemicals and Drugs,'E'-->Analytical, Diagnostic and Therapeutic Techniques and Equipment,'F'-->Psychiatry and Psychology,'G'-->Phenomena and Processes,'J'-->Technology and Food and Beverages,'N'-->Health Care 
  ##download_KEGG {clusterProfiler} 
  ##Mdotplot or Mbarplot 
  require(ReactomePA) 
  require(KEGG.db) 
  require(fgsea) 
  require(ggplot2) 
  require(stringr) 
  require(clusterProfiler) 
  require(meshes) 
  require(DOSE) 
  require(data.table) 
  require(plyr) 
  require(igraph) 
  require(corrplot) 
  require("org.Hs.eg.db") 
  #require(meshr) 
  if (type=='KEGG'){ 
    #geneClusterID=list(BXXXT=BXXXT_geneID,Colitis=Yan_geneID,DM=Diabet_geneID,GC=Cancer_geneID) 
    ckegg=compareCluster(geneCluster = geneClusterID, fun = "enrichKEGG", pAdjustMethod='fdr',qvalueCutoff=qvalueCutoff,,use_internal_data=T) 
    aa=attributes(ckegg) 
    temp_result=aa$compareClusterResult 
    temp_result=as.data.table(temp_result) 
    if (filtKEGGdisease){ 
      BXXXT_common_KeggID=unique(temp_result$ID) 
      BXXXT_common_KeggIDfilt=BXXXT_common_KeggID[!BXXXT_common_KeggID%like%'hsa05'] 
      BXXXT_common_KeggIDfilt=BXXXT_common_KeggIDfilt[!BXXXT_common_KeggIDfilt%in%c('hsa04930','hsa04940','hsa04950','hsa04931','hsa04932','hsa04933')] 
      temp_result=temp_result[ID%in%BXXXT_common_KeggIDfilt,] 
    } 
    temp_result=temp_result[p.adjust<qvalueCutoff,] 
    ckegg2=new("compareClusterResult", compareClusterResult = temp_result,geneClusters = geneClusterID, fun = 'Enrichment') 
    ### 
      #png('KeggEnrich.png',width=1024,height=800) 
      #Mdotplot(ckegg2,colorBy='qvalue',showCategory = showCategory,font.size=18) 
      #dev.off(which = dev.cur()) 
 
    return(list(Enrichs=temp_result,OBject=ckegg2)) 
  }else if (type=='GO'){ 
    cgo1=compareCluster(geneCluster = geneClusterID, OrgDb=org.Hs.eg.db,fun = "enrichGO",ont = ontGO,pAdjustMethod="fdr",qvalueCutoff=qvalueCutoff,readable=T) 
    #### 
      #png(paste('GO_Enrich_',ontGO,'.png',sep=''),width=1368,height=800) 
      #Mdotplot(cgo1,colorBy='qvalue',showCategory = showCategory,font.size=18) 
      #dev.off() 
    #### 
    return(list(Enrichs=cgo1@compareClusterResult,OBject=cgo1)) 
  }else if (type=='MESH'){ 
    BXXXT_geneID=geneClusterID[[ClusterNameToMESHenrich]] 
    meshEnrich=meshes::enrichMeSH(BXXXT_geneID, MeSHDb = "MeSH.Hsa.eg.db", database=meshDatabase, category = meshcategory,pAdjustMethod='fdr',qvalueCutoff=qvalueCutoff) 
    meshEnrichs=as.data.frame(meshEnrich) 
    meshObject<-new("enrichResult", result = meshEnrichs, qvalueCutoff = qvalueCutoff, 
                    pAdjustMethod = 'fdr', organism = 'Homo sapiens', 
                    ontology = 'MeSH', gene = BXXXT_geneID, 
                    universe = 'extID', 
                    readable = FALSE) 
   return(list(Enrichs=meshEnrichs,OBject=meshObject)) 
  }else if (type=='DO') { 
    data(EG2DOLite) 
    data(DOLiteTerm) 
    data(DOLite2EG) 
    cdo=compareCluster(geneCluster = geneClusterID, fun = "enrichDO",ont=ontDO,pAdjustMethod="fdr",qvalueCutoff=qvalueCutoff) 
    #### 
    #png('DO_Enrich.png',width=1368,height=800) 
    #Mdotplot(cdo,colorBy='qvalue',showCategory = showCategory,font.size=18) 
    #dev.off() 
    #### 
    return(list(Enrichs=cdo@compareClusterResult,OBject=cdo)) 
  }else if (type=='ReactomePA'){ 
    ReactomePA=compareCluster(geneCluster = geneClusterID, fun = "enrichPathway",pAdjustMethod="fdr",qvalueCutoff=qvalueCutoff,readable=T) 
    #### 
    #png(paste('GO_Enrich_',ontGO,'.png',sep=''),width=1368,height=800) 
    #Mdotplot(cgo1,colorBy='qvalue',showCategory = showCategory,font.size=18) 
    #dev.off() 
    #### 
    return(list(Enrichs=ReactomePA@compareClusterResult,OBject=ReactomePA)) 
 
  }else{ 
    print('ERROR:type not found!') 
  } 
} 
 
