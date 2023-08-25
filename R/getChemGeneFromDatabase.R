getChemGeneFromDatabase = function(ChemData,ChemDataType='inchikey',QEDset=0.2,geneScore=400,geneSelectPv=0.05,targetDatabase=c('HIT','TCMID','STITCH','TCMSP','CUSTOM'),TarDB='db/Tar.db'){ 
  #ChemData:data.frame:colnames-->'cid','chemical_name','key'.'key' is smiles or inchikey 
  #ChemDataType:'inchikey' or 'smiles'  #result:list(combinedGeneID=selectGeneID,GeneID_split=tempR_split,subD3=subD3):$GeneID_split-->data.frame(cid,chemical_name,inchikey,geneID).$subD3-->data.frame(cid,geneID,score,database,QED_DES,chemical_name,inchikey) 
  require(data.table) 
  require(stringr) 
  require(RSQLite) 
  require(plyr) 
  require(ggplot2) 
  require(pracma) 
  options(stringsAsFactors = F) 
  ChemData$key=as.character(ChemData$key) 
  ChemData$cid=as.character(ChemData$cid) 
  ChemData$chemical_name=as.character(ChemData$chemical_name) 
  #setwd('D:/PA2.1Plus') 
  ChemData2=ChemData 
  if (ChemDataType!='inchikey'){ 
    tempSmiles=ChemData$key 
    names(tempSmiles)=ChemData$cid 
    ChemData2$inchikey=unlist(TransInchikey(tempSmiles,molecularType='SMItxt')) 
  }else{ 
    ChemData2=plyr::rename(ChemData2,c('key'='inchikey')) 
  } 
  ############################ 
  ChemData2=ChemData2[!is.na(ChemData2$inchikey),] 
  ChemData2=as.data.table(ChemData2) 
  ChemData2=unique(ChemData2,by='inchikey') 
  QueryKey=unique(ChemData2$inchikey) 
  if (nrow(ChemData2)>0){ 
    paID2=paste("'",QueryKey,"'",sep='',collapse = ',') 
    paID3=paste("(",paID2,")",sep='') 
    chemIDOrder=paste("inchikey IN ",paID3,sep='') 
    OrderchemIDOrder=paste("SELECT * FROM STITCH_HIT_TCMID_TCMSP_Tar where",chemIDOrder,sep=' ') 
    tmp <- dbConnect(SQLite(), TarDB) 
    res <- dbSendQuery(tmp, OrderchemIDOrder) 
    subD0 <- fetch(res, n =-1) 
    dbDisconnect(tmp) 
    subD=join(subD0[,c('geneID','score','inchikey','database','QED_DES')],ChemData2,by='inchikey') 
    subD=as.data.table(subD) 
    subD=subD[QED_DES>=QEDset&score>=geneScore&database%in%targetDatabase,] 
    rm(subD0) 
    gc() 
    subD$geneID=as.character(subD$geneID) 
    subD2=unique(subD[,.(id=cid,geneID=geneID,score,database,QED_DES=QED_DES[1]),by=c('cid','geneID')]) 
    tempR=geneScoreCal(subD2,Pcutoff = geneSelectPv) 
    tempR_split=subD2[,.(geneID=toString(geneID),QED_DES=QED_DES[1]),by='id'] 
    tempR_split=plyr::rename(tempR_split,c('id'='cid')) 
    tempR_split=join(tempR_split,ChemData2,by='cid') 
    tempR_split=tempR_split[,.(cid,chemical_name,inchikey,geneID,QED_DES)] 
    selectGeneID=tempR$selectGeneId###selectGeneID 
    subD3=subD2[,.(cid,geneID,score,database,QED_DES)] 
    subD3=join(subD3,ChemData2,by='cid') 
  }else{ 
    selectGeneID=NA 
    tempR_split=NA 
    subD3=NA 
  } 
    return(list(combinedGeneID=selectGeneID,GeneID_split=tempR_split,subD3=subD3)) 
  } 
 
