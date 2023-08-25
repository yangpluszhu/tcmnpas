getHerbFromdatabase = function(ChemData,ChemDataType='inchikey',targetDatabase=c('HIT','TCMID','STITCH','TCMSP','CUSTOM'),IFFromtargetData=F,AllHerbListData=AllHerbListData){ 
  #ChemData:data.frame:colnames-->'cid','chemical_name','key'.'key' is smiles or inchikey 
  #IFFromtargetData:F-->From db/Tar.db. 
  #IFFromtargetData:T-->From herb_ChemDataAll 
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
  ChemData2=ChemData 
  if (ChemDataType!='inchikey'){ 
    tempSmiles=ChemData$key 
    names(tempSmiles)=ChemData$cid 
    ChemData2$inchikey=unlist(TransInchikey(tempSmiles,molecularType='SMItxt')) 
  }else{ 
    ChemData2=plyr::rename(ChemData2,c('key'='inchikey')) 
  } 
  if (!exists('AllHerbListData')){ 
    #ALL=fread('db/AllHerbListData.csv') 
    load('db/AllHerbListData.RData') 
  } 
  ############################ 
  ChemData2=ChemData2[!is.na(ChemData2$inchikey),] 
  ChemData2=as.data.table(ChemData2) 
  ChemData2=unique(ChemData2,by='inchikey') 
  QueryKey=unique(ChemData2$inchikey) 
if (!IFFromtargetData){ 
  if (nrow(ChemData2)>0){ 
    paID2=paste("'",QueryKey,"'",sep='',collapse = ',') 
    paID3=paste("(",paID2,")",sep='') 
    chemIDOrder=paste("inchikey IN ",paID3,sep='') 
    OrderchemIDOrder=paste("SELECT * FROM STITCH_HIT_TCMID_TCMSP_Tar where",chemIDOrder,sep=' ') 
    tmp <- dbConnect(SQLite(), 'db/Tar.db') 
    res <- dbSendQuery(tmp, OrderchemIDOrder) 
    subD0 <- fetch(res, n =-1) 
    #dbDisconnect(tmp) 
    ########## 
    subD0_herbID=unique(subD0$herbID) 
    paID2=paste("'",subD0_herbID,"'",sep='',collapse = ',') 
    paID3=paste("(",paID2,")",sep='') 
    herbIDOrder=paste("herbID IN ",paID3,sep='') 
    OrderherbIDOrder=paste("SELECT * FROM STITCH_HIT_TCMID_TCMSP_Tar where",herbIDOrder,sep=' ') 
    #tmp <- dbConnect(SQLite(), 'db/Tar.db') 
    res_herbID <- dbSendQuery(tmp, OrderherbIDOrder) 
    subD0_herb <- fetch(res_herbID, n =-1) 
    dbDisconnect(tmp) 
    ##################### 
    subD0_herb=join(subD0_herb,AllHerbListData,by='herbID') 
    subD0_herb=as.data.table(subD0_herb) 
    subD0_herb=subD0_herb[,.(TotalChem=length(unique(inchikey))),by='herb'] 
    ##################### 
    subD=join(subD0[,c('herbID','inchikey')],ChemData2,by='inchikey') 
    subD=join(subD,AllHerbListData,by='herbID') 
    subD=as.data.table(subD) 
    subD=subD[database%in%targetDatabase,] 
    subD=subD[,.(cid,chemical_name,herb,inchikey,QED_DES)] 
    subD=unique(subD,by=c('herb','inchikey')) 
    rm(subD0) 
    gc() 
    subDSummary=subD[,.(chemical=toString(unique(cid)),NumChem=length(unique(cid))),by='herb'] 
    subDSummary=join(subDSummary,subD0_herb,by='herb') 
    subDSummary$BackgroundRatio=subDSummary$NumChem/subDSummary$TotalChem 
    subDSummary$QueryRatio=subDSummary$NumChem/length(QueryKey) 
    subDSummary=subDSummary[order(-NumChem)] 
  }else{ 
    subD=NA 
    subDSummary=NA 
  } 
}else{ 
  if (nrow(ChemData2)>0){ 
    subD0=herb_ChemDataAll[herb_ChemDataAll$inchikey%in%QueryKey&herb_ChemDataAll$database%in%targetDatabase,] 
    subD0_herb=unique(subD0$herb) 
    subD0_herbData=herb_ChemDataAll[herb_ChemDataAll$herb%in%subD0_herb,] 
    subD0_herbData=as.data.table(subD0_herbData) 
    subD0_herbData=subD0_herbData[,.(TotalChem=length(unique(inchikey))),by='herb'] 
    subD0=join(subD0,ChemData2,by='inchikey') 
    subD=as.data.table(subD0) 
    subD=subD[,.(cid,chemical_name,herb,inchikey,QED_DES)] 
    subD=unique(subD,by=c('herb','inchikey')) 
    subDSummary=subD[,.(chemical=toString(unique(cid)),NumChem=length(unique(cid))),by='herb'] 
    subDSummary=join(subDSummary,subD0_herbData,by='herb') 
    subDSummary$BackgroundRatio=subDSummary$NumChem/subDSummary$TotalChem 
    subDSummary$QueryRatio=subDSummary$NumChem/length(QueryKey) 
    subDSummary=subDSummary[order(-NumChem)] 
  }else{ 
    subD=NA 
    subDSummary=NA 
  } 
} 
  return(list(subD=subD,subDSummary=subDSummary)) 
} 
 
