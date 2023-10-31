getHerbGene = function(herbName,QEDset=0.2,geneScore=800,geneSelectPv=0.01,targetDatabase=c('HIT','TCMID','STITCH','TCMSP','CUSTOM'),Lipinski_Viosset=0,Veber_Viosset=0){
  require(data.table)
  require(stringr)
  require(RSQLite)
  require(plyr)
  require(ggplot2)
  require(pracma)
  if (!exists('AllHerbListData')){
    #ALL=fread('db/AllHerbListData.csv')
    AllHerbListData=fread('db/AllHerbListData.csv',encoding='UTF-8')
  }
  herbList=herbName
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
  }else{
    selectGeneID=character()
  }
  if (length(selectGeneID)<1)selectGeneID='NA'
  return(selectGeneID)
}
 
