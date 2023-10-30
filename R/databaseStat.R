databaseStat = function(databasefile=NULL){
  require(RSQLite)
  require(lubridate)
  require(stringr)
  require(data.table)
  require(foreign)
  path='prescription'
  if (is.null(databasefile)){
    databaseName=paste(path,'database.db',sep='/')
  }else{
    databaseName=databasefile
  }
  ##########
  tmp <- dbConnect(SQLite(), databaseName)
  res1=dbSendQuery(tmp,'SELECT COUNT(DISTINCT drug) FROM prescription')
  a1=fetch(res1, n =-1)
  NumYao=a1[1,1]
  ##
  res2=dbSendQuery(tmp,'SELECT COUNT(DISTINCT id) FROM prescription')
  a2=fetch(res2, n =-1)
  NumF=a2[1,1]
  ###
  res3=dbSendQuery(tmp,'SELECT COUNT(DISTINCT doctor) FROM prescription')
  a3=fetch(res3, n =-1)
  NumDoctor=a3[1,1]
  ###
  res4=dbSendQuery(tmp,'SELECT COUNT(DISTINCT patient) FROM prescription')
  a4=fetch(res4, n =-1)
  Numpatient=a4[1,1]
  result=data.frame(处方数=NumF,药物数=NumYao,患者数=Numpatient,开方医生部门数=NumDoctor)
  return(result)
}

