queryHerbChem = function(herbList=NULL,Qinchikey=NULL,database=c('HIT','TCMSP','TCMID'),HerbChemDB='../npa/db/herb_ChemDataAll.db'){ 
  require(stringr) 
  require(RSQLite) 
  if (!is.null(herbList)&&!is.na(herbList)){ 
    herbList2=paste("'",herbList,"'",sep='',collapse = ',') 
    herbList3=paste("(",herbList2,")",sep='') 
    herbOrder00=paste('herbID IN',herbList3,sep=' ') 
  }else{ 
    herbOrder00=NULL 
  } 
  ################### 
  if (!is.null(Qinchikey)&&!is.na(Qinchikey)){ 
    inchikey2=paste("'",Qinchikey,"'",sep='',collapse = ',') 
    inchikey3=paste("(",inchikey2,")",sep='') 
    QinchikeyOrder00=paste('inchikey IN',inchikey3,sep=' ') 
  }else{ 
    QinchikeyOrder00=NULL 
  } 
  ################### 
  databaseList2=paste("'",database,"'",sep='',collapse = ',') 
  databaseList3=paste("(",databaseList2,")",sep='') 
  databaseOrder00=paste('database IN',databaseList3,sep=' ') 
  ################### 
  if (!is.null(herbOrder00)|!is.null(QinchikeyOrder00)){ 
    OrderTemp=c(herbOrder00,QinchikeyOrder00,databaseOrder00) 
    OrderTemp=OrderTemp[OrderTemp!=''] 
    if (length(OrderTemp)==1){ 
      Order=paste("SELECT * FROM herb_ChemDataAll where",OrderTemp,sep=' ') 
    }else{ 
      OrderTemp2=paste(OrderTemp,collapse = ' AND ') 
      Order=paste("SELECT * FROM herb_ChemDataAll where",OrderTemp2,sep=' ') 
    } 
  }else{ 
    Order=NULL 
  } 
  ################### 
  if (!is.null(Order)){ 
    con <- dbConnect(SQLite(),dbname=HerbChemDB) 
    res <- dbSendQuery(con, Order) 
    result <- fetch(res, n =-1) 
    dbDisconnect(con) 
  }else{ 
    result=data.frame() 
  } 
  return(result) 
} 
 
