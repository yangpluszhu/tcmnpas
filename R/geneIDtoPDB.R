geneIDtoPDB = function(geneID,multi=F,geneTopdbDB='../npa/db/GeneIDtoPDB.db',up='../npa/db/UniProt.ws_up9606.RData'){ 
  require(RSQLite) 
  require(UniProt.ws) 
  require(data.table) 
  options(stringsAsFactors = F) 
  load(up,baseenv()) 
  #load(geneNApdbListdb) 
  con <- dbConnect(SQLite(),dbname=geneTopdbDB) 
  afcOS1='SELECT * FROM geneNApdbListD' 
  eval(substitute({rsOS1<-dbSendQuery(con, afcOS1)},list(afcOS1=afcOS1))) 
  d1afcOS1 <- dbFetch(rsOS1) 
  geneNApdbListdb=as.character(d1afcOS1$value) 
  ###### 
  afcOS2='SELECT * FROM geneHavingPDBListD' 
  eval(substitute({rsOS2<-dbSendQuery(con, afcOS2)},list(afcOS2=afcOS2))) 
  d1afcOS2 <- dbFetch(rsOS2) 
  geneHavingPDBList=as.character(d1afcOS2$value) 
  ###### 
  afcOS3='SELECT * FROM genenoHavingPDBListD' 
  eval(substitute({rsOS3<-dbSendQuery(con, afcOS3)},list(afcOS3=afcOS3))) 
  d1afcOS3<- dbFetch(rsOS3) 
  genenoHavingPDBList=as.character(d1afcOS3$value) 
  ###### 
  afcOS4='SELECT * FROM geneInmypdbListD' 
  eval(substitute({rsOS4<-dbSendQuery(con, afcOS4)},list(afcOS4=afcOS4))) 
  d1afcOS4<- dbFetch(rsOS4) 
  geneInmypdbList=as.character(d1afcOS4$value) 
  ###### 
  dbDisconnect(con) 
  geneID=as.character(geneID) 
  geneID=geneID[!geneID%in%geneNApdbListdb] 
  if (length(geneID)>0){ 
    if (multi==F){ 
    geneIDIN=geneID[geneID%in%geneHavingPDBList] 
    geneIDnoIN=geneID[geneID%in%genenoHavingPDBList] 
    geneIDWeb=setdiff(geneID,c(geneIDIN,geneIDnoIN)) 
    con <- dbConnect(SQLite(),dbname=geneTopdbDB) 
    if (length(geneIDIN)>0){ 
      geneID2=paste("'",geneIDIN,"'",sep='',collapse = ',') 
      geneID3=paste("(",geneID2,")",sep='') 
      os1=paste0('SELECT * FROM havingPDB where geneID in',geneID3) 
      eval(substitute({rs1<-dbSendQuery(con, os1)},list(os1=os1))) 
      d1 <- dbFetch(rs1) 
      d1=d1[,c('geneID','PDB')] 
    }else{ 
      d1=data.frame() 
    } 
    ## 
    if (length(geneIDnoIN)){ 
      geneID4=paste("'",geneIDnoIN,"'",sep='',collapse = ',') 
      geneID5=paste("(",geneID4,")",sep='') 
      os2=paste0('SELECT * FROM nohavingPDB where geneID in',geneID5) 
      eval(substitute({rs2<-dbSendQuery(con, os2)},list(os2=os2))) 
      d2 <- dbFetch(rs2) 
      d2=d2[,c('geneID','PDB')] 
    }else{ 
      d2=data.frame() 
    } 
    ## 
    if (length(geneIDWeb)>0){ 
      if (class(try({res3 <- UniProt.ws::select(up, keys = geneIDWeb,columns = c("PDB"),keytype = "ENTREZ_GENE")},silent = T))!="try-error"){ 
        colnames(res3)=c('geneID','PDB') 
        res3=as.data.table(res3) 
        res3=unique(res3,by='geneID') 
        d3=as.data.frame(res3) 
      }else{ 
        d3=data.frame() 
      } 
    }else{ 
      d3=data.frame() 
    } 
    d_total=rbind(d1,d2,d3) 
    if (nrow(d_total)>0){ 
      result=unique(d_total) 
    }else{ 
      result='ERROR:geneID not found in PDB database!' 
    } 
    }else if (multi){ 
      #######################multi##################### 
      geneIDINList=geneID[geneID%in%geneInmypdbList] 
      geneIDnoINList=geneID[!geneID%in%geneInmypdbList] 
      con <- dbConnect(SQLite(),dbname=geneTopdbDB) 
      if (length(geneIDINList)>0){ 
        geneID22=paste("'",geneIDINList,"'",sep='',collapse = ',') 
        geneID23=paste("(",geneID22,")",sep='') 
        os21=paste0('SELECT * FROM geneTopdb where geneID in',geneID23) 
        eval(substitute({rs1<-dbSendQuery(con, os21)},list(os21=os21))) 
        d21 <- dbFetch(rs1) 
        d21=d21[,c('geneID','PDB')] 
      }else{ 
        d21=data.frame() 
      } 
      ############### 
        if (length(geneIDnoINList)>0){ 
        if (class(try({res23 <- UniProt.ws::select(up, keys = geneIDnoINList,columns = c("PDB"),keytype = "ENTREZ_GENE")},silent = T))!="try-error"){ 
          colnames(res23)=c('geneID','PDB') 
          res23=as.data.table(res23) 
          res23=unique(res23,by='geneID') 
          d23=as.data.frame(res23) 
        }else{ 
          d23=data.frame() 
        } 
        }else{ 
          d23=data.frame() 
      } 
      ################################################ 
      d_total=rbind(d21,d23) 
      if (nrow(d_total)>0){ 
        result=unique(d_total) 
      }else{ 
        result='ERROR:geneID not found in PDB database!' 
      } 
    } 
    dbDisconnect(con) 
  }else{ 
    result='ERROR:geneID not found in PDB database!' 
  } 
  return(result) 
} 
 
