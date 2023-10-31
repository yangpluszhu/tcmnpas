findStitchTarget2 = function(chem_inchikey,scoreSet=0,STITCH_CTdb='db/STITCH5_CTdb.db'){ 
  #chem_inchikey:inchikey 
  #STITCH_CTdb:'STITCH5_CTdb.db',colnames=c('cid','ENSP','score','geneID','inchikey') 
  #STITCH_CKEYdb:'STITCH_CKEYdb.db',colnames=c('inchikey','source_cid','stereo_chemical_id','flat_chemical_id') 
  require(RSQLite) 
  require(stringr) 
  require(data.table) 
  require(foreign) 
  tmp <- dbConnect(SQLite(), STITCH_CTdb) 
  Order0=paste('inchikey IN ',"('",chem_inchikey,"')",sep='') 
  Order1=paste('select * from STITCHdb where',Order0,sep=' ') 
  res1 <- dbSendQuery(tmp, Order1) 
  data1=fetch(res1) 
  dbDisconnect(tmp) 
  if (nrow(data1)!=0){ 
    data2=as.data.table(data1) 
    data2=unique(data2) 
    if (nrow(data2)!=0){ 
      data2=data2[score>=scoreSet,] 
    } 
	result=data2 
    result$inchikey=chem_inchikey     
  }else{ 
    result=data.table()    
  }   
  return(result) 
} 
 
