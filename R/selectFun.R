selectFun = function(protein,up=up){ 
  require(UniProt.ws) 
  if (!exists('up')){ 
    load('db/UniProtwsOB.RData') 
  } 
  columns='ENTREZ_GENE' 
  kt='ENSEMBL_PROTEIN' 
  if (class(try({b=select(up,protein,columns,kt)},silent = T))!='try-error'){ 
    a=b 
    return(a$ENTREZ_GENE[1]) 
  }else{ 
    a='NA' 
    return(a) 
  } 
  #a=select(up,protein,columns,kt) 
 
} 
 
