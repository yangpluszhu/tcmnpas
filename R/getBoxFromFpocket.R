getBoxFromFpocket = function(pdbID,offset=5){ 
  require(bio3d) 
  require(BiocGenerics) 
  require(stringr) 
  if (!str_detect(pdbID,'\\.')){ 
    a=bio3d::read.pdb(pdbID) 
    pdbName=paste0(pdbID,'.pdb') 
    if (file.exists(pdbName))file.remove(pdbName) 
    bio3d::write.pdb(a,pdbName) 
    pdbIDOK=pdbID 
    os1='fpocket -f' 
    os2=pdbName 
    system(command = paste(os1,os2),intern=T) 
    pqrFileLocat=paste0(pdbIDOK,'_out/',pdbIDOK,'_pockets.pqr') 
    result=ExtrFpocketvalue(pqrFile=pqrFileLocat,offsetvalue=offset) 
    file.remove(pdbName) 
    unlink(paste0(pdbIDOK,'_out'), recursive = T) 
  }else{ 
    pdbFIle=basename(pdbID) 
    pdbIDOK=unlist(str_split(pdbFIle,'\\.'))[1] 
    pdbFIlePATH=str_replace(pdbID,pdbFIle,'') 
    os1='fpocket -f' 
    os2=pdbID 
    system(command = paste(os1,os2),intern=T) 
    pqrFileLocat=paste0(pdbFIlePATH,pdbIDOK,'_out/',pdbIDOK,'_pockets.pqr') 
    result=ExtrFpocketvalue(pqrFile=pqrFileLocat,offsetvalue=offset) 
    #file.remove(pdbID) 
    unlink(paste0(pdbFIlePATH,pdbIDOK,'_out'), recursive = T) 
  } 
  return(result) 
} 
 
