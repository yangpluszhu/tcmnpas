##' SDF to SMILES
##' note:obabel must be installed
##' @title SDFtoSMILES
##' @param SDFfile molecular file *.sdf
##' @param outtype 'txt' or 'file'
##' @param outName a csv file with transformed SMILES character
##' @return a character vector or a dataframe
##' @export
##' @author Yang Ming
SDFtoSMILES = function(SDFfile,outtype='txt',outName='smiles.csv'){ 
  require(ChemmineR) 
  sdfdata=read.SDFset(SDFfile) 
  valid <- validSDF(sdfdata) 
  sdfdata <- sdfdata[valid] 
  NumChem=length(sdfdata) 
  if (class(try(sdfdata@ID,silent = T))=='try-error'){ 
    ChemID<-paste0('Chem',1:length(NumChem)) 
  }else{ 
    ChemID=sdfdata@ID 
  } 
    Tempsmiles <- as.character(sdf2smiles(sdfdata)) 
    if (outtype=='txt'){ 
      smiles=Tempsmiles 
      names(smiles)=ChemID 
    }else{ 
      smiles=data.frame(id=ChemID,smiles=Tempsmiles) 
      write.csv(smiles,outName,row.names = F) 
    } 
  return(smiles) 
} 
 
