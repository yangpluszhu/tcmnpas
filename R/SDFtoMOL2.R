##' SDF to MOL2
##' note:obabel must be installed
##' @title SDFtoMOL2
##' @param molFile molecular file *.sdf
##' @param molName molecular file name
##' @return a mol2 file
##' @export
##' @author Yang Ming
SDFtoMOL2 = function(molFile,molName='mole.mol2'){ 
  library(stringr) 
  ##molFile:SDFfileName 
  os1='obabel -isdf' 
  os2=molFile 
  os3='-omol2' 
  os4='-O' 
  os5=molName 
  os6='--gen3d' 
  Dos=paste(os1,os2,os3,os4,os5,os6) 
  system(command = Dos,intern=T) 
} 
 
