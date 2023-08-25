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
 
