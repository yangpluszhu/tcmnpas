MOL2toSMILES = function(molFile,molName='mole.smi'){ 
  library(stringr) 
  ##molFile:SDFfileName 
  os1='obabel -imol2' 
  os2=molFile 
  os3='-osmiles' 
  os4='-O' 
  os5=molName 
  #os6='--gen3d' 
  Dos=paste(os1,os2,os3,os4,os5) 
  system(command = Dos,intern=T) 
} 
 
