PDBQTtoPDB = function(molFile,molName='mole.pdb'){ 
  library(stringr) 
  ##molFile:PDBQTfileName 
  os1='obabel -ipdbqt' 
  os2=molFile 
  os3='-opdb' 
  os4='-O' 
  os5=molName 
  Dos=paste(os1,os2,os3,os4,os5) 
  system(command = Dos,intern=T) 
} 
 
