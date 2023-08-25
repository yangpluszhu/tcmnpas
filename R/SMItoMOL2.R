SMItoMOL2 = function(moltxt,molName='mole.mol2'){ 
  library(stringr) 
  library(rcdk) 
  m <- parse.smiles(moltxt)[[1]] 
  moltxt2=get.smiles(m) 
  ##moltxt:smiles 
  os1='obabel -:"' 
  os2=moltxt2 
  os3=paste0(os1,os2,'"') 
  os4='-omol2 -O' 
  os5=molName 
  os6='--gen3d' 
  Dos=paste(os3,os4,os5,os6) 
  system(command = Dos,intern=T) 
} 
 
