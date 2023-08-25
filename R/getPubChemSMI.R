getPubChemSMI = function(PubChemID){ 
  ##PubChemID:PubChemID 
  ##output:smiles 
  require(XML) 
  require(stringr) 
  require(RCurl) 
  require(ChemmineR) 
  URL1='https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/' 
  URL2='/record/SDF/?record_type=2d&response_type=save&response_basename=Structure2D_CID_' 
  URL=paste(URL1,PubChemID,URL2,PubChemID,sep='') 
  #URL='https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/100982/record/SDF/?record_type=2d&response_type=save&response_basename=Structure2D_CID_100982' 
  SdfFile=read.SDFset(URL) 
  SMIob=sdf2smiles(SdfFile) 
  SMI=as.character(SMIob) 
  return(SMI) 
} 
 
