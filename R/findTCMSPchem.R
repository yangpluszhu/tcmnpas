findTCMSPchem = function(ChemName){
  #ChemName:ChemName-->Enlish
  #TCMSPdatabase:TCMSPherbChem and TCMSPChemTarget
  require(data.table)
  require(stringr)
  if (!exists('TCMSPherbChem')){
    load('db/TCMSPdatabase.RData')
  }
 data=TCMSPherbChem[MOL_name%like%ChemName,.(MOL_ID,MOL_name,inchikey,smiles)]
 data=unique(data)
 data$InputName=ChemName
 return(data)
}

