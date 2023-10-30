findCTDtarget = function(ChemName){
  require(data.table)
  require(stringr)
  #CTDdatabase:CTDchemTarget--->colnames:c('ChemicalName','CTDid','geneName','geneID')
  if (!exists('CTDchemTarget')){
    load('db/CTDdatabase.RData')
  }
  data=CTDchemTarget[ChemicalName%like%ChemName,.(ChemicalName,CTDid,geneName,geneID,score)]
  data=unique(data)
  data$InputName=ChemName
  return(data)
}

