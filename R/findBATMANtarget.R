findBATMANtarget = function(pubMedid){
  require(data.table)
  require(stringr)
  #BATMANdatabase:BATMANChemTarget--->colnames:c('Pubchem_CID','geneID')
  if (!exists('BATMANChemTarget')){
    load('db/BATMANdatabase.RData')
  }
  data=BATMANChemTarget[Pubchem_CID%in%pubMedid,.(Pubchem_CID,geneID)]
  data$InputName=pubMedid
  return(data)
}

