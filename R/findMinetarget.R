findMinetarget = function(Name,NameType){
  #Name
  #NameType:chemicalName(CN) or herbChineseName(HN)
  #MineCTD:'herb','chemical_name','geneID','smi'
  require(data.table)
  require(stringr)
  require(plyr)
  if (!exists('MineCTD')){
    MineCTD=fread('db/Mine_herb_chem_gene_smi.csv',colClasses = 'character')
  }
 if (NameType=='CN'){
   data=MineCTD[chemical_name%like%Name,.(chemical_name,geneID,smiles)]
   data$InputName=Name
   return(data)
 }else if (NameType=='HN'){
   data=MineCTD[herb%like%Name,.(herb,chemical_name,geneID,smiles)]
   data$InputName=Name
   return(data)
 }else{
   print('ERROR:NameType not found!')
 }
}

