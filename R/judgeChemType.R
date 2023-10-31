judgeChemType = function(ChemTXT){ 
  require(stringr) 
  temp=ChemTXT 
  if (nchar(temp)==27&nchar(unlist(str_split(temp,'-')))[1]==14){ 
    ChemType='inchikey' 
  }else{ 
    ChemType='smiles' 
  } 
  return(ChemType) 
} 
 
