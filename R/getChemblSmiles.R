getChemblSmiles = function(chembl,type){ 
  #chembl:chemblID('CHEMBL150') or inchikey 
  #type:'ID' or 'KEY' 
  #output:c('Name','chemblID','smiles','InChiKey') 
  require(chemblr) 
  result=character() 
  length(result)=5 
  names(result)=c('Name','chemblID','smiles','InChiKey','InputName') 
  result['InputName']=chembl 
  if (type=='ID'){ 
    data=get.compounds(chembl,type='chemblid') 
    if (ncol(data)>1){ 
      result['Name']=data$synonyms 
      result['chemblID']=data$chemblId 
      result['smiles']=data$smiles 
      result['InChiKey']=data$stdInChiKey 
      #result['InputName']=chembl 
    } 
  }else if (type=='KEY'){ 
    data=get.compounds(chembl,type='stdinchi') 
    if (ncol(data)>1){ 
      result['Name']=data$synonyms 
      result['chemblID']=data$chemblId 
      result['smiles']=data$smiles 
      result['InChiKey']=data$stdInChiKey 
      #result['InputName']=chembl 
    } 
  } 
  return(result) 
} 
 
