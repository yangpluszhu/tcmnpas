findTCMSPChemFromHerb = function(HerbName){ 
    #HerbName:HerbName-->chinese 
    #TCMSPdatabase:TCMSPherbChem and TCMSPChemTarget 
    require(data.table) 
    require(stringr) 
    if (!exists('TCMSPherbChem')){ 
      load('db/TCMSPdatabase.RData') 
    } 
    data=TCMSPherbChem[herb_name%like%HerbName,.(herb_name,MOL_name,inchikey,smiles)] 
    data=unique(data) 
    data$InputName=HerbName 
    return(data) 
} 
 
