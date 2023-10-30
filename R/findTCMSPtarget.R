findTCMSPtarget = function(Name,NameType){
  #Name
  #NameType:chemicalName(CN) or herbChineseName(HN)
  #TCMSPdatabase:TCMSPherbChem and TCMSPChemTarget
  require(data.table)
  require(stringr)
  require(plyr)
  if (!exists('TCMSPherbChem')){
    load('db/TCMSPdatabase.RData')
  }
  if (NameType=='CN'){
    MOLidD=TCMSPherbChem[MOL_name%like%Name,.(MOL_name,MOL_ID)]
    if (nrow(MOLidD)!=0){
      MOLidD=unique(MOLidD)
      MOLid=unique(MOLid$MOL_ID)
      CT=TCMSPChemTarget[MOL_ID%in%MOLid,.(MOL_ID,geneName,geneID)]
      if (nrow(CT)!=0){
        CT2=join(CT,MOLidD,by='MOL_ID')
        CT2$InputName=Name
      }else{
        CT2=data.table()
      }
    }else{
      CT2=data.table()
    }
    return(CT2)
  }else if (NameType=='HN') {
    TargetD=TCMSPChemTarget[herb_name%like%Name,.(herb_name,MOL_ID,geneName,geneID)]
    temp=TCMSPherbChem[,.(MOL_ID,MOL_name)]
    if (nrow(TargetD)!=0){
      CT=join(TargetD,temp,by='MOL_ID')
      CT$InputName=Name
    }else{
      CT=data.table()
    }
    return(CT)
  }else{
    print('ERROR:NameType not found!')
  }
}

