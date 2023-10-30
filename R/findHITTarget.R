findHITTarget = function(Name,NameType){
  #Name
  #NameType:chemicalName(CN) or herbChineseName(HCN) or herbPingYinName(HPN) or Chem_pubmedID(PID)
  #HITdatabase:HitChemTarget and HitHerbChem
  require(data.table)
  require(stringr)
  require(plyr)
  if (!exists('HitHerbChem')){
    load('db/HITdatabase.RData')
  }
  if (NameType=='CN'){
    CNidD=HitHerbChem[name%like%Name,.(HITid=cid,name)]
    CNidD=unique(CNidD)
    CNid=unique(CNidD$HITid)
    CT=HitChemTarget[HITid%in%CNid,.(HITid,geneName,geneID)]
    if (nrow(CT)!=0){
      CT=join(CT,CNidD,by='HITid')
      CT$InputName=Name
      }
    return(CT)
  }else if (NameType=='HCN'){
    HCNid=HitHerbChem[herbCH%like%Name,.(cid,name,herbCH)]
    if (nrow(HCNid)!=0){
      HCNCid=unique(HCNid$cid)
      CT=HitChemTarget[HITid%in%HCNCid,.(cid=HITid,geneName,geneID)]
      CT=as.data.frame(CT)
      HCNCid=as.data.frame(HCNCid)
      CT2=join(CT,HCNid,by='cid')
      CT2$InputName=Name
    }else{
      CT2=data.table()
    }
    return(CT2)
  }else if (NameType=='HPN'){
    HCNid=HitHerbChem[herbPY==Name,.(cid,name,herbPY)]
    if (nrow(HCNid)!=0){
      HCNCid=unique(HCNid$cid)
      CT=HitChemTarget[HITid%in%HCNCid,.(cid=HITid,geneName,geneID)]
      CT=as.data.frame(CT)
      HCNCid=as.data.frame(HCNCid)
      CT2=join(CT,HCNid,by='cid')
      CT2$InputName=Name
    }else{
      CT2=data.table()
    }
    return(CT2)
  }else if (NameType=='PID'){
    Name2=paste('CID:',Name,sep='')
    temp=HitHerbChem[,.(HITid=cid,name)]
    CT=HitChemTarget[pub_id%in%Name2,.(HITid,PID=pub_id,geneName,geneID)]
    if (nrow(CT)!=0){CT=join(CT,temp,by='HITid');CT$InputName=Name}
    return(CT)
  }
}

