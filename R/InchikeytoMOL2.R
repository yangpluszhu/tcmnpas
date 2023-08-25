InchikeytoMOL2 = function(Inchikey,MOLdataDB='../npa/db/inchikey_smi_db.Rdata'){ 
  library(RCurl) 
  library(rvest) 
  library(data.table) 
  options(stringsAsFactors = F) 
  #MOLdata=fread('E:/DATABASE/MYpro/BanXiaXieXin/result/ToolBox/db/2020-01-15_16_27_34_STITCH_HIT_TCMID_TCMSP_Tar.csv') 
  #MOLdata=unique(MOLdata[,.(inchikey,smiles)]) 
  #save(lsit='MOLdata',file='inchikey_smi_db.Rdata') 
  load(MOLdataDB) 
  smi=MOLdata$smiles[chmatch(Inchikey,MOLdata$inchikey)] 
  return(smi) 
} 
 
