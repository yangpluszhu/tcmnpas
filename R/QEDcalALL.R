QEDcalALL = function(druglikeDesFile,sdfFile){ 
  require(ChemmineR) 
  require(ChemmineOB) 
  require(data.table) 
  require(plyr) 
  require(rio) 
  options(stringsAsFactors = F) 
  druglikeData=fread(druglikeDesFile) 
  druglikeData=as.data.frame(druglikeData) 
  RenamesChar=c(ALogP='ALOGP',Molecular_Mass='MW',Num_RotatableBonds='ROTB',Num_AromaticRings='AROM',Num_H_Acceptors_Lipinski='HBA',Num_H_Donors_Lipinski='HBD',Molecular_PolarSurfaceArea='PSA') 
  druglikeData=plyr::rename(druglikeData,RenamesChar) 
  druglikeData$id=as.character(druglikeData$id) 
  druglikeData$ALOGP=as.numeric(druglikeData$ALOGP) 
  druglikeData$MW=as.numeric(druglikeData$MW) 
  druglikeData$ROTB=as.numeric(druglikeData$ROTB) 
  druglikeData$AROM=as.numeric(druglikeData$AROM) 
  druglikeData$HBA=as.numeric(druglikeData$HBA) 
  druglikeData$HBD=as.numeric(druglikeData$HBD) 
  druglikeData$PSA=as.numeric(druglikeData$PSA) 
#ALogP-->ALOGP 
#Molecular_Mass-->MW 
#Num_RotatableBonds-->ROTB 
#Num_AromaticRings-->AROM 
#Num_H_Acceptors_Lipinski-->HBA 
#Num_H_Donors_Lipinski-->HBD 
#Molecular_PolarSurfaceArea-->PSA 
  alert=alertCal(sdfFile) 
  druglikeData=join(druglikeData,alert,by='id') 
  druglike_des=QEDcal(druglikeData) 
  return(druglike_des) 
} 
 
