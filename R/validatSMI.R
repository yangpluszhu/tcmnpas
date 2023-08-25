validatSMI = function(smiTXT){ 
  library(stringr) 
  tempOutName='TEMPvalidateSMI.mol2' 
  validat=logical() 
  for (i in 1:length(smiTXT)){ 
    print(i) 
    if (file.exists(tempOutName))file.remove(tempOutName) 
    SMItoMOL2(moltxt=smiTXT[i],molName=tempOutName) 
    if (file.exists(tempOutName)){ 
      if (file.info(tempOutName)['size']!=0){ 
        tempV=T 
      }else{ 
        tempV=F 
      } 
    }else{ 
      tempV=F 
    } 
    validat=c(validat,tempV) 
  } 
  return(validat) 
} 
 
