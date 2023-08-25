##' Molecular Docking. 
##' Note: vina, PSOvina,and mgltools must be installed and configured correctly. 
##' 
##' @title vinaDock 
##' @param preparedLigand prepared Ligand 
##' @param preparedProtein prepared Protein 
##' @param configfile a txt file 
##' @param configfileTYPE "file" OR "input" 
##' @param center_x input parameter 
##' @param center_y input parameter 
##' @param center_z input parameter 
##' @param size_x input parameter 
##' @param size_y input parameter 
##' @param size_z input parameter 
##' @param exhaustiveness input parameter 
##' @param num_modes input parameter 
##' @param energy_range input parameter 
##' @param cpu input parameter 
##' @param outPDBQT output docked file name 
##' @param PSO whether to use PSOvina 
##' @return a list Object 
##' @export 
##' @author Yang Ming 
vinaDock = function(preparedLigand,preparedProtein,configfile='966c_BOXconfig.txt',configfileTYPE='file',center_x=NULL,center_y=NULL,center_z=NULL,size_x=NULL,size_y=NULL,size_z=NULL,exhaustiveness=8,num_modes=9,energy_range=3,cpu=2,outPDBQT='DockedPDBQT.pdbqt',PSO=FALSE){ 
  ####configfileTYPE:file-->Filetxt;input-->center_x,center_y,center_z,size_x,size_y,size_z 
  ####output:dockResult(dataFrame)colname-->c('mode','affinity','dist1','dist2') 
  ####outPDBQT:Docked File 
  if (file.exists('tempDock.log'))file.remove('tempDock.log') 
  if (PSO){ 
    o1='~/psovina2ls/build/linux/release/psovina2ls --ligand' 
  }else{ 
    o1='vina --ligand' 
  } 
  o2=preparedLigand 
  o3='--receptor' 
  o4=preparedProtein 
  o5='--config' 
  o6=configfile 
  o7='--out' 
  o8=outPDBQT 
  o9='--log tempDock.log' 
  oset1='--center_x' 
  oset2='--center_y' 
  oset3='--center_z' 
  oset4='--size_x' 
  oset5='--size_y' 
  oset6='--size_z' 
  if (configfileTYPE=='file'){ 
    Dos=paste(o1,o2,o3,o4,o5,o6,o7,o8,o9,paste0('--exhaustiveness ',exhaustiveness),paste0('--num_modes ',num_modes),paste0('--energy_range ',energy_range),paste0('--cpu ',cpu),paste0('--seed 123456')) 
  }else if (configfileTYPE=='input'){ 
    Dos=paste(o1,o2,o3,o4,oset1,center_x,oset2,center_y,oset3,center_z,oset4,size_x,oset5,size_y,oset6,size_z,o7,o8,o9,paste0('--exhaustiveness ',exhaustiveness),paste0('--num_modes ',num_modes),paste0('--energy_range ',energy_range),paste0('--cpu ',cpu),paste0('--seed 123456')) 
  } 
  system(command = Dos,intern=T) 
  dockResult=ExtrVinaLogs(LogFile='tempDock.log') 
  return(dockResult) 
} 
 
