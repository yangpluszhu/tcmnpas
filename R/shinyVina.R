shinyVina = function(chemTXT,chemtype='smiles',pdb,pdbtype='id',getBoxtype='FromLigand',center_x=NULL,center_y=NULL,center_z=NULL,size_x=NULL,size_y=NULL,size_z=NULL,outDockedPDBQT='DockedPDBQT',python2PTAHs='c:/Progra~2/MGLTools-1.5.6rc3/python2',prepare_ligand4PATHs='E:/ADtools',prepare_receptor4PATHs='E:/ADtools',PDBboxDbs='../npa/db/PDBboxDB.db',exhaustivenessSet=8,num_modesSet=9,energy_rangeSet=3,cpuSet=2,PSOvina=F){ 
  ###chemtype:smiles;inchikey;mol2(fileName) 
  library(stringr) 
  library(bio3d) 
  library(RSQLite) 
  tmpPDB=pdb 
  temDIR0=tempdir() 
  if (system('pwd')!=0){ 
    #pythonDir='c:/Progra~2/MGLTools-1.5.6rc3/python2' 
    #ADTpath='E:/ADtools' 
    temDIR=unlist(str_split(temDIR0,'\\\\')) 
  }else{ 
    #pythonDir='~/mgltools_x86_64Linux2_1.5.6/bin/pythonsh' 
    #ADTpath='~/ADtools' 
    temDIR=unlist(str_split(temDIR0,'/')) 
  } 
  temDIR=temDIR[length(temDIR)] 
  outDockedResultPDBQT=paste0(temDIR,outDockedPDBQT,'.pdbqt') 
  outDockedResultPDBQTpath=paste0(temDIR0,'/',temDIR,outDockedPDBQT,'.pdbqt') 
  temppreparedMOLfile=paste0(temDIR0,'/',temDIR,'_Ligand.pdbqt') 
  temppreparedProfile=paste0(temDIR0,'/',temDIR,'_Protein.pdbqt') 
  if (file.exists(outDockedResultPDBQT))file.remove(outDockedResultPDBQT) 
  if (file.exists(temppreparedMOLfile))file.remove(temppreparedMOLfile) 
  if (file.exists(temppreparedProfile))file.remove(temppreparedProfile) 
  prepareLigand(molfileTXT=chemTXT,moltype=chemtype,preparedMOLfile=temppreparedMOLfile,python2PTAH=python2PTAHs,prepare_ligand4PATH=prepare_ligand4PATHs) 
  BOXinfo=getPDBbox(PDBid=tmpPDB,type=getBoxtype,center_x=center_x,center_y=center_y,center_z=center_z,size_x=size_x,size_y=size_y,size_z=size_z,python2PTAH=python2PTAHs,prepare_ligand4PATH=prepare_ligand4PATHs,PDBboxDb=PDBboxDbs) 
  prepareProtein(pdbFile=tmpPDB,pdbtype=pdbtype,preparedProfile=temppreparedProfile,python2PTAH=python2PTAHs,prepare_receptor4PATH=prepare_receptor4PATHs,getLigand=F) 
  if (nrow(BOXinfo)>0&file.exists(temppreparedProfile)&file.exists(temppreparedMOLfile)){ 
    if (class(try({DockResult<-vinaDock(preparedLigand=temppreparedMOLfile,preparedProtein=temppreparedProfile,configfileTYPE='input',center_x=BOXinfo$x,center_y=BOXinfo$y,center_z=BOXinfo$z,size_x=BOXinfo$sizex,size_y=BOXinfo$sizey,size_z=BOXinfo$sizez,exhaustiveness=exhaustivenessSet,num_modes=num_modesSet,energy_range=energy_rangeSet,cpu=cpuSet,outPDBQT=outDockedResultPDBQTpath,PSO=PSOvina)},silent = T))=='try-error'){ 
      prepareProtein(pdbFile=tmpPDB,pdbtype=pdbtype,preparedProfile=temppreparedProfile,python2PTAH=python2PTAHs,prepare_receptor4PATH=prepare_receptor4PATHs,getLigand=F,repair = F) 
      DockResult<-vinaDock(preparedLigand=temppreparedMOLfile,preparedProtein=temppreparedProfile,configfileTYPE='input',center_x=BOXinfo$x,center_y=BOXinfo$y,center_z=BOXinfo$z,size_x=BOXinfo$sizex,size_y=BOXinfo$sizey,size_z=BOXinfo$sizez,exhaustiveness=exhaustivenessSet,num_modes=num_modesSet,energy_range=energy_rangeSet,cpu=cpuSet,outPDBQT=outDockedResultPDBQTpath,PSO=PSOvina) 
    } 
    #DockResult=vinaDock(preparedLigand=temppreparedMOLfile,preparedProtein=temppreparedProfile,configfileTYPE='input',center_x=BOXinfo$x,center_y=BOXinfo$y,center_z=BOXinfo$z,size_x=BOXinfo$sizex,size_y=BOXinfo$sizey,size_z=BOXinfo$sizez,exhaustiveness=exhaustivenessSet,num_modes=num_modesSet,energy_range=energy_rangeSet,cpu=cpuSet,outPDBQT=outDockedResultPDBQT,PSO=PSOvina) 
  }else{ 
    DockResult=data.frame() 
  } 
  return(list(DockResult=DockResult,BOXinfo=BOXinfo,outDockedPDBQT=outDockedResultPDBQT,preparedProtein=temppreparedProfile,preparedLigand=temppreparedMOLfile,outDockedResultPDBQTpath=outDockedResultPDBQTpath)) 
} 
 
