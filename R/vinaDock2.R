vinaDock2 = function(chemTXT,chemtype='smiles',pdb,pdbtype='id',configfileTYPE='file',configfile='966c_BOXconfig.txt',center_x=NULL,center_y=NULL,center_z=NULL,size_x=NULL,size_y=NULL,size_z=NULL,outPDBQT='DockedPDBQT.pdbqt',python2PTAH='c:/Progra~2/MGLTools-1.5.6rc3/python2',prepare_ligand4PATH='E:/ADtools',prepare_receptor4PATH='E:/ADtools',getLigand=T,RefLigandFile='Ligand.pdb',IFProLigandPrepared=T,PSOvina=F){ 
  ###chemtype:smiles;inchikey;mol2(fileName) 
  ###pdbtype:id(fetch from web) or not id(PdbFileName) 
  ###configfileTYPE:file-->Filetxt;input-->center_x,center_y,center_z,size_x,size_y,size_z;getFromPDB-->getboxFromPDB 
  ####output:result--dataFrame[mode affinity dist1 dist2];outPDBQTFile[Docked_pdbqt];getLigand=T-->RefLigandFile['Ligand.pdb'] 
  if (file.exists('DockedPDBQT.pdbqt'))file.remove('DockedPDBQT.pdbqt') 
  if (file.exists('TempOutLigandPrepared.pdbqt'))file.remove('TempOutLigandPrepared.pdbqt') 
  if (file.exists('TempOutProteinPrepared.pdbqt'))file.remove('TempOutProteinPrepared.pdbqt') 
  prepareLigand(molfileTXT=chemTXT,moltype=chemtype,preparedMOLfile='TempOutLigandPrepared.pdbqt',python2PTAH=python2PTAH,prepare_ligand4PATH=prepare_ligand4PATH) 
  prepareProtein(pdbFile=pdb,pdbtype=pdbtype,preparedProfile='TempOutProteinPrepared.pdbqt',python2PTAH=python2PTAH,prepare_receptor4PATH=prepare_receptor4PATH,getLigand=getLigand,RefLigandFile=RefLigandFile,IFProLigandPrepared=IFProLigandPrepared) 
  if (configfileTYPE=='getFromPDB'){ 
    aaa=getBoxFromPdb(pdbFile=pdb,outLigand=RefLigandFile,outVINAconfigFile='config.txt',IFprepare=T,python2PTAH=python2PTAH,prepare_ligand4PATH=prepare_ligand4PATH) 
    if (!is.na(aaa)){ 
      result=vinaDock(preparedLigand='TempOutLigandPrepared.pdbqt',preparedProtein='TempOutProteinPrepared.pdbqt',configfileTYPE='input',center_x=aaa['xmean'],center_y=aaa['ymean'],center_z=aaa['zmean'],size_x=aaa['xmax'],size_y=aaa['ymax'],size_z=aaa['zmax'],outPDBQT=outPDBQT,PSO=PSOvina) 
    }else{ 
      result='ERROR:no ligand found in PDB!' 
    } 
  }else{ 
    result=vinaDock(preparedLigand='TempOutLigandPrepared.pdbqt',preparedProtein='TempOutProteinPrepared.pdbqt',configfile=configfile,configfileTYPE=configfileTYPE,center_x=center_x,center_y=center_y,center_z=center_z,size_x=size_x,size_y=size_y,size_z=size_z,outPDBQT=outPDBQT,PSO=PSOvina) 
  } 
  return(result) 
} 
 
