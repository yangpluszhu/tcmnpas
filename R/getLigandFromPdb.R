getLigandFromPdb = function(pdbFile,outLigand='RefoutLigand',IFprepare=T,python2PTAH='c:/Progra~2/MGLTools-1.5.6rc3/python2',prepare_ligand4PATH='E:/ADtools'){ 
  require(bio3d) 
  ###pdbFile:pdb or pdb id 
  ###IFprepare: if output parpared ligand(*.pdb) 
  ###if no ligand then output NA 
  ####### 
  temDIR=tempdir() 
  if (system('pwd')!=0){ 
    #pythonDir='c:/Progra~2/MGLTools-1.5.6rc3/python2' 
    #ADTpath='E:/ADtools' 
    temDIR=unlist(str_split(temDIR,'\\\\')) 
  }else{ 
    #pythonDir='~/mgltools_x86_64Linux2_1.5.6/bin/pythonsh' 
    #ADTpath='~/ADtools' 
    temDIR=unlist(str_split(temDIR,'/')) 
  } 
  temDIR=temDIR[length(temDIR)] 
  tempoutLigand=paste0(temDIR,outLigand,'.pdb') 
  ####### 
  tempPDB=bio3d::read.pdb(pdbFile) 
  indsLigand <- atom.select(tempPDB, c("ligand")) 
  #indsMetal=atom.select(tempPDB, c("metal")) 
  if (length(indsLigand$atom)>2){ 
    #indsMetal=atom.select(tempPDB, c("metal")) 
    #RemoveMetal=list(atom=setdiff(indsLigand$atom,indsMetal$atom),xyz=setdiff(indsLigand$xyz,indsMetal$xyz),call=indsLigand$call) 
    #LigandPdb=trim.pdb(tempPDB, inds=RemoveMetal) 
    LigandPdb=trim.pdb(tempPDB, inds=indsLigand) 
    bio3d::write.pdb(LigandPdb,file=tempoutLigand) 
    selectLargeLigand(tempoutLigand) 
    if (class(try({bio3d::read.pdb(tempoutLigand)},silent = T))!='try-error'){ 
      if (IFprepare){ 
        if (file.exists('TEMPPPPLoutLigand.pdbqt'))file.remove('TEMPPPPLoutLigand.pdbqt') 
        prepareLigand(molfileTXT=tempoutLigand,moltype='pdb',preparedMOLfile='TEMPPPPLoutLigand.pdbqt',python2PTAH=python2PTAH,prepare_ligand4PATH=prepare_ligand4PATH) 
        PDBQTtoPDB(molFile='TEMPPPPLoutLigand.pdbqt',molName=tempoutLigand) 
      } 
      return(tempoutLigand) 
    }else{ 
      return(NA) 
    } 
  }else{ 
    return(NA) 
  } 
} 
 
