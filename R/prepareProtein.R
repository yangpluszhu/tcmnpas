##' PrepareReceptor for Molecular Docking. 
##' Note: vina, PSOvina,and mgltools must be installed and configured correctly. 
##' 
##' @title prepareProtein 
##' @param pdbFile Receptor file
##' @param pdbtype Receptor type: id pdb or pdbqt 
##' @param preparedProfile prepared file name for output: *.pdbqt
##' @param python2PTAH python2.exe PTAH 
##' @param prepare_receptor4PATH "prepare_receptor4.py" PTAH. The default Vina receptor protein preparation program (addh–>true for adding hydrogen; nphs–>true for merging charges and removing non-polar hydrogen; lps–>true for merging charges and removing lone pairs; waters–>true for removing water molecules; nonstdres–>true for removing chains composed of less than 20 standard amino acid residues) is used.
##' @return a prepared receptor file
##' @export 
##' @author Yang Ming 
prepareProtein = function(pdbFile,pdbtype='id',preparedProfile='outPro.pdbqt',python2PTAH='c:/Progra~2/MGLTools-1.5.6rc3/python2',prepare_receptor4PATH=system.file('tools',package='tcmnpas'),getLigand=T,RefLigandFile='Ligand',IFProLigandPrepared=T,repair=T){ 
  require(bio3d) 
  ##pdbtype:id-->fromWeb,pdb,pdbqt 
  ### 
  if (file.exists('TempProtein.pdb'))file.remove('TempProtein.pdb') 
  if (pdbtype=='id'){ 
  tempPDB=bio3d::read.pdb(pdbFile) 
  indsLigand <- atom.select(tempPDB, c("ligand")) 
  indsWater=atom.select(tempPDB, c("water")) 
  indsALL=atom.select(tempPDB) 
  indsRemove=list(atom=setdiff(setdiff(indsALL$atom,indsLigand$atom),indsWater$atom),xyz=setdiff(setdiff(indsALL$xyz,indsLigand$xyz),indsWater$xyz),call=indsALL$call) 
  RemoveligandWatertempPDB <- trim.pdb(tempPDB, inds=indsRemove) 
  bio3d::write.pdb(RemoveligandWatertempPDB,file='TempProtein.pdb') 
  o0=paste0(prepare_receptor4PATH,'/prepare_receptor4.py') 
    o1='-r TempProtein.pdb' 
    o2='-o' 
    if (repair){ 
      Dos=paste(python2PTAH,o0,o1,o2,preparedProfile,'-A hydrogens') 
    }else{ 
      Dos=paste(python2PTAH,o0,o1,o2,preparedProfile) 
    } 
    system(command = Dos,intern=T) 
    ####### 
    if (getLigand&length(indsLigand$atom)!=0){ 
      #indsMetal=atom.select(tempPDB, c("metal")) 
      #RemoveMetal=list(atom=setdiff(indsLigand$atom,indsMetal$atom),xyz=setdiff(indsLigand$xyz,indsMetal$xyz),call=indsLigand$call) 
      #LigandPdb=trim.pdb(tempPDB, inds=RemoveMetal) 
      LigandPdb=trim.pdb(tempPDB, inds=indsLigand) 
      bio3d::write.pdb(LigandPdb,file=paste0(RefLigandFile,'.pdb')) 
      selectLargeLigand(paste0(RefLigandFile,'.pdb')) 
      if (IFProLigandPrepared){ 
        #prepareLigand(molfileTXT=RefLigandFile,moltype='pdb',preparedMOLfile='TEMPPoutLigand.pdbqt',python2PTAH=python2PTAH,prepare_ligand4PATH=prepare_receptor4PATH) 
        #PDBQTtoPDB(molFile='TEMPPoutLigand.pdbqt',molName=RefLigandFile) 
        prepareLigand(molfileTXT=paste0(RefLigandFile,'.pdb'),moltype='pdb',preparedMOLfile=paste0(RefLigandFile,'.pdbqt'),python2PTAH=python2PTAH,prepare_ligand4PATH=prepare_receptor4PATH) 
      } 
    } 
    ####### 
  }else if (pdbtype=='pdb'){ 
    pdb <- bio3d::read.pdb(pdbFile) 
    indsLigand <- atom.select(pdb, c("ligand")) 
    indsWater=atom.select(pdb, c("water")) 
    indsALL=atom.select(pdb) 
    indsRemove=list(atom=setdiff(setdiff(indsALL$atom,indsLigand$atom),indsWater$atom),xyz=setdiff(setdiff(indsALL$xyz,indsLigand$xyz),indsWater$xyz),call=indsALL$call) 
    RemoveligandWaterpdb <- trim.pdb(pdb, inds=indsRemove) 
    bio3d::write.pdb(RemoveligandWaterpdb,'TempProtein.pdb') 
    o0=paste0(prepare_receptor4PATH,'/prepare_receptor4.py') 
    o1='-r' 
    o2='TempProtein.pdb' 
    o3='-o' 
    if (repair){ 
      Dos=paste(python2PTAH,o0,o1,o2,o3,preparedProfile,'-A hydrogens') 
    }else{ 
      Dos=paste(python2PTAH,o0,o1,o2,o3,preparedProfile) 
    } 
    system(command = Dos,intern=T) 
    ###### 
    if (getLigand&length(indsLigand$atom)!=0){ 
      #indsMetal=atom.select(pdb, c("metal")) 
      #RemoveMetal=list(atom=setdiff(indsLigand$atom,indsMetal$atom),xyz=setdiff(indsLigand$xyz,indsMetal$xyz),call=indsLigand$call) 
      #LigandPdb=trim.pdb(pdb, inds=RemoveMetal) 
      LigandPdb=trim.pdb(pdb, inds=indsLigand) 
      bio3d::write.pdb(LigandPdb,file=paste0(RefLigandFile,'.pdb')) 
      selectLargeLigand(paste0(RefLigandFile,'.pdb')) 
      if (IFProLigandPrepared){ 
        prepareLigand(molfileTXT=paste0(RefLigandFile,'.pdb'),moltype='pdb',preparedMOLfile=paste0(RefLigandFile,'.pdbqt'),python2PTAH=python2PTAH,prepare_ligand4PATH=prepare_receptor4PATH) 
        #PDBQTtoPDB(molFile='TEMPPoutLigand.pdbqt',molName=RefLigandFile) 
      } 
    } 
    ####### 
  }else if (pdbtype=='pdbqt'){ 
    o0=paste0(prepare_receptor4PATH,'/prepare_receptor4.py') 
    o1='-r' 
    o2=pdbFile 
    o3='-o' 
    if (repair){ 
      Dos=paste(python2PTAH,o0,o1,o2,o3,preparedProfile,'-A hydrogens') 
    }else{ 
      Dos=paste(python2PTAH,o0,o1,o2,o3,preparedProfile) 
    } 
    system(command = Dos,intern=T) 
  } 
} 
 
