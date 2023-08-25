geneDockNoRefligand = function(chemTXT,chemtype='smiles',OutpreparedMOLfile='outPreparedLigand.pdbqt',geneID,geneToPDBmulti=F,LocalPDBconfigList=LocalPDBList,LocalPDBconfigDIR='dock_data/box_config',python2PTAH='c:/Progra~2/MGLTools-1.5.6rc3/python2',prepare_ligand4PATH='E:/ADtools',prepare_receptor4PATH='E:/ADtools',PSOvina=F){ 
  ###chemtype:smiles;inchikey;mol2(fileName) 
  ###pdbtype:id(fetch from web) or not id(PdbFileName) 
  ####RefLigandFile(*.pdb):NULL--No RMSDcal; 
  if (file.exists('DockedPDBQT.pdbqt'))file.remove('DockedPDBQT.pdbqt') 
  if (file.exists(OutpreparedMOLfile))file.remove(OutpreparedMOLfile) 
  if (file.exists('TempOutProteinPrepared.pdbqt'))file.remove('TempOutProteinPrepared.pdbqt') 
  if (file.exists('TempRefOutLigandPrepared.pdbqt'))file.remove('TempRefOutLigandPrepared.pdbqt') 
  #if (file.exists(outPDBQT))file.remove(outPDBQT) 
  prepareLigand(molfileTXT=chemTXT,moltype=chemtype,preparedMOLfile=OutpreparedMOLfile,python2PTAH=python2PTAH,prepare_ligand4PATH=prepare_ligand4PATH) 
  pdbID=geneIDtoPDB(geneID = geneID,multi=geneToPDBmulti) 
  if (class(pdbID)=='data.frame'){ 
    pdbID[,'Bestaffinity']=numeric() 
    for (i in 1:nrow(pdbID)){ 
      TempPdb=tolower(pdbID$PDB[i]) 
      TempConfigfile=paste0(LocalPDBconfigDIR,'/',TempPdb,'_BOXconfig.txt') 
      TempBox=getBoxFromPdb(pdbFile=TempPdb,IFprepare=T,python2PTAH=python2PTAH,prepare_ligand4PATH=prepare_ligand4PATH) 
      if (TempPdb%in%LocalPDBconfigList){ 
        TempDockResult=vinaDock3(preparedLigandPDBQT=OutpreparedMOLfile,pdb=TempPdb,pdbtype='id',configfileTYPE='file',configfile=TempConfigfile,python2PTAH=python2PTAH,prepare_ligand4PATH=prepare_ligand4PATH,prepare_receptor4PATH=prepare_receptor4PATH,getLigand=F,PSOvina=PSOvina) 
        if (class(try({Tempaffinity=TempDockResult$affinity},silent = T))!='try-error'){ 
          TempBest=min(Tempaffinity) 
        }else{ 
         TempBest=NA 
        } 
      }else if (sum(!is.na(TempBox))>0){ 
        TempDockResult=vinaDock3(preparedLigandPDBQT=OutpreparedMOLfile,pdb=TempPdb,pdbtype='id',configfileTYPE='getFromPDB',python2PTAH=python2PTAH,prepare_ligand4PATH=prepare_ligand4PATH,prepare_receptor4PATH=prepare_receptor4PATH,getLigand=F,PSOvina=PSOvina) 
        if (class(try({Tempaffinity=TempDockResult$affinity},silent = T))!='try-error'){ 
          TempBest=min(Tempaffinity) 
        }else{ 
          TempBest=NA 
        } 
      }else{ 
        TempBest=NA 
      } 
      pdbID$Bestaffinity[i]=TempBest 
    } 
  } 
  result=pdbID 
  return(result) 
} 
 
