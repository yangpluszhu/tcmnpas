##' PrepareLigand for Molecular Docking. 
##' Note: vina, PSOvina,and mgltools must be installed and configured correctly. 
##' 
##' @title prepareLigand 
##' @param molfileTXT Molecular txt character
##' @param moltype Molecular type: smiles inchikey mol2 or pdb 
##' @param preparedMOLfile prepared file name for output: *.pdbqt
##' @param python2PTAH python2.exe PTAH 
##' @param prepare_ligand4PATH "prepare_ligand4.py" PTAH.The default Vina compound ligand preparation program (nphs–>: true to merge charges and remove non-polar hydrogens; lps–> true to merge charges and remove lone pairs of atoms) is used.
##' @return a prepared molecular file
##' @export 
##' @author Yang Ming 
prepareLigand = function(molfileTXT,moltype='smiles',preparedMOLfile='outLigand.pdbqt',python2PTAH='c:/Progra~2/MGLTools-1.5.6rc3/python2',prepare_ligand4PATH=system.file('tools',package='tcmnpas')){ 
  ##moltype:smiles;inchikey;mol2;pdb 
  if (file.exists('Temp.mol2'))file.remove('Temp.mol2') 
  if (moltype=='smiles'){ 
    SMItoMOL2(moltxt=molfileTXT,molName='Temp.mol2') 
    o0=paste0(prepare_ligand4PATH,'/prepare_ligand4.py') 
    o1='-l Temp.mol2' 
    o2='-o' 
    o3='-F do' 
    Dos=paste(python2PTAH,o0,o1,o2,preparedMOLfile,o3) 
    #Dos=paste(python2PTAH,o0,o1,o2,preparedMOLfile) 
    system(command = Dos,intern=T) 
  }else if (moltype=='inchikey'){ 
    smitxt=InchikeytoMOL2(Inchikey=molfileTXT) 
    if (!is.na(smitxt)){ 
      SMItoMOL2(moltxt=smitxt,molName='Temp.mol2') 
      o0=paste0(prepare_ligand4PATH,'/prepare_ligand4.py') 
      o1='-l Temp.mol2' 
      o2='-o' 
      o3='-F do' 
      Dos=paste(python2PTAH,o0,o1,o2,preparedMOLfile,o3) 
      #Dos=paste(python2PTAH,o0,o1,o2,preparedMOLfile) 
      system(command = Dos,intern=T) 
    }else{ 
      print('Error:no smiles found!Please retrive it from web!') 
      stop() 
    } 
  }else if (moltype%in%c('mol2','pdb')){ 
    o0=paste0(prepare_ligand4PATH,'/prepare_ligand4.py') 
    o1='-l' 
    o2=molfileTXT 
    o3='-o' 
    o4='-F do' 
    Dos=paste(python2PTAH,o0,o1,o2,o3,preparedMOLfile,o4) 
    #Dos=paste(python2PTAH,o0,o1,o2,o3,preparedMOLfile) 
    system(command = Dos,intern=T) 
  }else if (moltype=='sdf'){ 
    SDFtoMOL2(molFile=molfileTXT,molName='Temp.mol2') 
    o0=paste0(prepare_ligand4PATH,'/prepare_ligand4.py') 
    o1='-l Temp.mol2' 
    o2='-o' 
    o3='-F do' 
    Dos=paste(python2PTAH,o0,o1,o2,preparedMOLfile,o3) 
    #Dos=paste(python2PTAH,o0,o1,o2,preparedMOLfile) 
    system(command = Dos,intern=T) 
  } 
} 
 
