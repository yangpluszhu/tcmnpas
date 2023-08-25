##' Batch Molecular Docking. 
##' Note: mgltools must be installed and configured correctly. 
##' 
##' @title BatchDock 
##' @param chemList chemList[named character] OR file-->chemList:*.csv[colname:id,smiles] 
##' @param chemListtype "smiles" OR "file" 
##' @param PDBList PDBList:*.csv[colname:PDB,x,y,z,sizex,sizey,sizez] OR dataFrame-->PDBList[colname:PDB,x,y,z,sizex,sizey,sizez] 
##' @param PDBListtype "file" OR "dataFrame" 
##' @param dockType "cross" OR "parallel" 
##' @param python2PTAH pythonPath for mgltools 
##' @param prepare_ligand4PATH prepare_ligandPath for prepare functions 
##' @return a list Object 
##' @export 
##' @author Yang Ming 
BatchDock = function(chemList,chemListtype='smiles',PDBList,PDBListtype='file',dockType='cross',python2PTAH='~/mgltools_x86_64Linux2_1.5.6/bin/pythonsh',prepare_ligand4PATH='~/ADtools',PSOvina=F){ 
  #####chemListtype:smiles-->chemList[named character] OR file-->chemList:*.csv[colname:id,smiles] 
  #####PDBListtype:file-->PDBList:*.csv[colname:PDB,x,y,z,sizex,sizey,sizez] OR dataFrame-->PDBList[colname:PDB,x,y,z,sizex,sizey,sizez] 
  #####dockType:cross-->each chem*pdb OR parallel-->each row 
  #####output:list(DockingResult=DockingResult,outDockFilePath=outDockFilePath);DockingResult-->data.frame(Chem=tempChemID,PDB=tempPDB,Bestaffinity=tempResult$affinity[1],DockedFile=outDockFileName) 
  require(rio) 
  if (dockType=='cross'){ 
    BatchDockResult=BatchDock1(chemList=chemList,chemListtype=chemListtype,PDBList=PDBList,PDBListtype=PDBListtype,python2PTAH=python2PTAH,prepare_ligand4PATH=prepare_ligand4PATH,PSOvina=PSOvina) 
  }else if (dockType=='parallel'){ 
    BatchDockResult=BatchDockRef(chemList=chemList,chemListtype=chemListtype,PDBList=PDBList,PDBListtype=PDBListtype,python2PTAH=python2PTAH,prepare_ligand4PATH=prepare_ligand4PATH,PSOvina=PSOvina) 
  } 
  return(BatchDockResult) 
} 
 
