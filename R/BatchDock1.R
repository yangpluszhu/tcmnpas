BatchDock1 = function(chemList,chemListtype='smiles',PDBList,PDBListtype='file',python2PTAH='~/mgltools_x86_64Linux2_1.5.6/bin/pythonsh',prepare_ligand4PATH='~/ADtools',PSOvina=F){ 
  #####chemListtype:smiles-->chemList[named character] OR file-->chemList:*.csv[colname:id,smiles] 
  #####PDBListtype:file-->PDBList:*.csv[colname:PDB,x,y,z,sizex,sizey,sizez] OR dataFrame-->PDBList[colname:PDB,x,y,z,sizex,sizey,sizez] 
  require(rio) 
  if (chemListtype=='smiles'){ 
    chemID=names(chemList) 
    chemSmiles=chemList 
  }else{ 
    chemFile=import(chemList) 
    chemID=chemFile$id 
    chemSmiles=chemFile$smiles 
  } 
  ############ 
  if (PDBListtype=='file'){ 
    PDBListdata=import(PDBList) 
  }else{ 
    PDBListdata=PDBList 
  } 
  ############# 
  NumChem=length(chemID) 
  NumPDB=length(PDBListdata$PDB) 
  DockingResult=data.frame() 
  outDockFilePath=character() 
  if (NumChem*NumPDB<51){ 
  for (i in 1:NumChem){ 
    tempChem=chemSmiles[i] 
    tempChemID=chemID[i] 
    for (k in 1:NumPDB){ 
      tempPDB=PDBListdata$PDB[k] 
      Tcenter_x=PDBListdata$x[k] 
      Tcenter_y=PDBListdata$y[k] 
      Tcenter_z=PDBListdata$z[k] 
      Tsize_x=PDBListdata$sizex[k] 
      Tsize_y=PDBListdata$sizey[k] 
      Tsize_z=PDBListdata$sizez[k] 
      tmpdir <- tempdir() 
      outDockFileName=paste0(tempChemID,'_',tempPDB,'_DockedPDBQT') 
      #outDockFilePath=paste0(tmpdir,'/',tempChemID,'_',tempPDB,'_DockedPDBQT') 
      if (class(try({tempResult<-shinyVina(chemTXT=tempChem,chemtype='smiles',pdb=tempPDB,pdbtype='id',getBoxtype='Input',center_x=Tcenter_x,center_y=Tcenter_y,center_z=Tcenter_z,size_x=Tsize_x,size_y=Tsize_y,size_z=Tsize_z,outDockedPDBQT=outDockFileName,python2PTAHs=python2PTAH,prepare_ligand4PATHs=prepare_ligand4PATH,prepare_receptor4PATHs=prepare_ligand4PATH,PDBboxDbs='PDBboxDB.db',exhaustivenessSet=8,num_modesSet=1,energy_rangeSet=1,cpuSet=1,PSOvina=PSOvina)},silent = T))!='try-error'){ 
        if (nrow(tempResult$DockResult)>0){ 
          tempOutdockFile=tempResult$outDockedPDBQT 
		  tempOutdockFilepath=tempResult$outDockedResultPDBQTpath 
          outDockFilePath=c(outDockFilePath,tempOutdockFilepath) 
          output=data.frame(Chem=tempChemID,PDB=tempPDB,Bestaffinity=tempResult$DockResult$affinity[1],DockedFile=outDockFileName) 
        }else{ 
          tempOutdockFile=NULL 
		  tempOutdockFilepath=NULL 
          outDockFilePath=c(outDockFilePath,tempOutdockFile) 
          output=data.frame(Chem=tempChemID,PDB=tempPDB,Bestaffinity=NA,DockedFile=NA) 
        } 
      }else{ 
        tempOutdockFile=NULL 
		tempOutdockFilepath=NULL 
        outDockFilePath=c(outDockFilePath,tempOutdockFile) 
        output=data.frame(Chem=tempChemID,PDB=tempPDB,Bestaffinity=NA,DockedFile=NA) 
      } 
      #tempResult=shinyVina(chemTXT=tempChem,chemtype='smiles',pdb=tempPDB,pdbtype='id',getBoxtype='Input',center_x=Tcenter_x,center_y=Tcenter_y,center_z=Tcenter_z,size_x=Tsize_x,size_y=Tsize_y,size_z=Tsize_z,outDockedPDBQT=outDockFileName,python2PTAHs=python2PTAH,prepare_ligand4PATHs=prepare_ligand4PATH,prepare_receptor4PATHs=prepare_ligand4PATH,PDBboxDbs='PDBboxDB.db',exhaustivenessSet=8,num_modesSet=1,energy_rangeSet=1,cpuSet=1,PSOvina=PSOvina) 
      #tempOutdockFile=tempResult$outDockedPDBQT 
      #outDockFilePath=c(outDockFilePath,tempOutdockFile) 
      #output=data.frame(Chem=tempChemID,PDB=tempPDB,Bestaffinity=tempResult$DockResult$affinity[1],DockedFile=outDockFileName) 
      DockingResult=rbind(DockingResult,output) 
    } 
  } 
  }else{ 
    DockingResult=data.frame() 
    outDockFilePath=character() 
  } 
  return(list(DockingResult=DockingResult,outDockFilePath=outDockFilePath)) 
} 
 
