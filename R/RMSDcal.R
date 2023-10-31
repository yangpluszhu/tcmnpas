RMSDcal = function(DockedPDBQTfile,LigandPdbFile,IFfit=F){ 
  require(bio3d) 
  require(stringr) 
  ##LigandPdbFile:Refmolecular 
  ##DockedPDBQTfile:Vina output 
  #if (file.exists('TemppDockedPDBQT.pdb'))file.remove('TemppDockedPDBQT.pdb') 
  #PDBQTtoPDB(molFile=DockedPDBQTfile,molName='TemppDockedPDBQT.pdb') 
  if (str_detect(LigandPdbFile,'\\.cif')){ 
    Refbb=bio3d::read.cif(LigandPdbFile) 
  }else{ 
    Refbb=bio3d::read.pdb(LigandPdbFile) 
  } 
  aa=bio3d::read.pdb(DockedPDBQTfile,multi=T)       
  Mininds=min(dim(aa$xyz)[2],dim(Refbb$xyz)[2]) 
  inds <- bio3d::gap.inspect(Refbb$xyz) 
  RMSDresult=bio3d::rmsd(Refbb,aa$xyz,a.inds=inds$f.inds[1:Mininds], b.inds=inds$f.inds[1:Mininds],fit=IFfit) 
  return(data.frame(mode=1:dim(aa$xyz)[1],RMSD=RMSDresult)) 
} 
 
