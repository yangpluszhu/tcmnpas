selectLargeLigand = function(pdbFile,type='file'){ 
  require(bio3d) 
  options(warn =-1) 
  if (type=='file'){ 
    pdb <- bio3d::read.pdb(pdbFile) 
  }else{ 
    pdb=pdbFile 
  } 
  #lid=table(pdb$atom$resid) 
  indsMetal=atom.select(pdb, c("metal")) 
  lid=table(pdb$atom$chain) 
  maxChain=names(lid)[lid==max(lid)][1] 
  ATOMinds=which(pdb$atom$chain==maxChain) 
  ATOMinds=setdiff(ATOMinds,indsMetal$atom) 
  #maxChain=table(pdb$atom$chain[ATOMinds]) 
  #maxChainName=names(maxChain)[maxChain==max(maxChain)[1]] 
  #selectedATOMinds=which((pdb$atom$resno==maxID)&(pdb$atom$chain==maxChainName)) 
  #selectedXYZinds=atom2xyz(which((pdb$atom$resno==maxID)&(pdb$atom$chain==maxChainName))) 
  selectedATOMinds=ATOMinds 
  selectedXYZinds=atom2xyz(ATOMinds) 
  selectedinds=list(atom=selectedATOMinds,xyz=selectedXYZinds,call=pdb$call) 
  subPDB=trim.pdb(pdb, inds=selectedinds) 
  #selectedATOM=pdb$atom[(pdb$atom$resno==maxID)&(pdb$atom$chain==names(table(pdb$atom$chain))[1]),] 
  #selectedXYZ=pdb$xyz[atom2xyz(which((pdb$atom$resno==maxID)&(pdb$atom$chain==names(table(pdb$atom$chain))[1])))] 
  #subPDB<-class(pdb) 
  #subPDB$atom<-selectedATOM 
  #subPDB$xyz<-selectedXYZ 
  #subPDB$call<-pdb$call 
  if (type=='file'){ 
    bio3d::write.pdb(subPDB,pdbFile) 
    os1='obabel -ipdb' 
    os2=pdbFile 
    os3='-opdb -O' 
    os4=pdbFile 
    Dos=paste(os1,os2,os3,os4) 
    system(command = Dos,intern=T) 
  }else{ 
    bio3d::write.pdb(subPDB,'TempNewLargeLigand.pdb') 
    os1='obabel -ipdb' 
    os2='TempNewLargeLigand.pdb' 
    os3='-opdb -O' 
    os4='TempNewLargeLigand.pdb' 
    Dos=paste(os1,os2,os3,os4) 
    system(command = Dos,intern=T) 
  } 
} 
 
