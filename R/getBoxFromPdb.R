getBoxFromPdb = function(pdbFile,outLigand='outLigand.pdb',outVINAconfigFile='config.txt',IFprepare=T,python2PTAH='c:/Progra~2/MGLTools-1.5.6rc3/python2',prepare_ligand4PATH='E:/ADtools',offvalue=5){ 
  ###pdbFile:pdb or pdb id 
  ###IFprepare: if output parpared ligand(*.pdb) 
  ###if no ligand then output NA 
  require(bio3d) 
  aa=getLigandFromPdb(pdbFile=pdbFile,outLigand=outLigand,IFprepare=IFprepare,python2PTAH=python2PTAH,prepare_ligand4PATH=prepare_ligand4PATH) 
  if (!is.na(aa)){ 
    configF=character() 
    TempPdb <- bio3d::read.pdb(aa) 
    xyz=rep(c('x','y','z'),length(TempPdb$xyz)/3) 
    xmean=mean(TempPdb$xyz[xyz=='x']) 
    ymean=mean(TempPdb$xyz[xyz=='y']) 
    zmean=mean(TempPdb$xyz[xyz=='z']) 
    #xmax=max(abs(TempPdb$xyz[xyz=='x'])) 
    #ymax=max(abs(TempPdb$xyz[xyz=='y'])) 
    #zmax=max(abs(TempPdb$xyz[xyz=='z'])) 
    xmax=max(TempPdb$xyz[xyz=='x'])-min(TempPdb$xyz[xyz=='x'])+offvalue 
    ymax=max(TempPdb$xyz[xyz=='y'])-min(TempPdb$xyz[xyz=='y'])+offvalue 
    zmax=max(TempPdb$xyz[xyz=='z'])-min(TempPdb$xyz[xyz=='z'])+offvalue 
    configF[1]=paste('center_x =',xmean) 
    configF[2]=paste('center_y =',ymean) 
    configF[3]=paste('center_z =',zmean) 
    configF[4]=paste('size_x =',xmax) 
    configF[5]=paste('size_y =',ymax) 
    configF[6]=paste('size_z =',zmax) 
    writeLines(configF,outVINAconfigFile) 
    return(c(xmean=xmean,ymean=ymean,zmean=zmean,xmax=xmax,ymax=ymax,zmax=zmax)) 
  }else{ 
    return(NA) 
  } 
} 
 
