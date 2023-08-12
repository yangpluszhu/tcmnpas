#Sys.setlocale("LC_ALL", "us")
geneDockNoRefligand=function(chemTXT,chemtype='smiles',OutpreparedMOLfile='outPreparedLigand.pdbqt',geneID,geneToPDBmulti=F,LocalPDBconfigList=LocalPDBList,LocalPDBconfigDIR='dock_data/box_config',python2PTAH='c:/Progra~2/MGLTools-1.5.6rc3/python2',prepare_ligand4PATH='E:/ADtools',prepare_receptor4PATH='E:/ADtools',PSOvina=F){
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
#####################################################################################################
##########################################################################################################
vinaDock2=function(chemTXT,chemtype='smiles',pdb,pdbtype='id',configfileTYPE='file',configfile='966c_BOXconfig.txt',center_x=NULL,center_y=NULL,center_z=NULL,size_x=NULL,size_y=NULL,size_z=NULL,outPDBQT='DockedPDBQT.pdbqt',python2PTAH='c:/Progra~2/MGLTools-1.5.6rc3/python2',prepare_ligand4PATH='E:/ADtools',prepare_receptor4PATH='E:/ADtools',getLigand=T,RefLigandFile='Ligand.pdb',IFProLigandPrepared=T,PSOvina=F){
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
################################################################################################################
vinaDock3=function(preparedLigandPDBQT='TempOutLigandPrepared.pdbqt',pdb,pdbtype='id',configfileTYPE='file',configfile='966c_BOXconfig.txt',center_x=NULL,center_y=NULL,center_z=NULL,size_x=NULL,size_y=NULL,size_z=NULL,outPDBQT='DockedPDBQT.pdbqt',python2PTAH='c:/Progra~2/MGLTools-1.5.6rc3/python2',prepare_ligand4PATH='E:/ADtools',prepare_receptor4PATH='E:/ADtools',getLigand=T,RefLigandFile='Ligand.pdb',IFProLigandPrepared=T,PSOvina=F){
  ###chemtype:smiles;inchikey;mol2(fileName)
  ###pdbtype:id(fetch from web) or not id(PdbFileName)
  ###configfileTYPE:file-->Filetxt;input-->center_x,center_y,center_z,size_x,size_y,size_z;getFromPDB-->getboxFromPDB
  ####output:result--dataFrame[mode affinity dist1 dist2];outPDBQTFile[Docked_pdbqt];getLigand=T-->RefLigandFile['Ligand.pdb']
  if (file.exists('DockedPDBQT.pdbqt'))file.remove('DockedPDBQT.pdbqt')
  if (file.exists('TempOutProteinPrepared.pdbqt'))file.remove('TempOutProteinPrepared.pdbqt')
  prepareProtein(pdbFile=pdb,pdbtype=pdbtype,preparedProfile='TempOutProteinPrepared.pdbqt',python2PTAH=python2PTAH,prepare_receptor4PATH=prepare_receptor4PATH,getLigand=getLigand,RefLigandFile=RefLigandFile,IFProLigandPrepared=IFProLigandPrepared)
  if (configfileTYPE=='getFromPDB'){
    aaa=getBoxFromPdb(pdbFile=pdb,outLigand=RefLigandFile,outVINAconfigFile='config.txt',IFprepare=T,python2PTAH=python2PTAH,prepare_ligand4PATH=prepare_ligand4PATH)
    if (!is.na(aaa)){
      result=vinaDock(preparedLigand=preparedLigandPDBQT,preparedProtein='TempOutProteinPrepared.pdbqt',configfileTYPE='input',center_x=aaa['xmean'],center_y=aaa['ymean'],center_z=aaa['zmean'],size_x=aaa['xmax'],size_y=aaa['ymax'],size_z=aaa['zmax'],outPDBQT=outPDBQT,PSO=PSOvina)
    }else{
      result='ERROR:no ligand found in PDB!'
    }
  }else{
    result=vinaDock(preparedLigand=preparedLigandPDBQT,preparedProtein='TempOutProteinPrepared.pdbqt',configfile=configfile,configfileTYPE=configfileTYPE,center_x=center_x,center_y=center_y,center_z=center_z,size_x=size_x,size_y=size_y,size_z=size_z,outPDBQT=outPDBQT,PSO=PSOvina)
  }
  return(result)
}
#############################################################################################
RMSDcal=function(DockedPDBQTfile,LigandPdbFile,IFfit=F){
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
####################################################################
vinaDock=function(preparedLigand,preparedProtein,configfile='966c_BOXconfig.txt',configfileTYPE='file',center_x=NULL,center_y=NULL,center_z=NULL,size_x=NULL,size_y=NULL,size_z=NULL,exhaustiveness=8,num_modes=9,energy_range=3,cpu=2,outPDBQT='DockedPDBQT.pdbqt',PSO=FALSE){
  ####configfileTYPE:file-->Filetxt;input-->center_x,center_y,center_z,size_x,size_y,size_z
  ####output:dockResult(dataFrame)colname-->c('mode','affinity','dist1','dist2')
  ####outPDBQT:Docked File
  if (file.exists('tempDock.log'))file.remove('tempDock.log')
  if (PSO){
    o1='~/psovina2ls/build/linux/release/psovina2ls --ligand'
  }else{
    o1='vina --ligand'
  }
  o2=preparedLigand
  o3='--receptor'
  o4=preparedProtein
  o5='--config'
  o6=configfile
  o7='--out'
  o8=outPDBQT
  o9='--log tempDock.log'
  oset1='--center_x'
  oset2='--center_y'
  oset3='--center_z'
  oset4='--size_x'
  oset5='--size_y'
  oset6='--size_z'
  if (configfileTYPE=='file'){
    Dos=paste(o1,o2,o3,o4,o5,o6,o7,o8,o9,paste0('--exhaustiveness ',exhaustiveness),paste0('--num_modes ',num_modes),paste0('--energy_range ',energy_range),paste0('--cpu ',cpu),paste0('--seed 123456'))
  }else if (configfileTYPE=='input'){
    Dos=paste(o1,o2,o3,o4,oset1,center_x,oset2,center_y,oset3,center_z,oset4,size_x,oset5,size_y,oset6,size_z,o7,o8,o9,paste0('--exhaustiveness ',exhaustiveness),paste0('--num_modes ',num_modes),paste0('--energy_range ',energy_range),paste0('--cpu ',cpu),paste0('--seed 123456'))
  }
  system(command = Dos,intern=T)
  dockResult=ExtrVinaLogs(LogFile='tempDock.log')
  return(dockResult)
}
#######################################################################################################
#######################################################################################################
ExtrVinaLogs=function(LogFile='result.log'){
  require(stringr)
  options(stringsAsFactors = F)
  ExtractAffvalue=function(x){
    x=unlist(str_split(x,' '))
    x=x[x!='']
    x=as.numeric(x)
    return(x)
  }
  ####
  a=readLines(LogFile,skipNul=T)
  a=a[(which(str_detect(a,'mode'))+3):(length(a)-1)]
  Affvalue=numeric()
  for (i in 1:length(a)){
    temp=ExtractAffvalue(a[i])
    Affvalue=rbind(Affvalue,temp)
  }
  colnames(Affvalue)=c('mode','affinity','dist1','dist2')
  Affvalue=as.data.frame(Affvalue)
  rownames(Affvalue)<-NULL
  ####
  return(Affvalue)
}
######################################################################################
getBoxFromFpocket=function(pdbID,offset=5){
  require(bio3d)
  require(BiocGenerics)
  require(stringr)
  if (!str_detect(pdbID,'\\.')){
    a=bio3d::read.pdb(pdbID)
    pdbName=paste0(pdbID,'.pdb')
    if (file.exists(pdbName))file.remove(pdbName)
    bio3d::write.pdb(a,pdbName)
    pdbIDOK=pdbID
    os1='fpocket -f'
    os2=pdbName
    system(command = paste(os1,os2),intern=T)
    pqrFileLocat=paste0(pdbIDOK,'_out/',pdbIDOK,'_pockets.pqr')
    result=ExtrFpocketvalue(pqrFile=pqrFileLocat,offsetvalue=offset)
    file.remove(pdbName)
    unlink(paste0(pdbIDOK,'_out'), recursive = T)
  }else{
    pdbFIle=basename(pdbID)
    pdbIDOK=unlist(str_split(pdbFIle,'\\.'))[1]
    pdbFIlePATH=str_replace(pdbID,pdbFIle,'')
    os1='fpocket -f'
    os2=pdbID
    system(command = paste(os1,os2),intern=T)
    pqrFileLocat=paste0(pdbFIlePATH,pdbIDOK,'_out/',pdbIDOK,'_pockets.pqr')
    result=ExtrFpocketvalue(pqrFile=pqrFileLocat,offsetvalue=offset)
    #file.remove(pdbID)
    unlink(paste0(pdbFIlePATH,pdbIDOK,'_out'), recursive = T)
  }
  return(result)
}
##########################################
ExtrFpocketvalue=function(pqrFile='966c_out/966c_pockets.pqr',offsetvalue=5){
  require(stringr)
  require(data.table)
  options(stringsAsFactors = F)
  a=readLines(pqrFile)
  a=a[which(str_detect(a,'ATOM'))]
  ExtrFpocket=function(x){
    x=unlist(str_split(x,' '))
    x=x[x!='']
    return(x)
  }
  Fpocketvalue=character()
  for (i in 1:length(a)){
    temp=ExtrFpocket(a[i])
    Fpocketvalue=rbind(Fpocketvalue,temp)
  }
  colnames(Fpocketvalue)=c('atom','n1','n2','n3','pockID','x','y','z','d','r')
  Fpocketvalue=as.data.frame(Fpocketvalue)
  Fpocketvalue$x=as.numeric(Fpocketvalue$x)
  Fpocketvalue$y=as.numeric(Fpocketvalue$y)
  Fpocketvalue$z=as.numeric(Fpocketvalue$z)
  Fpocketvalue$pockID=as.numeric(Fpocketvalue$pockID)
  rownames(Fpocketvalue)<-NULL
  Fpocketvalue=as.data.table(Fpocketvalue)
  Fpocketvalue$pockID=as.character(Fpocketvalue$pockID)
  MaxPocketIDtable=table(Fpocketvalue$pockID)
  MaxPocketID=names(MaxPocketIDtable)[MaxPocketIDtable==max(MaxPocketIDtable)]
  FpocketvalueStat=Fpocketvalue[,.(x=mean(x),y=mean(y),z=mean(z),sizex=offsetvalue+max(x)-min(x),sizey=offsetvalue+max(y)-min(y),sizez=offsetvalue+max(z)-min(z)),by='pockID']
  FpocketvalueStat=as.data.frame(FpocketvalueStat)
  firstPid=FpocketvalueStat[FpocketvalueStat$pockID==MaxPocketID,]
  return(list(firstPid=firstPid,FpocketvalueStat=FpocketvalueStat))
}
###########################################################################
prepareLigand=function(molfileTXT,moltype='smiles',preparedMOLfile='outLigand.pdbqt',python2PTAH='c:/Progra~2/MGLTools-1.5.6rc3/python2',prepare_ligand4PATH='E:/ADtools'){
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
######################
prepareProtein=function(pdbFile,pdbtype='id',preparedProfile='outPro.pdbqt',python2PTAH='c:/Progra~2/MGLTools-1.5.6rc3/python2',prepare_receptor4PATH='E:/ADtools',getLigand=T,RefLigandFile='Ligand',IFProLigandPrepared=T,repair=T){
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
########################################################################
getBoxFromPdb=function(pdbFile,outLigand='outLigand.pdb',outVINAconfigFile='config.txt',IFprepare=T,python2PTAH='c:/Progra~2/MGLTools-1.5.6rc3/python2',prepare_ligand4PATH='E:/ADtools',offvalue=5){
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
####
getLigandFromPdb=function(pdbFile,outLigand='RefoutLigand',IFprepare=T,python2PTAH='c:/Progra~2/MGLTools-1.5.6rc3/python2',prepare_ligand4PATH='E:/ADtools'){
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
######################################################################################################################
selectLargeLigand=function(pdbFile,type='file'){
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
#######################################################
InchikeytoMOL2=function(Inchikey,MOLdataDB='../npa/db/inchikey_smi_db.Rdata'){
  library(RCurl)
  library(rvest)
  library(data.table)
  options(stringsAsFactors = F)
  #MOLdata=fread('E:/DATABASE/MYpro/BanXiaXieXin/result/ToolBox/db/2020-01-15_16_27_34_STITCH_HIT_TCMID_TCMSP_Tar.csv')
  #MOLdata=unique(MOLdata[,.(inchikey,smiles)])
  #save(lsit='MOLdata',file='inchikey_smi_db.Rdata')
  load(MOLdataDB)
  smi=MOLdata$smiles[chmatch(Inchikey,MOLdata$inchikey)]
  return(smi)
}
#######################################################
validatSMI=function(smiTXT){
  library(stringr)
  tempOutName='TEMPvalidateSMI.mol2'
  validat=logical()
  for (i in 1:length(smiTXT)){
    print(i)
    if (file.exists(tempOutName))file.remove(tempOutName)
    SMItoMOL2(moltxt=smiTXT[i],molName=tempOutName)
    if (file.exists(tempOutName)){
      if (file.info(tempOutName)['size']!=0){
        tempV=T
      }else{
        tempV=F
      }
    }else{
      tempV=F
    }
    validat=c(validat,tempV)
  }
  return(validat)
}
##########################################################
SMItoMOL2=function(moltxt,molName='mole.mol2'){
  library(stringr)
  library(rcdk)
  m <- parse.smiles(moltxt)[[1]]
  moltxt2=get.smiles(m)
  ##moltxt:smiles
  os1='obabel -:"'
  os2=moltxt2
  os3=paste0(os1,os2,'"')
  os4='-omol2 -O'
  os5=molName
  os6='--gen3d'
  Dos=paste(os3,os4,os5,os6)
  system(command = Dos,intern=T)
}
####################################################
SDFtoSMILES=function(SDFfile,outtype='txt',outName='smiles.csv'){
  require(ChemmineR)
  sdfdata=read.SDFset(SDFfile)
  valid <- validSDF(sdfdata)
  sdfdata <- sdfdata[valid]
  NumChem=length(sdfdata)
  if (class(try(sdfdata@ID,silent = T))=='try-error'){
    ChemID<-paste0('Chem',1:length(NumChem))
  }else{
    ChemID=sdfdata@ID
  }
    Tempsmiles <- as.character(sdf2smiles(sdfdata))
    if (outtype=='txt'){
      smiles=Tempsmiles
      names(smiles)=ChemID
    }else{
      smiles=data.frame(id=ChemID,smiles=Tempsmiles)
      write.csv(smiles,outName,row.names = F)
    }
  return(smiles)
}
##############################################
SDFtoMOL2=function(molFile,molName='mole.mol2'){
  library(stringr)
  ##molFile:SDFfileName
  os1='obabel -isdf'
  os2=molFile
  os3='-omol2'
  os4='-O'
  os5=molName
  os6='--gen3d'
  Dos=paste(os1,os2,os3,os4,os5,os6)
  system(command = Dos,intern=T)
}
################################################
MOL2toSMILES=function(molFile,molName='mole.smi'){
  library(stringr)
  ##molFile:SDFfileName
  os1='obabel -imol2'
  os2=molFile
  os3='-osmiles'
  os4='-O'
  os5=molName
  #os6='--gen3d'
  Dos=paste(os1,os2,os3,os4,os5)
  system(command = Dos,intern=T)
}
###############################################
PDBQTtoPDB=function(molFile,molName='mole.pdb'){
  library(stringr)
  ##molFile:PDBQTfileName
  os1='obabel -ipdbqt'
  os2=molFile
  os3='-opdb'
  os4='-O'
  os5=molName
  Dos=paste(os1,os2,os3,os4,os5)
  system(command = Dos,intern=T)
}
################################################
judgeChemType=function(ChemTXT){
  require(stringr)
  temp=ChemTXT
  if (nchar(temp)==27&nchar(unlist(str_split(temp,'-')))[1]==14){
    ChemType='inchikey'
  }else{
    ChemType='smiles'
  }
  return(ChemType)
}
##################################################
geneIDtoPDB=function(geneID,multi=F,geneTopdbDB='../npa/db/GeneIDtoPDB.db',up='../npa/db/UniProt.ws_up9606.RData'){
  require(RSQLite)
  require(UniProt.ws)
  require(data.table)
  options(stringsAsFactors = F)
  load(up,baseenv())
  #load(geneNApdbListdb)
  con <- dbConnect(SQLite(),dbname=geneTopdbDB)
  afcOS1='SELECT * FROM geneNApdbListD'
  eval(substitute({rsOS1<-dbSendQuery(con, afcOS1)},list(afcOS1=afcOS1)))
  d1afcOS1 <- dbFetch(rsOS1)
  geneNApdbListdb=as.character(d1afcOS1$value)
  ######
  afcOS2='SELECT * FROM geneHavingPDBListD'
  eval(substitute({rsOS2<-dbSendQuery(con, afcOS2)},list(afcOS2=afcOS2)))
  d1afcOS2 <- dbFetch(rsOS2)
  geneHavingPDBList=as.character(d1afcOS2$value)
  ######
  afcOS3='SELECT * FROM genenoHavingPDBListD'
  eval(substitute({rsOS3<-dbSendQuery(con, afcOS3)},list(afcOS3=afcOS3)))
  d1afcOS3<- dbFetch(rsOS3)
  genenoHavingPDBList=as.character(d1afcOS3$value)
  ######
  afcOS4='SELECT * FROM geneInmypdbListD'
  eval(substitute({rsOS4<-dbSendQuery(con, afcOS4)},list(afcOS4=afcOS4)))
  d1afcOS4<- dbFetch(rsOS4)
  geneInmypdbList=as.character(d1afcOS4$value)
  ######
  dbDisconnect(con)
  geneID=as.character(geneID)
  geneID=geneID[!geneID%in%geneNApdbListdb]
  if (length(geneID)>0){
    if (multi==F){
    geneIDIN=geneID[geneID%in%geneHavingPDBList]
    geneIDnoIN=geneID[geneID%in%genenoHavingPDBList]
    geneIDWeb=setdiff(geneID,c(geneIDIN,geneIDnoIN))
    con <- dbConnect(SQLite(),dbname=geneTopdbDB)
    if (length(geneIDIN)>0){
      geneID2=paste("'",geneIDIN,"'",sep='',collapse = ',')
      geneID3=paste("(",geneID2,")",sep='')
      os1=paste0('SELECT * FROM havingPDB where geneID in',geneID3)
      eval(substitute({rs1<-dbSendQuery(con, os1)},list(os1=os1)))
      d1 <- dbFetch(rs1)
      d1=d1[,c('geneID','PDB')]
    }else{
      d1=data.frame()
    }
    ##
    if (length(geneIDnoIN)){
      geneID4=paste("'",geneIDnoIN,"'",sep='',collapse = ',')
      geneID5=paste("(",geneID4,")",sep='')
      os2=paste0('SELECT * FROM nohavingPDB where geneID in',geneID5)
      eval(substitute({rs2<-dbSendQuery(con, os2)},list(os2=os2)))
      d2 <- dbFetch(rs2)
      d2=d2[,c('geneID','PDB')]
    }else{
      d2=data.frame()
    }
    ##
    if (length(geneIDWeb)>0){
      if (class(try({res3 <- UniProt.ws::select(up, keys = geneIDWeb,columns = c("PDB"),keytype = "ENTREZ_GENE")},silent = T))!="try-error"){
        colnames(res3)=c('geneID','PDB')
        res3=as.data.table(res3)
        res3=unique(res3,by='geneID')
        d3=as.data.frame(res3)
      }else{
        d3=data.frame()
      }
    }else{
      d3=data.frame()
    }
    d_total=rbind(d1,d2,d3)
    if (nrow(d_total)>0){
      result=unique(d_total)
    }else{
      result='ERROR:geneID not found in PDB database!'
    }
    }else if (multi){
      #######################multi#####################
      geneIDINList=geneID[geneID%in%geneInmypdbList]
      geneIDnoINList=geneID[!geneID%in%geneInmypdbList]
      con <- dbConnect(SQLite(),dbname=geneTopdbDB)
      if (length(geneIDINList)>0){
        geneID22=paste("'",geneIDINList,"'",sep='',collapse = ',')
        geneID23=paste("(",geneID22,")",sep='')
        os21=paste0('SELECT * FROM geneTopdb where geneID in',geneID23)
        eval(substitute({rs1<-dbSendQuery(con, os21)},list(os21=os21)))
        d21 <- dbFetch(rs1)
        d21=d21[,c('geneID','PDB')]
      }else{
        d21=data.frame()
      }
      ###############
        if (length(geneIDnoINList)>0){
        if (class(try({res23 <- UniProt.ws::select(up, keys = geneIDnoINList,columns = c("PDB"),keytype = "ENTREZ_GENE")},silent = T))!="try-error"){
          colnames(res23)=c('geneID','PDB')
          res23=as.data.table(res23)
          res23=unique(res23,by='geneID')
          d23=as.data.frame(res23)
        }else{
          d23=data.frame()
        }
        }else{
          d23=data.frame()
      }
      ################################################
      d_total=rbind(d21,d23)
      if (nrow(d_total)>0){
        result=unique(d_total)
      }else{
        result='ERROR:geneID not found in PDB database!'
      }
    }
    dbDisconnect(con)
  }else{
    result='ERROR:geneID not found in PDB database!'
  }
  return(result)
}
##############################################################################################################################################################################################################################
getPDBbox=function(PDBid,type='FromLigand',center_x=NULL,center_y=NULL,center_z=NULL,size_x=NULL,size_y=NULL,size_z=NULL,python2PTAH='~/mgltools_x86_64Linux2_1.5.6/bin/pythonsh',prepare_ligand4PATH='~/ADtools',PDBboxDb='PDBboxDB.db'){
  ##ubuntu
  ##type:FromLigand;FromFpocket;Input
  #PDBboxDb first!
  require(bio3d)
  require(stringr)
  tempa=PDBid
  if (str_detect(tempa,'\\.')){
    #tempa=unlist(str_split(tempa,'\\.'))
    tempa2=basename(normalizePath(tempa))
    tempa2=str_replace_all(tempa2,'.pdb','')
    #tempa2=tempa2[length(tempa2)-1]
    #tempb=str_sub(tempa,start = nchar(tempa)-3,end = nchar(tempa))
    tempPDBid=tolower(tempa2)
  }else{
    tempPDBid=tolower(tempa)
  }
  if (type=='FromLigand'){
    if (file.exists(PDBboxDb)){
      os=paste0('SELECT * FROM PDBboxDB where PDB=','"',tempPDBid,'"')
      con <- dbConnect(SQLite(),dbname=PDBboxDb)
      rs <- dbSendQuery(con, os)
      d1 <- dbFetch(rs)
      if (nrow(d1)>0){
        PDBbox=d1
      }else{
        if (class(try({d2=getBoxFromPdb(pdbFile=tempa,IFprepare=T,python2PTAH=python2PTAH,prepare_ligand4PATH=prepare_ligand4PATH)},silent =T))!='try-error'){#c(xmean=xmean,ymean=ymean,zmean=zmean,xmax=xmax,ymax=ymax,zmax=zmax)
          if (sum(!is.na(d2))>0){
            d2=as.data.frame(t(d2))
            colnames(d2)=c('x','y','z','sizex','sizey','sizez')
            PDBbox=cbind(PDB=tempPDBid,d2,method='FromLigand')
          }else{
            PDBbox=data.frame()
          }
        }else{
          PDBbox=data.frame()
        }
      }
    }else{
      if (class(try({d2=getBoxFromPdb(pdbFile=tempa,IFprepare=T,python2PTAH=python2PTAH,prepare_ligand4PATH=prepare_ligand4PATH)},silent =T))!='try-error'){#c(xmean=xmean,ymean=ymean,zmean=zmean,xmax=xmax,ymax=ymax,zmax=zmax)
        if (sum(!is.na(d2))>0){
          d2=as.data.frame(t(d2))
          colnames(d2)=c('x','y','z','sizex','sizey','sizez')
          PDBbox=cbind(PDB=tempPDBid,d2,method='FromLigand')
        }else{
          PDBbox=data.frame()
        }
      }else{
        PDBbox=data.frame()
      }
    }
  }else if (type=='FromFpocket'){
    if (file.exists(PDBboxDb)){
      os=paste0('SELECT * FROM PDBboxDB where PDB=','"',tempPDBid,'"')
      con <- dbConnect(SQLite(),dbname=PDBboxDb)
      rs <- dbSendQuery(con, os)
      d1 <- dbFetch(rs)
      if (nrow(d1)>0){
        PDBbox=d1
      }else{
        if (class(try({d3=getBoxFromFpocket(pdbID=tempa)},silent = T))!='try-error'){
          PDBbox=cbind(PDB=tempPDBid,d3$firstPid[,c('x','y','z','sizex','sizey','sizez')],method='FromFpocket')
        }else{
          PDBbox=data.frame()
        }
      }
    }else{
      if (class(try({d3=getBoxFromFpocket(pdbID=tempa)},silent = T))!='try-error'){
        PDBbox=cbind(PDB=tempPDBid,d3$firstPid[,c('x','y','z','sizex','sizey','sizez')],method='FromFpocket')
      }else{
        PDBbox=data.frame()
      }
    }
  }else if (type=='Input'){
    PDBbox=data.frame(PDB=tempPDBid,x=center_x,y=center_y,z=center_z,sizex=size_x,sizey=size_y,sizez=size_z,method='input')
  }
  return(PDBbox)
}
################################################################################################################
shinyVina=function(chemTXT,chemtype='smiles',pdb,pdbtype='id',getBoxtype='FromLigand',center_x=NULL,center_y=NULL,center_z=NULL,size_x=NULL,size_y=NULL,size_z=NULL,outDockedPDBQT='DockedPDBQT',python2PTAHs='c:/Progra~2/MGLTools-1.5.6rc3/python2',prepare_ligand4PATHs='E:/ADtools',prepare_receptor4PATHs='E:/ADtools',PDBboxDbs='../npa/db/PDBboxDB.db',exhaustivenessSet=8,num_modesSet=9,energy_rangeSet=3,cpuSet=2,PSOvina=F){
  ###chemtype:smiles;inchikey;mol2(fileName)
  library(stringr)
  library(bio3d)
  library(RSQLite)
  tmpPDB=pdb
  temDIR0=tempdir()
  if (system('pwd')!=0){
    #pythonDir='c:/Progra~2/MGLTools-1.5.6rc3/python2'
    #ADTpath='E:/ADtools'
    temDIR=unlist(str_split(temDIR0,'\\\\'))
  }else{
    #pythonDir='~/mgltools_x86_64Linux2_1.5.6/bin/pythonsh'
    #ADTpath='~/ADtools'
    temDIR=unlist(str_split(temDIR0,'/'))
  }
  temDIR=temDIR[length(temDIR)]
  outDockedResultPDBQT=paste0(temDIR,outDockedPDBQT,'.pdbqt')
  outDockedResultPDBQTpath=paste0(temDIR0,'/',temDIR,outDockedPDBQT,'.pdbqt')
  temppreparedMOLfile=paste0(temDIR0,'/',temDIR,'_Ligand.pdbqt')
  temppreparedProfile=paste0(temDIR0,'/',temDIR,'_Protein.pdbqt')
  if (file.exists(outDockedResultPDBQT))file.remove(outDockedResultPDBQT)
  if (file.exists(temppreparedMOLfile))file.remove(temppreparedMOLfile)
  if (file.exists(temppreparedProfile))file.remove(temppreparedProfile)
  prepareLigand(molfileTXT=chemTXT,moltype=chemtype,preparedMOLfile=temppreparedMOLfile,python2PTAH=python2PTAHs,prepare_ligand4PATH=prepare_ligand4PATHs)
  BOXinfo=getPDBbox(PDBid=tmpPDB,type=getBoxtype,center_x=center_x,center_y=center_y,center_z=center_z,size_x=size_x,size_y=size_y,size_z=size_z,python2PTAH=python2PTAHs,prepare_ligand4PATH=prepare_ligand4PATHs,PDBboxDb=PDBboxDbs)
  prepareProtein(pdbFile=tmpPDB,pdbtype=pdbtype,preparedProfile=temppreparedProfile,python2PTAH=python2PTAHs,prepare_receptor4PATH=prepare_receptor4PATHs,getLigand=F)
  if (nrow(BOXinfo)>0&file.exists(temppreparedProfile)&file.exists(temppreparedMOLfile)){
    if (class(try({DockResult<-vinaDock(preparedLigand=temppreparedMOLfile,preparedProtein=temppreparedProfile,configfileTYPE='input',center_x=BOXinfo$x,center_y=BOXinfo$y,center_z=BOXinfo$z,size_x=BOXinfo$sizex,size_y=BOXinfo$sizey,size_z=BOXinfo$sizez,exhaustiveness=exhaustivenessSet,num_modes=num_modesSet,energy_range=energy_rangeSet,cpu=cpuSet,outPDBQT=outDockedResultPDBQTpath,PSO=PSOvina)},silent = T))=='try-error'){
      prepareProtein(pdbFile=tmpPDB,pdbtype=pdbtype,preparedProfile=temppreparedProfile,python2PTAH=python2PTAHs,prepare_receptor4PATH=prepare_receptor4PATHs,getLigand=F,repair = F)
      DockResult<-vinaDock(preparedLigand=temppreparedMOLfile,preparedProtein=temppreparedProfile,configfileTYPE='input',center_x=BOXinfo$x,center_y=BOXinfo$y,center_z=BOXinfo$z,size_x=BOXinfo$sizex,size_y=BOXinfo$sizey,size_z=BOXinfo$sizez,exhaustiveness=exhaustivenessSet,num_modes=num_modesSet,energy_range=energy_rangeSet,cpu=cpuSet,outPDBQT=outDockedResultPDBQTpath,PSO=PSOvina)
    }
    #DockResult=vinaDock(preparedLigand=temppreparedMOLfile,preparedProtein=temppreparedProfile,configfileTYPE='input',center_x=BOXinfo$x,center_y=BOXinfo$y,center_z=BOXinfo$z,size_x=BOXinfo$sizex,size_y=BOXinfo$sizey,size_z=BOXinfo$sizez,exhaustiveness=exhaustivenessSet,num_modes=num_modesSet,energy_range=energy_rangeSet,cpu=cpuSet,outPDBQT=outDockedResultPDBQT,PSO=PSOvina)
  }else{
    DockResult=data.frame()
  }
  return(list(DockResult=DockResult,BOXinfo=BOXinfo,outDockedPDBQT=outDockedResultPDBQT,preparedProtein=temppreparedProfile,preparedLigand=temppreparedMOLfile,outDockedResultPDBQTpath=outDockedResultPDBQTpath))
}
#####################################################################################################################################################################################################################################################################################HighLevel
#source('E:/ADtools/ADtoolBox.R')
#########################################################################################################################
BatchDock=function(chemList,chemListtype='smiles',PDBList,PDBListtype='file',dockType='cross',python2PTAH='~/mgltools_x86_64Linux2_1.5.6/bin/pythonsh',prepare_ligand4PATH='~/ADtools',PSOvina=F){
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
##############################################################################################
getPDBListbox=function(PDBList,PDBListtype='file',python2PTAH='~/mgltools_x86_64Linux2_1.5.6/bin/pythonsh',prepare_ligand4PATH='~/ADtools',PDBboxDb='PDBboxDB.db'){
  ###DB-->FromLigand-->Fpocket
  ####PDBList:character or *.txt
  ####python2PTAH='~/mgltools_x86_64Linux2_1.5.6/bin/pythonsh',prepare_ligand4PATH='~/ADtools',PDBboxDb='PDBboxDB.db'
  require(bio3d)
  require(stringr)
  if (PDBListtype=='file'){
    PDBquerylist=readLines(PDBList)
  }else{
    if (str_detect(PDBList,',|;|，|、')){
      PDBquerylist=unlist(str_split(PDBList,',|;|，|、'))
      PDBquerylist=str_trim(PDBquerylist,'both')
    }else{
      PDBquerylist=str_trim(PDBList,'both')
    }
  }
  PDBquerylist=unique(tolower(PDBquerylist))
  NumPdb=length(PDBquerylist)
  result=data.frame()
  for (i in 1:NumPdb){
    tempPDB=PDBquerylist[i]
    r1=getPDBbox(PDBid=tempPDB,type='FromLigand',python2PTAH=python2PTAH,prepare_ligand4PATH=prepare_ligand4PATH,PDBboxDb=PDBboxDb)
    if (nrow(r1)==0&system('pwd')==0){
      r1=getPDBbox(PDBid=tempPDB,type='FromFpocket',python2PTAH=python2PTAH,prepare_ligand4PATH=prepare_ligand4PATH,PDBboxDb=PDBboxDb)
    }
    if (nrow(r1)==0){
      r1=data.frame(PDB=tempPDB,x=NA,y=NA,z=NA,sizex=NA,sizey=NA,sizez=NA,method=NA)
    }
    result=rbind(result,r1)
  }
  return(result)
}
####################################################################################
BatchDock1=function(chemList,chemListtype='smiles',PDBList,PDBListtype='file',python2PTAH='~/mgltools_x86_64Linux2_1.5.6/bin/pythonsh',prepare_ligand4PATH='~/ADtools',PSOvina=F){
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
#############################################################################################################
BatchDockRef=function(chemList,chemListtype='smiles',PDBList,PDBListtype='file',python2PTAH='~/mgltools_x86_64Linux2_1.5.6/bin/pythonsh',prepare_ligand4PATH='~/ADtools',PSOvina=F){
  #####chemListtype:smiles-->chemList[named character] OR file-->chemList:*.csv[colname:id,smiles]
  #####PDBListtype:file-->PDBList:*.csv[colname:PDB,x,y,z,sizex,sizey,sizez] OR dataFrame-->PDBList[colname:PDB,x,y,z,sizex,sizey,sizez]
  require(rio)
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
  if (NumChem==NumPDB&NumChem<51){
    DockingResult=data.frame()
    outDockFilePath=character()
    for (i in 1:NumChem){
      tempChem=chemSmiles[i]
      tempChemID=chemID[i]
      tempPDB=PDBListdata$PDB[i]
      Tcenter_x=PDBListdata$x[i]
      Tcenter_y=PDBListdata$y[i]
      Tcenter_z=PDBListdata$z[i]
      Tsize_x=PDBListdata$sizex[i]
      Tsize_y=PDBListdata$sizey[i]
      Tsize_z=PDBListdata$sizez[i]
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
          outDockFilePath=c(outDockFilePath,tempOutdockFilepath)
          output=data.frame(Chem=tempChemID,PDB=tempPDB,Bestaffinity=NA,DockedFile=NA)
        }
      }else{
        tempOutdockFile=NULL
		tempOutdockFilepath=NULL
        outDockFilePath=c(outDockFilePath,tempOutdockFilepath)
        output=data.frame(Chem=tempChemID,PDB=tempPDB,Bestaffinity=NA,DockedFile=NA)
      }
      DockingResult=rbind(DockingResult,output)
    }
  }else{
    DockingResult=data.frame()
    outDockFilePath=character()
  }
  return(list(DockingResult=DockingResult,outDockFilePath=outDockFilePath))
}
###############################################################################################################
queryHerbChem=function(herbList=NULL,Qinchikey=NULL,database=c('HIT','TCMSP','TCMID'),HerbChemDB='../npa/db/herb_ChemDataAll.db'){
  require(stringr)
  require(RSQLite)
  if (!is.null(herbList)&&!is.na(herbList)){
    herbList2=paste("'",herbList,"'",sep='',collapse = ',')
    herbList3=paste("(",herbList2,")",sep='')
    herbOrder00=paste('herbID IN',herbList3,sep=' ')
  }else{
    herbOrder00=NULL
  }
  ###################
  if (!is.null(Qinchikey)&&!is.na(Qinchikey)){
    inchikey2=paste("'",Qinchikey,"'",sep='',collapse = ',')
    inchikey3=paste("(",inchikey2,")",sep='')
    QinchikeyOrder00=paste('inchikey IN',inchikey3,sep=' ')
  }else{
    QinchikeyOrder00=NULL
  }
  ###################
  databaseList2=paste("'",database,"'",sep='',collapse = ',')
  databaseList3=paste("(",databaseList2,")",sep='')
  databaseOrder00=paste('database IN',databaseList3,sep=' ')
  ###################
  if (!is.null(herbOrder00)|!is.null(QinchikeyOrder00)){
    OrderTemp=c(herbOrder00,QinchikeyOrder00,databaseOrder00)
    OrderTemp=OrderTemp[OrderTemp!='']
    if (length(OrderTemp)==1){
      Order=paste("SELECT * FROM herb_ChemDataAll where",OrderTemp,sep=' ')
    }else{
      OrderTemp2=paste(OrderTemp,collapse = ' AND ')
      Order=paste("SELECT * FROM herb_ChemDataAll where",OrderTemp2,sep=' ')
    }
  }else{
    Order=NULL
  }
  ###################
  if (!is.null(Order)){
    con <- dbConnect(SQLite(),dbname=HerbChemDB)
    res <- dbSendQuery(con, Order)
    result <- fetch(res, n =-1)
    dbDisconnect(con)
  }else{
    result=data.frame()
  }
  return(result)
}
##########################################################################################################################



