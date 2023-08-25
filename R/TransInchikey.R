TransInchikey = function(molecularFile,molecularType='SMI'){ 
  #molecularType='SMI':#molecularFile:SMI.csv-->colnames(molecularFile)=c('id','smiles') 
  #molecularType='SDF':#molecularFile:sdf-->datablocktag(sdf,'id') 
  #molecularType='SMItxt':SMI character text 
  require(data.table) 
  require(plyr) 
  require(stringr) 
  require(clipr) 
  require(igraph) 
  require(readr) 
  require(ggplot2) 
  require(ChemmineR) 
  if (molecularType=='SMI'){ 
    SMI=fread(molecularFile) 
    SMI2=SMI[smiles!='',] 
    smiles=SMI2$smiles 
    names(smiles)=SMI2$id 
    SMISdfSets=ChemmineR::smiles2sdf(smiles) 
    SMISdfSets=SMISdfSets[validSDF(SMISdfSets)] 
    #sum(validSDF(SMISdfSets))==length(SMISdfSets) 
    datablock(SMISdfSets)<-data.frame(id=SMISdfSets@ID) 
    CID=datablocktag(SMISdfSets,'id') 
    Nmol=length(CID) 
    inchikey=character() 
    length(inchikey)=Nmol 
    names(inchikey)=CID 
    for (i in 1:Nmol){ 
      TempSdf=write.SDF(SMISdfSets[i],'TempSdf.sdf',cid=T) 
      TempName=CID[i] 
      #Outputfilename=paste(TempName,'mol2',sep='.') 
      Dos=paste('obabel -isdf TempSdf.sdf -oinchikey','-Oinchikey.txt',sep=' ') 
      shell(Dos,intern=T) 
      if (file.exists('inchikey.txt')){ 
        key=read_lines('inchikey.txt') 
        inchikey[TempName]=key 
      }else{ 
        inchikey[TempName]=NA 
      } 
    } 
  }else if (molecularType=='SDF'){ 
    SMISdfSets=read.SDFset(molecularFile) 
    Nmol=length(SMISdfSets) 
    CID=datablocktag(SMISdfSets,'id') 
    inchikey=character() 
    length(inchikey)=Nmol 
    names(inchikey)=CID 
    for (i in 1:Nmol){ 
      TempSdf=write.SDF(SMISdfSets[i],'TempSdf.sdf',cid=T) 
      TempName=CID[i] 
      #Outputfilename=paste(TempName,'mol2',sep='.') 
      Dos=paste('obabel -isdf TempSdf.sdf -oinchikey','-Oinchikey.txt',sep=' ') 
      shell(Dos,intern=T) 
      if (file.exists('inchikey.txt')){ 
        key=read_lines('inchikey.txt') 
        inchikey[TempName]=key 
      }else{ 
        inchikey[TempName]=NA 
      } 
    } 
  }else if (molecularType=='SMItxt'){ 
    SMISdfSets=ChemmineR::smiles2sdf(molecularFile) 
    #SMISdfSets=SMISdfSets[validSDF(SMISdfSets)] 
    #sum(validSDF(SMISdfSets))==length(SMISdfSets) 
    datablock(SMISdfSets)<-data.frame(id=paste('temp',1:length(SMISdfSets),sep='')) 
    CID=datablocktag(SMISdfSets,'id') 
    Nmol=length(CID) 
    inchikey=character() 
    length(inchikey)=Nmol 
    names(inchikey)=CID 
    for (i in 1:Nmol){ 
      TempSdf=write.SDF(SMISdfSets[i],'TempSdf.sdf',cid=T) 
      TempName=CID[i] 
      #Outputfilename=paste(TempName,'mol2',sep='.') 
      Dos=paste('obabel -isdf TempSdf.sdf -oinchikey','-Oinchikey.txt',sep=' ') 
      shell(Dos,intern=T) 
      if (file.exists('inchikey.txt')){ 
        key=read_lines('inchikey.txt') 
        inchikey[TempName]=key 
      }else{ 
        inchikey[TempName]=NA 
      } 
    } 
  }else{ 
    print('ERROR:molecularType not found!') 
    inchikey=character() 
  } 
return(as.data.frame(inchikey,stringsAsFactors = FALSE)) 
} 
 
