TransMol = function(sdfFile){ 
  require(data.table) 
  require(plyr) 
  require(stringr) 
  require(igraph) 
  require(ggplot2) 
  require(ChemmineR) 
  require(PaDEL) 
  SDF=read.SDFset(sdfFile) 
  Nmol=length(SDF) 
  CID=datablocktag(SDF,'id') 
  for (i in 1:Nmol){ 
    TempSdf=write.SDF(SDF[i],'TempSdf.sdf',cid=T) 
    TempName=CID[i] 
    Outputfilename=paste(TempName,'mol2',sep='.') 
    Dos=paste('obabel -isdf TempSdf.sdf -oMOL2','-O',Outputfilename,'--gen3D',sep=' ') 
    shell(Dos,intern=T) 
  } 
} 
 
