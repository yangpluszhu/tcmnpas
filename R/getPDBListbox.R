getPDBListbox = function(PDBList,PDBListtype='file',python2PTAH='~/mgltools_x86_64Linux2_1.5.6/bin/pythonsh',prepare_ligand4PATH='~/ADtools',PDBboxDb='PDBboxDB.db'){ 
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
 
