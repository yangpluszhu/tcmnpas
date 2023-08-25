StitchChemTargetMake = function(StitchChemInchikey,Stitchprotein_chemical,up=up,STITCH_CTdb='STITCH_CTdb.db',STITCH_CKEYdb='STITCH_CKEYdb.db'){ 
#columns='ENTREZ_GENE' 
#kt='ENSEMBL_PROTEIN' 
#select(up,'ENSP00000257254',columns,kt) 
#StitchChemInchikey='E:\\DATABASE\\MYpro\\chemicals.inchikeys.v5.0.tsv' 
#Stitchprotein_chemical='E:\\DATABASE\\STITCH\\5.0\\9606.protein_chemical.links.v5.0.tsv' 
#output:STITCH_CTdb.db-->colnames(STITCHdb.db)=c('cid','ENSP','score','geneID','inchikey') 
  #STITCH_CKEYdb.db-->colnames=c('inchikey','source_cid','stereo_chemical_id','flat_chemical_id') 
  require(data.table) 
  require(stringr) 
  require(UniProt.ws) 
  require(RSQLite) 
  require(lubridate) 
  require(foreign) 
  if (!exists('up')){ 
    load('db/UniProtwsOB.RData') 
  } 
  ChemKey=fread(StitchChemInchikey,colClasses='character') 
  protein_chemical=fread(Stitchprotein_chemical) 
  protein_chemical$chemical=str_replace_all(protein_chemical$chemical,'[a-zA-Z]','') 
  protein_chemical$chemical=as.numeric(protein_chemical$chemical) 
  protein_chemical$protein=str_replace_all(protein_chemical$protein,'9606\\.','') 
  protein=unique(protein_chemical$protein) 
  tempD=data.table(id=as.character(1:length(protein)),protein=protein) 
  ENTREZ_GENE=tempD[,.(protein=protein,geneID=selectFun(protein)),by='id'] 
  ENTREZ_GENE=ENTREZ_GENE[geneID!='NA',] 
  proteinIn=unique(ENTREZ_GENE$protein) 
protein_chemicalOK=protein_chemical[protein%in%proteinIn,] 
setkey(ENTREZ_GENE,protein) 
protein_chemicalOK[,c('geneID'):=list(ENTREZ_GENE[.(protein_chemicalOK$protein),geneID])] 
protein_chemicalOK$chemical=as.character(protein_chemicalOK$chemical) 
setkey(ChemKey,source_cid) 
protein_chemicalOK[,c('inchikey'):=list(ChemKey[.(protein_chemicalOK$chemical),inchikey])] 
colnames(protein_chemicalOK)=c('cid','ENSP','score','geneID','inchikey') 
protein_chemicalOK=unique(protein_chemicalOK) 
# 
tmp <- dbConnect(SQLite(), STITCH_CTdb) 
dbWriteTable(tmp,'STITCHdb',protein_chemicalOK,append=F) 
dbSendQuery(tmp,'create index index_key on STITCHdb (inchikey)') 
dbSendQuery(tmp,'create index index_cid on STITCHdb (cid)') 
dbSendQuery(tmp,'create index index_ENSP on STITCHdb (ENSP)') 
dbSendQuery(tmp,'create index index_geneID on STITCHdb (geneID)') 
dbDisconnect(tmp) 
tmp2 <- dbConnect(SQLite(), STITCH_CKEYdb) 
dbWriteTable(tmp2,'STITCH_CKEYdb',ChemKey,append=F) 
dbSendQuery(tmp2,'create index index_key on STITCH_CKEYdb (inchikey)') 
dbSendQuery(tmp2,'create index index_cid on STITCH_CKEYdb (source_cid)') 
dbSendQuery(tmp2,'create index index_cid2 on STITCH_CKEYdb (stereo_chemical_id)') 
dbSendQuery(tmp2,'create index index_cid3 on STITCH_CKEYdb (flat_chemical_id)') 
dbDisconnect(tmp2) 
} 
 
