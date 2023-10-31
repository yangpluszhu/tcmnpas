findStitchTargetAll = function(chem_inchikey,scoreSet=0,STITCH_CTdb='db/STITCH5_CTdb.db'){ 
#findStitchTargetAll=function(chem_inchikey,scoreSet=0,STITCH_CTdb='db/STITCH5_CTdb.db',STITCH_CKEYdb='db/STITCH5_CKEYdb.db') 
#chem_inchikey:inchikey 
  #STITCH_CTdb:'STITCH5_CTdb.db',colnames=c('cid','ENSP','score','geneID','inchikey') 
  #STITCH_CKEYdb:'STITCH_CKEYdb.db',colnames=c('inchikey','source_cid','stereo_chemical_id','flat_chemical_id') 
  #out：data.frame(cid,ENSP,geneID,score,inchikey) 
  require(RSQLite) 
  require(stringr) 
  require(data.table) 
  require(foreign) 
  result=data.table() 
  for (i in 1:length(chem_inchikey)){ 
  tempKey=chem_inchikey[i] 
  #tempKey_target=findStitchTarget(chem_inchikey=tempKey,scoreSet=scoreSet,STITCH_CTdb=STITCH_CTdb,STITCH_CKEYdb=STITCH_CKEYdb) 
  tempKey_target=findStitchTarget2(chem_inchikey=tempKey,scoreSet=scoreSet,STITCH_CTdb=STITCH_CTdb) 
  result=rbind(result,tempKey_target) 
  }   
  return(result) 
} 
 
