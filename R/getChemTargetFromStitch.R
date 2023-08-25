getChemTargetFromStitch = function(ChemData,ChemDataType='inchikey',geneScore=400,geneSelectPv=0.05,STITCH_CTdb='db/STITCH5_CTdb.db'){ 
  #ChemData:data.frame:colnames-->'cid','chemical_name','key'.'key' is smiles or inchikey 
  #STITCH_CTdb:'STITCH5_CTdb.db',colnames=c('cid','ENSP','score','geneID','inchikey') 
  #out:list(selectGeneID=selectGeneID,stitchTargetDataResult=stitchTargetDataResult) 
  require(data.table) 
  require(stringr) 
  require(RSQLite) 
  require(plyr) 
  require(ggplot2) 
  require(pracma) 
  require(readr) 
  options(stringsAsFactors = F) 
  ChemData$key=as.character(ChemData$key) 
  ChemData$cid=as.character(ChemData$cid) 
  ChemData$chemical_name=as.character(ChemData$chemical_name) 
  ChemData2=ChemData 
  if (ChemDataType!='inchikey'){ 
    tempSmiles=ChemData$key 
    names(tempSmiles)=ChemData$cid 
    ChemData2$inchikey=unlist(TransInchikey(tempSmiles,molecularType='SMItxt')) 
  }else{ 
    ChemData2=plyr::rename(ChemData2,c('key'='inchikey')) 
  } 
  ChemData2=ChemData2[!is.na(ChemData2$inchikey),] 
  ChemData2=as.data.table(ChemData2) 
  ChemData2=unique(ChemData2,by='inchikey') 
  QueryKey=unique(ChemData2$inchikey) 
  if (nrow(ChemData2)>0){ 
    stitchTargetData=findStitchTargetAll(chem_inchikey=QueryKey,scoreSet=geneScore,STITCH_CTdb=STITCH_CTdb) 
    stitchTargetData=as.data.table(stitchTargetData) 
    stitchTargetData=stitchTargetData[,.(id=cid,inchikey,geneID,score)] 
    stitchTargetDatasub=unique(stitchTargetData[,.(id,geneID)]) 
    tempR=geneScoreCal(stitchTargetDatasub,Pcutoff = geneSelectPv) 
    selectGeneID=tempR$selectGeneId###selectGeneID 
    stitchTargetDataResult=join(stitchTargetData[,.(inchikey,geneID,score)],ChemData2,by='inchikey') 
    stitchTargetDataResult=stitchTargetDataResult[,.(cid,chemical_name,inchikey,geneID,score)] 
  }else{ 
    selectGeneID=NA 
    stitchTargetDataResult=NA 
  } 
  #write_lines(selectGeneID,path='selectGeneID.txt') 
  return(list(selectGeneID=selectGeneID,stitchTargetDataResult=stitchTargetDataResult)) 
} 
 
