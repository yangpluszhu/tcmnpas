getChemGeneFromDatabaseGUI = function(化合物数据=NULL,化合物结构类型='inchikey',检索药靶数据库HIT='是',检索药靶数据库TCMID='是',检索药靶数据库STITCH='是',检索药靶数据库TCMSP='是',检索药靶数据库CUSTOM='是',类药性指数阈值=0.2,基因关联得分=400,基因关联显著性阈值=0.05,Savefilename=NULL,Lipinski_Viosset=0,Veber_Viosset=0){
  #getChemGeneFromDatabase(ChemData,ChemDataType='inchikey',QEDset=0.2,geneScore=400,geneSelectPv=0.05,targetDatabase=c('HIT','TCMID','STITCH','TCMSP','CUSTOM'),Lipinski_Viosset=0,Veber_Viosset=0){
    #ChemData:data.frame:colnames-->'cid','chemical_name','key'.'key' is smiles or inchikey
    #ChemDataType:'inchikey' or 'smiles'
    #result:list(combinedGeneID=selectGeneID,GeneID_split=tempR_split):$GeneID_split-->data.frame(cid,chemical_name,inchikey,geneID)
  require(data.table)
  require(stringr)
  require(RSQLite)
  require(plyr)
  require(ggplot2)
  require(pracma)
  require(rio)
  require(r2excel)
  require(readr)
  options(stringsAsFactors = F)
  ChemData=化合物数据
  ChemData=import(ChemData)
  ChemDataType=化合物结构类型
  HIT=检索药靶数据库HIT
  TCMID=检索药靶数据库TCMID
  STITCH=检索药靶数据库STITCH
  TCMSP=检索药靶数据库TCMSP
  CUSTOM=检索药靶数据库CUSTOM
  QEDset=类药性指数阈值
  geneScore=基因关联得分
  geneSelectPv=基因关联显著性阈值
  targetDatabase=c('HIT','TCMID','STITCH','TCMSP','CUSTOM')
  if (HIT=='是'){
    HIT=T
  }else{
    HIT=F
  }
  #######
  if (TCMID=='是'){
    TCMID=T
  }else{
    TCMID=F
  }
  #######
  #######
  if (STITCH=='是'){
    STITCH=T
  }else{
    STITCH=F
  }
  #######
  if (TCMSP=='是'){
    TCMSP=T
  }else{
    TCMSP=F
  }
  #######
  if (CUSTOM=='是'){
    CUSTOM=T
  }else{
    CUSTOM=F
  }
  #######
  #######
  if (is.null(Savefilename))Savefilename='化合物靶标库靶标检索结果'
  targetDatabase=targetDatabase[c(HIT,TCMID,STITCH,TCMSP,CUSTOM)]
  ######################
  result=getChemGeneFromDatabase(ChemData=ChemData,ChemDataType=ChemDataType,QEDset=QEDset,geneScore=geneScore,geneSelectPv=geneSelectPv,targetDatabase=targetDatabase,Lipinski_Viosset=Lipinski_Viosset,Veber_Viosset=Veber_Viosset)
  write_lines(result$combinedGeneID,path='Report/ChemselectGeneID合并结果.txt')
  write.csv(result$GeneID_split, paste('Report/',Savefilename,'.csv',sep=''),row.names = F)
  write.csv(result$subD3, paste('Report/',Savefilename,'_供合并.csv',sep=''),row.names = F)
  ##########
  txt=paste('化合物靶标库靶标检索结果已保存在Report文件夹中的',Savefilename,'.csv',sep='')
  outputfile=file('OutPut_message.txt','wb')
  writeChar(txt,outputfile)
  close(outputfile)
  shell('OutPut_message.txt',wait=FALSE)
}

