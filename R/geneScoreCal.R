geneScoreCal = function(chem_target,Pcutoff){
  #chem_target:data.frame-->colnames(chem_target)=c('id','geneID')
  #Pcutoff:筛选的P阈值
  library(data.table)
  library(plyr)
  library(stringr)
  options(stringsAsFactors = F)
  chem_target=as.data.frame(chem_target)
  chem_target$geneID=as.character(chem_target$geneID)
  TC_num=mean(table(chem_target$geneID))####gene被靶平均化合物数
  #C_num=length(unique(chem_target$id))
  C_num=nrow(chem_target)
  G_C_num=table(chem_target$geneID)
  geneID=names(G_C_num)
  geneP=numeric()
  length(geneP)=length(geneID)
  names(geneP)=geneID
  for (i in geneID){
    tempP=binom.test(G_C_num[i],C_num,TC_num/C_num,alternative = 'greater')$p.value
    geneP[i]=tempP
  }
  geneRank=rank(geneP)
  geneP2=ifelse(geneP<Pcutoff,geneP,1)
  geneScore=(-1*log10(geneP2))/geneRank
  names(geneScore)=names(geneP)
  geneScore=sort(geneScore,decreasing=T)
  geneScore=geneScore[geneScore>0]
  tempID=sort(geneP)
  selectGeneId=names(tempID)[tempID<Pcutoff]###length(selectGeneId)
  return(list(geneScore=geneScore,selectGeneId=selectGeneId))
}

