QEDcal = function(druglikeDesData){ 
  #colnames(druglikeDesData)=c('id',colnames(QED_par)) 
  require(ChemmineR) 
  require(ChemmineOB) 
  require(data.table) 
  require(plyr) 
  para=fread('db/QED_par.csv') 
  para=as.data.frame(para) 
#druglike=read.csv(druglikeDesFile,head=T,stringsAsFactors=F,colClasses=c('character',rep('numeric',8))) 
  druglike=druglikeDesData 
  #if (!setequal(setdiff(colnames(druglike),'id'),colnames(para))){ 
    #NamesTrans=c('ALOGP','MW','ROTB','AROM','HBA','HBD','PSA') 
  #names(NamesTrans)=c('ALogP','Molecular_Mass','Num_RotatableBonds','Num_AromaticRings','Num_H_Acceptors_Lipinski','Num_H_Donors_Lipinski','Molecular_PolarSurfaceArea') 
  #druglike=rename(druglike,NamesTrans) 
  #print('ERROR:colnames') 
  #stop() 
  #} 
  druglike=druglike[,c('id',colnames(para))] 
  rownames(druglike)=druglike$id 
  #para=read.table('clipboard',head=T) 
  rownames(para)=letters[1:7] 
  druglike_des=matrix(0,nrow=nrow(druglike),ncol=ncol(para)) 
  colnames(druglike_des)=colnames(para) 
  rownames(druglike_des)=rownames(druglike) 
  for (i in colnames(para)){ 
    for (j in rownames(druglike)){ 
      druglike_des[j,i]=(para['a',i]+para['b',i]/(1+exp(-(druglike[j,i]-para['c',i]+para['d',i]/2)/para['e',i]))*(1-1/(1+exp(-(druglike[j,i]-para['c',i]-para['d',i]/2)/para['f',i]))))/para['g',i] 
    } 
  } 
  druglike_des=as.data.frame(druglike_des) 
  druglike_des$QED=exp(rowSums(log(druglike_des[,colnames(para)]))/8) 
  colnames(druglike_des)=paste(colnames(druglike_des),'_DES',sep='') 
  return(druglike_des) 
} 
 
