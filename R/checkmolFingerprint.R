checkmolFingerprint = function(SDFfile){
  library("ChemmineR")
  library(stringr)
  library(plyr)
  library("reshape")
  #setwd('d:/R/checkmol/')
  if (!dir.exists('tempSDF')){
    dir.create('tempSDF')
  }
  #s=read.SDFset('combine_three_final_filt_prepared.sdf')#####sdf文件名
  s=read.SDFset(SDFfile)
  #################check#########
  # index=which(validSDF(s)!=TRUE)
  # idx=numeric()
  # for (i in 1:9){
  #   idx[i]=datablock(s)[[index[i]]]['id']
  # }
  ###########################
  odlwd=getwd()
  #setwd(paste(odlwd,'tempSDF',sep='/'))
  filename=vector(mode='character',length=length(s))
  a=NULL
  for (i in 1:length(s)){
    tempname=datablock(s)[[i]]['id']######id标识列名作为文件名
    filename[i]=paste(tempname,'.sdf',sep='')
    write.SDF(s[i], file=filename[i])
    cmd=paste('checkmol -c',filename[i],sep=' ');
    a[[tempname]]=shell(cmd,inter=T)
    file.remove(filename[i])
  }
  aa=sapply(a,str_split,pattern=';')
  melt_aa=melt(aa)
  melt_aa=melt_aa[-which(melt_aa[,'value']==''),]
  melt_aa$value=as.character(melt_aa$value)
  check_finger=cast(melt_aa,L1~value,length)####L1列为id标识列length用于计数
  id=str_replace_all(filename,pattern='.sdf',replacement='')
  zero_id=id[which(!id%in%check_finger$L1)]
  zero_finger=data.frame(L1=zero_id,matrix(0,nrow=length(zero_id),ncol=ncol(check_finger)-1))
  colnames(zero_finger)=colnames(check_finger)
  check_finger_todal=rbind(check_finger,zero_finger)#L1列为id标识列length用于计数
  check_finger_todal=rename(check_finger_todal,c('L1'='id'))
  setwd(odlwd)
  return(check_finger_todal)
}

