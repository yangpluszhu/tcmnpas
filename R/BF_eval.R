BF_eval = function(BF_result,yao){ 
   
  require(stringr) 
  require(data.table) 
  require(plyr) 
  yao=as.data.frame(yao) 
  BF=BF_result$BF_result 
  countBF=nrow(BF) 
  patientStat=matrix(0,nrow=countBF,ncol=4) 
  colnames(patientStat)=c('id','Rs0.8','Rs0.9','Rs1') 
  for (i in 1:countBF){ 
    herb=BF$makeup[i] 
    herbV=str_split(herb,',')[[1]] 
    yao_sub=yao[,c('Pid',herbV)] 
    yao_sub$Rsupport=rowSums(yao_sub[,herbV])/length(herbV) 
    yao_sub=as.data.table(yao_sub) 
    yao_sub2=yao_sub[,.(Rsmax=max(Rsupport)),by='Pid'] 
    Rs0.8=sum(yao_sub2$Rsmax>=0.8)/nrow(yao_sub2) 
    Rs0.9=sum(yao_sub2$Rsmax>=0.9)/nrow(yao_sub2) 
    Rs1=sum(yao_sub2$Rsmax==1)/nrow(yao_sub2) 
    patientStat[i,'id']=BF$id[i] 
    #patientStat[i,'N']=BF$N[i] 
    patientStat[i,'Rs0.8']=Rs0.8 
    patientStat[i,'Rs0.9']=Rs0.9 
    patientStat[i,'Rs1']=Rs1 
  } 
  patientStat=as.data.table(patientStat) 
  patientStatAll=plyr::join(BF,patientStat,by='id') 
  patientStatAll=as.data.table(patientStatAll) 
  patientStatAll=patientStatAll[,.(id,N,CBWN,S0.8,Rs0.8,S0.9,Rs0.9,S1,Rs1,makeup)] 
  patientStatAll$CBWN=as.numeric(patientStatAll$CBWN) 
  patientStatAll$S0.8=as.numeric(patientStatAll$S0.8) 
  patientStatAll$Rs0.8=as.numeric(patientStatAll$Rs0.8) 
  patientStatAll$S0.9=as.numeric(patientStatAll$S0.9) 
  patientStatAll$Rs0.9=as.numeric(patientStatAll$Rs0.9) 
  patientStatAll$S1=as.numeric(patientStatAll$S1) 
  patientStatAll$Rs1=as.numeric(patientStatAll$Rs1) 
  patientStatAll$CBWN=round(patientStatAll$CBWN,4) 
  patientStatAll$S0.8=round(patientStatAll$S0.8,4) 
  patientStatAll$Rs0.8=round(patientStatAll$Rs0.8,4) 
  patientStatAll$S0.9=round(patientStatAll$S0.9,4) 
  patientStatAll$Rs0.9=round(patientStatAll$Rs0.9,4) 
  patientStatAll$S1=round(patientStatAll$S1,4) 
  patientStatAll$Rs1=round(patientStatAll$Rs1,4) 
  return(patientStatAll) 
} 
 
