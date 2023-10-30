basic_trans = function(prescription,mine_data,local=TRUE){
  ##prescription-->"id","time", "patient","diagnosis", "doctor", "drug" +"drug_trans"[innComponentName of mine_data]
  ##mine_data--->"num","drug_id","mine_id","mine_Cname", "drug_commonName", "mine_name" ,"class","fufang", "jixing","guige""生产企业"         "规格", "jixing2" , "innComponentName","type"合并了<--herb_data
  ##output:prescription-->"id","time", "patient","diagnosis", "doctor","drug"[binary]
  require(stringr)
  require(data.table)
  require(lubridate)
  require(Matrix)
  if (local){
  hf_name=unique(mine_data$mine_name[mine_data$class=='hf'])
  herb_name=unique(mine_data$mine_name[mine_data$class=='h'])
  presc_data=numeric()
  length(presc_data)=6
  names(presc_data)=c("id","time", "patient","diagnosis", "doctor", "drug")
  for (i in 1:nrow(prescription)){
    temp1=numeric()
    length(temp1)=6
    names(temp1)=c("id","time", "patient","diagnosis", "doctor", "drug")
    tDrug=prescription$drug[i]
    tDrug=str_split(tDrug,',')[[1]]
    if (sum(tDrug%in%hf_name)>0){
      tDrug_hf=tDrug[tDrug%in%hf_name]
      for (j in 1:length(tDrug_hf)){
        tDrug=str_replace_all(tDrug,tDrug_hf[j],mine_data$innComponentName[match(tDrug_hf[j],mine_data$mine_name)])
      }
    }
    tDrug=paste(tDrug,collapse = ',')
    temp1["id"]=prescription$id[i]
    temp1["time"]=prescription$time[i]
    temp1["patient"]=prescription$patient[i]
    temp1["diagnosis"]=prescription$diagnosis[i]
    temp1["doctor"]=prescription$doctor[i]
    temp1["drug"]=tDrug
    presc_data=rbind(presc_data,temp1)
  }
  presc_data=presc_data[-1,,drop=F]
  rownames(presc_data)=NULL
  presc_data=as.data.frame(presc_data,stringsAsFactors = F)
  yao=unique(str_split(paste(presc_data$drug,collapse = ','),',')[[1]])
  tData=Matrix(0,nrow=nrow(presc_data),ncol=length(yao),sparse=T)
  colnames(tData)=yao
  #rownames(tData)=presc_data$id
  for (k in 1:nrow(presc_data)){
    kt=unique(str_split(presc_data$drug[k],',')[[1]])
    #tData[k,'id']=presc_data$id[k]
    tData[k,kt]=1
  }
  tData=as.data.frame(as.matrix(tData),stringsAsFactors = F)
  tData$id=presc_data$id
  presc_data_sub=subset(prescription,select=c("id","time", "patient","diagnosis", "doctor"))
  prescription_basic=merge(presc_data_sub,tData,by='id')
#   if (length(yao)==1){
#     prescription_basic[,yao]=as.numeric(prescription_basic[,yao])
#   }else{
#     prescription_basic[,yao]=apply(prescription_basic[,yao],2,as.numeric)
#   }
  ##########not_manual####
  }else{
    hf_name=unique(mine_data$mine_Cname[mine_data$class=='hf'])
    herb_name=unique(mine_data$mine_Cname[mine_data$class=='h'])
    presc_data=numeric()
    length(presc_data)=6
    names(presc_data)=c("id","time", "patient","diagnosis", "doctor", "drug")
    for (i in 1:nrow(prescription)){
      temp1=numeric()
      length(temp1)=6
      names(temp1)=c("id","time", "patient","diagnosis", "doctor", "drug")
      tDrug=prescription$drug[i]
      tDrug=str_split(tDrug,',')[[1]]
      if (sum(tDrug%in%hf_name)>0){
        tDrug_hf=tDrug[tDrug%in%hf_name]
        for (j in 1:length(tDrug_hf)){
          tDrug=str_replace_all(tDrug,tDrug_hf[j],mine_data$innComponentName[match(tDrug_hf[j],mine_data$mine_Cname)])
        }
      }
      tDrug=paste(tDrug,collapse = ',')
      temp1["id"]=prescription$id[i]
      temp1["time"]=prescription$time[i]
      temp1["patient"]=prescription$patient[i]
      temp1["diagnosis"]=prescription$diagnosis[i]
      temp1["doctor"]=prescription$doctor[i]
      temp1["drug"]=tDrug
      presc_data=rbind(presc_data,temp1)
    }
    presc_data=presc_data[-1,,drop=F]
    rownames(presc_data)=NULL
    presc_data=as.data.frame(presc_data,stringsAsFactors = F)
    yao=unique(str_split(paste(presc_data$drug,collapse = ','),',')[[1]])
    tData=Matrix(0,nrow=nrow(presc_data),ncol=length(yao),sparse = T)
    colnames(tData)=yao
    for (k in 1:nrow(presc_data)){
      kt=unique(str_split(presc_data$drug[k],',')[[1]])
      #tData[k,'id']=presc_data$id[k]
      tData[k,kt]=1
    }
    tData=as.data.frame(as.matrix(tData),stringsAsFactors = F)
    tData$id=presc_data$id
    presc_data_sub=subset(prescription,select=c("id","time", "patient","diagnosis", "doctor"))
    prescription_basic=merge(presc_data_sub,tData,by='id')
#     if (length(yao)==1){
#       prescription_basic[,yao]=as.numeric(prescription_basic[,yao])
#     }else{
#       prescription_basic[,yao]=apply(prescription_basic[,yao],2,as.numeric)
#     }
  }
  #prescription_basic=Matrix(prescription_basic,sparse=T)
  return(list(prescription_basic=prescription_basic,yao=yao))
}

