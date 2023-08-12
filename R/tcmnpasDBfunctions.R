AdjNetScore = function(TotalScore,NumSeedGene,NumScoreGene,multiple=1000){
#ScoreGeneData$AdjTotalScore=KATZresult$TotalScore/KATZresult$NumSeedGeneID*log10(1+KATZresult$NumSeedGeneID/KATZresult$NumScoreGeneID/2)*1000
Score=TotalScore/NumSeedGene*log10(1+NumSeedGene/NumScoreGene/2)*multiple
return(Score)
}

alertCal = function(sdfFile){
##sdfFile有id字段！
  library(ChemmineR)
  library(ChemmineOB)
  unwant=read.csv('db/unwanted_structure.csv',head=T,stringsAsFactors=F)
  sdf=read.SDFset(sdfFile)
  valid <- validSDF(sdf)
  #sum(valid)
  # datablocktag(sdf[which(!valid)],'id')
  sdf <- sdf[valid]
  cid(sdf)=datablocktag(sdf,'id')
  alert=matrix(0,nrow=length(sdf),ncol=nrow(unwant))
  rownames(alert)=datablocktag(sdf,'id')
  colnames(alert)=unwant$Name
  for (i in cid(sdf)){
    for (j in 1:nrow(unwant)){
      alert[i,j]=smartsSearchOB(sdf[i],unwant$smarts[j])
    }
  }
  alert_value=rowSums(alert)
  alert=as.data.frame(alert)
  alert$id=datablocktag(sdf,'id')
  alert$ALERTS=alert_value
  return(alert)
}

asSymmetric = function(x,rule='upper'){
  require(sna)
  y=symmetrize(x,rule=rule)
  detach('package:sna',unload=TRUE)
  rownames(y)=rownames(x)
  colnames(y)=colnames(x)
  return(y)
}

BF_analysis = function(Pdata,tS0.9=0.03,min_yao=8,threshold_leastYao=10,exact=FALSE,combineSim=FALSE,combine_method=2,combine_threshold=0.7,csvInput=T){
  ##P--->matrix,colnames=yao,rownames=prescription$id OR csv[colname:id,herb]
  ##S0.9筛选阈值
  ##threshold_leastYao-->自适应筛选时BF至少含药数
  ##exact-->自适应筛选时BF含药数是否恰好等于threshold_leastYao
  ##Pdata-->含Pid，Vid，herb列 OR P--->matrix
  ##output:list$BF_result:makeup+statstatistic;list$yao_fre
  require(igraph)
  require(plyr)
  require(data.table)
  options(stringsAsFactors = F)
  if (csvInput){
    if ('character'%in%class(Pdata)){
      tempP0=fread(Pdata,encoding='UTF-8')
    }else if ('data.frame'%in%class(Pdata)){
      tempP0=Pdata
      tempP0=as.data.table(tempP0)
    }
    if (!'Vid'%in%colnames(tempP0))tempP0$Vid=0
    tempP0=tempP0[!(is.na(tempP0$Pid)|is.na(tempP0$herb)|is.na(tempP0$Vid)),]
    tempHerbName=unique(tempP0$herb)
    tempP0=data.table::dcast(Pid+Vid~herb,data=tempP0,length)
    tempP=tempP0
    tempP$id=paste(tempP$Pid,tempP$Vid,sep='_')
    P=as.matrix(tempP[,tempHerbName,with=F])
    P[P>0]<-1
    rownames(P)=as.character(tempP$id)
  }else{
    P=Pdata
    P[P>0]<-1
  }
  #P[!is.na(P)]<-as.numerc(1)
  #P[P>0]<-1
  na_table=data.frame(id='ALL',result='None')
  if (nrow(P)==1){
    BF_result=colnames(P)
    yao_fre=colSums(P)
    names(yao_fre)=colnames(P)
    k_opt=NULL
    Degree_yao2=data.frame(drug=colnames(P),Degree=1)
  }else{
    yao=colnames(P)
    yao_yao=t(P)%*%P
    yao_ya0_Binary=yao_yao
    diag(yao_ya0_Binary)=0
    yao_ya0_Binary[yao_ya0_Binary>=1]=1
    Degree_yao=colSums(yao_ya0_Binary)
    Degree_yao2=data.frame(drug=names(Degree_yao),Degree=Degree_yao)
    Degree_yao2=arrange(Degree_yao2,desc(Degree))
    yao_fre=colSums(P)#####药物频率
    names(yao_fre)=colnames(P)
    diag(yao_yao)=0
    yao_yaoNetwork=graph_from_adjacency_matrix(yao_yao,mode='undirected',weighted=T)
    rowS=rowSums(P)
    colS=colSums(P)
    ###########筛选二值化阈值
    kseq=seq(0.01,0.2,by=0.005)
    num_BF=numeric(length(kseq))
    num_yaoBF=numeric(length(kseq))
    for (k in 1:length(kseq)){
      yao_yao_g=yao_yao#####yao_yao_g
      yao_yao_g[which(yao_yao_g<nrow(P)*kseq[k])]=0#####yao_yao_g
      yao_yao_g[which(yao_yao_g>=nrow(P)*kseq[k])]=1#####yao_yao_g
      g2=graph.adjacency(yao_yao_g,mode='undirected',add.rownames=T,diag=F)
      g2_clique=maximal.cliques(g2,min=min_yao)
      g2_clique2=lapply(g2_clique,function(x) sort(as.numeric(unlist(x))))
      if (length(g2_clique2)!=0){
        identical(yao,V(g2)$name)
        #source('combine_sim.R')
        if (combineSim){
          combine=combine_sim(g2_clique2,method=combine_method,threshold=combine_threshold)
          Bf=combine$combine######
        }else{
          Bf=g2_clique2
        }
        ####
        #Bf=g2_clique2
        #####
        names(Bf)=paste('BF',1:length(Bf),sep='')
        bf_name=list()
        for (i in names(Bf)){
          bf_name[[i]]=c(yao[Bf[[i]]])###############bf_name
        }
        ##########
        bf_stat=list()
        for (i in names(Bf)){
          tem=P[,Bf[[i]]]
          tem_rowS=rowSums(tem)/length(Bf[[i]])
          CBWN=mean(tem_rowS)
          S0.8=sum(tem_rowS>=0.8)/nrow(P)
          S0.9=sum(tem_rowS>=0.9)/nrow(P)
          S1=sum(tem_rowS==1)/nrow(P)
          bf_stat[[i]]=c(CBWN=CBWN,S0.8=S0.8,S0.9=S0.9,S1=S1) ######bf_stat
        }
        ###按S1阈值筛选BF
        tS0.9=tS0.9
        for (i in names(Bf)){
          if (bf_stat[[i]]['S0.9']<tS0.9){
            bf_name[[i]]=NULL
            bf_stat[[i]]=NULL
          }
        }
        num_BF[k]=length(bf_name)
		if (exact){
		num_yaoBF[k]=sum(sapply(bf_name,length)==threshold_leastYao)###统计大于等于10个药的基本方数
		}else{
		num_yaoBF[k]=sum(sapply(bf_name,length)>=threshold_leastYao)###统计大于等于10个药的基本方数
		}        
      }else{
        num_BF[k]=0
        num_yaoBF[k]=0
      }
    }
    if (sum(num_BF)==0){
      k_opt='none'
      BF_result=na_table
    }else{
      k_opt=tail(kseq[num_yaoBF==max(num_yaoBF)],1)
      k=k_opt
      yao_yao_g=yao_yao#####yao_yao_g
      yao_yao_g[which(yao_yao_g<nrow(P)*k)]=0#####yao_yao_g
      yao_yao_g[which(yao_yao_g>=nrow(P)*k)]=1#####yao_yao_g
      g2=graph.adjacency(yao_yao_g,mode='undirected',add.rownames=T,diag=F)
      g2_clique=maximal.cliques(g2,min=min_yao)
      g2_clique2=lapply(g2_clique,function(x) sort(as.numeric(unlist(x))))
      #identical(yaoname,V(g2)$name)
      if (combineSim){
        combine=combine_sim(g2_clique2,method=combine_method,threshold=combine_threshold)
        Bf=combine$combine######
      }else{
        Bf=g2_clique2
      }
      #source('combine_sim.R')
      #combine=combine_sim(g2_clique2,method=2,threshold=0.6)
      #Bf=combine$combine######Bf
      ####
      #Bf=g2_clique2
      #####
      if (length(Bf)==0){
        k_opt='none'
        BF_result=na_table
      }else{
        names(Bf)=paste('BF',1:length(Bf),sep='')
        bf_name=list()
        for (i in names(Bf)){
          bf_name[[i]]=c(yao[Bf[[i]]])###############bf_name
        }
        ##########
        bf_stat=list()
        for (i in names(Bf)){
          tem=P[,Bf[[i]]]
          tem_rowS=rowSums(tem)/length(Bf[[i]])
          CBWN=mean(tem_rowS)
          S0.8=sum(tem_rowS>=0.8)/nrow(P)
          S0.9=sum(tem_rowS>=0.9)/nrow(P)
          S1=sum(tem_rowS==1)/nrow(P)
          bf_stat[[i]]=c(CBWN=CBWN,S0.8=S0.8,S0.9=S0.9,S1=S1) ######bf_stat
        }
        ###按S1阈值筛选BF
        tS0.9=tS0.9
        for (i in names(Bf)){
          if (bf_stat[[i]]['S0.9']<tS0.9){
            bf_name[[i]]=NULL
            bf_stat[[i]]=NULL
          }
        }
        ###rename##
        BF_name=paste('BF',1:length(bf_name),sep='')
        names(bf_name)=BF_name
        names(bf_stat)=BF_name
        BF_result=ldply(bf_stat)
        colnames(BF_result)=c('id','CBWN','S0.8','S0.9','S1')
        for (i in 1:nrow(BF_result)){
          tempName=BF_result$id[i]
          BF_result$makeup[i]=paste(bf_name[[tempName]],collapse = ',')
          BF_result$N[i]=length(unique(bf_name[[tempName]]))
        }
      BF_result=BF_result[,c('id','N','CBWN','S0.8','S0.9','S1','makeup')]
      BF_result[,'CBWN']=round(BF_result[,'CBWN'],digits = 4)
      BF_result[,'S0.8']=round(BF_result[,'S0.8'],digits = 4)
      BF_result[,'S0.9']=round(BF_result[,'S0.9'],digits = 4)
      BF_result[,'S1']=round(BF_result[,'S1'],digits = 4)
      }

    }
  }
  #BF_result=BF_result[,c('id','N','CBWN','S0.8','S0.9','S1','makeup')]
  return(list(BF_result=BF_result,yao_fre=data.frame(Name=names(yao_fre),Fre=yao_fre),opt_a=k_opt,Degree_yao=Degree_yao2,BFnetwork=g2,yao_yaoNetwork=yao_yaoNetwork,dataP0=tempP0))
}

BF_analysis_shiny = function(Pdata,tS0.9=0.03,min_yao=8,threshold_leastYao=10,exact=FALSE,combineSim=FALSE,combine_method=2,combine_threshold=0.7,csvInput=T,calRs=F){
  #method：1一起合并，2逐个合并
  ##P--->matrix,colnames=yao,rownames=prescription$id OR csv[colname:id,herb]
  ##S0.9筛选阈值
  ##threshold_leastYao-->自适应筛选时BF至少含药数
  ##exact-->自适应筛选时BF含药数是否恰好等于threshold_leastYao
  ##Pdata-->含Pid，Vid，herb列 OR P--->matrix
  ##output:list$BF_result:makeup+statstatistic;list$yao_fre
  require(igraph)
  require(plyr)
  require(data.table)
  require(stringr)
  options(stringsAsFactors = F)
  BF_analysisOB=BF_analysis(Pdata=Pdata,tS0.9=tS0.9,min_yao=min_yao,threshold_leastYao=threshold_leastYao,exact=exact,combineSim=combineSim,combine_method=combine_method,combine_threshold=combine_threshold,csvInput=csvInput)
  ##return(list(BF_result=BF_result,yao_fre=yao_fre,opt_a=k_opt,Degree_yao=Degree_yao2,network=g2,dataP0=tempP0))
  if (nrow(BF_analysisOB$BF_result)>1){
    coMemberMatrixResult=coMemberMatrix(BF_analysisOB)
  }else{
    coMemberMatrixResult=NULL
  }
  #####
  if (calRs){
    BF_evalResult=BF_eval(BF_result=BF_analysisOB,yao=BF_analysisOB$dataP0)
  }else{
    BF_evalResult=NULL
  }
  return(list(BF_result=BF_analysisOB$BF_result,yao_fre=BF_analysisOB$yao_fre,opt_a=BF_analysisOB$opt_a,Degree_yao=BF_analysisOB$Degree_yao,BFnetwork=BF_analysisOB$BFnetwork,dataP0=BF_analysisOB$dataP0,BF_evalResult=BF_evalResult,coMemberMatrixResult=coMemberMatrixResult,yao_yaoNetwork=BF_analysisOB$yao_yaoNetwork))
}

BF_eval = function(BF_result,yao){
  #BF_result-->由BF_analysis得到，yao-->data.frame,列名包含'Pid','Vid'与药名
  ##计算基于'patient'的统计量
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

checkNodeType = function(id){
  library(stringr)
  idTrans=as.numeric(id)
  idTrans2=str_detect(id,'c')
  geneType=rep('gene',length(id))
  chemType=rep('chem',length(id))
  herbType=rep('herb',length(id))
  nodeType=ifelse(is.na(idTrans),ifelse(idTrans2,chemType,herbType),geneType)
  return(nodeType)
}

Chem2DSim = function(SMI){
  #SMI:data.frame-->colnames(SMI)=c('id','smiles')
  #calMolecularDes--->set padel_types.xml!
  library(data.table)
  library(plyr)
  library(stringr)
  library(igraph)
  library(ggplot2)
  library(ChemmineR)
  library(PaDEL)
  library(fingerprint)
  if (!dir.exists('MOL')){
    dir.create('MOL')
  }
  SMI=as.data.table(SMI)
  SMI2=SMI[smiles!='',]
  smiles=SMI2$smiles
  names(smiles)=SMI2$id
  SMISdfSets=ChemmineR::smiles2sdf(smiles)
  SMISdfSets=SMISdfSets[validSDF(SMISdfSets)]
  #sum(validSDF(SMISdfSets))==length(SMISdfSets)
  write.SDF(SMISdfSets,'MOL/SMISdfSets.sdf',cid=T)
  FP=CalDescriptors('MOL/SMISdfSets.sdf')
  FPv=setdiff(colnames(FP),'Name')
  FPlist=list()
  CID=FP$Name
  for (i in 1:nrow(FP)){
    temp=FP[i,FPv,with=F]
    tempObject=new("fingerprint", nbit=length(FPv), bits=which(temp!=0))
    FPlist=c(FPlist,tempObject)
  }
  length(FPlist)==length(CID)
  ChemSimMatrix=fp.sim.matrix(FPlist,method='tanimoto')
  rownames(ChemSimMatrix)=CID
  colnames(ChemSimMatrix)=CID
  return(list(ChemSimMatrix=ChemSimMatrix,FP=FP))
}

Chem3DSim = function(SMI){
  #SMI:data.frame-->colnames(SMI)=c('id','smiles')
  library(data.table)
  library(plyr)
  library(stringr)
  library(clipr)
  library(igraph)
  library(ggplot2)
  library(ChemmineR)
  if (!dir.exists('MOL'))dir.create('MOL')
  SMI=as.data.table(SMI)
  SMI2=SMI[smiles!='',]
  smiles=SMI2$smiles
  names(smiles)=SMI2$id
  SMISdfSets=ChemmineR::smiles2sdf(smiles)
  SMISdfSets=SMISdfSets[validSDF(SMISdfSets)]
  #sum(validSDF(SMISdfSets))==length(SMISdfSets)
  datablock(SMISdfSets)<-data.frame(id=SMISdfSets@ID)
  write.SDF(SMISdfSets,'MOL/SMISdfSets.sdf',cid=T)
  ###
  setwd(paste(getwd(),'/MOL',sep=''))
  TransMol('SMISdfSets.sdf')
  DrugId=unique(SMISdfSets@ID)
  ChemShapeSimMatrix=matrix(0,nrow=length(DrugId),ncol=length(DrugId))
  colnames(ChemShapeSimMatrix)=DrugId
  rownames(ChemShapeSimMatrix)=DrugId
  ChemFeatureSimMatrix=matrix(0,nrow=length(DrugId),ncol=length(DrugId))
  colnames(ChemFeatureSimMatrix)=DrugId
  rownames(ChemFeatureSimMatrix)=DrugId
  ChemHybridSimMatrix=matrix(0,nrow=length(DrugId),ncol=length(DrugId))
  colnames(ChemHybridSimMatrix)=DrugId
  rownames(ChemHybridSimMatrix)=DrugId
  for (i in 1:length(DrugId)){
    cid1=DrugId[i]
    molFilename=paste(cid1,'mol2',sep='.')
    print(i)
    j=i+1
    while (j<=length(DrugId)){
      cid2=DrugId[j]
      molFilename2=paste(cid2,'mol2',sep='.')
      dos=paste('Cynthia -q',molFilename,'-t',molFilename2)
      shell(dos,intern=T)
      if (class(try(fread('Result.list'),silent = T))=='try-error'){
        ChemHybridSimMatrix[cid1,cid2]=NA
        ChemShapeSimMatrix[cid1,cid2]=NA
        ChemFeatureSimMatrix[cid1,cid2]=NA
      }else{
        tempResult=fread('Result.list')
        if (nrow(tempResult)==0){
          ChemHybridSimMatrix[cid1,cid2]=NA
          ChemShapeSimMatrix[cid1,cid2]=NA
          ChemFeatureSimMatrix[cid1,cid2]=NA
        }else{
          ChemHybridSimMatrix[cid1,cid2]=tempResult$HybridScore
          ChemShapeSimMatrix[cid1,cid2]=tempResult$ShapeScore
          ChemFeatureSimMatrix[cid1,cid2]=tempResult$FeatureScore
        }
      }
      j=j+1
    }
  }
  setwd(str_replace(getwd(),'/MOL',''))
  ChemHybridSimMatrix=asSymmetric(ChemHybridSimMatrix)
  ChemShapeSimMatrix=asSymmetric(ChemShapeSimMatrix)
  ChemFeatureSimMatrix=asSymmetric(ChemFeatureSimMatrix)
 return(list(ChemHybridSimMatrix=ChemHybridSimMatrix,ChemShapeSimMatrix=ChemShapeSimMatrix,ChemFeatureSimMatrix=ChemFeatureSimMatrix))
}

chemScoreCal = function(chem_target,geneScore,selectGeneId,method=1,cutoffValue=0){
  #chem_target:data.frame-->colnames(chem_target)=c('id','geneID')
  #geneScore,selectGeneId-->由geneScoreCal得到
  #method筛选方法：1--》使用cutoffValue；2--》按照最小100%覆盖selectGeneId
  #cutoffValue：使用method1时的阈值,使用method2时的覆盖百分比
  library(data.table)
  library(plyr)
  library(stringr)
  options(stringsAsFactors = F)
  chem_target=as.data.frame(chem_target)
  chem_target$geneID=as.character(chem_target$geneID)
  chemID=unique(chem_target$id)
  geneID=unique(chem_target$geneID)
  chemScore=numeric()
  length(chemScore)=length(chemID)
  names(chemScore)=chemID
  IdentitySelectGene=rep(0,length=length(geneID))
  names(IdentitySelectGene)=geneID
  IdentitySelectGene[selectGeneId]=1
  for (i in chemID){
    tempD=chem_target[chem_target$id==i,]
    tempGene=as.character(tempD$geneID)
    if (sum(tempGene%in%selectGeneId)>0){
      tempGene=tempGene[tempGene%in%selectGeneId]
      chemScore[i]=mean(geneScore[tempGene]*IdentitySelectGene[tempGene])
    }else{
      chemScore[i]=0
    }
  }
  chemScore=sort(chemScore,decreasing = T)
  if (method==1){
    selectChemId=names(chemScore)[chemScore>cutoffValue]
    return(list(chemScore=chemScore,selectChemId=selectChemId))
  }else{
    n=1
    percent=0
    while (percent<cutoffValue){
      percent=sum(unique(chem_target$geneID[chem_target$id%in%names(chemScore)[1:n]])%in%selectGeneId)/length(selectGeneId)
      n=n+1
    }
    selectChemId=names(chemScore)[1:n]
    return(list(chemScore=chemScore,selectChemId=selectChemId))
  }
}

ChemTargetNetScore = function(ChemData,diseaseGeneID,Methodnet='KATZ',IFcombine=F){
  #ChemData:data.frame:colnames-->'cid','geneID'
  #Methodnet:'KATZ' or 'RW'
  #IFcombine:F-->split;T-->combined
  #output:ScoreGeneData-->data.frame:
  #(1)IFcombine=T:data.frame(SeedGeneIDPPI=0,NumSeedGeneID=0,ScoreGeneIDPPI=0,NumScoreGeneID=0,overlapScore=0,Path1Score=0,Path2Score=0,Path3Score=0,TotalScore=0,AdjTotalScore=0)
  #(2)IFcombine=F:data.frame(cid=0,SeedGeneIDPPI=0,NumSeedGeneID=0,ScoreGeneIDPPI=0,NumScoreGeneID=0,overlapScore=0,Path1Score=0,Path2Score=0,Path3Score=0,TotalScore=0,AdjTotalScore=0)
  require(data.table)
  require(stringr)
  require(RSQLite)
  require(plyr)
  require(ggplot2)
  require(pracma)
  options(stringsAsFactors = F)
  ChemData$geneID=as.character(ChemData$geneID)
  if (IFcombine){
    queryGeneID=unique(ChemData$geneID)
    NetCorrResult=NetCorr(SeedGeneID=diseaseGeneID,ScoreGeneID=queryGeneID,Methodnet=Methodnet,randomtimes=1,gamma=0.7,testM=2)
    ScoreGeneData=data.frame(SeedGeneIDPPI=0,NumSeedGeneID=0,ScoreGeneIDPPI=0,NumScoreGeneID=0,overlapScore=0,Path1Score=0,Path2Score=0,Path3Score=0,TotalScore=0,AdjTotalScore=0)
    if (Methodnet=='KATZ'){
      KATZresult=NetCorrResult$resultKATZ
      #KATZresult=data.frame(SeedGeneIDPPI=toString(SeedGeneID),NumSeedGeneID=length(SeedGeneID),ScoreGeneIDPPI=toString(ScoreGeneID),NumScoreGeneID=length(ScoreGeneID),overlapScore=0,Path1Score=0,Path2Score=0,Path3Score=0,TotalScore=0,RandomScoreMedian=0,Pvalue=1)
      ScoreGeneData$SeedGeneIDPPI=KATZresult$SeedGeneIDPPI
      ScoreGeneData$NumSeedGeneID=KATZresult$NumSeedGeneID
      ScoreGeneData$ScoreGeneIDPPI=KATZresult$ScoreGeneIDPPI
      ScoreGeneData$NumScoreGeneID=KATZresult$NumScoreGeneID
      ScoreGeneData$overlapScore=KATZresult$overlapScore
      ScoreGeneData$Path1Score=KATZresult$Path1Score
      ScoreGeneData$Path2Score=KATZresult$Path2Score
      ScoreGeneData$Path3Score=KATZresult$Path3Score
      ScoreGeneData$TotalScore=KATZresult$TotalScore
      ##score/n0*exp(score/n1-1)*1000
      #ScoreGeneData$AdjTotalScore=KATZresult$TotalScore/KATZresult$NumSeedGeneID*exp(KATZresult$TotalScore/KATZresult$NumScoreGeneID-1)*1000
      #score/n0*log10(1+n0/n1/2)
      ScoreGeneData$AdjTotalScore=KATZresult$TotalScore/KATZresult$NumSeedGeneID*log10(1+KATZresult$NumSeedGeneID/KATZresult$NumScoreGeneID/2)*1000
    }else if (Methodnet=='RW'){
      RWresult=NetCorrResult$resultRWTable
      #resultRWTable:data.frame(SeedGeneIDPPI=toString(resultRW$SeedGeneID),NumSeedGeneID=resultRW$NumSeedGeneIDPPI,ScoreGeneIDPPI=toString(resultRW$ScoreGeneID),NumScoreGeneID=resultRW$NumScoreGeneIDPPI,TotalScore=resultRW$TotalScore,RandomScoreMean=resultRW$RandomScoreMedian,Pvalue=resultRW$Pvalue)
      ScoreGeneData$SeedGeneIDPPI=RWresult$SeedGeneIDPPI
      ScoreGeneData$NumSeedGeneID=RWresult$NumSeedGeneID
      ScoreGeneData$ScoreGeneIDPPI=RWresult$ScoreGeneIDPPI
      ScoreGeneData$NumScoreGeneID=RWresult$NumScoreGeneID
      ScoreGeneData$TotalScore=RWresult$TotalScore
      ##score/n0*exp(score/n1-1)*1000
      #ScoreGeneData$AdjTotalScore=RWresult$TotalScore/RWresult$NumSeedGeneID*exp(RWresult$TotalScore/RWresult$NumScoreGeneID-1)
      #score/n0*log10(1+n0/n1/2)
      ScoreGeneData$AdjTotalScore=RWresult$TotalScore/RWresult$NumSeedGeneID*log10(1+RWresult$NumSeedGeneID/RWresult$NumScoreGeneID/2)
    }
  }else{
    queryGeneIDdata=ChemData
    queryCID=unique(queryGeneIDdata$cid)
    ScoreGeneData=data.frame()
    for (i in 1:length(queryCID)){
      tempqueryGeneID=queryGeneIDdata$geneID[queryGeneIDdata$cid==queryCID[i]]
      tempNetCorrResult=NetCorr(SeedGeneID=diseaseGeneID,ScoreGeneID=tempqueryGeneID,Methodnet=Methodnet,randomtimes=1,gamma=0.7,testM=1)
      tempScoreGeneData=data.frame(cid=0,SeedGeneIDPPI=0,NumSeedGeneID=0,ScoreGeneIDPPI=0,NumScoreGeneID=0,overlapScore=0,Path1Score=0,Path2Score=0,Path3Score=0,TotalScore=0,AdjTotalScore=0)
      if (Methodnet=='KATZ'){
        tempKATZresult=tempNetCorrResult$resultKATZ
        tempScoreGeneData$cid=queryCID[i]
        #tempScoreGeneData$chemical_name=queryGeneIDdata$chemical_name[queryGeneIDdata$cid==queryCID[i]][1]
        #tempScoreGeneData$inchikey=queryGeneIDdata$inchikey[queryGeneIDdata$cid==queryCID[i]][1]
        tempScoreGeneData$SeedGeneIDPPI=tempKATZresult$SeedGeneIDPPI
        tempScoreGeneData$NumSeedGeneID=tempKATZresult$NumSeedGeneID
        tempScoreGeneData$ScoreGeneIDPPI=tempKATZresult$ScoreGeneIDPPI
        tempScoreGeneData$NumScoreGeneID=tempKATZresult$NumScoreGeneID
        tempScoreGeneData$overlapScore=tempKATZresult$overlapScore
        tempScoreGeneData$Path1Score=tempKATZresult$Path1Score
        tempScoreGeneData$Path2Score=tempKATZresult$Path2Score
        tempScoreGeneData$Path3Score=tempKATZresult$Path3Score
        tempScoreGeneData$TotalScore=tempKATZresult$TotalScore
        ##score/n0*exp(score/n1-1)*1000
        ###score/n0*log10(1+n0/n1/2)
        #tempScoreGeneData$AdjTotalScore=tempKATZresult$TotalScore/tempKATZresult$NumSeedGeneID*exp(tempKATZresult$TotalScore/tempKATZresult$NumScoreGeneID-1)*1000
        tempScoreGeneData$AdjTotalScore=tempKATZresult$TotalScore/tempKATZresult$NumSeedGeneID*log10(1+tempKATZresult$NumSeedGeneID/tempKATZresult$NumScoreGeneID/2)*1000
        ScoreGeneData=rbind(ScoreGeneData,tempScoreGeneData)
      }else if (Methodnet=='RW'){
        tempKATZresult=tempNetCorrResult$resultRWTable
        tempScoreGeneData$cid=queryCID[i]
        #tempScoreGeneData$chemical_name=queryGeneIDdata$chemical_name[queryGeneIDdata$cid==queryCID[i]][1]
        #tempScoreGeneData$inchikey=queryGeneIDdata$inchikey[queryGeneIDdata$cid==queryCID[i]][1]
        tempScoreGeneData$SeedGeneIDPPI=tempKATZresult$SeedGeneIDPPI
        tempScoreGeneData$NumSeedGeneID=tempKATZresult$NumSeedGeneID
        tempScoreGeneData$ScoreGeneIDPPI=tempKATZresult$ScoreGeneIDPPI
        tempScoreGeneData$NumScoreGeneID=tempKATZresult$NumScoreGeneID
        ##score/n0*exp(score/n1-1)*1000
        ###score/n0*log10(1+n0/n1/2)
        #tempScoreGeneData$AdjTotalScore=tempKATZresult$TotalScore/tempKATZresult$NumSeedGeneID*exp(tempKATZresult$TotalScore/tempKATZresult$NumScoreGeneID-1)*1000
        tempScoreGeneData$AdjTotalScore=tempKATZresult$TotalScore/tempKATZresult$NumSeedGeneID*log10(1+tempKATZresult$NumSeedGeneID/tempKATZresult$NumScoreGeneID/2)
        ScoreGeneData=rbind(ScoreGeneData,tempScoreGeneData)
      }
    }
  }
  return(ScoreGeneData)
}

combine_sim = function(g_fre_sub_cli,method=1,threshold=0.7){####g_fre_sub_cli是maxiclique得到的list对象;method：1一起合并，2逐个合并
  similarity=matrix(0,nrow=length(g_fre_sub_cli),ncol=length(g_fre_sub_cli))
  for (i in 1:length(g_fre_sub_cli)){
    for (j in 1:length(g_fre_sub_cli)){
      similarity[i,j]=length(intersect(g_fre_sub_cli[[i]],g_fre_sub_cli[[j]]))/length(union(g_fre_sub_cli[[i]],g_fre_sub_cli[[j]]))
      #similarity[i,j]=length(intersect(g_fre_sub_cli[[i]],g_fre_sub_cli[[j]]))/min(length(g_fre_sub_cli[[i]]),length(g_fre_sub_cli[[j]]))
    }
  }
  diag(similarity)=0
  g_fre_sub_cli_combine=g_fre_sub_cli
  maxsim=max(similarity)
  while (maxsim>=threshold){
    maxind=which(similarity==maxsim,arr.ind=T)
    if (method==1){
      maxind=sort(union(maxind[,1],maxind[,2]))
    }else{
      if (nrow(maxind)==1){
        maxind=sort(union(maxind[,1],maxind[,2]))
      }else{
        maxind=sort(union(maxind[1,1],maxind[1,2]))
      }
    }
    maxind=setdiff(maxind,nrow(similarity))
    g_fre_sub_cli_combine[maxind]=NULL
    g_fre_sub_cli_combine[[paste(maxsim)]]=sort(unique(unlist(g_fre_sub_cli[maxind])))


    similarity=matrix(0,nrow=length(g_fre_sub_cli_combine),ncol=length(g_fre_sub_cli_combine))
    for (i in 1:length(g_fre_sub_cli_combine)){
      for (j in 1:length(g_fre_sub_cli_combine)){
        similarity[i,j]=length(intersect(g_fre_sub_cli_combine[[i]],g_fre_sub_cli_combine[[j]]))/length(union(g_fre_sub_cli_combine[[i]],g_fre_sub_cli_combine[[j]]))
        #similarity[i,j]=length(intersect(g_fre_sub_cli_combine[[i]],g_fre_sub_cli_combine[[j]]))/min(length(g_fre_sub_cli_combine[[i]]),length(g_fre_sub_cli_combine[[j]]))
      }
    }
    diag(similarity)=0
    ##############合并结束输出g_fre_sub_cli_combine######################
    maxsim=max(similarity)
  }
  combine=lapply(g_fre_sub_cli_combine,sort)
  return(list(similarity=similarity,combine=combine))
}

combineDataFrame = function(D1,D2){
  name1=colnames(D1)
  name2=colnames(D2)
  name=unique(c(name1,name2))
  D1out=setdiff(name,name1)
  if (length(D1out)>0){
    D1outData=matrix(0,nrow=nrow(D1),ncol=length(D1out))
    D1outData=as.data.frame(D1outData)
    colnames(D1outData)=D1out
    D1C=cbind(D1,D1outData)
  }else{
    D1C=D1
  }
  D1C=D1C[,name]
  ######
  D2out=setdiff(name,name2)
  if (length(D1out)>0){
    D2outData=matrix(0,nrow=nrow(D2),ncol=length(D2out))
    D2outData=as.data.frame(D2outData)
    colnames(D2outData)=D2out
    D2C=cbind(D2,D2outData)
  }else{
    D2C=D2
  }
  D2C=D2C[,name]
  result=rbind(D1C,D2C)
}

combineIgraphOB = function(g1,g2){
  require(igraph)
  Node1_attr=vertex_attr_names(g1)
  Node2_attr=vertex_attr_names(g2)
  common_attr=intersect(Node1_attr,Node2_attr)
  common_attr=setdiff(common_attr,'name')
  Node1_attr_1=paste(common_attr,'1',sep='_')
  Node1_attr_2=paste(common_attr,'2',sep='_')
  ##########
  E1_attr=edge_attr_names(g1)
  E2_attr=edge_attr_names(g2)
  Ecommon_attr=intersect(E1_attr,E2_attr)
  E1_attr_1=paste(Ecommon_attr,'1',sep='_')
  E1_attr_2=paste(Ecommon_attr,'2',sep='_')
  #########
  g=g1+g2
  ##########

  if (length(common_attr)!=0){
    for (i in 1:length(common_attr)){
      tempText1=paste('V(g)$',Node1_attr_2[i],'[is.na(V(g)$',Node1_attr_2[i],')]','<-','V(g)$',Node1_attr_1[i],'[is.na(V(g)$',Node1_attr_2[i],')]',sep='')
      eval(parse(text=tempText1))
      order1=paste('V(g)$',common_attr[i],'=','V(g)$',Node1_attr_2[i],sep='')
      eval(parse(text=order1))
      ######################
      TTorder1=paste('g=delete_vertex_attr(g,','"',Node1_attr_1[i],'"',')',sep='')
      TTorder2=paste('g=delete_vertex_attr(g,','"',Node1_attr_2[i],'"',')',sep='')
      eval(parse(text=TTorder1))
      eval(parse(text=TTorder2))
    }
    ###
  }
  ########################
  if (length(Ecommon_attr)!=0){
    for (i in 1:length(Ecommon_attr)){
      tempText2=paste('E(g)$',E1_attr_2[i],'[is.na(E(g)$',E1_attr_2[i],')]','<-','E(g)$',E1_attr_1[i],'[is.na(E(g)$',E1_attr_2[i],')]',sep='')
      eval(parse(text=tempText2))
      ###
      order2=paste('E(g)$',Ecommon_attr[i],'=','E(g)$',E1_attr_2[i],sep='')
      eval(parse(text=order2))
      ###
      ETTorder1=paste('g=delete_edge_attr(g,','"',E1_attr_1[i],'"',')',sep='')
      ETTorder2=paste('g=delete_edge_attr(g,','"',E1_attr_2[i],'"',')',sep='')
      eval(parse(text=ETTorder1))
      eval(parse(text=ETTorder2))
    }
    #######
  }
  return(g)
}

combineLabDataFrame = function(labData1,labData2){
  Name1=unlist(labData1[1,])
  Name2=unlist(labData2[1,])
  NameAll=unique(c(Name1,Name2))
  #First6Names=Name2[1:6]
  #NameAllDel6name=setdiff(NameAll,First6Names)
  NameAllD=data.frame(ID=paste('V',1:length(NameAll),sep=''),Name=NameAll,stringsAsFactors=F)
  Name1Out=setdiff(NameAll,Name1)
  if (length(Name1Out)>0){
    Name1OutD=as.data.frame(matrix('*',nrow=nrow(labData1),ncol=length(Name1Out)),stringsAsFactors=F)
    Name1OutD[1,]=Name1Out
    Name1D=cbind(labData1,Name1OutD)
  }else{
    Name1D=labData1
  }
  Name1List=unlist(Name1D[1,])
  colName1=NameAllD$ID[match(Name1List,NameAllD$Name)]
  colnames(Name1D)=colName1
  ###
  Name2Out=setdiff(NameAll,Name2)
  if (length(Name2Out)>0){
    Name2OutD=as.data.frame(matrix('*',nrow=nrow(labData2),ncol=length(Name2Out)),stringsAsFactors=F)
    Name2OutD[1,]=Name2Out
    Name2D=cbind(labData2,Name2OutD)
  }else{
    Name2D=labData2
  }
  Name2List=unlist(Name2D[1,])
  colName2=NameAllD$ID[match(Name2List,NameAllD$Name)]
  colnames(Name2D)=colName2
  Name2D=Name2D[,colName1]
  Name2D=Name2D[-1,]
  ###
  dataAll=rbind(Name1D,Name2D)
  #dataAll=dataAll[,c(First6Names,NameAllDel6name)]
  return(dataAll)
}

coMemberMatrix = function(BF_result){
  require(data.table)
  require(igraph)
  BF_result2=BF_result$BF_result
  NBF=nrow(BF_result2)
  CoMatrix=matrix(0,nrow=NBF,ncol=NBF)
  colnames(CoMatrix)=BF_result2$id
  rownames(CoMatrix)=BF_result2$id
  for (i in 1:NBF){
    for (j in 1:NBF){
      tempi=BF_result2$id[i]
      tempj=BF_result2$id[j]
      makeupi=BF_result2$makeup[BF_result2$id==tempi]
      makeupi2=unlist(str_split(makeupi,','))
      makeupj=BF_result2$makeup[BF_result2$id==tempj]
      makeupj2=unlist(str_split(makeupj,','))
      CoMatrix[tempi,tempj]=length(intersect(makeupi2,makeupj2))
    }
  }
  CoMatrix=as.data.frame(CoMatrix)
  return(CoMatrix)
}

DisAccess = function(SeedGeneID,ScoreGeneID,randomtimes=1000,testM=1){
  require(data.table)
  require(plyr)
  require(stringr)
  require(igraph)
  #testM:1-->wilcox;testM:2-->perm
  #require("org.Hs.eg.db")
  require(pracma)
  if (!exists('ppiBinaryNet')){
   load('db/ppiNetData.db')
	 #load('db/ppiNetData.db')
  }
  ##
  if (sum(str_detect(SeedGeneID,',| |;'))>0){
    SeedGeneID=as.character(unlist(str_split(SeedGeneID,',| |;')))
    SeedGeneID=str_replace_all(SeedGeneID,'"','')
    SeedGeneID=str_replace_all(SeedGeneID,"'",'')
    SeedGeneID=SeedGeneID[SeedGeneID!='']
    SeedGeneID=unique(str_trim(SeedGeneID,side='both'))
  }
  SeedGeneID=as.character(unique(SeedGeneID))
  SeedGeneID=SeedGeneID[SeedGeneID%in%V(ppiBinaryNet)$name]
  ##
  if (sum(str_detect(ScoreGeneID,',| |;'))>0){
    ScoreGeneID=as.character(unlist(str_split(ScoreGeneID,',| |;')))
    ScoreGeneID=str_replace_all(ScoreGeneID,'"','')
    ScoreGeneID=str_replace_all(ScoreGeneID,"'",'')
    ScoreGeneID=ScoreGeneID[ScoreGeneID!='']
    ScoreGeneID=unique(str_trim(ScoreGeneID,side='both'))
  }
  ScoreGeneID=as.character(unique(ScoreGeneID))
  ScoreGeneID=ScoreGeneID[ScoreGeneID%in%V(ppiBinaryNet)$name]
  ##
  ppiDiameter=diameter(ppiBinaryNet)
  ##############
  Yan_inPPI_geneID=SeedGeneID
  BXXXT_inPPI_geneID=ScoreGeneID
  Dis_Yan=distances(ppiBinaryNet,BXXXT_inPPI_geneID,Yan_inPPI_geneID)
  dtimes=randomtimes
  #BXXXT_DRandomYan=matrix(0,nrow=dtimes,ncol=length(BXXXT_inPPI_geneID))
  ##
  BXXXT_DRandomYan=numeric()
  seedNum=12345678
  for (i in 1:dtimes){
    set.seed(seedNum+i)
    RandomNodes=sample(V(ppiBinaryNet)$name,length(BXXXT_inPPI_geneID))
    ##
    RdisY=distances(ppiBinaryNet,RandomNodes,Yan_inPPI_geneID)
    RdisY[RdisY==Inf]<-ppiDiameter
    #RdisY=rowMeans(RdisY)
    RdisY=mean(RdisY)
    #BXXXT_DRandomYan[i,]=RdisY
	BXXXT_DRandomYan[i]=RdisY
    ##
  }
  #########
  rm('ppiBinaryNet')
  gc()
  BXXXT_DRandomYan2=BXXXT_DRandomYan
  #BXXXT_DRandomYan2[which(BXXXT_DRandomYan2==Inf)]=NA
  #Dis_YanScore=rowMeans(Dis_Yan)
  Dis_YanScore=mean(Dis_Yan)
  #BXXXT_DRandomYan2Score=rowMeans(BXXXT_DRandomYan2)
  BXXXT_DRandomYan2Score=BXXXT_DRandomYan2
  if (testM==1){
  Pvalue=wilcox.test(Dis_YanScore,BXXXT_DRandomYan2Score)$p.value
  #Pvalue=sum(BXXXT_DRandomYan2Score<median(Dis_YanScore))/randomtimes
  }else{
  Pvalue=sum(BXXXT_DRandomYan2Score<median(Dis_YanScore))/randomtimes
  }
  
  return(list(SeedGeneID=SeedGeneID,NumSeedGeneID=length(SeedGeneID),ScoreGeneID=ScoreGeneID,NumScoreGeneID=length(ScoreGeneID),MedianDistance=paste(Dis_YanScore,' [',quantile(Dis_Yan,0.025),',',quantile(Dis_Yan,0.975),']',sep=''),RandomScoreMedian=paste(mean(BXXXT_DRandomYan2Score),' [',quantile(BXXXT_DRandomYan2Score,0.025),',',quantile(BXXXT_DRandomYan2Score,0.975),']',sep=''),Pvalue=Pvalue))
}

DotPlot = function(Data,cluster,item,colorBy,sizeby,showCategory=5,ascending=T,xlabSize=8,ylabSize=8,plotWith=20,plotHeight=18,filename=NULL){
  require(data.table)
  require(ggplot2)
  if (is.null(filename))filename='DotPlot'
  Data=as.data.table(Data)
  Data=Data[,c(cluster,item,colorBy,sizeby),with=F]
  if (ascending){
    DataSub=Data[order(get(cluster),get(colorBy))]
  }else{
    DataSub=Data[order(get(cluster),-get(colorBy))]
  }
  DataSub=DataSub[,head(.SD, showCategory),.SDcols=setdiff(colnames(DataSub),cluster),by=cluster]
  #DataSub[,c('sizebyVar'):=list((get(sizeby)-mean(get(sizeby)))/mean(get(sizeby))+1.5)]
  #p=ggplot(DataSub)+aes(x=get(cluster),y=get(item),colour=get(colorBy))+geom_point(aes(size=get(sizeby)))+labs(size=sizeby,colour=colorBy)+scale_size_continuous(range(1,4))+scale_colour_gradient(low="green", high="red")
  p=ggplot(DataSub)+aes(x=get(cluster),y=get(item),colour=get(colorBy))+geom_point(aes(size=get(sizeby)))+labs(size=sizeby,colour=colorBy)+scale_colour_gradient(low="red", high="green")
  p=p+xlab('')+ylab('')
  p=p+theme(axis.text.y=element_text(size=ylabSize,colour='black'),axis.text.x=element_text(size=xlabSize,colour='black'))
  print(p)
  if (!dir.exists('Report'))dir.create('Report')
  ggsave(paste('Report/',filename,'.png',sep=''),plot=p,width=plotWith,height=plotHeight,units='cm',dpi=600)
}

downLoadKEGGmap = function(geneID=NULL,keggid){
  require(pathview)
  require(stringr)
  #############
  pathWayIDtoPlot=unlist(str_split(keggid,',|;| '))
  pathWayIDtoPlot=str_replace_all(pathWayIDtoPlot,'"','')
  pathWayIDtoPlot=str_replace_all(pathWayIDtoPlot,"'",'')
  pathWayIDtoPlot=pathWayIDtoPlot[pathWayIDtoPlot!='']
  pathWayIDtoPlot=str_trim(pathWayIDtoPlot,side='both')
  ############
  if (str_detect(geneID,',|;| ')){
    geneID=unlist(str_split(geneID,',|;| '))
    geneID=str_replace_all(geneID,'"','')
    geneID=str_replace_all(geneID,"'",'')
    geneID=geneID[geneID!='']
    geneID=str_trim(geneID,side='both')
  }
  geneID=unique(geneID)
  ############
  NumPath=length(unique(pathWayIDtoPlot))
  pathview(gene.data=geneID,pathway.id=pathWayIDtoPlot,species="hsa",out.suffix='PA')
  for (i in 1:NumPath){
    tempSheetName=pathWayIDtoPlot[i]
    #sheetTemp=createSheet(wb, sheetName =tempSheetName)
    picfile=paste(tempSheetName,'.PA.png',sep='')
    #addPicture(picfile,sheetTemp,scale=5)
    file.copy(picfile,to='Report/',overwrite=T)
    file.remove(picfile)
    file.remove(paste(tempSheetName,'.png',sep=''))
    file.remove(paste(tempSheetName,'.xml',sep=''))
  }
}

EnrichMentAnalysis = function(geneClusterID,type,qvalueCutoff=0.05,filtKEGGdisease=F,ontGO='BP',ontDO='DOLite',ClusterNameToMESHenrich=names(geneClusterID)[1],meshDatabase='gene2pubmed',meshcategory='C'){
  #geneClusterID:list-->names(geneClusterID)
  #type:'KEGG' or 'GO' 'DO' or 'MESH' or 'ReactomePA'
  #ont:'BP','MF','CC'
  #MeshEnrich:
  ##ClusterNameToMESHenrich:names of Input List to analysis
  #            meshDatabase:'gendoo','gene2pubmed','RBBH'
  #            meshcategory:'C'-->disease,'D'-->Chemicals and Drugs,'E'-->Analytical, Diagnostic and Therapeutic Techniques and Equipment,'F'-->Psychiatry and Psychology,'G'-->Phenomena and Processes,'J'-->Technology and Food and Beverages,'N'-->Health Care
  ##download_KEGG {clusterProfiler}
  ##Mdotplot or Mbarplot
  require(ReactomePA)
  require(KEGG.db)
  require(fgsea)
  require(ggplot2)
  require(stringr)
  require(clusterProfiler)
  require(meshes)
  require(DOSE)
  require(data.table)
  require(plyr)
  require(igraph)
  require(corrplot)
  require("org.Hs.eg.db")
  #require(meshr)
  if (type=='KEGG'){
    #geneClusterID=list(BXXXT=BXXXT_geneID,Colitis=Yan_geneID,DM=Diabet_geneID,GC=Cancer_geneID)
    ckegg=compareCluster(geneCluster = geneClusterID, fun = "enrichKEGG", pAdjustMethod='fdr',qvalueCutoff=qvalueCutoff,,use_internal_data=T)
    aa=attributes(ckegg)
    temp_result=aa$compareClusterResult
    temp_result=as.data.table(temp_result)
    if (filtKEGGdisease){
      BXXXT_common_KeggID=unique(temp_result$ID)
      BXXXT_common_KeggIDfilt=BXXXT_common_KeggID[!BXXXT_common_KeggID%like%'hsa05']
      BXXXT_common_KeggIDfilt=BXXXT_common_KeggIDfilt[!BXXXT_common_KeggIDfilt%in%c('hsa04930','hsa04940','hsa04950','hsa04931','hsa04932','hsa04933')]
      temp_result=temp_result[ID%in%BXXXT_common_KeggIDfilt,]
    }
    temp_result=temp_result[p.adjust<qvalueCutoff,]
    ckegg2=new("compareClusterResult", compareClusterResult = temp_result,geneClusters = geneClusterID, fun = 'Enrichment')
    ###
      #png('KeggEnrich.png',width=1024,height=800)
      #Mdotplot(ckegg2,colorBy='qvalue',showCategory = showCategory,font.size=18)
      #dev.off(which = dev.cur())

    return(list(Enrichs=temp_result,OBject=ckegg2))
  }else if (type=='GO'){
    cgo1=compareCluster(geneCluster = geneClusterID, OrgDb=org.Hs.eg.db,fun = "enrichGO",ont = ontGO,pAdjustMethod="fdr",qvalueCutoff=qvalueCutoff,readable=T)
    ####
      #png(paste('GO_Enrich_',ontGO,'.png',sep=''),width=1368,height=800)
      #Mdotplot(cgo1,colorBy='qvalue',showCategory = showCategory,font.size=18)
      #dev.off()
    ####
    return(list(Enrichs=cgo1@compareClusterResult,OBject=cgo1))
  }else if (type=='MESH'){
    BXXXT_geneID=geneClusterID[[ClusterNameToMESHenrich]]
    meshEnrich=meshes::enrichMeSH(BXXXT_geneID, MeSHDb = "MeSH.Hsa.eg.db", database=meshDatabase, category = meshcategory,pAdjustMethod='fdr',qvalueCutoff=qvalueCutoff)
    meshEnrichs=as.data.frame(meshEnrich)
    # meshParams=new("MeSHHyperGParams",geneIds=BXXXT_geneID,universeGeneIds=universeGeneIDs,annotation="MeSH.Hsa.eg.db",category="C",database="gendoo",pvalueCutoff=pvalueCutoff,pAdjust="lFDR")
    # BXXXT_mesh_gendoo=meshHyperGTest(meshParams)
    # BXXXT_mesh_gendoos=summary(BXXXT_mesh_gendoo)
    # BXXXT_mesh_gendoos=as.data.table(BXXXT_mesh_gendoos)
    # database(meshParams)<-"gene2pubmed"
    # BXXXT_mesh_pubMed=meshHyperGTest(meshParams)
    # BXXXT_mesh_pubMeds=summary(BXXXT_mesh_pubMed)
    # BXXXT_mesh_pubMeds=as.data.table(BXXXT_mesh_pubMeds)
    # BXXXT_meshTotal=rbind(BXXXT_mesh_gendoos,BXXXT_mesh_pubMeds)
    # BXXXT_meshTotal=unique(BXXXT_meshTotal,by='MESHID')
    # BXXXT_meshTotal_Includ=BXXXT_meshTotal[Pvalue<pvalueCutoff,.(ID=MESHID,Description=MESHTERM,GeneRatio=OddsRatio,pvalue=Pvalue,p.adjust=lFDR,qvalue=lFDR,geneID=GENEID,Count=Count)]
    # BXXXT_meshTotal_Includ=BXXXT_meshTotal_Includ[order(-Count)]
    meshObject<-new("enrichResult", result = meshEnrichs, qvalueCutoff = qvalueCutoff,
                    pAdjustMethod = 'fdr', organism = 'Homo sapiens',
                    ontology = 'MeSH', gene = BXXXT_geneID,
                    universe = 'extID',
                    readable = FALSE)
    ####
      #png('Mesh_plot.png',width=1024,height=800)
      #Mbarplot(meshObject,drop=TRUE, colorBy='qvalue',showCategory=showCategory,font.size=18)
      #dev.off(which = dev.cur())
   #####
    return(list(Enrichs=meshEnrichs,OBject=meshObject))
  }else if (type=='DO') {
    data(EG2DOLite)
    data(DOLiteTerm)
    data(DOLite2EG)
    cdo=compareCluster(geneCluster = geneClusterID, fun = "enrichDO",ont=ontDO,pAdjustMethod="fdr",qvalueCutoff=qvalueCutoff)
    ####
    #png('DO_Enrich.png',width=1368,height=800)
    #Mdotplot(cdo,colorBy='qvalue',showCategory = showCategory,font.size=18)
    #dev.off()
    ####
    return(list(Enrichs=cdo@compareClusterResult,OBject=cdo))
  }else if (type=='ReactomePA'){
    ReactomePA=compareCluster(geneCluster = geneClusterID, fun = "enrichPathway",pAdjustMethod="fdr",qvalueCutoff=qvalueCutoff,readable=T)
    ####
    #png(paste('GO_Enrich_',ontGO,'.png',sep=''),width=1368,height=800)
    #Mdotplot(cgo1,colorBy='qvalue',showCategory = showCategory,font.size=18)
    #dev.off()
    ####
    return(list(Enrichs=ReactomePA@compareClusterResult,OBject=ReactomePA))

  }else{
    print('ERROR:type not found!')
  }
}

ExpandNet = function(geneID,method=1,cutoffValue=2){
  #method:扩展方法：0-->核心网不进行扩展；1-->最短路径法；2-->最近邻法
  #cutoffValue:method1时的阈值--BXXXT_Dis[BXXXT_Dis>cutoffValue]<-0
  #ppiBinaryNet
  library(data.table)
  library(plyr)
  library(stringr)
  library(igraph)
  library("org.Hs.eg.db")
  if (!exists('ppiBinaryNet')){
    load('db/ppiNetData.db')
  }
  geneID=as.character(geneID)
  geneID=geneID[geneID%in%V(ppiBinaryNet)$name]
  if (method==1){
    BXXXT_Dis=distances(ppiBinaryNet,geneID,geneID)
    BXXXT_Dis[BXXXT_Dis>cutoffValue]<-0
    Index_dis=which(BXXXT_Dis!=0,arr.ind=T)
    kpath_BXXXT_Dis=list()
    length(kpath_BXXXT_Dis)=nrow(Index_dis)
    for (i in 1:nrow(Index_dis)){
      templocation1=Index_dis[i,1]
      templocation2=Index_dis[i,2]
      tempPath=shortest_paths(ppiBinaryNet,rownames(BXXXT_Dis)[templocation1],colnames(BXXXT_Dis)[templocation2],output='vpath')
      kpath_BXXXT_Dis[[i]]=as_ids(tempPath$vpath[[1]])
    }
    BXXXT_kpath_nodes=unique(unlist(kpath_BXXXT_Dis))
    BXXXT_kpath_Net=induced_subgraph(ppiBinaryNet,vids=BXXXT_kpath_nodes)
    BXXXT_kpath_Net=igraph::simplify(BXXXT_kpath_Net)
    return(BXXXT_kpath_Net)
  }else if (method==2){
    BXXXT_Neighbor=neighborhood(ppiBinaryNet,1,nodes=geneID[1])
    for (i in 2:length(geneID)){
      tempID=geneID[i]
      temp=neighborhood(ppiBinaryNet,1,nodes=tempID)
      if (i==2){
        temp2=c(BXXXT_Neighbor[[1]],temp[[1]])
      }else{
        temp2=c(BXXXT_Neighbor,temp[[1]])
      }
      BXXXT_Neighbor=unique(temp2)
    }
    ####
    BXXXT_Net=induced_subgraph(ppiBinaryNet,vids=BXXXT_Neighbor)
    return(BXXXT_Net)
  }else{
    BXXXT_CoreNet=induced_subgraph(ppiBinaryNet,vids=geneID)
    return(BXXXT_CoreNet)
  }
}

ExportChemTarNet = function(chemTar,geneListName='geneList',diseaseGeneID=NULL,diseaseName='Disease',expandMethod=c('Core','EP','EN'),cutoffValue=2,IncludeHerb=T,IncludeChem=T,IncludePPI=T,Savefilename='EpChemTarNet'){
  ##chemTar:colnames:herb,cid,chemical_name,geneID
  library(data.table)
  library(plyr)
  library(igraph)
  if (!dir.exists('Report')) dir.create('Report')
  chemTar=as.data.table(chemTar)
  if (!'geneName'%in%colnames(chemTar))chemTar=transform(chemTar,geneName=geneTrans(geneID,type='ID'))
  if (!'chemical_name'%in%colnames(chemTar))chemTar=transform(chemTar,chemical_name=cid)
  if (!'herb'%in%colnames(chemTar))chemTar=transform(chemTar,herb='drug')
  chemTar=chemTar[!is.na(geneName),]
  chemTar=unique(chemTar,by=c('herb','cid','geneID'))
  #################chemHerbNet
  chemHerb=chemTar[,.N,by=c('herb','cid')]
  chemHerb=rename(chemHerb,c('N'='weight'))
  chemHerbNodeName1=unique(data.table(name=unique(chemHerb$herb),symbol=unique(chemHerb$herb),class='herb'))
  chemHerbNodeName2=unique(chemTar[,.(name=cid,symbol=chemical_name,class='chemical')],by='name')
  chemHerbNodeName=rbind(chemHerbNodeName1,chemHerbNodeName2)
  chemHerbNet=graph_from_data_frame(chemHerb,directed = F,vertices =chemHerbNodeName )
  ###############chemGeneNet
  chemGene=chemTar[,.N,by=c('cid','geneID')]
  chemGene=rename(chemGene,c('N'='weight'))
  chemGeneNodeName1=unique(chemTar[,.(name=cid,symbol=chemical_name,class='chemical')],by='name')
  chemGeneNodeName2=unique(chemTar[,.(name=geneID,symbol=geneName,class='gene')],by='name')
  chemGeneNodeName=rbind(chemGeneNodeName1,chemGeneNodeName2)
  chemGeneNet=graph_from_data_frame(chemGene,directed = F,vertices =chemGeneNodeName )
  ###############herbGeneNet###
  herbGene=chemTar[,.N,by=c('herb','geneID')]
  herbGene=rename(herbGene,c('N'='weight'))
  herbGeneNodeName1=unique(data.table(name=unique(herbGene$herb),symbol=unique(herbGene$herb),class='herb'))
  herbGeneNodeName2=unique(chemTar[,.(name=geneID,symbol=geneName,class='gene')],by='name')
  herbGeneNodeName=rbind(herbGeneNodeName1,herbGeneNodeName2)
  herbGeneNet=graph_from_data_frame(herbGene,directed = F,vertices =herbGeneNodeName )
  #############PPI
  #if ('Core'%in%expandMethod)method=0
  #if ('EP'%in%expandMethod)method=1
  #if ('EN'%in%expandMethod)method=2
  if (is.null(Savefilename))Savefilename='EpChemTarNet'
  ExNet=ExportNetwithin(geneID=unique(chemTar$geneID),geneIDName=geneListName,diseaseGeneID=diseaseGeneID,diseaseName=diseaseName,expandMethod=expandMethod,cutoffValue=cutoffValue)
  E(ExNet)$weight<-1
  if (IncludeHerb&IncludeChem&IncludePPI){
    NetALL=combineIgraphOB(chemHerbNet,chemGeneNet)
    NetALL=combineIgraphOB(NetALL,ExNet)
  }
  if (!IncludeHerb&IncludeChem&IncludePPI)NetALL=combineIgraphOB(chemGeneNet,ExNet)
  if (IncludeHerb&!IncludeChem&IncludePPI)NetALL=combineIgraphOB(herbGeneNet,ExNet)
  if (IncludeHerb&IncludeChem&!IncludePPI)NetALL=combineIgraphOB(chemHerbNet,chemGeneNet)
  if (!IncludeHerb&!IncludeChem&IncludePPI)NetALL=ExNet
  if (!IncludeHerb&IncludeChem&!IncludePPI)NetALL=chemGeneNet
  if (IncludeHerb&!IncludeChem&!IncludePPI)NetALL=herbGeneNet
  if (!IncludeHerb&!IncludeChem&!IncludePPI)NetALL=make_empty_graph()
  file=paste('Report/',Savefilename,'.gml',sep='')
  write_graph(NetALL, file=file, format='gml')
  return(NetALL)
}

ExportNet = function(geneID,geneIDName='geneList',diseaseGeneID,diseaseName='Disease',expandMethod=c('Core','EP','EN'),cutoffValue=2,Savefilename='EpNet'){
  #return(list(overlapGene=EXoverlapNodes,overlapGeneCount=length(EXoverlapNodes),ExGene=V(ExNet)$name,Exdisease=V(ExDiseaseNet)$name,OverlapResult=OverlapResult,NetStat=NetStat))
  library(data.table)
  library(plyr)
  library(stringr)
  library(igraph)
  #require(NetPathMiner)
  if (sum(str_detect(geneID,',|;'))){
    geneID=as.character(unlist(str_split(geneID,',|;')))
    geneID=str_replace_all(geneID,'"','')
    geneID=str_replace_all(geneID,"'",'')
    geneID=geneID[geneID!='']
    geneID=unique(str_trim(geneID,side='both'))
  }
  #######
  if (is.null(geneIDName))geneIDName='geneList'
  if (is.null(diseaseName))diseaseName='Disease'
  if ('Core'%in%expandMethod)method=0
  if ('EP'%in%expandMethod)method=1
  if ('EN'%in%expandMethod)method=2
  ####
  if (!is.null(diseaseGeneID)){
  if (sum(str_detect(diseaseGeneID,',|;'))){
    diseaseGeneID=as.character(unlist(str_split(diseaseGeneID,',|;')))
    diseaseGeneID=str_replace_all(diseaseGeneID,'"','')
    diseaseGeneID=str_replace_all(diseaseGeneID,"'",'')
    diseaseGeneID=diseaseGeneID[diseaseGeneID!='']
    diseaseGeneID=unique(str_trim(diseaseGeneID,side='both'))
  }
  ####
  CoreOverlap=intersect(geneID,diseaseGeneID)
  ExNet=ExpandNet(geneID=geneID,method=method,cutoffValue=cutoffValue)
  ExDiseaseNet=ExpandNet(geneID=diseaseGeneID,method=method,cutoffValue=cutoffValue)
  NetAll=ExNet+ExDiseaseNet
  overlapNodes=intersect(V(ExNet)$name,V(ExDiseaseNet)$name)
  EXoverlapNodes=setdiff(intersect(V(ExNet)$name,V(ExDiseaseNet)$name),CoreOverlap)
  NonlapDiseaseNodes=setdiff(V(ExDiseaseNet)$name,c(EXoverlapNodes,diseaseGeneID))
  NongeneIDNods=setdiff(V(ExNet)$name,c(EXoverlapNodes,geneID))
  netGeneID=V(NetAll)$name
  netGeneName=geneTrans(netGeneID,type='ID')
  V(NetAll)$symbol=netGeneName
  V(NetAll)$type<-'NA'
  V(NetAll)$type[V(NetAll)$name%in%geneID]<-geneIDName
  V(NetAll)$type[V(NetAll)$name%in%diseaseGeneID]<-diseaseName
  V(NetAll)$type[V(NetAll)$name%in%CoreOverlap]<-'CoreOverlap'
  V(NetAll)$type[V(NetAll)$name%in%EXoverlapNodes]<-'ExpandOverlap'
  V(NetAll)$type[V(NetAll)$name%in%NonlapDiseaseNodes]<-paste('Expand',diseaseName,sep='_')
  V(NetAll)$type[V(NetAll)$name%in%NongeneIDNods]<-paste('Expand',geneIDName,sep='_')
  V(NetAll)$degree=igraph::degree(NetAll)
  if (is.null(Savefilename))Savefilename='EpNet'
  file=paste('Report/',Savefilename,'.gml',sep='')
  write_graph(NetAll, file=file, format='gml')
  #plotCytoscapeGML(NetAll, file=file)
  order=paste('NetList=list(',geneIDName,'=ExNet',',',diseaseName,'=ExDiseaseNet',')',sep='')
  eval(parse(text=order))
  OverlapResult=as.data.frame(NodeOverlapAnalysis(NetList=NetList,Plot=TRUE))
  NetStat=NetStatCal(NetList)
  NetStat=as.data.frame(NetStat,stringsAsFactors=F)
  NetStat=cbind(Item=rownames(NetStat),NetStat)
  return(list(overlapGene=overlapNodes,overlapGeneCount=length(overlapNodes),ExGene=V(ExNet)$name,Exdisease=V(ExDiseaseNet)$name,OverlapResult=OverlapResult,NetStat=NetStat))
  }else{
  ###################
    ExNet=ExpandNet(geneID=geneID,method=method,cutoffValue=cutoffValue)
    NetAll=ExNet
    netGeneID=V(NetAll)$name
    netGeneName=geneTrans(netGeneID,type='ID')
    V(NetAll)$symbol=netGeneName
    V(NetAll)$degree=igraph::degree(NetAll)
    if (is.null(Savefilename))Savefilename='EpNet'
    file=paste('Report/',Savefilename,'.gml',sep='')
    write_graph(NetAll, file=file, format='gml')
    #plotCytoscapeGML(NetAll, file=file)
    order=paste('NetList=list(',geneIDName,'=ExNet',')',sep='')
    eval(parse(text=order))
    #OverlapResult=as.data.frame(NodeOverlapAnalysis(NetList=NetList,Plot=TRUE))
    NetStat=NetStatCal(NetList)
    NetStat=as.data.frame(NetStat,stringsAsFactors=F)
    NetStat=cbind(Item=rownames(NetStat),NetStat)
##############################
    return(list(overlapGene=NULL,overlapGeneCount=0,ExGene=V(ExNet)$name,Exdisease=NULL,OverlapResult=data.frame(),NetStat=NetStat))
    }
}

ExportNetwithin = function(geneID,geneIDName='geneList',diseaseGeneID,diseaseName='Disease',expandMethod=c('Core','EP','EN'),cutoffValue=2){
  #return(list(overlapGene=EXoverlapNodes,overlapGeneCount=length(EXoverlapNodes),ExGene=V(ExNet)$name,Exdisease=V(ExDiseaseNet)$name,OverlapResult=OverlapResult,NetStat=NetStat))
  library(data.table)
  library(plyr)
  library(stringr)
  library(igraph)
  #require(NetPathMiner)
  if (sum(str_detect(geneID,',|;'))){
    geneID=as.character(unlist(str_split(geneID,',|;')))
    geneID=str_replace_all(geneID,'"','')
    geneID=str_replace_all(geneID,"'",'')
    geneID=geneID[geneID!='']
    geneID=unique(str_trim(geneID,side='both'))
  }
  #######
  if (is.null(geneIDName))geneIDName='geneList'
  if (is.null(diseaseName))diseaseName='Disease'
  if ('Core'%in%expandMethod)method=0
  if ('EP'%in%expandMethod)method=1
  if ('EN'%in%expandMethod)method=2
  ####
  if (!is.null(diseaseGeneID)){
    if (sum(str_detect(diseaseGeneID,',|;'))){
      diseaseGeneID=as.character(unlist(str_split(diseaseGeneID,',|;')))
      diseaseGeneID=str_replace_all(diseaseGeneID,'"','')
      diseaseGeneID=str_replace_all(diseaseGeneID,"'",'')
      diseaseGeneID=diseaseGeneID[diseaseGeneID!='']
      diseaseGeneID=unique(str_trim(diseaseGeneID,side='both'))
    }
    ####
    CoreOverlap=intersect(geneID,diseaseGeneID)
    ExNet=ExpandNet(geneID=geneID,method=method,cutoffValue=cutoffValue)
    ExDiseaseNet=ExpandNet(geneID=diseaseGeneID,method=method,cutoffValue=cutoffValue)
    NetAll=ExNet+ExDiseaseNet
    overlapNodes=intersect(V(ExNet)$name,V(ExDiseaseNet)$name)
    EXoverlapNodes=setdiff(intersect(V(ExNet)$name,V(ExDiseaseNet)$name),CoreOverlap)
    NonlapDiseaseNodes=setdiff(V(ExDiseaseNet)$name,c(EXoverlapNodes,diseaseGeneID))
    NongeneIDNods=setdiff(V(ExNet)$name,c(EXoverlapNodes,geneID))
    netGeneID=V(NetAll)$name
    netGeneName=geneTrans(netGeneID,type='ID')
    V(NetAll)$symbol=netGeneName
    V(NetAll)$type<-'NA'
    V(NetAll)$type[V(NetAll)$name%in%geneID]<-geneIDName
    V(NetAll)$type[V(NetAll)$name%in%diseaseGeneID]<-diseaseName
    V(NetAll)$type[V(NetAll)$name%in%CoreOverlap]<-'CoreOverlap'
    V(NetAll)$type[V(NetAll)$name%in%EXoverlapNodes]<-'ExpandOverlap'
    V(NetAll)$type[V(NetAll)$name%in%NonlapDiseaseNodes]<-paste('Expand',diseaseName,sep='_')
    V(NetAll)$type[V(NetAll)$name%in%NongeneIDNods]<-paste('Expand',geneIDName,sep='_')
    #V(NetAll)$degree=igraph::degree(NetAll)
    return(NetAll)
  }else{
    ###################
    ExNet=ExpandNet(geneID=geneID,method=method,cutoffValue=cutoffValue)
    NetAll=ExNet
    netGeneID=V(NetAll)$name
    netGeneName=geneTrans(netGeneID,type='ID')
    V(NetAll)$symbol=netGeneName
    #V(NetAll)$degree=igraph::degree(NetAll)
    return(NetAll)
  }
}

findBATMANtarget = function(pubMedid){
  require(data.table)
  require(stringr)
  #BATMANdatabase:BATMANChemTarget--->colnames:c('Pubchem_CID','geneID')
  if (!exists('BATMANChemTarget')){
    load('db/BATMANdatabase.RData')
  }
  data=BATMANChemTarget[Pubchem_CID%in%pubMedid,.(Pubchem_CID,geneID)]
  data$InputName=pubMedid
  return(data)
}

findCTDtarget = function(ChemName){
  require(data.table)
  require(stringr)
  #CTDdatabase:CTDchemTarget--->colnames:c('ChemicalName','CTDid','geneName','geneID')
  if (!exists('CTDchemTarget')){
    load('db/CTDdatabase.RData')
  }
  data=CTDchemTarget[ChemicalName%like%ChemName,.(ChemicalName,CTDid,geneName,geneID,score)]
  data=unique(data)
  data$InputName=ChemName
  return(data)
}

findHITChemFromHerb = function(HerbName,type='CH'){
  #HerbName
  #type:type of HerbName-->chinese(CH) or PingYin(PY)
  require(data.table)
  require(stringr)
  if (!exists('HitHerbChem')){
    load('db/HITdatabase.RData')
  }
  if (type=='CH'){
    data=HitHerbChem[herbCH%like%HerbName,.(name,herbCH)]
  }else if (type=='PY'){
    data=HitHerbChem[herbCH==HerbName,.(name,herbCH)]
  }
  return(data)
}

findHITTarget = function(Name,NameType){
  #Name
  #NameType:chemicalName(CN) or herbChineseName(HCN) or herbPingYinName(HPN) or Chem_pubmedID(PID)
  #HITdatabase:HitChemTarget and HitHerbChem
  require(data.table)
  require(stringr)
  require(plyr)
  if (!exists('HitHerbChem')){
    load('db/HITdatabase.RData')
  }
  if (NameType=='CN'){
    CNidD=HitHerbChem[name%like%Name,.(HITid=cid,name)]
    CNidD=unique(CNidD)
    CNid=unique(CNidD$HITid)
    CT=HitChemTarget[HITid%in%CNid,.(HITid,geneName,geneID)]
    if (nrow(CT)!=0){
      CT=join(CT,CNidD,by='HITid')
      CT$InputName=Name
      }
    return(CT)
  }else if (NameType=='HCN'){
    HCNid=HitHerbChem[herbCH%like%Name,.(cid,name,herbCH)]
    if (nrow(HCNid)!=0){
      HCNCid=unique(HCNid$cid)
      CT=HitChemTarget[HITid%in%HCNCid,.(cid=HITid,geneName,geneID)]
      CT=as.data.frame(CT)
      HCNCid=as.data.frame(HCNCid)
      CT2=join(CT,HCNid,by='cid')
      CT2$InputName=Name
    }else{
      CT2=data.table()
    }
    return(CT2)
  }else if (NameType=='HPN'){
    HCNid=HitHerbChem[herbPY==Name,.(cid,name,herbPY)]
    if (nrow(HCNid)!=0){
      HCNCid=unique(HCNid$cid)
      CT=HitChemTarget[HITid%in%HCNCid,.(cid=HITid,geneName,geneID)]
      CT=as.data.frame(CT)
      HCNCid=as.data.frame(HCNCid)
      CT2=join(CT,HCNid,by='cid')
      CT2$InputName=Name
    }else{
      CT2=data.table()
    }
    return(CT2)
  }else if (NameType=='PID'){
    Name2=paste('CID:',Name,sep='')
    temp=HitHerbChem[,.(HITid=cid,name)]
    CT=HitChemTarget[pub_id%in%Name2,.(HITid,PID=pub_id,geneName,geneID)]
    if (nrow(CT)!=0){CT=join(CT,temp,by='HITid');CT$InputName=Name}
    return(CT)
  }
}

findMinetarget = function(Name,NameType){
  #Name
  #NameType:chemicalName(CN) or herbChineseName(HN)
  #MineCTD:'herb','chemical_name','geneID','smi'
  require(data.table)
  require(stringr)
  require(plyr)
  if (!exists('MineCTD')){
    MineCTD=fread('db/Mine_herb_chem_gene_smi.csv',colClasses = 'character')
  }
 if (NameType=='CN'){
   data=MineCTD[chemical_name%like%Name,.(chemical_name,geneID,smiles)]
   data$InputName=Name
   return(data)
 }else if (NameType=='HN'){
   data=MineCTD[herb%like%Name,.(herb,chemical_name,geneID,smiles)]
   data$InputName=Name
   return(data)
 }else{
   print('ERROR:NameType not found!')
 }
}

findStitchTarget = function(chem_inchikey,scoreSet=0,STITCH_CTdb='db/STITCH5_CTdb.db',STITCH_CKEYdb='db/STITCH5_CKEYdb.db'){
  #chem_inchikey:inchikey
  #STITCH_CTdb:'STITCH5_CTdb.db',colnames=c('cid','ENSP','score','geneID','inchikey')
  #STITCH_CKEYdb:'STITCH_CKEYdb.db',colnames=c('inchikey','source_cid','stereo_chemical_id','flat_chemical_id')
  require(RSQLite)
  require(stringr)
  require(data.table)
  require(foreign)
  tmp <- dbConnect(SQLite(), STITCH_CKEYdb)
  Order0=paste('inchikey IN ',"('",chem_inchikey,"')",sep='')
  Order1=paste('select source_cid from STITCH_CKEYdb where',Order0,sep=' ')
  res1 <- dbSendQuery(tmp, Order1)
  data1=fetch(res1)
  dbDisconnect(tmp)
  if (nrow(data1)!=0){
    fethCID=data1$source_cid[1]
    Order00=paste('cid IN ',"('",fethCID,"')",sep='')
    Order2=paste('select cid,ENSP,geneID,score from STITCHdb where',Order00,sep=' ')
    tmp2 <- dbConnect(SQLite(), STITCH_CTdb)
    res2 <- dbSendQuery(tmp2, Order2)
    data2 <- fetch(res2)
    dbDisconnect(tmp2)
    data2=as.data.table(data2)
    data2=unique(data2)
    if (nrow(data2)!=0){
      data2=data2[score>=scoreSet,]
    }
	result=data2
    result$inchikey=chem_inchikey    
  }else{
    result=data.table()   
  }  
  return(result)
}

findStitchTarget2 = function(chem_inchikey,scoreSet=0,STITCH_CTdb='db/STITCH5_CTdb.db'){
  #chem_inchikey:inchikey
  #STITCH_CTdb:'STITCH5_CTdb.db',colnames=c('cid','ENSP','score','geneID','inchikey')
  #STITCH_CKEYdb:'STITCH_CKEYdb.db',colnames=c('inchikey','source_cid','stereo_chemical_id','flat_chemical_id')
  require(RSQLite)
  require(stringr)
  require(data.table)
  require(foreign)
  tmp <- dbConnect(SQLite(), STITCH_CTdb)
  Order0=paste('inchikey IN ',"('",chem_inchikey,"')",sep='')
  Order1=paste('select * from STITCHdb where',Order0,sep=' ')
  res1 <- dbSendQuery(tmp, Order1)
  data1=fetch(res1)
  dbDisconnect(tmp)
  if (nrow(data1)!=0){
    data2=as.data.table(data1)
    data2=unique(data2)
    if (nrow(data2)!=0){
      data2=data2[score>=scoreSet,]
    }
	result=data2
    result$inchikey=chem_inchikey    
  }else{
    result=data.table()   
  }  
  return(result)
}

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

findTCMSPchem = function(ChemName){
  #ChemName:ChemName-->Enlish
  #TCMSPdatabase:TCMSPherbChem and TCMSPChemTarget
  require(data.table)
  require(stringr)
  if (!exists('TCMSPherbChem')){
    load('db/TCMSPdatabase.RData')
  }
 data=TCMSPherbChem[MOL_name%like%ChemName,.(MOL_ID,MOL_name,inchikey,smiles)]
 data=unique(data)
 data$InputName=ChemName
 return(data)
}

findTCMSPChemFromHerb = function(HerbName){
    #HerbName:HerbName-->chinese
    #TCMSPdatabase:TCMSPherbChem and TCMSPChemTarget
    require(data.table)
    require(stringr)
    if (!exists('TCMSPherbChem')){
      load('db/TCMSPdatabase.RData')
    }
    data=TCMSPherbChem[herb_name%like%HerbName,.(herb_name,MOL_name,inchikey,smiles)]
    data=unique(data)
    data$InputName=HerbName
    return(data)
}

findTCMSPtarget = function(Name,NameType){
  #Name
  #NameType:chemicalName(CN) or herbChineseName(HN)
  #TCMSPdatabase:TCMSPherbChem and TCMSPChemTarget
  require(data.table)
  require(stringr)
  require(plyr)
  if (!exists('TCMSPherbChem')){
    load('db/TCMSPdatabase.RData')
  }
  if (NameType=='CN'){
    MOLidD=TCMSPherbChem[MOL_name%like%Name,.(MOL_name,MOL_ID)]
    if (nrow(MOLidD)!=0){
      MOLidD=unique(MOLidD)
      MOLid=unique(MOLid$MOL_ID)
      CT=TCMSPChemTarget[MOL_ID%in%MOLid,.(MOL_ID,geneName,geneID)]
      if (nrow(CT)!=0){
        CT2=join(CT,MOLidD,by='MOL_ID')
        CT2$InputName=Name
      }else{
        CT2=data.table()
      }
    }else{
      CT2=data.table()
    }
    return(CT2)
  }else if (NameType=='HN') {
    TargetD=TCMSPChemTarget[herb_name%like%Name,.(herb_name,MOL_ID,geneName,geneID)]
    temp=TCMSPherbChem[,.(MOL_ID,MOL_name)]
    if (nrow(TargetD)!=0){
      CT=join(TargetD,temp,by='MOL_ID')
      CT$InputName=Name
    }else{
      CT=data.table()
    }
    return(CT)
  }else{
    print('ERROR:NameType not found!')
  }
}

fortify__compareClusterResult = function(model, data, showCategory=5, by="geneRatio", includeAll=TRUE){
  clProf.df <- as.data.frame(model)

  ## get top 5 (default) categories of each gene cluster.
  if (is.null(showCategory)) {
    result <- clProf.df
  } else {
    Cluster <- NULL # to satisfy codetools
    result <- ddply(.data = clProf.df,
                    .variables = .(Cluster),
                    .fun = function(df, N) {
                      if (length(df$Count) > N) {
                        if (any(colnames(df) == "pvalue")) {
                          idx <- order(df$pvalue, decreasing=FALSE)[1:N]
                        } else {
                          ## for groupGO
                          idx <- order(df$Count, decreasing=T)[1:N]
                        }
                        return(df[idx,])
                      } else {
                        return(df)
                      }
                    },
                    N=showCategory
    )

  }
  ID <- NULL
  if (includeAll == TRUE) {
    result = subset(clProf.df, ID %in% result$ID)
  }

  ## remove zero count
  result$Description <- as.character(result$Description) ## un-factor
  GOlevel <- result[,c("ID", "Description")] ## GO ID and Term
  GOlevel <- unique(GOlevel)

  result <- result[result$Count != 0, ]
  result$Description <- factor(result$Description,
                               levels=rev(GOlevel[,2]))


  if (by=="rowPercentage") {
    Description <- Count <- NULL # to satisfy codetools
    result <- ddply(result,
                    .(Description),
                    transform,
                    Percentage = Count/sum(Count),
                    Total = sum(Count))

    ## label GO Description with gene counts.
    x <- mdply(result[, c("Description", "Total")], paste, sep=" (")
    y <- sapply(x[,3], paste, ")", sep="")
    result$Description <- y

    ## restore the original order of GO Description
    xx <- result[,c(2,3)]
    xx <- unique(xx)
    rownames(xx) <- xx[,1]
    Termlevel <- xx[as.character(GOlevel[,1]),2]

    ##drop the *Total* column
    result <- result[, colnames(result) != "Total"]

    result$Description <- factor(result$Description,
                                 levels=rev(Termlevel))

  } else if (by == "count") {
    ## nothing
  } else if (by == "geneRatio") {
    gsize <- as.numeric(sub("/\\d+$", "", as.character(result$GeneRatio)))
    gcsize <- as.numeric(sub("^\\d+/", "", as.character(result$GeneRatio)))
    result$GeneRatio = gsize/gcsize
    result$Cluster <- paste(as.character(result$Cluster),"\n", "(", gcsize, ")", sep="")
  } else {
    ## nothing
  }
  return(result)
}

fortify_enrichResult = function(model, data, showCategory=5, by = "Count", order=FALSE, drop=FALSE, ...){
  res <- as.data.frame(model)
  if (drop) {
    res <- res[res$Count != 0, ]
  }
  ####
  gsize <- as.numeric(sub("/\\d+$", "", as.character(res$GeneRatio)))
  gcsize <- as.numeric(sub("^\\d+/", "", as.character(res$GeneRatio)))
  res$GeneRatio = gsize/gcsize
  ####
  #res$GeneRatio <- parse_ratio(res$GeneRatio)

  if (order) {
    if (by == "Count") {
      idx <- order(res$Count, decreasing=TRUE)
    } else {
      idx <- order(res$GeneRatio, decreasing=TRUE)
    }
    res <- res[idx,]
  }

  if ( is.numeric(showCategory) ) {
    if ( showCategory <= nrow(res) ) {
      res <- res[1:showCategory,]
    }
  } else { ## selected categories
    res <- res[res$ID %in% showCategory,]
  }

  res$Description <- factor(res$Description,
                            levels=rev(res$Description))

  return(res)
}

GeneListNP = function(GeneList,GeneListName='CUSTOM',diseaseGeneID=NULL,GOMFenrich=T,GOBPenrich=T,GOCCenrich=T,KEGGenrich=T,PAenrich=T,DOenrich=T,MESHenrich=T,meshcategory='C',qvalueCutoff=0.05,TopNIDCurve=10){
  #diseaseGeneID:data.frame('class','geneID') or character
  require(data.table)
  require(stringr)
  require(RSQLite)
  require(plyr)
  require(ggplot2)
  require(pracma)
  ##############################################
  if (!is.null(diseaseGeneID)){
    if (class(diseaseGeneID)%in%c('data.table','data.frame')){
      diseaseGeneID$geneID=as.character(diseaseGeneID$geneID)
      DGeneID=diseaseGeneID$geneID
      names(DGeneID)=diseaseGeneID$class
      DiseaseName=unique(names(DGeneID))
      DiseaseNum=length(DiseaseName)
    }else{
      if (is.null(names(diseaseGeneID))){
        diseaseGeneID=as.character(unlist(str_split(diseaseGeneID,',| |;')))
        diseaseGeneID=str_replace_all(diseaseGeneID,'"','')
        diseaseGeneID=str_replace_all(diseaseGeneID,"'",'')
        diseaseGeneID=diseaseGeneID[diseaseGeneID!='']
        diseaseGeneID=unique(str_trim(diseaseGeneID,side='both'))
        DGeneID=diseaseGeneID
        names(DGeneID)=rep('Disease',length(diseaseGeneID))
        DiseaseName='Disease'
        DiseaseNum=1
      }else{
        diseaseGeneID=as.character(unlist(str_split(diseaseGeneID,',| |;')))
        diseaseGeneID=str_replace_all(diseaseGeneID,'"','')
        diseaseGeneID=str_replace_all(diseaseGeneID,"'",'')
        diseaseGeneID=diseaseGeneID[diseaseGeneID!='']
        diseaseGeneID=unique(str_trim(diseaseGeneID,side='both'))
        DGeneID=diseaseGeneID
        DiseaseName=unique(names(DGeneID))
        DiseaseNum=length(DiseaseName)
      }
    }
  }
#############################################################################################
#########################################################################################################
    selectGeneID=GeneList
    herbListName=GeneListName
    selectGeneName=geneTrans(selectGeneID,type='ID')
    ##############
        if (is.null(diseaseGeneID)){
          Order=paste('geneCluster=list(',herbListName,'=selectGeneID)',sep='')
        }else{
          if (DiseaseNum==1){
            Order=paste('geneCluster=list(',herbListName,'=selectGeneID,',DiseaseName,'=DGeneID[names(DGeneID)==DiseaseName])',sep='')
          }
          ##
          if (DiseaseNum==2){
            Order=paste('geneCluster=list(',herbListName,'=selectGeneID,',DiseaseName[1],'=DGeneID[names(DGeneID)==DiseaseName[1]],',DiseaseName[2],'=DGeneID[names(DGeneID)==DiseaseName[2]])',sep='')
          }
          if (DiseaseNum==3){
            Order=paste('geneCluster=list(',herbListName,'=selectGeneID,',DiseaseName[1],'=DGeneID[names(DGeneID)==DiseaseName[1]],',DiseaseName[2],'=DGeneID[names(DGeneID)==DiseaseName[2]],',DiseaseName[3],'=DGeneID[names(DGeneID)==DiseaseName[3]])',sep='')
          }
        }
        eval(parse(text=Order))
        ###################################
        if (KEGGenrich){
          if (class(try({EnrichKEGG=EnrichMentAnalysis(geneCluster,type='KEGG',qvalueCutoff=qvalueCutoff)},silent=T))=='try-error')EnrichKEGG=NULL
        }else{
          EnrichKEGG=NULL
        }
        ##############
        if (PAenrich){
          if (class(try({EnrichPA=EnrichMentAnalysis(geneCluster,type='ReactomePA',qvalueCutoff=qvalueCutoff)},silent=T))=='try-error')EnrichPA=NULL
        }else{
          EnrichPA=NULL
        }
        ############
        if (GOBPenrich){
          if (class(try({enrichGObp=EnrichMentAnalysis(geneCluster,type='GO',ontGO='BP',qvalueCutoff=qvalueCutoff)},silent=T))=='try-error')enrichGObp=NULL
        }else{
          enrichGObp=NULL
        }
        #############
        if (GOMFenrich){
          if (class(try({enrichGOmf=EnrichMentAnalysis(geneCluster,type='GO',ontGO='MF',qvalueCutoff=qvalueCutoff)},silent=T))=='try-error')enrichGOmf=NULL
        }else{
          enrichGOmf=NULL
        }
        ############
        if (GOCCenrich){
          if (class(try({enrichGOcc=EnrichMentAnalysis(geneCluster,type='GO',ontGO='CC',qvalueCutoff=qvalueCutoff)},silent=T))=='try-error')enrichGOcc=NULL
        }else{
          enrichGOcc=NULL
        }
        ############
        if (DOenrich){
          data(EG2DOLite)
          data(DOLiteTerm)
          data(DOLite2EG)
          if (class(try({enrichmentDO=EnrichMentAnalysis(geneCluster,type='DO',qvalueCutoff=qvalueCutoff)},silent=T))=='try-error')enrichmentDO=NULL
        }else{
          enrichmentDO=NULL
        }
        ############
        if (MESHenrich){
          if (class(try({enrichMESH=EnrichMentAnalysis(geneCluster,type='MESH',meshcategory=meshcategory,qvalueCutoff=qvalueCutoff)},silent=T))=='try-error')enrichMESH=NULL
        }else{
          enrichMESH=NULL
        }
        ##############################################################################
        if (!is.null(diseaseGeneID)&!is.null(EnrichKEGG)&(herbListName%in%EnrichKEGG$Enrichs$Cluster)){
          tempD=EnrichKEGG$Enrichs
          tempD$Cluster=as.character(tempD$Cluster)
          tempD$ID=as.character(tempD$ID)
          item=unique(tempD$Cluster)
          diseaseItem=item[!item%in%herbListName]
          CommonDiseaseKEGGSim=numeric()
          herbListNameItem=tempD$ID[tempD$Cluster==herbListName]
          if (length(diseaseItem)>0){
            length(CommonDiseaseKEGGSim)=length(diseaseItem)
            names(CommonDiseaseKEGGSim)=diseaseItem
            for (i in 1:length(diseaseItem)){
              diseaseItemtemp=tempD$ID[tempD$Cluster==diseaseItem[i]]
              CommonDiseaseKEGGSim[i]=length(intersect(herbListNameItem,diseaseItemtemp))/min(length(herbListNameItem),length(diseaseItemtemp))
            }
            diseaseToCurve=diseaseItem[1]
            diseasetNameItem=tempD$ID[tempD$Cluster==diseaseToCurve]
            TopNIDCurve=min(TopNIDCurve,length(herbListNameItem))
            common_KeggDisaeae=numeric()
            for (k in 1:TopNIDCurve){
              tempC=herbListNameItem[1:k]
              common_KeggDisaeae[k]=length(intersect(tempC,diseasetNameItem))/k
            }
            common_Keggdata=data.frame(ID=1:TopNIDCurve,Disease=common_KeggDisaeae)
            common_Keggdata=as.data.table(common_Keggdata)
            common_Keggdata=melt(common_Keggdata,id='ID',variable='Disease',value='Percentage')
            common_Keggdata=as.data.frame(common_Keggdata)
            CommonKEGGauc=trapz(1:TopNIDCurve,common_KeggDisaeae)
            PlotCommonKEGG=ggplot(common_Keggdata)+aes(x=ID,y=Percentage,group=Disease,colour=Disease,shape=Disease,fill=Disease)+geom_line()+geom_point(size=2)+xlab(paste('Top N enriched-Kegg pathway in ',herbListName))+ylab('Percentage of common enriched-Kegg pathway')+theme(legend.text=element_text(size=10),legend.title=element_text(size=10))+theme(axis.title.x=element_text(size=8),axis.title.y=element_text(size=8))+theme(axis.text.x=element_text(size=6),axis.text.y=element_text(size=6))
          }
          else{
            CommonDiseaseKEGGSim=NULL
            CommonKEGGauc=NULL
            PlotCommonKEGG=NULL
          }
          ######
        }else{
          CommonDiseaseKEGGSim=NULL
          CommonKEGGauc=NULL
          PlotCommonKEGG=NULL
        }
        ###############################################################################
        ###########ReactomePA#############################
        ##############################################################################
        if (!is.null(diseaseGeneID)&!is.null(EnrichPA)&(herbListName%in%EnrichPA$Enrichs$Cluster)){
          tempD=EnrichPA$Enrichs
          tempD$Cluster=as.character(tempD$Cluster)
          tempD$ID=as.character(tempD$ID)
          item=unique(tempD$Cluster)
          diseaseItem=item[!item%in%herbListName]
          CommonDiseasePASim=numeric()
          herbListNameItem=tempD$ID[tempD$Cluster==herbListName]
          if (length(diseaseItem)>0){
            length(CommonDiseasePASim)=length(diseaseItem)
            names(CommonDiseasePASim)=diseaseItem
            for (i in 1:length(diseaseItem)){
              diseaseItemtemp=tempD$ID[tempD$Cluster==diseaseItem[i]]
              CommonDiseasePASim[i]=length(intersect(herbListNameItem,diseaseItemtemp))/min(length(herbListNameItem),length(diseaseItemtemp))
            }
            diseaseToCurve=diseaseItem[1]
            diseasetNameItem=tempD$ID[tempD$Cluster==diseaseToCurve]
            TopNIDCurve=min(TopNIDCurve,length(herbListNameItem))
            common_PADisaeae=numeric()
            for (k in 1:TopNIDCurve){
              tempC=herbListNameItem[1:k]
              common_PADisaeae[k]=length(intersect(tempC,diseasetNameItem))/k
            }
            common_PAdata=data.frame(ID=1:TopNIDCurve,Disease=common_PADisaeae)
            common_PAdata=as.data.table(common_PAdata)
            common_PAdata=melt(common_PAdata,id='ID',variable='Disease',value='Percentage')
            common_PAdata=as.data.frame(common_PAdata)
            CommonPAauc=trapz(1:TopNIDCurve,common_PADisaeae)
            PlotCommonPA=ggplot(common_PAdata)+aes(x=ID,y=Percentage,group=Disease,colour=Disease,shape=Disease,fill=Disease)+geom_line()+geom_point(size=2)+xlab(paste('Top N enriched-Reactome pathway in ',herbListName))+ylab('Percentage of common enriched-Reactome pathway')+theme(legend.text=element_text(size=10),legend.title=element_text(size=10))+theme(axis.title.x=element_text(size=8),axis.title.y=element_text(size=8))+theme(axis.text.x=element_text(size=6),axis.text.y=element_text(size=6))
          }
          else{
            CommonDiseasePASim=NULL
            CommonPAauc=NULL
            PlotCommonPA=NULL
          }
          ######
        }else{
          CommonDiseasePASim=NULL
          CommonPAauc=NULL
          PlotCommonPA=NULL
        }
        ################PA END
        ################################################################################
        ################
        BlankGO=(!is.null(enrichGObp))|(!is.null(enrichGOcc))|(!is.null(enrichGOmf))
        BlankherbListName=(herbListName%in%enrichGObp$Enrichs$Cluster)|(herbListName%in%enrichGOcc$Enrichs$Cluster)|(herbListName%in%enrichGOmf$Enrichs$Cluster)
        if (!is.null(diseaseGeneID)&BlankGO&BlankherbListName){
          tempD=rbind(enrichGObp$Enrichs,enrichGOcc$Enrichs,enrichGOmf$Enrichs)
          tempD=arrange(tempD,qvalue)
          tempD$Cluster=as.character(tempD$Cluster)
          tempD$ID=as.character(tempD$ID)
          item=unique(tempD$Cluster)
          diseaseItem=item[!item%in%herbListName]
          CommonDiseaseGOSim=numeric()
          herbListNameItem=unique(tempD$ID[tempD$Cluster==herbListName])
          if (length(diseaseItem)>0){
            length(CommonDiseaseGOSim)=length(diseaseItem)
            names(CommonDiseaseGOSim)=diseaseItem
            for (i in 1:length(diseaseItem)){
              diseaseItemtemp=unique(tempD$ID[tempD$Cluster==diseaseItem[i]])
              CommonDiseaseGOSim[i]=length(intersect(herbListNameItem,diseaseItemtemp))/min(length(herbListNameItem),length(diseaseItemtemp))
            }
            diseaseToCurve=diseaseItem[1]
            diseasetNameItem=unique(tempD$ID[tempD$Cluster==diseaseToCurve])
            TopNIDCurve=min(TopNIDCurve,length(herbListNameItem))
            common_GODisaeae=numeric()
            for (k in 1:TopNIDCurve){
              tempC=herbListNameItem[1:k]
              common_GODisaeae[k]=length(intersect(tempC,diseasetNameItem))/k
            }
            common_GOdata=data.frame(ID=1:TopNIDCurve,Disease=common_GODisaeae)
            common_GOdata=as.data.table(common_GOdata)
            common_GOdata=melt(common_GOdata,id='ID',variable='Disease',value='Percentage')
            common_GOdata=as.data.frame(common_GOdata)
            CommonGOauc=trapz(1:TopNIDCurve,common_GODisaeae)
            PlotCommonGO=ggplot(common_GOdata)+aes(x=ID,y=Percentage,group=Disease,colour=Disease,shape=Disease,fill=Disease)+geom_line()+geom_point(size=2)+xlab(paste('Top N enriched-GO in ',herbListName))+ylab('Percentage of common enriched-GO')+theme(legend.text=element_text(size=10),legend.title=element_text(size=10))+theme(axis.title.x=element_text(size=8),axis.title.y=element_text(size=8))+theme(axis.text.x=element_text(size=6),axis.text.y=element_text(size=6))
          }else{
            CommonDiseaseGOSim=NULL
            CommonGOauc=NULL
            PlotCommonGO=NULL
          }
          ##############
        }else{
          CommonDiseaseGOSim=NULL
          CommonGOauc=NULL
          PlotCommonGO=NULL
        }

        result=list(selectGene=data.table(GeneID=selectGeneID,GeneName=selectGeneName),enrich=list(KEGG=EnrichKEGG,PA=EnrichPA,GObp=enrichGObp,GOmf=enrichGOmf,GOcc=enrichGOcc,DO=enrichmentDO,MESH=enrichMESH),Disease=list(CommonDiseaseKEGGSim=CommonDiseaseKEGGSim,CommonKEGGauc=CommonKEGGauc,PlotCommonKEGG=PlotCommonKEGG,CommonDiseasePASim=CommonDiseasePASim,CommonPAauc=CommonPAauc,PlotCommonPA=PlotCommonPA,CommonDiseaseGOSim=CommonDiseaseGOSim,CommonGOauc=CommonGOauc,PlotCommonGO=PlotCommonGO))
  #######################################
  return(result)
}

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

geneTrans = function(geneTerms,type){
  #geneTerms:geneID or geneName
  #type:'ID'-->ID to Name;'NAME'-->Name to ID
  require(clusterProfiler)
  require(data.table)
  require(plyr)
  geneTerms=as.character(geneTerms)
  if (type=='ID'){
    eg = bitr(geneTerms, fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Hs.eg.db")
    eg=as.data.table(eg)
    setkey(eg,ENTREZID)
    SYMBOL=eg[.(geneTerms),SYMBOL,mult='first']
    result=SYMBOL
	#return(SYMBOL)
  }else if (type=='NAME'){
    eg = bitr(geneTerms, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
    eg=as.data.table(eg)
    setkey(eg,SYMBOL)
    ENTREZID=eg[.(geneTerms),ENTREZID,mult='first']
	result=ENTREZID
    #return(ENTREZID)
  }
  return(result)
}

geneTransFUN = function(geneTerms,type){
#geneTerms:geneID or geneName
  #type:'ID'-->ID to Name;'NAME'-->Name to ID
  require(clusterProfiler)
  require(data.table)
  require(plyr)
  require(stringr)
  if (sum(str_detect(geneTerms,',| |;|/'),na.rm=T)>0){
    geneTerms=as.character(unlist(str_split(geneTerms,',| |;|/')))
    geneTerms=str_replace_all(geneTerms,'"','')
    geneTerms=str_replace_all(geneTerms,"'",'')
    geneTerms=geneTerms[geneTerms!='']
    geneTerms=unique(str_trim(geneTerms,side='both'))
	TransName=geneTrans(geneTerms=geneTerms,type=type)
    result=toString(TransName)
  }else{
    TransName=geneTrans(geneTerms=geneTerms,type=type)
	result=TransName
  }
  #TransName=geneTrans(geneTerms=geneTerms,type=type)
  #result=toString(TransName)
  return(result)
}

getChemblSmiles = function(chembl,type){
  #chembl:chemblID('CHEMBL150') or inchikey
  #type:'ID' or 'KEY'
  #output:c('Name','chemblID','smiles','InChiKey')
  require(chemblr)
  result=character()
  length(result)=5
  names(result)=c('Name','chemblID','smiles','InChiKey','InputName')
  result['InputName']=chembl
  if (type=='ID'){
    data=get.compounds(chembl,type='chemblid')
    if (ncol(data)>1){
      result['Name']=data$synonyms
      result['chemblID']=data$chemblId
      result['smiles']=data$smiles
      result['InChiKey']=data$stdInChiKey
      #result['InputName']=chembl
    }
  }else if (type=='KEY'){
    data=get.compounds(chembl,type='stdinchi')
    if (ncol(data)>1){
      result['Name']=data$synonyms
      result['chemblID']=data$chemblId
      result['smiles']=data$smiles
      result['InChiKey']=data$stdInChiKey
      #result['InputName']=chembl
    }
  }
  return(result)
}

getChemFromHerb = function(HerbList,QEDset=0.2,Database=c('HIT','TCMID','TCMSP','CUSTOM'),herb_ChemDataAll=herb_ChemDataAll){
  require(data.table)
  require(stringr)
  require(RSQLite)
  require(plyr)
  require(ggplot2)
  require(pracma)
  options(stringsAsFactors = F)
  if (str_detect(HerbList,',| |;')){
    HerbList=as.character(unlist(str_split(HerbList,',| |;')))
    HerbList=str_replace_all(HerbList,'"','')
    HerbList=str_replace_all(HerbList,"'",'')
    HerbList=HerbList[HerbList!='']
    HerbList=unique(str_trim(HerbList,side='both'))
  }
  HerbListData=herb_ChemDataAll[herb_ChemDataAll$herb%in%HerbList&herb_ChemDataAll$QED_DES>=QEDset&herb_ChemDataAll$database%in%Database,]
  return(HerbListData)
}

getChemGeneFromDatabase = function(ChemData,ChemDataType='inchikey',QEDset=0.2,geneScore=400,geneSelectPv=0.05,targetDatabase=c('HIT','TCMID','STITCH','TCMSP','CUSTOM'),TarDB='db/Tar.db'){
  #ChemData:data.frame:colnames-->'cid','chemical_name','key'.'key' is smiles or inchikey
  #ChemDataType:'inchikey' or 'smiles'
  #result:list(combinedGeneID=selectGeneID,GeneID_split=tempR_split,subD3=subD3):$GeneID_split-->data.frame(cid,chemical_name,inchikey,geneID).$subD3-->data.frame(cid,geneID,score,database,QED_DES,chemical_name,inchikey)
  require(data.table)
  require(stringr)
  require(RSQLite)
  require(plyr)
  require(ggplot2)
  require(pracma)
  options(stringsAsFactors = F)
  ChemData$key=as.character(ChemData$key)
  ChemData$cid=as.character(ChemData$cid)
  ChemData$chemical_name=as.character(ChemData$chemical_name)
  #setwd('D:/PA2.1Plus')
  ChemData2=ChemData
  if (ChemDataType!='inchikey'){
    tempSmiles=ChemData$key
    names(tempSmiles)=ChemData$cid
    ChemData2$inchikey=unlist(TransInchikey(tempSmiles,molecularType='SMItxt'))
  }else{
    ChemData2=plyr::rename(ChemData2,c('key'='inchikey'))
  }
  ############################
  ChemData2=ChemData2[!is.na(ChemData2$inchikey),]
  ChemData2=as.data.table(ChemData2)
  ChemData2=unique(ChemData2,by='inchikey')
  QueryKey=unique(ChemData2$inchikey)
  if (nrow(ChemData2)>0){
    paID2=paste("'",QueryKey,"'",sep='',collapse = ',')
    paID3=paste("(",paID2,")",sep='')
    chemIDOrder=paste("inchikey IN ",paID3,sep='')
    OrderchemIDOrder=paste("SELECT * FROM STITCH_HIT_TCMID_TCMSP_Tar where",chemIDOrder,sep=' ')
    tmp <- dbConnect(SQLite(), TarDB)
    res <- dbSendQuery(tmp, OrderchemIDOrder)
    subD0 <- fetch(res, n =-1)
    dbDisconnect(tmp)
    subD=join(subD0[,c('geneID','score','inchikey','database','QED_DES')],ChemData2,by='inchikey')
    subD=as.data.table(subD)
    subD=subD[QED_DES>=QEDset&score>=geneScore&database%in%targetDatabase,]
    rm(subD0)
    gc()
    subD$geneID=as.character(subD$geneID)
    subD2=unique(subD[,.(id=cid,geneID=geneID,score,database,QED_DES=QED_DES[1]),by=c('cid','geneID')])
    tempR=geneScoreCal(subD2,Pcutoff = geneSelectPv)
    tempR_split=subD2[,.(geneID=toString(geneID),QED_DES=QED_DES[1]),by='id']
    tempR_split=plyr::rename(tempR_split,c('id'='cid'))
    tempR_split=join(tempR_split,ChemData2,by='cid')
    tempR_split=tempR_split[,.(cid,chemical_name,inchikey,geneID,QED_DES)]
    selectGeneID=tempR$selectGeneId###selectGeneID
    subD3=subD2[,.(cid,geneID,score,database,QED_DES)]
    subD3=join(subD3,ChemData2,by='cid')
  }else{
    selectGeneID=NA
    tempR_split=NA
    subD3=NA
  }
    return(list(combinedGeneID=selectGeneID,GeneID_split=tempR_split,subD3=subD3))
  }

getChemGeneTargetScore = function(ChemData,diseaseGeneID,ChemDataType='inchikey',QEDset=0.2,geneScore=400,geneSelectPv=0.05,targetDatabase=c('HIT','TCMID','STITCH','TCMSP','CUSTOM'),Methodnet='KATZ',IFcombine=F){
  #ChemData:data.frame:colnames-->'cid','chemical_name','key'.'key' is smiles or inchikey
  #ChemDataType:'inchikey' or 'smiles'
  #Methodnet:'KATZ' or 'RW'
  #IFcombine:F-->split;T-->combined
  #output:ScoreGeneData-->data.frame:
  #(1)IFcombine=T:data.frame(SeedGeneIDPPI=0,NumSeedGeneID=0,ScoreGeneIDPPI=0,NumScoreGeneID=0,overlapScore=0,Path1Score=0,Path2Score=0,Path3Score=0,TotalScore=0,AdjTotalScore=0)
  #(2)IFcombine=F:data.frame(cid=0,chemical_name=0,inchikey=0,SeedGeneIDPPI=0,NumSeedGeneID=0,ScoreGeneIDPPI=0,NumScoreGeneID=0,overlapScore=0,Path1Score=0,Path2Score=0,Path3Score=0,TotalScore=0,AdjTotalScore=0)
  require(data.table)
  require(stringr)
  require(RSQLite)
  require(plyr)
  require(ggplot2)
  require(pracma)
  options(stringsAsFactors = F)
  getChemGeneFromDatabaseResult=getChemGeneFromDatabase(ChemData=ChemData,ChemDataType=ChemDataType,QEDset=QEDset,geneScore=geneScore,geneSelectPv=geneSelectPv,targetDatabase=targetDatabase)
  if (IFcombine){
    queryGeneID=getChemGeneFromDatabaseResult$combinedGeneID
    NetCorrResult=NetCorr(SeedGeneID=diseaseGeneID,ScoreGeneID=queryGeneID,Methodnet=Methodnet,randomtimes=1,gamma=0.7,testM=1)
    ScoreGeneData=data.frame(SeedGeneIDPPI=0,NumSeedGeneID=0,ScoreGeneIDPPI=0,NumScoreGeneID=0,overlapScore=0,Path1Score=0,Path2Score=0,Path3Score=0,TotalScore=0,AdjTotalScore=0)
    if (Methodnet=='KATZ'){
      KATZresult=NetCorrResult$resultKATZ
      #KATZresult=data.frame(SeedGeneIDPPI=toString(SeedGeneID),NumSeedGeneID=length(SeedGeneID),ScoreGeneIDPPI=toString(ScoreGeneID),NumScoreGeneID=length(ScoreGeneID),overlapScore=0,Path1Score=0,Path2Score=0,Path3Score=0,TotalScore=0,RandomScoreMedian=0,Pvalue=1)
      ScoreGeneData$SeedGeneIDPPI=KATZresult$SeedGeneIDPPI
      ScoreGeneData$NumSeedGeneID=KATZresult$NumSeedGeneID
      ScoreGeneData$ScoreGeneIDPPI=KATZresult$ScoreGeneIDPPI
      ScoreGeneData$NumScoreGeneID=KATZresult$NumScoreGeneID
      ScoreGeneData$overlapScore=KATZresult$overlapScore
      ScoreGeneData$Path1Score=KATZresult$Path1Score
      ScoreGeneData$Path2Score=KATZresult$Path2Score
      ScoreGeneData$Path3Score=KATZresult$Path3Score
      ScoreGeneData$TotalScore=KATZresult$TotalScore
      ##score/n0*exp(score/n1-1)*1000
      #ScoreGeneData$AdjTotalScore=KATZresult$TotalScore/KATZresult$NumSeedGeneID*exp(KATZresult$TotalScore/KATZresult$NumScoreGeneID-1)*1000
      #score/n0*log10(1+n0/n1/2)
      ScoreGeneData$AdjTotalScore=KATZresult$TotalScore/KATZresult$NumSeedGeneID*log10(1+KATZresult$NumSeedGeneID/KATZresult$NumScoreGeneID/2)*1000
    }else if (Methodnet=='RW'){
      RWresult=NetCorrResult$resultRWTable
      #resultRWTable:data.frame(SeedGeneIDPPI=toString(resultRW$SeedGeneID),NumSeedGeneID=resultRW$NumSeedGeneIDPPI,ScoreGeneIDPPI=toString(resultRW$ScoreGeneID),NumScoreGeneID=resultRW$NumScoreGeneIDPPI,TotalScore=resultRW$TotalScore,RandomScoreMean=resultRW$RandomScoreMedian,Pvalue=resultRW$Pvalue)
      ScoreGeneData$SeedGeneIDPPI=RWresult$SeedGeneIDPPI
      ScoreGeneData$NumSeedGeneID=RWresult$NumSeedGeneID
      ScoreGeneData$ScoreGeneIDPPI=RWresult$ScoreGeneIDPPI
      ScoreGeneData$NumScoreGeneID=RWresult$NumScoreGeneID
      ScoreGeneData$TotalScore=RWresult$TotalScore
      ##score/n0*exp(score/n1-1)*1000
      #ScoreGeneData$AdjTotalScore=RWresult$TotalScore/RWresult$NumSeedGeneID*exp(RWresult$TotalScore/RWresult$NumScoreGeneID-1)
      #score/n0*log10(1+n0/n1/2)
      ScoreGeneData$AdjTotalScore=RWresult$TotalScore/RWresult$NumSeedGeneID*log10(1+RWresult$NumSeedGeneID/RWresult$NumScoreGeneID/2)
    }
  }else{
    queryGeneIDdata=getChemGeneFromDatabaseResult$GeneID_split
    queryCID=unique(queryGeneIDdata$cid)
    ScoreGeneData=data.frame()
    for (i in 1:length(queryCID)){
      tempqueryGeneID=queryGeneIDdata$geneID[queryGeneIDdata$cid==queryCID[i]]
      tempNetCorrResult=NetCorr(SeedGeneID=diseaseGeneID,ScoreGeneID=tempqueryGeneID,Methodnet=Methodnet,randomtimes=1,gamma=0.7,testM=1)
      tempScoreGeneData=data.frame(cid=0,chemical_name=0,inchikey=0,SeedGeneIDPPI=0,NumSeedGeneID=0,ScoreGeneIDPPI=0,NumScoreGeneID=0,overlapScore=0,Path1Score=0,Path2Score=0,Path3Score=0,TotalScore=0,AdjTotalScore=0)
      if (Methodnet=='KATZ'){
        tempKATZresult=tempNetCorrResult$resultKATZ
        tempScoreGeneData$cid=queryCID[i]
        tempScoreGeneData$chemical_name=queryGeneIDdata$chemical_name[queryGeneIDdata$cid==queryCID[i]][1]
        tempScoreGeneData$inchikey=queryGeneIDdata$inchikey[queryGeneIDdata$cid==queryCID[i]][1]
        tempScoreGeneData$SeedGeneIDPPI=tempKATZresult$SeedGeneIDPPI
        tempScoreGeneData$NumSeedGeneID=tempKATZresult$NumSeedGeneID
        tempScoreGeneData$ScoreGeneIDPPI=tempKATZresult$ScoreGeneIDPPI
        tempScoreGeneData$NumScoreGeneID=tempKATZresult$NumScoreGeneID
        tempScoreGeneData$overlapScore=tempKATZresult$overlapScore
        tempScoreGeneData$Path1Score=tempKATZresult$Path1Score
        tempScoreGeneData$Path2Score=tempKATZresult$Path2Score
        tempScoreGeneData$Path3Score=tempKATZresult$Path3Score
        tempScoreGeneData$TotalScore=tempKATZresult$TotalScore
        ##score/n0*exp(score/n1-1)*1000
        ###score/n0*log10(1+n0/n1/2)
        #tempScoreGeneData$AdjTotalScore=tempKATZresult$TotalScore/tempKATZresult$NumSeedGeneID*exp(tempKATZresult$TotalScore/tempKATZresult$NumScoreGeneID-1)*1000
        tempScoreGeneData$AdjTotalScore=tempKATZresult$TotalScore/tempKATZresult$NumSeedGeneID*log10(1+tempKATZresult$NumSeedGeneID/tempKATZresult$NumScoreGeneID/2)*1000
        ScoreGeneData=rbind(ScoreGeneData,tempScoreGeneData)
      }else if (Methodnet=='RW'){
        tempKATZresult=tempNetCorrResult$resultRWTable
        tempScoreGeneData$cid=queryCID[i]
        tempScoreGeneData$chemical_name=queryGeneIDdata$chemical_name[queryGeneIDdata$cid==queryCID[i]][1]
        tempScoreGeneData$inchikey=queryGeneIDdata$inchikey[queryGeneIDdata$cid==queryCID[i]][1]
        tempScoreGeneData$SeedGeneIDPPI=tempKATZresult$SeedGeneIDPPI
        tempScoreGeneData$NumSeedGeneID=tempKATZresult$NumSeedGeneID
        tempScoreGeneData$ScoreGeneIDPPI=tempKATZresult$ScoreGeneIDPPI
        tempScoreGeneData$NumScoreGeneID=tempKATZresult$NumScoreGeneID
        ##score/n0*exp(score/n1-1)*1000
        ###score/n0*log10(1+n0/n1/2)
        #tempScoreGeneData$AdjTotalScore=tempKATZresult$TotalScore/tempKATZresult$NumSeedGeneID*exp(tempKATZresult$TotalScore/tempKATZresult$NumScoreGeneID-1)*1000
        tempScoreGeneData$AdjTotalScore=tempKATZresult$TotalScore/tempKATZresult$NumSeedGeneID*log10(1+tempKATZresult$NumSeedGeneID/tempKATZresult$NumScoreGeneID/2)
        ScoreGeneData=rbind(ScoreGeneData,tempScoreGeneData)
      }
    }
  }
  return(ScoreGeneData)
}

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
  library(readr)
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

getCTDId_and_pubID = function(name,cookiefileS=NULL){####chemical_name--->CTDMeshID&pubID
  library(RCurl)
  library(stringr)
  library(rvest)
  library(httr)
  library(data.table)
  #name='DDFG'
  name=tolower(name)
  name2=str_replace_all(name,' ','%20')
  name3=paste('%22',name2,'%22',sep='')
  myheader=c('User-Agent'='Mozilla/5.0 (Windows NT 6.1; WOW64)','Accept'='text/html,application/xhyml+xml,application/xml;q=0.9,*/*;q=0.8','Accept-Language'='zh-CN,zh;q=0.8,en;q=0.6','Connection'='keep-alive','Accept-Charset'='GB2312,utf-8;q=0.7,*;q=0.7')
  chemCSV_url=paste('http://ctdbase.org/basicQuery.go?bq=',name3,'&6578706f7274=1&d-4029212-e=1&bqCat=chem',sep='')
  if (c(class(try({chemTable=suppressWarnings(fread(chemCSV_url))},silent=T))=='try-error')[1]){
    MeshID=NA
    pubID=NA
    sure=NA
  }else{
    #chemTable=suppressWarnings(fread(chemCSV_url))
    if (sum(str_detect(chemTable,'<meta'))>0){
      d2 =debugGatherer()
      cHandle2<- getCurlHandle(httpheader=myheader,followlocation=1,
                               debugfunction=d2$update,verbose=TRUE,
                               cookiefile=cookiefileS)
      au=getURL(chemCSV_url,curl=cHandle2,.encoding="UTF-8")
      MeshID0=str_extract(au,'term=[A-Z0-9]+')
      MeshID=str_replace(MeshID0,'term=','')
    }else{
      chemTable[,Chemical:=tolower(Chemical)]
      setkey(chemTable,Chemical)
      if (is.na(as.character(chemTable[name,'Accession ID',with=F]))){
        MeshID=as.character(chemTable[1,'Accession ID',with=FALSE])
      }else{
        MeshID=as.character(chemTable[name,'Accession ID',with=F])
      }
    }
    pubSurl=paste('http://www.ncbi.nlm.nih.gov/pcsubstance/?term=',MeshID,sep='')
    d2 =debugGatherer()
    cHandle2<- getCurlHandle(httpheader=myheader,followlocation=1,
                             debugfunction=d2$update,verbose=TRUE,
                             cookiefile=cookiefileS)
    au<- getURL(pubSurl,curl=cHandle2,.encoding="UTF-8")
    #au=getURL(pubSurl,httpheader=myheader,followlocation=T)
    if (str_detect(au,'No items found')){
      pubID=NA
      sure=NA
    }else{
      sid0=str_extract(au,'SID:[0-9]+')
      sid=str_extract(sid0,'[0-9]+')
      url_sid=paste('pubchem.ncbi.nlm.nih.gov/rest/rdf/substance/SID',sid,'.html',sep='')
      a=getURL(url_sid,curl=cHandle2,.encoding="UTF-8")
      #       lo=str_locate(au,'data-pubchem-title=')
      #       lo2=str_locate(au,'data-pubchem-record')
      #       temp=str_sub(au,lo[1,2]+1,lo2[1,1]-3)
      #       tempName=str_replace_all(temp,'\"','')
      #       result=getPubID(tempName)
      tempcid=str_extract(a,'CID[0-9]+')
      pubID=str_extract(tempcid,'[0-9]+')
      sure='ok'
    }

  }
  Re=c(MeshID=MeshID,pubID=pubID,sure=sure)
  names(Re)=c('MeshID','pubID','sure')
  return(Re)
}

getDCdiseaseGene = function(keywords,Type=c('DEG','OMIM','GWAS')){
  require(data.table)
  require(stringr)
  require(clusterProfiler)
  ##get "http://disease-connect.org/" disease Genes
#DCdatabase=fread('E:/DATABASE/Disease-Gene/Disease-Gene.csv',colClasses = 'character')
#DCdatabase$Diseas=tolower(DCdatabase$Diseas)
if (!exists('DCdatabase')){
  load('db/DCdatabase.RData')
}
keywords=tolower(keywords)
result=DCdatabase[(DCdatabase$Diseas%like%keywords)&(DCdatabase$Type%in%Type),]
if (nrow(result)!=0){
  genes=unique(result$GeneName)
  #genesTran=bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  #geneID=genesTran$ENTREZID
  geneID=unique(result$geneID)
}else{
  genes=''
  geneID=''
}
return(list(geneNames=genes,geneID=geneID,result=result))
}

getDiseaseGene = function(keyWords,DGdatabase=c('OMIM','TTD','DC','GeneCard'),GeneCardScore=10){
  require(data.table)
  library(RCurl)
  library(stringr)
  library(XML)
  library(readr)
  require(clusterProfiler)
  GeneCardScore=as.numeric(GeneCardScore)
  if ('OMIM'%in%DGdatabase){
    result1=getOMIMdiseaseGene(keyWords)
	resultFromDC=getDCdiseaseGene(keywords=keyWords,Type='OMIM')
	result1$geneNames=unique(c(result1$geneNames,resultFromDC$geneNames))
	result1$geneID=unique(c(result1$geneID,resultFromDC$geneID))
  }else{
    result1=list(geneNames=character(),geneID=character(),result=data.frame())
  }
  ##
  if ('TTD'%in%DGdatabase){
    result2=getTTDdiseaseGene(keyWords)
  }else{
    result2=list(geneNames=character(),geneID=character(),result=data.frame())
  }
  ##
  if ('DC'%in%DGdatabase){
    result3=getDCdiseaseGene(keyWords)
  }else{
    result3=list(geneNames=character(),geneID=character(),result=data.frame())
  }
  #
  if ('GeneCard'%in%DGdatabase){
    result4=getGeneCardGene(keywords=keyWords,score=GeneCardScore)
  }else{
    result4=list(geneNames=character(),geneID=character(),result=data.frame())
  }
  ##
  result=list(geneNames=unique(c(result1$geneNames,result2$geneNames,result3$geneNames,result4$geneNames)),geneID=unique(c(result1$geneID,result2$geneID,result3$geneID,result4$geneID)),result=list(OMIM=result1$result,TTD=result2$result,DC=result3$result,GeneCard=result4$result))
  return(result)
}

getGeneCardGene = function(keywords,score=10){
  require(data.table)
  require(stringr)
  require(RCurl)
  require(clusterProfiler)
  require(readr)
  #http://www.genecards.org/Search/Export?queryString=%27fat%20liver%27&searchType=Keywords
  url1='http://www.genecards.org/Search/Export?queryString=%27'
  url2='%27&searchType=Keywords'
  keywords=str_replace_all(keywords,' ','%20')
  url=paste(url1,keywords,url2,sep='')
  ##
  myHttpheader=c('User-Agent'='Mozilla/5.0(Windows;U;Windows NT 6.1;zh-CN;rv:1.9.1.6)','Accept'='text/html,application/xhyml+xml,application/xml;q=0.9,*/*;q=0.8','Accept-Language'='en-us','Connection'='keep-alive','Accept-Charset'='GB2312,utf-8;q=0.7,*;q=0.7')
  d2 =debugGatherer()
  cHandle2<- getCurlHandle(httpheader=myHttpheader,followlocation=1,
                           debugfunction=d2$update,verbose=TRUE)
  temp2<- getURLContent(url,curl=cHandle2,.encoding="UTF-8")
  result=read_csv(temp2)
  result=result[result$`Relevance score`>score,]
  if (nrow(result)!=0){
    geneName0=unique(result$`Gene Symbol`)
    geneID=geneTrans(geneName0,type = 'NAME')
    geneID=unique(geneID[!is.na(geneID)])
    geneName=geneTrans(geneID,type = 'ID')
    #genesTran=bitr(geneID, fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Hs.eg.db")
    #genes=genesTran$SYMBOL
  }else{
    geneID=''
    geneName=''
  }
  return(list(geneNames=geneName,geneID=geneID,result=result))
}

getGeneCardGeneStep1 = function(keywords){
  require(data.table)
  require(stringr)
  require(RCurl)
  require(clusterProfiler)
  require(readr)
  require(rvest)
  require(stringdist)
  #library(reticulate)
  #library(jsonlite)
  #library(xml2)
  library("httr")
  #setwd('D:/PA2.1Plus/')
  url1='https://www.malacards.org/search/results?query='
  #url2='%27&searchType=Keywords'
  #keywords='fatty liver'
  keywords2=str_replace_all(keywords,' ','+')
  url=paste(url1,keywords2,sep='')
  UA <- "Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:93.0) Gecko/20100101 Firefox/93.0"###update 20230212
  a=GET(url,add_headers(`Connection` = "keep-alive", `User-Agent` = UA),config(sslversion=6,ssl_verifypeer=1))###update 20230212
  b=httr::content(a,'text')
  bb=read_html(b)
  tablesBB=html_table(bb)
  #tablesBB=tablesBB[[1]]###Before 20230212
  tablesBB=tablesBB[[2]]###update 20230212
  tablesBB=tablesBB[!is.na(tablesBB$MIFTS),]
  #TempsimScore=stringsim(keywords,tablesBB$Name,method='soundex')
  #tablesBB$simScore=TempsimScore
  tablesBB=as.data.table(tablesBB)
  tablesBB=tablesBB[order(-Score)]
  tablesBB=tablesBB[,-2]
  #write.csv(tablesBB,'Report/GeneCardSearchStep1.csv',row.names = F)
  return(tablesBB)
}

getGeneCardGeneStep2 = function(MCID){
  require(data.table)
  require(stringr)
  require(RCurl)
  require(clusterProfiler)
  require(readr)
  require(rvest)
  require(stringdist)
  #library(reticulate)
  #library(jsonlite)
  #library(xml2)
  library("httr")
  options(stringsAsFactors = F)
  Symbol_detect=function(x){
    require(stringr)
    y=numeric()
    temp=unlist(x)
    y=sum(str_detect(unlist(temp),'Symbol'),na.rm =T)
    return(y)
  }
  #setwd('D:/PA2.1Plus/')
  MCID2=unlist(str_split(MCID,',|，|;'))
  MCID2=str_trim(MCID2,'both')
  MCID2=unique(MCID2)
  url1='https://www.malacards.org/card/'
  urlEND='?limit[RelatedGenes]=0#RelatedGenes-table'
  UA <- "Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:93.0) Gecko/20100101 Firefox/93.0"###update 20230212
  searchResult=data.frame()
  for (i in MCID2){
    url=paste0(url1,i,urlEND)
    a=GET(url,add_headers(`Connection` = "keep-alive", `User-Agent` = UA),config(sslversion=6,ssl_verifypeer=1))###update 20230212
    b=httr::content(a,'text')
    bb=read_html(b)
    tableNOde=html_node(x=bb,xpath='//*[@id="RelatedGenes-table"]')
    if (!is.na(tableNOde)){
      tablesBB=html_table(tableNOde)
      #GeneID=geneTrans(geneTerms=tablesBB$Symbol,type='NAME')
      if (class(try({GeneID=geneTrans(geneTerms=tablesBB$Symbol,type='NAME')},silent = T))!="try-error"){
        tablesBB$GeneID=GeneID
        tablesBB$MCID=i
        #tablesBB=tablesBB[,c('MCID','Symbol','GeneID','Description','Category','Score','Evidence','PubmedIds')]###Before 20230212
        tablesBB=tablesBB[,c('MCID','Symbol','GeneID','Description','Category','Score','Evidence','PubMed IDs')]###update 20230212
      }else{
        tablesBB$GeneID=NA
        tablesBB$MCID=i
        #tablesBB=tablesBB[,c('MCID','Symbol','GeneID','Description','Category','Score','Evidence','PubmedIds')]###Before 20230212
        tablesBB=tablesBB[,c('MCID','Symbol','GeneID','Description','Category','Score','Evidence','PubMed IDs')]###update 20230212
      }
    }else{
      tablesBB=data.frame()
    }
    searchResult=rbind(searchResult,tablesBB)
  }
  if (nrow(searchResult)>0){
    CombineResult=searchResult[,c('Symbol','GeneID','Description','Category','Score','Evidence')]
    CombineResult=as.data.table(CombineResult)
    CombineResult=unique(CombineResult,by='GeneID')
  }else{
    CombineResult=searchResult
  }
  return(list(searchResult=searchResult,CombineResult=CombineResult))
}

getHerbFromdatabase = function(ChemData,ChemDataType='inchikey',targetDatabase=c('HIT','TCMID','STITCH','TCMSP','CUSTOM'),IFFromtargetData=F,AllHerbListData=AllHerbListData){
  #ChemData:data.frame:colnames-->'cid','chemical_name','key'.'key' is smiles or inchikey
  #IFFromtargetData:F-->From db/Tar.db.
  #IFFromtargetData:T-->From herb_ChemDataAll
  require(data.table)
  require(stringr)
  require(RSQLite)
  require(plyr)
  require(ggplot2)
  require(pracma)
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
  if (!exists('AllHerbListData')){
    #ALL=fread('db/AllHerbListData.csv')
    load('db/AllHerbListData.RData')
  }
  ############################
  ChemData2=ChemData2[!is.na(ChemData2$inchikey),]
  ChemData2=as.data.table(ChemData2)
  ChemData2=unique(ChemData2,by='inchikey')
  QueryKey=unique(ChemData2$inchikey)
if (!IFFromtargetData){
  if (nrow(ChemData2)>0){
    paID2=paste("'",QueryKey,"'",sep='',collapse = ',')
    paID3=paste("(",paID2,")",sep='')
    chemIDOrder=paste("inchikey IN ",paID3,sep='')
    OrderchemIDOrder=paste("SELECT * FROM STITCH_HIT_TCMID_TCMSP_Tar where",chemIDOrder,sep=' ')
    tmp <- dbConnect(SQLite(), 'db/Tar.db')
    res <- dbSendQuery(tmp, OrderchemIDOrder)
    subD0 <- fetch(res, n =-1)
    #dbDisconnect(tmp)
    ##########
    subD0_herbID=unique(subD0$herbID)
    paID2=paste("'",subD0_herbID,"'",sep='',collapse = ',')
    paID3=paste("(",paID2,")",sep='')
    herbIDOrder=paste("herbID IN ",paID3,sep='')
    OrderherbIDOrder=paste("SELECT * FROM STITCH_HIT_TCMID_TCMSP_Tar where",herbIDOrder,sep=' ')
    #tmp <- dbConnect(SQLite(), 'db/Tar.db')
    res_herbID <- dbSendQuery(tmp, OrderherbIDOrder)
    subD0_herb <- fetch(res_herbID, n =-1)
    dbDisconnect(tmp)
    #####################
    subD0_herb=join(subD0_herb,AllHerbListData,by='herbID')
    subD0_herb=as.data.table(subD0_herb)
    subD0_herb=subD0_herb[,.(TotalChem=length(unique(inchikey))),by='herb']
    #####################
    subD=join(subD0[,c('herbID','inchikey')],ChemData2,by='inchikey')
    subD=join(subD,AllHerbListData,by='herbID')
    subD=as.data.table(subD)
    subD=subD[database%in%targetDatabase,]
    subD=subD[,.(cid,chemical_name,herb,inchikey,QED_DES)]
    subD=unique(subD,by=c('herb','inchikey'))
    rm(subD0)
    gc()
    subDSummary=subD[,.(chemical=toString(unique(cid)),NumChem=length(unique(cid))),by='herb']
    subDSummary=join(subDSummary,subD0_herb,by='herb')
    subDSummary$BackgroundRatio=subDSummary$NumChem/subDSummary$TotalChem
    subDSummary$QueryRatio=subDSummary$NumChem/length(QueryKey)
    subDSummary=subDSummary[order(-NumChem)]
  }else{
    subD=NA
    subDSummary=NA
  }
}else{
  if (nrow(ChemData2)>0){
    subD0=herb_ChemDataAll[herb_ChemDataAll$inchikey%in%QueryKey&herb_ChemDataAll$database%in%targetDatabase,]
    subD0_herb=unique(subD0$herb)
    subD0_herbData=herb_ChemDataAll[herb_ChemDataAll$herb%in%subD0_herb,]
    subD0_herbData=as.data.table(subD0_herbData)
    subD0_herbData=subD0_herbData[,.(TotalChem=length(unique(inchikey))),by='herb']
    subD0=join(subD0,ChemData2,by='inchikey')
    subD=as.data.table(subD0)
    subD=subD[,.(cid,chemical_name,herb,inchikey,QED_DES)]
    subD=unique(subD,by=c('herb','inchikey'))
    subDSummary=subD[,.(chemical=toString(unique(cid)),NumChem=length(unique(cid))),by='herb']
    subDSummary=join(subDSummary,subD0_herbData,by='herb')
    subDSummary$BackgroundRatio=subDSummary$NumChem/subDSummary$TotalChem
    subDSummary$QueryRatio=subDSummary$NumChem/length(QueryKey)
    subDSummary=subDSummary[order(-NumChem)]
  }else{
    subD=NA
    subDSummary=NA
  }
}
  return(list(subD=subD,subDSummary=subDSummary))
}

getHerbGene = function(herbName,QEDset=0.2,geneScore=800,geneSelectPv=0.01,targetDatabase=c('HIT','TCMID','STITCH','TCMSP','CUSTOM'),AllHerbListData=AllHerbListData){
  require(data.table)
  require(stringr)
  require(RSQLite)
  require(plyr)
  require(ggplot2)
  require(pracma)
  if (!exists('AllHerbListData')){
    #ALL=fread('db/AllHerbListData.csv')
    AllHerbListData=fread('db/AllHerbListData.csv',encoding='UTF-8')
  }
  herbList=herbName
  herbListALL=AllHerbListData$herb
  herbList2=unlist(str_split(herbList,',| |;'))
  actHerb=herbList2[herbList2%in%herbListALL]
  actHerbID=AllHerbListData[AllHerbListData$herb%in%actHerb,herbID]
  if (length(actHerbID)!=0){
    paID2=paste("'",actHerbID,"'",sep='',collapse = ',')
    paID3=paste("(",paID2,")",sep='')
    herbIDOrder=paste("herbID IN ",paID3,sep='')
    OrderherbID=paste("SELECT * FROM STITCH_HIT_TCMID_TCMSP_Tar where",herbIDOrder,sep=' ')
    tmp <- dbConnect(SQLite(), 'db/Tar.db')
    res <- dbSendQuery(tmp, OrderherbID)
    subD0 <- fetch(res, n =-1)
    dbDisconnect(tmp)
    subD=join(subD0,AllHerbListData,by='herbID')
    subD=as.data.table(subD)
    subD=subD[,!'herbID',with=F]
    subD=subD[,c('herb',setdiff(colnames(subD),'herb')),with=F]
    subD=subD[QED_DES>=QEDset&score>=geneScore&database%in%targetDatabase,]
    rm(subD0)
    gc()
    subD$geneID=as.character(subD$geneID)
    tempR=geneScoreCal(unique(subD[,.(herb,id=cid,geneID=geneID)]),Pcutoff = geneSelectPv)
    selectGeneID=tempR$selectGeneId###selectGeneID
  }else{
    selectGeneID=character()
  }
  if (length(selectGeneID)<1)selectGeneID='NA'
  return(selectGeneID)
}

getOMIMdiseaseGene = function(keyWords){
  library(RCurl)
  library(stringr)
  library(XML)
  library(readr)
  require(clusterProfiler)
  keyWords=str_replace_all(keyWords,' ','+')
  url=paste('http://www.omim.org/search/?index=geneMap&search=','"',keyWords,'"','&start=1&limit=20000&format=tsv',sep='')
  myHttpheader=c('User-Agent'='Mozilla/5.0(Windows;U;Windows NT 6.1;zh-CN;rv:1.9.1.6)','Accept'='text/html,application/xhyml+xml,application/xml;q=0.9,*/*;q=0.8','Accept-Language'='en-us','Connection'='keep-alive','Accept-Charset'='GB2312,utf-8;q=0.7,*;q=0.7')
  d2 =debugGatherer()
  cHandle2<- getCurlHandle(httpheader=myHttpheader,followlocation=1,
                           debugfunction=d2$update,verbose=TRUE)
  temp2<- getURL(url,curl=cHandle2,.encoding="UTF-8")
  content=htmlTreeParse(temp2, asText = TRUE)
  a=content$children$html[['body']]
  b=a[['p']][['text']]
  b=as.character(b)
  write_lines(b,'temp.txt')
  bb=read_lines('temp.txt')
  blank=which(bb=='')
  Copyright=which(str_detect(bb,'Copyright'))
  lo1=blank[which(blank>=Copyright)[1]]
  lo2=blank[which(blank>=Copyright)[2]]
  result=character()
  for (i in (lo1+1):(lo2-1)){
    tempD=unlist(str_split_fixed(bb[i],'\t',n=11))
    result=rbind(result,tempD)
  }
  colnames(result)=unlist(result[1,])
  result=result[-1,,drop=F]
  #result=apply(result,2,as.character)
  result=as.data.frame(result,stringsAsFactors = F)
  if (nrow(result)==0){
    genes=''
    geneID=''
  }else{
    genes=paste(result$`Gene/Locus`,collapse = ',')
    genes=unique(str_trim(unlist(str_split(genes,',')),side='both'))
    #genesTran=bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
    #geneID=genesTran$ENTREZID
    if (class(try({geneID=geneTrans(genes,type = 'ID')},silent = T))=='try-error')geneID=NA
  }

  return(list(geneNames=genes,geneID=geneID,result=result))
}

getPubChemSMI = function(PubChemID){
  ##PubChemID:PubChemID
  ##output:smiles
  library(XML)
  library(stringr)
  library(RCurl)
  library(ChemmineR)
  URL1='https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/'
  URL2='/record/SDF/?record_type=2d&response_type=save&response_basename=Structure2D_CID_'
  URL=paste(URL1,PubChemID,URL2,PubChemID,sep='')
  #URL='https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/100982/record/SDF/?record_type=2d&response_type=save&response_basename=Structure2D_CID_100982'
  SdfFile=read.SDFset(URL)
  SMIob=sdf2smiles(SdfFile)
  SMI=as.character(SMIob)
  return(SMI)
}

getPubID = function(name,cookiefileS=NULL){###chemical_name-->get PubmedChemID
  library(RCurl)
  library(stringr)
  library(rvest)
  library(httr)
  library(XML)
  name=tolower(name)
  #name='protocatechualdehyde'
  if (str_detect(name,'\\s')){
    name2=str_replace_all(name,' ','+')
    #name3=paste('%22',name2,'%22',sep='')
    name3=name2
  }else{
    name3=name
  }
  url=paste('http://www.ncbi.nlm.nih.gov/pccompound/?term=',name3,sep='')
  myheader=c('User-Agent'='Mozilla/5.0 (Windows NT 6.1; WOW64)','Accept'='text/html,application/xhyml+xml,application/xml;q=0.9,*/*;q=0.8','Accept-Language'='zh-CN,zh;q=0.8,en;q=0.6','Connection'='keep-alive','Accept-Charset'='GB2312,utf-8;q=0.7,*;q=0.7')
  d2 =debugGatherer()
  cHandle2<- getCurlHandle(httpheader=myheader,followlocation=1,
                           debugfunction=d2$update,verbose=TRUE,
                           cookiefile=cookiefileS)
  au<- getURL(url,curl=cHandle2,.encoding="UTF-8")
  #au=getURL(url,httpheader=myheader,followlocation=T)
  if (str_detect(au,'title=\"CID:')){
    #txt=str_sub(au,str_locate(au,'title=\"CID:')[1,2]+1,str_locate(au,'title=\"CID:')[1,2]+15)
    txt=str_extract(au,'CID[0-9]+')
    cid=str_extract(txt,pattern='[0-9]+')
    sure='ok'
  }else if (str_detect(au,'No items found')){
    cid=NA
    sure=NA
  }else if (str_detect(au,'Pccompound_ResultsPanel')){
    #   txt=str_sub(au,str_locate(au,'CID:')[1,2]+1,str_locate(au,'CID:')[1,2]+25)
    #   cid=str_extract(txt,pattern='[0-9]+')
    #   sure='multi'
    #aurl=getURL(url,httpheader=myheader,followlocation=T)
    bparse=htmlParse(au)
    bparse2=xpathSApply(bparse,"//p[@class='title']",fun=xmlValue)
    bparse0=xpathSApply(bparse,"//p[@class='title']")
    attrtemp=lapply(bparse0,xmlToList)
    href=sapply(attrtemp,tempcatch)
    tempsplit=lapply(bparse2,str_split,pattern=';')
    tempsplit2=lapply(tempsplit,unlist)
    exact_sum=unlist(lapply(tempsplit2,function(x)sum(tolower(x)%in%name)))
    if (sum(exact_sum)>0){
      cid=href[which(exact_sum==1)[1]]
      sure='ok'
    }else{
      cid=href[1]
      sure='multi'
    }
  }
  return(c(cid=cid,sure=sure))##cid-->PubmedChemID;sure-->ok:only one match;multi-->multi match;NA-->not available
}

getSpiderChem = function(name,cookiefileS=NULL){####chemical_name-->get SpederID&Smiles
  library(XML)
  library(stringr)
  library(RCurl)
  name=tolower(name)
  #name='Methyl palmitate'
  name=str_replace_all(name,' ','%20')
  if (str_detect(name,"'")){
    name2=str_replace(name,"'","")
    a0=paste('http://rdf.chemspider.com/search/',"%22",name2,"%22",sep='')
  }else{
    a0=paste('http://rdf.chemspider.com/search/',name,sep='')
  }
  myheader=c('User-Agent'='Mozilla/5.0 (Windows NT 6.1; WOW64)','Accept'='text/html,application/xhyml+xml,application/xml;q=0.9,*/*;q=0.8','Accept-Language'='zh-CN,zh;q=0.8,en;q=0.6','Connection'='keep-alive','Accept-Charset'='GB2312,utf-8;q=0.7,*;q=0.7')
  d2 =debugGatherer()
  cHandle2<- getCurlHandle(httpheader=myheader,followlocation=1,
                           debugfunction=d2$update,verbose=TRUE,
                           cookiefile=cookiefileS)
  a<- getURL(a0,curl=cHandle2,.encoding="UTF-8")
  #a=getURL(a,httpheader=myheader,followlocation=T)
  if (str_detect(a,'(no matches)|(Error)|(Bad Request)')){
    chemspider_id=NA
    smiles=NA
    sure=NA
  }else{
    if(str_detect(a,'ChemSpider ID:')){
      b1=htmlParse(a)
      c1=xpathSApply(b1,"//p/span",fun=xmlValue)
      chemspider_id=c1[[1]]
      smiles=c1[[4]]
      sure='multi'
    }else{
      tempa=str_extract(a,'Chemical-Structure\\.[0-9]+\\.rdf')
      chemspider_id=str_extract(tempa,'[0-9]+')
      b1=htmlParse(a)
      c1=xpathSApply(b1,"//smiles",fun=xmlValue)
      smiles=c1
      sure='ok'
    }
  }
  return(c(SpiderID=chemspider_id,smiles=smiles,sure=sure))
}

getTCMIDchemFromHerb = function(herbName){
  #herbName:herbChineseName
  #TCMIDdatabase:TCMID_data:colnames:c("chemical_name","pub_id","herb_name","herb_CHname")
  require(data.table)
  require(stringr)
  if (!exists('TCMID_data')){
    load('db/TCMIDdatabase.RData')
  }
  data=TCMID_data[herb_CHname%like%herbName,.(chemical_name,pub_id,herb_CHname)]
  if (nrow(data)!=0){
    data$InputName=herbName
  }
  return(data)
}

getTTDdiseaseGene = function(keywords){
  require(data.table)
  require(stringr)
  require(clusterProfiler)
  if (!exists('TTDdiseasedatabase')){
    load('db/TTDdiseasedatabase.RData')
  }
  keywords=tolower(keywords)
  result=TTDdiseasedatabase[TTDdiseasedatabase$disease%like%keywords,]
  if (nrow(result)!=0){
    geneID=unique(result$geneID)
    #genesTran=bitr(geneID, fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Hs.eg.db")
    #genes=genesTran$SYMBOL
	genes=unique(result$geneName)
  }else{
    geneID=''
    genes=''
  }
  return(list(geneNames=genes,geneID=geneID,result=result))
}

getZincID = function(name,cookiefileS=NULL){###chemical_name-->get ZincID
  library(XML)
  library(stringr)
  library(RCurl)
  name=tolower(name)
  name=str_replace_all(name,' ','%20')
  a=paste('http://zinc.docking.org/synonym/',name,sep='')
  myheader=c('User-Agent'='Mozilla/5.0 (Windows NT 6.1; WOW64)','Accept'='text/html,application/xhyml+xml,application/xml;q=0.9,*/*;q=0.8','Accept-Language'='zh-CN,zh;q=0.8,en;q=0.6','Connection'='keep-alive','Accept-Charset'='GB2312,utf-8;q=0.7,*;q=0.7')
  d2 =debugGatherer()
  cHandle2<- getCurlHandle(httpheader=myheader,followlocation=1,
                           debugfunction=d2$update,verbose=TRUE,
                           cookiefile=cookiefileS)
  a<- getURL(a,curl=cHandle2,.encoding="UTF-8")
  #a=getURL(a,httpheader=myheader,followlocation=T)
  b=htmlParse(a)
  c=getNodeSet(b,"//input[@type='checkbox']")
  if (class(try({d=xmlToList(c[[1]])['value']},silent =T))!='try-error'){
    #d=xmlToList(c[[1]])['value']
    chemical_Zid=paste('ZINC',d,sep='')
  }else{
    chemical_Zid=NA
  }
  return(ZincID=chemical_Zid)
}

getZincSmiles = function(zinc_id){
  library(XML)
  library(stringr)
  library(RCurl)
  myheader=c('User-Agent'='Mozilla/5.0 (Windows NT 6.1; WOW64)','Accept'='text/html,application/xhyml+xml,application/xml;q=0.9,*/*;q=0.8','Accept-Language'='zh-CN,zh;q=0.8,en;q=0.6','Connection'='keep-alive','Accept-Charset'='GB2312,utf-8;q=0.7,*;q=0.7')
  name=zinc_id
  a=paste('http://zinc.docking.org/substance/',name,sep='')
  a=htmlParse(a)
  b=getNodeSet(a,"//input[@id='item-smiles']")
  smiles=xmlToList(b[[1]])['value']
  return(smiles)
}

GODOsimilarity = function(genelist1,genelist2){
  require(DOSE)
  DOGOsimilarity=DOSE::clusterSim(genelist1,genelist2)
  return(DOGOsimilarity)
}

GOSimcal = function(genelist1,genelist2,method=c('GO','DO')){
  require(GOSemSim)
  require(DOSE)
  require(stringr)
  #####
  if (str_detect(genelist1,',| |;')){
    genelist1=as.character(unlist(str_split(genelist1,',| |;')))
    genelist1=str_replace_all(genelist1,'"','')
    genelist1=str_replace_all(genelist1,"'",'')
    genelist1=genelist1[genelist1!='']
    genelist1=unique(str_trim(genelist1,side='both'))
  }
  genelist1=as.character(unique(genelist1))
  #####
  if (str_detect(genelist2,',| |;')){
    genelist2=as.character(unlist(str_split(genelist2,',| |;')))
    genelist2=str_replace_all(genelist2,'"','')
    genelist2=str_replace_all(genelist2,"'",'')
    genelist2=genelist2[genelist2!='']
    genelist2=unique(str_trim(genelist2,side='both'))
  }
  genelist2=as.character(unique(genelist2))
  #####
  if ('GO'%in%method){
    resultGO=GOsimilarity(genelist1=genelist1,genelist2=genelist2)
  }else{
    resultGO=NULL
  }
  ######
  if ('DO'%in%method){
    resultDO=GODOsimilarity(genelist1=genelist1,genelist2=genelist2)
  }else{
    resultDO=NULL
  }
  #####
  return(list(GO=resultGO,DO=resultDO))
}

GOsimilarity = function(genelist1,genelist2){
  require(GOSemSim)
  d <- godata('org.Hs.eg.db', ont="BP", computeIC=FALSE)
  resultBP=GOSemSim::clusterSim(genelist1,genelist2,semData=d, measure="Wang")
  d <- godata('org.Hs.eg.db', ont="CC", computeIC=FALSE)
  resultCC=GOSemSim::clusterSim(genelist1,genelist2,semData=d, measure="Wang")
  d <- godata('org.Hs.eg.db', ont="MF", computeIC=FALSE)
  resultMF=GOSemSim::clusterSim(genelist1,genelist2,semData=d, measure="Wang")
  return(data.frame(BP=resultBP,MF=resultMF,CC=resultCC,MeanScore=mean(c(resultBP,resultMF,resultCC),na.rm=T)))
}

KATZ = function(SeedGeneID,ScoreGeneID,randomtimes=1000,testM=1,ppiNetMatrix){
  #SeedGeneID-->疾病靶标
  #ScoreGeneID-->药物靶标
  #testM:1-->wilcox;testM:2-->perm
  require(stringr)
  if (!exists('ppiNetMatrix')){
  load('db/PPIpath.db')
	#load('db/PPIpath.db')
  }
  ###
  if (sum(str_detect(SeedGeneID,',| |;'))>0){
    SeedGeneID=as.character(unlist(str_split(SeedGeneID,',| |;')))
    SeedGeneID=str_replace_all(SeedGeneID,'"','')
    SeedGeneID=str_replace_all(SeedGeneID,"'",'')
    SeedGeneID=SeedGeneID[SeedGeneID!='']
    SeedGeneID=unique(str_trim(SeedGeneID,side='both'))
  }
  SeedGeneID=SeedGeneID[SeedGeneID%in%colnames(ppiNetMatrix)]
  ##
  if (sum(str_detect(ScoreGeneID,',| |;'))>0){
  ScoreGeneID=as.character(unlist(str_split(ScoreGeneID,',| |;')))
  ScoreGeneID=str_replace_all(ScoreGeneID,'"','')
  ScoreGeneID=str_replace_all(ScoreGeneID,"'",'')
  ScoreGeneID=ScoreGeneID[ScoreGeneID!='']
  ScoreGeneID=unique(str_trim(ScoreGeneID,side='both'))
  }
  ScoreGeneID=ScoreGeneID[ScoreGeneID%in%colnames(ppiNetMatrix)]
  ##
  ScoreGeneData=data.frame(SeedGeneIDPPI=toString(SeedGeneID),NumSeedGeneID=length(SeedGeneID),ScoreGeneIDPPI=toString(ScoreGeneID),NumScoreGeneID=length(ScoreGeneID),overlapScore=0,Path1Score=0,Path2Score=0,Path3Score=0,TotalScore=0,RandomScoreMedian=0,Pvalue=1)
  if (length(SeedGeneID)!=0&length(ScoreGeneData)!=0){
  tempGene=SeedGeneID
  overlapScore=length(intersect(tempGene,ScoreGeneID))
  N1Score=sum(ppiNetMatrix[tempGene,ScoreGeneID])*0.001
  N2Score=sum(ppiNetMatrix2[tempGene,ScoreGeneID])*0.000001
  N3Score=sum(ppiNetMatrix3[tempGene,ScoreGeneID])*0.000000001
  ScoreGeneData[1,'overlapScore']=overlapScore
  ScoreGeneData[1,'Path1Score']=N1Score
  ScoreGeneData[1,'Path2Score']=N2Score
  ScoreGeneData[1,'Path3Score']=N3Score
  ScoreGeneData[1,'TotalScore']=overlapScore+N1Score+N2Score+N3Score
  }
  #####################
  #randomScore=numeric()
  randomScore2=numeric()
  seedNum=1234567
  for (i in 1:randomtimes){
    #randomGeneID=sample(colnames(ppiNetMatrix),length(SeedGeneID))
    #overlapScore=length(intersect(randomGeneID,ScoreGeneID))
    #N1Score=sum(ppiNetMatrix[randomGeneID,ScoreGeneID])*0.001
    #N2Score=sum(ppiNetMatrix2[randomGeneID,ScoreGeneID])*0.000001
    #N3Score=sum(ppiNetMatrix3[randomGeneID,ScoreGeneID])*0.000000001
    #randomScore[i]=overlapScore+N1Score+N2Score+N3Score
    #####
    set.seed(seedNum+i)
    randomGeneID2=sample(colnames(ppiNetMatrix),length(ScoreGeneID))
    overlapScore2=length(intersect(randomGeneID2,SeedGeneID))
    N1Score2=sum(ppiNetMatrix[randomGeneID2,SeedGeneID])*0.001
    N2Score2=sum(ppiNetMatrix2[randomGeneID2,SeedGeneID])*0.000001
    N3Score2=sum(ppiNetMatrix3[randomGeneID2,SeedGeneID])*0.000000001
    randomScore2[i]=overlapScore2+N1Score2+N2Score2+N3Score2
  }
  #Pvalue1=sum(randomScore>ScoreGeneData[1,'TotalScore'])/randomtimes
  if (testM==1){
  #Pvalue1=wilcox.test(randomScore,mu=ScoreGeneData[1,'TotalScore'],alternative='less')$p.value
  Pvalue2=wilcox.test(randomScore2,mu=ScoreGeneData[1,'TotalScore'],alternative='less')$p.value
  #ScoreGeneData[1,'Pvalue']=max(Pvalue1,Pvalue2)
  }else{
  Pvalue2=sum(randomScore2>ScoreGeneData[1,'TotalScore'])/randomtimes
  }
  ScoreGeneData[1,'Pvalue']=Pvalue2
  #Pvalue2=sum(randomScore2>ScoreGeneData[1,'TotalScore'])/randomtimes
  ScoreGeneData[1,'RandomScoreMedian']=paste(mean(randomScore2),' [',quantile(randomScore2,0.025),',',quantile(randomScore2,0.975),']',sep='')
  rm('ppiNetMatrix')
  rm('ppiNetMatrix2')
  rm('ppiNetMatrix3')
  gc()
  return(ScoreGeneData)
}

keggSimScore = function(chem_target,keggPath){
  #chem_target:data.frame-->colnames(chem_target)=c('id','geneID')
  library(data.table)
  library(plyr)
  library(stringr)
  library(clipr)
  library(igraph)
  library(ggplot2)
  library(ChemmineR)
  library(pracma)
  #chem_target$geneID=as.character(chem_target$geneID)
  if (!exists('keggPath')){
    #keggPath=fread('db/KEGG_data.csv')
    keggPath=fread('db/KEGG_data.db',encoding='UTF-8')
  }
  #keggPath=fread('KEGG_data160517.csv')
  keggPath$geneID=as.character(keggPath$geneID)
  DrugId=unique(chem_target$id)
  chem_target$geneID=as.character(chem_target$geneID)
  keggSimMatrix=matrix(0,nrow=length(DrugId),ncol=length(DrugId))
  colnames(keggSimMatrix)=DrugId
  rownames(keggSimMatrix)=DrugId
  for (i in 1:length(DrugId)){
    print(i)
    j=i+1
    while (j<=length(DrugId)){
      cid1=DrugId[i]
      cid1Genetemp=chem_target[id==cid1,geneID]
      keggid1=keggPath[geneID%in%cid1Genetemp,unique(PathwayID)]
      cid2=DrugId[j]
      cid2Genetemp=chem_target[id==cid2,geneID]
      keggid2=keggPath[geneID%in%cid2Genetemp,unique(PathwayID)]
      if (length(keggid1)==0|length(keggid2)==0){
        keggSimMatrix[cid1,cid2]=0
      }else{
        keggSimMatrix[cid1,cid2]=length(intersect(keggid1,keggid2))/length(union(keggid1,keggid2))
      }
      j=j+1
    }
  }
  keggSimMatrix=asSymmetric(keggSimMatrix)
  return(keggSimMatrix)
}

logScal = function(SDFfile){
  require(stringr)
  if (file.exists('tempReult.txt')){
    file.remove('tempReult.txt')
  }
  cmd=paste('xlogs',SDFfile,'tempReult.txt',sep=' ')
  shell(cmd,inter=T)
  a=readLines('tempReult.txt')
  bLogSValue=numeric()
  for (i in 1:length(a)){
    temp=a[i]
    b=unlist(str_split(temp,':'))
    b=str_trim(b,side='both')
    if (length(b)>1){
      bName=str_trim(str_replace_all(b[1],'XLOGS of MOL',''))
      if (str_detect(b[2],'WARNING')){
        b[2]=str_trim(str_replace_all(b[2],'\\(WARNING\\)',''))
      }
      bValue=as.numeric(b[2])
      names(bValue)=bName
      bLogSValue=c(bLogSValue,bValue)
    }
  }
  bLogSValueD=data.frame(id=names(bLogSValue),logS=bLogSValue)
  return(bLogSValueD)
}

makeHerb_ChemNet = function(NetworkData,NetworkType='HerbChem',INppi=T){
  library(igraph)
  library(rio)
  library(data.table)
  #library(visNetwork)
  ###NetworkData-->from TCMNPAS outPut[colnames:"herb","cid","chemical_name","geneID","score","inchikey", "smiles","database","QED_DES"]
  ###NetworkType:'HerbChem','ChemGene',or 'HerbChemGene'
  ###list(outputNetwork=outputNetwork,ChemCid_ChemName=ChemCid_ChemName):outputNetwork-->igraphObj;ChemCid_ChemName--dataFrame[colname->('cid','chemical_name')]
  options(stringsAsFactors = F)
  #dd=import('herb-chem-target.csv')
  dd=NetworkData
  dd$geneID=as.character(dd$geneID)
  ChemCid_ChemName=dd[,c('cid','chemical_name')]
  ChemCid_ChemName=as.data.table(ChemCid_ChemName)
  ChemCid_ChemName=unique(ChemCid_ChemName,by='cid')
  ChemCid_ChemName=as.data.frame(ChemCid_ChemName)
  ###
  Herb_Chem=dd[,c('herb','cid')]
  Herb_Chem=unique(Herb_Chem)
  Nodes_Herb_Chem=c(Herb_Chem$herb,Herb_Chem$cid)
  Nodes_Herb_Chem_Group=c(rep('Herb',length(Herb_Chem$herb)),rep('Chemical',length(Herb_Chem$cid)))
  Nodes_Herb_ChemData=data.frame(name=Nodes_Herb_Chem,group=Nodes_Herb_Chem_Group)
  Nodes_Herb_ChemData=unique(Nodes_Herb_ChemData)
  Herb_ChemNet=graph_from_data_frame(d=Herb_Chem, directed = F, vertices = Nodes_Herb_ChemData)#####Herb_ChemNet
  ####################
  Chem_Gene=dd[,c('cid','geneID')]
  #Chem_Gene=as.data.table(Chem_Gene)
  Chem_Gene=unique(Chem_Gene)
  #Chem_Gene$geneID=as.character(Chem_Gene$geneID)
  Nodes_Chem_Gene=c(Chem_Gene$cid,Chem_Gene$geneID)
  Nodes_Chem_Gene_Group=c(rep('Chemical',length(Chem_Gene$cid)),rep('Gene',length(Chem_Gene$geneID)))
  Nodes_Chem_GeneData=data.frame(name=Nodes_Chem_Gene,group=Nodes_Chem_Gene_Group)
  Nodes_Chem_GeneData=unique(Nodes_Chem_GeneData)
  Chem_GeneNet=graph_from_data_frame(d=Chem_Gene, directed = F, vertices = Nodes_Chem_GeneData)####Chem_GeneNet
  #####################
  #Herb_Chem_GeneNet=combineIgraphOB(Herb_ChemNet,Chem_GeneNet)
  Herb_Chem_GeneNet=Herb_ChemNet+Chem_GeneNet
  g1_group=V(Herb_Chem_GeneNet)$group_1
  g2_group=V(Herb_Chem_GeneNet)$group_2
  g_group=g1_group
  g_group[is.na(g1_group)]<-g2_group[is.na(g1_group)]
  V(Herb_Chem_GeneNet)$group<-g_group######Herb_Chem_GeneNet
  Herb_Chem_GeneNet=delete_vertex_attr(Herb_Chem_GeneNet,'group_1')
  Herb_Chem_GeneNet=delete_vertex_attr(Herb_Chem_GeneNet,'group_2')
  if (INppi){
    Chem_GeneNetCore=ExpandNet(geneID=unique(Chem_Gene$geneID),method=0)
    Chem_GeneNetCombinePPI=Chem_GeneNet+Chem_GeneNetCore
    Herb_Chem_GeneNetCombinePPI=Herb_Chem_GeneNet+Chem_GeneNetCore
    Herb_Chem_GeneNet=Herb_Chem_GeneNetCombinePPI###
    Chem_GeneNet=Chem_GeneNetCombinePPI#####
  }
  if (NetworkType=='HerbChem'){
    outputNetwork=Herb_ChemNet
  }else if (NetworkType=='ChemGene'){
    outputNetwork=Chem_GeneNet
  }else{
    outputNetwork=Herb_Chem_GeneNet
  }
  return(list(outputNetwork=outputNetwork,ChemCid_ChemName=ChemCid_ChemName))
}

Mbarplot = function(height, x="Count", colorBy='pvalue', showCategory=5, font.size=12, title="", ...){
  ## use *height* to satisy barplot generic definition
  ## actually here is an enrichResult object.
  object <- height

  colorBy <- match.arg(colorBy, c("pvalue", "p.adjust", "qvalue"))
  if (x == "geneRatio" || x == "GeneRatio") {
    x <- "GeneRatio"
  }
  else if (x == "count" || x == "Count") {
    x <- "Count"
  }

  Description <- Count <- NULL # to satisfy codetools
  df <- fortify_enrichResult(object, showCategory=showCategory, by=x, ...)

  p <- ggplot(df, aes_string(x = "Description", y = x))
  p <- p + geom_bar(stat = "identity") + coord_flip() + theme_dose(font.size)

  if("pvalue" %in% colnames(p$data)) {
    pvalue <- NULL # to satisfy codetools
    p <- p + aes_string(fill=colorBy) +
      scale_fill_continuous(low="red", high="blue")
  } else {
    p <- p+aes(fill=Description) +
      theme(legend.position="none")
  }
  p <- p + ggtitle(title) + xlab("") + ylab("")
  return(p)
}

Mdotplot = function(object, x=~Cluster, colorBy="p.adjust", showCategory=5, by="geneRatio", includeAll=TRUE, font.size=12, title=""){
  df <- fortify__compareClusterResult(object, showCategory=showCategory, by=by, includeAll=includeAll)
  Mplotting(df, x=x, type="dot", colorBy=colorBy, by=by, title=title, font.size=font.size)
}

membershipToList = function(Members){
  require(igraph)
  clusterName=unique(Members)
  Nnodes=length(Members)
  LM=list()
  length(LM)=length(clusterName)
  names(LM)=clusterName
  for (k in clusterName){
    LM[[k]]=names(Members)[which(Members==k)]
  }
  return(LM)
}

MolDesCal = function(SMI){
  #SMI:data.frame-->colnames(SMI)=c('id','smiles')
  #calMolecularDes--->set padel_types.xml!
  library(data.table)
  library(plyr)
  library(stringr)
  library(igraph)
  library(ggplot2)
  library(ChemmineR)
  library(PaDEL)
  library(fingerprint)
  if (!dir.exists('MOL')){
    dir.create('MOL')
  }
  SMI=as.data.table(SMI)
  SMI2=SMI[smiles!='',]
  smiles=SMI2$smiles
  names(smiles)=SMI2$id
  SMISdfSets=ChemmineR::smiles2sdf(smiles)
  SMISdfSets=SMISdfSets[validSDF(SMISdfSets)]
  #sum(validSDF(SMISdfSets))==length(SMISdfSets)
  write.SDF(SMISdfSets,'MOL/SMISdfSets.sdf',cid=T)
  FP=CalDescriptors('MOL/SMISdfSets.sdf')
  return(FP)
}

Mplotting = function(clProf.reshape.df,x = ~Cluster,type = "dot",colorBy = "p.adjust",by = "geneRatio",title="",font.size=12){
  Description <- Percentage <- Count <- Cluster <- GeneRatio <- p.adjust <- pvalue <- NULL # to satisfy codetools
  if (type == "bar") {
    if (by == "percentage") {
      p <- ggplot(clProf.reshape.df,
                  aes(x=Description, y = Percentage, fill=Cluster))
    } else if (by == "count") {
      p <- ggplot(clProf.reshape.df,
                  aes(x=Description, y = Count, fill=Cluster))
    } else {

    }
    p <- p +
      geom_bar() +
      coord_flip()
  }
  if (type == "dot") {
    if (by == "rowPercentage") {
      p <- ggplot(clProf.reshape.df,
                  aes_(x = x, y = ~Description, size = ~Percentage))
    } else if (by == "count") {
      p <- ggplot(clProf.reshape.df,
                  aes_(x = x, y = ~Description, size = ~Count))
    } else if (by == "geneRatio") {
      p <- ggplot(clProf.reshape.df,
                  aes_(x = x, y = ~Description, size = ~GeneRatio))
    } else {
      ## nothing here
    }
    if (any(colnames(clProf.reshape.df) == colorBy)) {
      p <- p +
        geom_point() +
        aes_string(color=colorBy) +
        scale_colour_gradient(low="red", high="blue")
    } else {
      p <- p + geom_point(colour="steelblue")
    }
  }
  p <- p + xlab("") + ylab("") + ggtitle(title) +
    theme_dose(font.size)
  ## theme(axis.text.x = element_text(colour="black", size=font.size, vjust = 1)) +
  ##     theme(axis.text.y = element_text(colour="black",
  ##           size=font.size, hjust = 1)) +
  ##               ggtitle(title)+theme_bw()
  ## p <- p + theme(axis.text.x = element_text(angle=angle.axis.x,
  ##                    hjust=hjust.axis.x,
  ##                    vjust=vjust.axis.x))
  return(p)
}

MyModulatity = function(g,LisMembership){
  ##g-->igraph object,Must have node name;LisMembership-->list,names(list):name of cluster
  require(igraph)
  require(clValid)
  Ncluster=length(LisMembership)
  mem=unique(unlist(LisMembership))
  if (length(mem)==vcount(g)&(sum(mem%in%V(g)$name)==length(mem))){
    membership_Mat=annotationListToMatrix(LisMembership,V(g)$name)
    Co_membership=membership_Mat%*%t(membership_Mat)
    adj=as_adj(g)
    D=degree(g)
    D1=replicate(length(D),as.matrix(D),simplify ='matrix' )
    D2=replicate(length(D),as.matrix(D),simplify ='matrix' )
    D2=t(D2)
    DD=D1*D2
    EC=ecount(g)
    MQ=(adj-DD/(2*EC))*Co_membership
    Q=sum(MQ)/(2*EC)
    return(Q)
  }else{
    g_sub=induced_subgraph(g,mem)
    membership_Mat=annotationListToMatrix(LisMembership,V(g_sub)$name)
    Co_membership=membership_Mat%*%t(membership_Mat)
    adj=as_adj(g_sub)
    D=degree(g_sub)
    D1=replicate(length(D),as.matrix(D),simplify ='matrix' )
    D2=replicate(length(D),as.matrix(D),simplify ='matrix' )
    D2=t(D2)
    DD=D1*D2
    EC=ecount(g_sub)
    MQ=(adj-DD/(2*EC))*Co_membership
    Q=sum(MQ)/(2*EC)
    return(Q)
  }
}

myscale = function(x){
  minx=min(x)
  maxx=max(x)
  y=9*(x-minx)/(maxx-minx)+1
}

NetCorr = function(SeedGeneID,ScoreGeneID,Methodnet,randomtimes=1000,gamma=0.7,testM=1){
  if (Methodnet=='RW'){
    resultRW=RWanalysis(SeedGeneID=SeedGeneID,ScoreGeneID=ScoreGeneID,randomtimes=randomtimes,gamma=gamma,testM=testM)
    resultRWTable=data.frame(SeedGeneIDPPI=toString(resultRW$SeedGeneID),NumSeedGeneID=resultRW$NumSeedGeneIDPPI,ScoreGeneIDPPI=toString(resultRW$ScoreGeneID),NumScoreGeneID=resultRW$NumScoreGeneIDPPI,TotalScore=resultRW$TotalScore,RandomScoreMean=resultRW$RandomScoreMedian,Pvalue=resultRW$Pvalue)
  }else{
    resultRW=NULL
    resultRWTable=NULL
  }
  #####
  if (Methodnet=='KATZ'){
    resultKATZ=KATZ(SeedGeneID,ScoreGeneID,randomtimes=randomtimes,testM=testM)
  }else{
    resultKATZ=NULL
  }
  #####
  if (Methodnet=='DIS'){
    resultDis=DisAccess(SeedGeneID,ScoreGeneID,randomtimes=randomtimes,testM=testM)
    resultDisTable=data.frame(SeedGeneIDPPI=toString(resultDis$SeedGeneID),NumSeedGeneID=resultDis$NumSeedGeneID,ScoreGeneIDPPI=toString(resultDis$ScoreGeneID),NumScoreGeneID=resultDis$NumScoreGeneID,MeanDistance=resultDis$MedianDistance,RandomScoreMean=resultDis$RandomScoreMedian,Pvalue=resultDis$Pvalue)
  }else{
    resultDis=NULL
    resultDisTable=NULL
  }
  #####
return(list(resultRW=resultRW,resultRWTable=resultRWTable,resultKATZ=resultKATZ,resultDis=resultDis,resultDisTable=resultDisTable))
}

NetStatCal = function(NetList){
  #NetList:list-->names(list)=NetName(igraphObject)
  library(data.table)
  library(plyr)
  library(stringr)
  library(igraph)
  library("org.Hs.eg.db")
  metric=c('Mean_degree','Density','Mean_betweenness','Count_nodes','Cluster_coefficient','Diameter','Shortest_path')
  #NetList=list(PPI=ppiBinaryNet,BXXXT=BXXXT_Net,Colitis=Yan_Net,Diabet=Diabet_Net,Cancer=Cancer_Net)
  statistics_matrix=matrix(0,nrow=length(metric),ncol=length(NetList))
  colnames(statistics_matrix)=names(NetList)
  rownames(statistics_matrix)=metric
  ###
  statistics_matrix['Mean_degree',]=sapply(NetList,function(x)paste(round(mean(igraph::degree(x)),2),round(sd(igraph::degree(x)),2),sep= '±'))
  statistics_matrix['Mean_betweenness',]=sapply(NetList,function(x)paste(round(mean(igraph::betweenness(x,directed=F)),2),round(sd(igraph::betweenness(x,directed=F)),2),sep= '±'))
  statistics_matrix['Density',]=sapply(NetList,igraph::edge_density)
  statistics_matrix['Count_nodes',]=sapply(NetList,igraph::vcount)
  statistics_matrix['Cluster_coefficient',]=sapply(NetList,igraph::transitivity,type='global')
  statistics_matrix['Diameter',]=sapply(NetList,igraph::diameter,directed=F)
  statistics_matrix['Shortest_path',]=sapply(NetList,igraph::mean_distance,directed=F)
  statistics_matrix=as.data.frame(statistics_matrix)
  rownames(statistics_matrix)=metric
  return(statistics_matrix)
}

NodeOverlapAnalysis = function(NetList,Plot=TRUE){
  #NetList:list-->names(list)=NetName(igraphObject)
  library(data.table)
  library(plyr)
  library(stringr)
  library(igraph)
  library(corrplot)
  library("org.Hs.eg.db")
  #NetList=list(BXXXT=BXXXT_kpath_Net,Colitis=Yan_kpath_Net,Diabet=Diabet_kpath_Net,Cancer=Cancer_kpath_Net)
  overlap_DISmatrix=matrix(0,nrow=length(NetList),ncol=length(NetList))
  colnames(overlap_DISmatrix)=names(NetList)
  rownames(overlap_DISmatrix)=names(NetList)
  for (i in names(NetList)){
    for (j in names(NetList)){
      tempNet_i=NetList[[i]]
      tempNet_j=NetList[[j]]
      overlap_DISmatrix[i,j]=length(intersect(V(tempNet_i)$name,V(tempNet_j)$name))/min(c(vcount(tempNet_i),vcount(tempNet_j)))
    }
  }
  diag(overlap_DISmatrix)=0
  colnames(overlap_DISmatrix)=names(NetList)
  rownames(overlap_DISmatrix)=names(NetList)
  if (Plot){
    png('NetNodeOverlap_matrix.png',width=1280,height=1024)
    corrplot(overlap_DISmatrix,method='pie',diag=F,cl.lim=c(0,1),tl.cex=3,cl.cex=3)
    dev.off()
  }
  return(overlap_DISmatrix)
}

P_analysis = function(prescription_basic,yao,mine_data,tS0.9=0.03,min_yao=8,threshold_leastYao=10,onlyHerb=FALSE){
  ##P_analysis(prescription_basic,yao,mine_data,TopDisease=5,TopDrug=10,tS0.9=0.03,min_yao=8,threshold_leastYao=10,onlyHerb=FALSE,Drug_pair_analysis=FALSE,local=TRUE)
  ##prescription_basic-->"id","time", "patient","diagnosis", "doctor","","drug"[binary]
  ##mine_data--->"num","drug_id","mine_id","mine_Cname", "drug_commonName", "mine_name" ,"class","fufang", "jixing","guige""生产企业"         "规格", "jixing2" , "innComponentName","type"合并了<--herb_data
  ##yao-->drug向量包含西药及成药，由basic_trans输出得到
  ##TopDisease-->统计前几位的疾病
  ##TopDrug--->统计关联度最强的前几位的药物
  ##min_yao-->基本方发现时需要指定的最小含药数
  ##threshold_leastYao-->自适应优化二值化阈值时采用的预期含药数
  ##herb_name-->由herb_data得到unique(herb_data$mine_name+herb_data$innComponentName)
  ##z_name-->表示成药向量
  ##x_name-->表示西药向量
  require(plyr)
  require(igraph)
  require(data.table)
  require(Matrix)
  prescription_basic=as.data.frame(prescription_basic)
  herb_name=unique(mine_data[class=='h',mine_Cname])
  ##########################################
  if (onlyHerb){
    noherb=yao[!yao%in%herb_name]
    Pyao=yao[yao%in%herb_name]
    P=subset(prescription_basic,select=Pyao)
    P=as.matrix(P)
    P=Matrix(P,sparse=T)
    rownames(P)=prescription_basic$id
    # if (Drug_pair_analysis){
    #   Drug_pair_result=drug_pair(P)
    # }else{
    #   Drug_pair_result=na_table
    # }
    #####
    BF_result=BF_analysis(P,tS0.9=tS0.9,min_yao=min_yao,threshold_leastYao=threshold_leastYao)
    BF=BF_result$BF_result####BF基本方信息
    }else{
    P=subset(prescription_basic,select=yao)
    P=as.matrix(P)
    P=Matrix(P,sparse=T)
    rownames(P)=prescription_basic$id
    # if (Drug_pair_analysis){
    #   Drug_pair_result=drug_pair(P)
    # }else{
    #   Drug_pair_result=na_table
    # }
    #####
    BF_result=BF_analysis(P,tS0.9=tS0.9,min_yao=min_yao,threshold_leastYao=threshold_leastYao)
    BF=BF_result$BF_result####BF基本方信息

  }
  ###########################################

  #############################################
  #############################################
  #return(list(BF=BF,herb_fre=herb_fre,opt_a=BF_result$opt_a,Degree_yao=BF_result$Degree_yao,Drug_pair_result=Drug_pair_result,noherb_fre=noherb_fre,disease_fre=disease_fre,P_stat=P_stat,disease_yao=disease_yao))
  return(list(BF=BF,opt_a=BF_result$opt_a,Degree_yao=BF_result$Degree_yao))
}

PathWayScore = function(geneID1,geneID1Name='Name1',geneID2,geneID2Name='Name2',method=c('WNS','keggSim')){
  library(data.table)
  library(plyr)
  library(stringr)
  ################################
  if (str_detect(geneID1,',|;| ')){
    geneID1=unlist(str_split(geneID1,',|;| '))
    geneID1=str_replace_all(geneID1,'"','')
    geneID1=str_replace_all(geneID1,"'",'')
    geneID1=geneID1[geneID1!='']
    geneID1=str_trim(geneID1,side='both')
  }
  geneID1=unique(geneID1)
  Data1=data.table(id=geneID1Name,geneID=geneID1)
  ##############################
  if (str_detect(geneID2,',|;| ')){
    geneID2=unlist(str_split(geneID2,',|;| '))
    geneID2=str_replace_all(geneID2,'"','')
    geneID2=str_replace_all(geneID2,"'",'')
    geneID2=geneID2[geneID2!='']
    geneID2=str_trim(geneID2,side='both')
  }
  geneID2=unique(geneID2)
  Data2=data.table(id=geneID2Name,geneID=geneID2)
  ##############################
  Data=rbind(Data1,Data2)
  #############################
  if (length(method)==2){
    result1=WNSscore(Data)
    result2=keggSimScore(Data)
    result=data.frame(geneID1Name=geneID1Name,geneID2Name=geneID2Name,WNSscore=result1[geneID1Name,geneID2Name],keggSimscore=result2[geneID1Name,geneID2Name])
  }else if ('WNS'%in%method){
    result1=WNSscore(Data)
    result=data.frame(geneID1Name=geneID1Name,geneID2Name=geneID2Name,WNSscore=result1[geneID1Name,geneID2Name])
  }else if ('keggSim'%in%method){
    result2=keggSimScore(Data)
    result=data.frame(geneID1Name=geneID1Name,geneID2Name=geneID2Name,keggSimscore=result2[geneID1Name,geneID2Name])
  }
  return(result)
}

pdftotxt = function(pdffile){
  output=str_replace(pdffile,'\\.pdf','.txt')
  dos=paste('pdftotext',pdffile,'-enc UTF-8')
  shell(dos)
  #Sys.sleep(1)
  return(output)
}

plotCommonItems = function(DataFrame,Ca,Cb,TopNIDCurve=30,pointSize=2,legendTextSize=8,axisTitleSize=8,axisTextSize=6,plotName='Commonplot.png',plotWidth=10,plotHeight=8,xLab='Top N enriched-Item',yLab='Percentage of common enriched-Items'){
library(ggplot2)
library(pracma)
library(data.table)
options(stringsAsFactors = F)
tempD=DataFrame
tempD$Cluster=as.character(tempD$Cluster)
tempD$ID=as.character(tempD$ID)
item=unique(tempD$Cluster)
tempD=as.data.table(tempD)
tempD=tempD[order(Cluster,qvalue)]
herbItem=tempD$ID[tempD$Cluster==Ca]
diseaseItem=tempD$ID[tempD$Cluster==Cb]
CommonItem=intersect(herbItem,diseaseItem)
if (length(diseaseItem)>0){
  TopNIDCurve=min(TopNIDCurve,length(herbItem))
  common_CommonItem=numeric()
  for (k in 1:TopNIDCurve){
    tempC=herbItem[1:k]
    common_CommonItem[k]=length(intersect(tempC,diseaseItem))/k
  }
  common_data=data.frame(ID=1:TopNIDCurve,Disease=common_CommonItem)
  common_data=as.data.table(common_data)
  common_data=melt(common_data,id='ID',variable='Disease',value='Percentage')
  common_data=as.data.frame(common_data)
  Commonauc=trapz(1:TopNIDCurve,common_CommonItem)
  PlotCommon=ggplot(common_data)+aes(x=ID,y=Percentage,group=Disease,colour=Disease,shape=Disease,fill=Disease)+geom_line()+geom_point(size=pointSize)+xlab(xLab)+ylab(yLab)+theme(legend.text=element_text(size=legendTextSize),legend.title=element_text(size=legendTextSize))+theme(axis.title.x=element_text(size=axisTitleSize),axis.title.y=element_text(size=axisTitleSize))+theme(axis.text.x=element_text(size=axisTextSize),axis.text.y=element_text(size=axisTextSize))+theme(legend.position = 'none')
ggsave(filename=plotName,plot=PlotCommon,width=plotWidth,height=plotHeight,dpi=600,units='in')
  }else{
  PlotCommon=NULL
  Commonauc=NULL
}
return(list(PlotObject=PlotCommon,AUC=Commonauc))
}

QEDcal = function(druglikeDesData){
  #colnames(druglikeDesData)=c('id',colnames(QED_par))
  library(ChemmineR)
  library(ChemmineOB)
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

QEDcalALL = function(druglikeDesFile,sdfFile){
  library(ChemmineR)
  library(ChemmineOB)
  require(data.table)
  require(plyr)
  require(rio)
  options(stringsAsFactors = F)
#druglikeDesFile直接由DS计算得到的csv，colnames(druglikeDesFile)=c('id',colnames(QED_par)[1:7]);##sdfFile是DS计算的文件有id字段！
  #druglikeData=read.csv(druglikeDesFile,head=T,stringsAsFactors=F,colClasses=c('character',rep('numeric',7)))
  druglikeData=fread(druglikeDesFile)
  druglikeData=as.data.frame(druglikeData)
  RenamesChar=c(ALogP='ALOGP',Molecular_Mass='MW',Num_RotatableBonds='ROTB',Num_AromaticRings='AROM',Num_H_Acceptors_Lipinski='HBA',Num_H_Donors_Lipinski='HBD',Molecular_PolarSurfaceArea='PSA')
  druglikeData=plyr::rename(druglikeData,RenamesChar)
  druglikeData$id=as.character(druglikeData$id)
  druglikeData$ALOGP=as.numeric(druglikeData$ALOGP)
  druglikeData$MW=as.numeric(druglikeData$MW)
  druglikeData$ROTB=as.numeric(druglikeData$ROTB)
  druglikeData$AROM=as.numeric(druglikeData$AROM)
  druglikeData$HBA=as.numeric(druglikeData$HBA)
  druglikeData$HBD=as.numeric(druglikeData$HBD)
  druglikeData$PSA=as.numeric(druglikeData$PSA)
#ALogP-->ALOGP
#Molecular_Mass-->MW
#Num_RotatableBonds-->ROTB
#Num_AromaticRings-->AROM
#Num_H_Acceptors_Lipinski-->HBA
#Num_H_Donors_Lipinski-->HBD
#Molecular_PolarSurfaceArea-->PSA
  alert=alertCal(sdfFile)
  druglikeData=join(druglikeData,alert,by='id')
  druglike_des=QEDcal(druglikeData)
  return(druglike_des)
}

randomWalk = function(igraphM,queryGenes,EdgeWeight = FALSE, gamma = 0.7){
  require(PAGI)
  require(igraph)
  ##igraphM-->igraph should has vertex name!
  ##queryGenes-->seedGeneNames_vector,vertex name
  ##return-score--->Named Vector
  #g=make_ring(10)
  #V(g)$name <- letters[1:10]
  vg=numeric(length(V(igraphM)))
  names(vg)=V(igraphM)$name
  index=queryGenes
  vg[index]=1
  score=RandomWalk2igraph(igraphM,vg,EdgeWeight,gamma)
  names(score)=names(vg)
  return(score)
}

read_pdf = function(pdffile,sep=FALSE){
  library(stringr)
  output=str_replace(pdffile,'\\.pdf','.txt')
  dos=paste('pdftotext',pdffile,'-enc UTF-8')
  shell(dos)
  Sys.sleep(1)
  txt=scan(output,what='character',encoding = 'UTF-8')
  if (sep){
    txt2=txt
  }else{
    txt2=paste(txt,collapse = '')
  }
  return(txt2)
  }

replaceDatabase = function(x){
  require(stringr)
  temp=str_replace_all(x,'c\\(|\\)|"','')
  return(temp)
}

ReportGeneListNP = function(GeneListFile=NULL,GeneList=NULL,GeneListName='CUSTOM',diseaseGeneID=NULL,GOMFenrich=T,GOBPenrich=T,GOCCenrich=T,KEGGenrich=T,PAenrich=T,DOenrich=T,MESHenrich=T,meshcategory='C',qvalueCutoff=0.05,TopNIDCurve=10,TOPNkeggPlot,showCategory,PLOTfontsize,PLOTwidth,PLOTheight,Savefilename=NULL){
  require(r2excel)
  require(xlsx)
  require(data.table)
  require(pathview)
  require(stringr)
  if (!is.null(GeneListFile)){
    GeneList=readLines(GeneListFile)
  }
  ######
  if (sum(str_detect(GeneList,',|;| '))){
    GeneList=unlist(str_split(GeneList,',|;| '))
    GeneList=str_replace_all(GeneList,'"','')
    GeneList=str_replace_all(GeneList,"'",'')
    GeneList=GeneList[GeneList!='']
    GeneList=str_trim(GeneList,side='both')
  }
  GeneList=unique(GeneList)
  ######
  ######
  TCMNPobj=GeneListNP(GeneList=GeneList,GeneListName=GeneListName,diseaseGeneID=diseaseGeneID,GOMFenrich=GOMFenrich,GOBPenrich=GOBPenrich,GOCCenrich=GOCCenrich,KEGGenrich=KEGGenrich,PAenrich=PAenrich,DOenrich=DOenrich,MESHenrich=MESHenrich,meshcategory=meshcategory,qvalueCutoff=qvalueCutoff,TopNIDCurve=TopNIDCurve)
  wb <- createWorkbook(type="xlsx")
  ##################
  if (GOMFenrich&(!is.null(TCMNPobj$enrich))){
    sheetGOMF=createSheet(wb, sheetName = "基因本体分子功能富集")
    xlsx.addTable(wb, sheetGOMF,TCMNPobj$enrich$GOmf$Enrichs)
    sheetGOMFPIC=createSheet(wb, sheetName = "基因本体分子功能富集PIC")
    plotFunction<-function(){
      p=Mdotplot(TCMNPobj$enrich$GOmf$OBject,showCategory=showCategory,font.size=PLOTfontsize)
      print(p)
    }
    xlsx.addPlot(wb, sheetGOMFPIC, plotFunction=plotFunction,width=PLOTwidth,height=PLOTheight)
  }
  ###
  if (GOBPenrich&(!is.null(TCMNPobj$enrich))){
    sheetGOBP=createSheet(wb, sheetName = "基因本体生物过程富集")
    xlsx.addTable(wb, sheetGOBP,TCMNPobj$enrich$GObp$Enrichs)
    sheetGOBPPIC=createSheet(wb, sheetName = "基因本体生物过程富集PIC")
    plotFunction<-function(){
      p=Mdotplot(TCMNPobj$enrich$GObp$OBject,showCategory=showCategory,font.size=PLOTfontsize)
      print(p)
    }
    xlsx.addPlot(wb, sheetGOBPPIC, plotFunction=plotFunction,width=PLOTwidth,height=PLOTheight)
  }
  ###
  if (GOCCenrich&(!is.null(TCMNPobj$enrich))){
    sheetGOCC=createSheet(wb, sheetName = "基因本体细胞组件富集")
    xlsx.addTable(wb, sheetGOCC,TCMNPobj$enrich$GOcc$Enrichs)
    sheetGOCCPIC=createSheet(wb, sheetName = "基因本体细胞组件富集PIC")
    plotFunction=function(){
      p=Mdotplot(TCMNPobj$enrich$GOcc$OBject,showCategory=showCategory,font.size=PLOTfontsize)
      print(p)
    }
    xlsx.addPlot(wb, sheetGOCCPIC, plotFunction=plotFunction,width=PLOTwidth,height=PLOTheight)
  }
  ####
  if ((GOMFenrich|GOBPenrich|GOCCenrich)&(!is.null(diseaseGeneID))){
    sheetGODisease=createSheet(wb, sheetName = "复方与疾病基因共关联GO")
    txt1=paste('复方与疾病基因共关联GO相似度：',TCMNPobj$Disease$CommonDiseaseGOSim,sep='')
    txt2=paste('复方与疾病基因共关联前',TopNIDCurve,'个GO本体曲线下面积AUC：',TCMNPobj$Disease$CommonGOauc,sep='')
    xlsx.addParagraph(wb, sheetGODisease, txt1, fontSize=12, isItalic=TRUE,
                      fontColor="darkred", backGroundColor="gray",startRow=1,
                      colSpan=10, rowSpan=1)
    xlsx.addParagraph(wb, sheetGODisease, txt2, fontSize=12, isItalic=TRUE,
                      fontColor="darkred", backGroundColor="gray",startRow=2,
                      colSpan=10, rowSpan=1)
    plotFunction<-function(){
      print(TCMNPobj$Disease$PlotCommonGO)
    }
    xlsx.addPlot(wb, sheetGODisease, plotFunction=plotFunction,startRow=5)
  }
  ####
  ###
  if (KEGGenrich&(!is.null(TCMNPobj$enrich))){
    sheetKEGG=createSheet(wb, sheetName = "KEGG信号通路富集")
    xlsx.addTable(wb, sheetKEGG,TCMNPobj$enrich$KEGG$Enrichs)
    sheetKEGGPIC=createSheet(wb, sheetName = "KEGG信号通路富集PIC")
    plotFunction<-function(){
      p=Mdotplot(TCMNPobj$enrich$KEGG$OBject,showCategory=showCategory,font.size=PLOTfontsize)
      print(p)
    }
    xlsx.addPlot(wb, sheetKEGGPIC, plotFunction=plotFunction,width=PLOTwidth,height=PLOTheight)
    ##
    if (TOPNkeggPlot!=0){
      TOPNkeggPlot=as.numeric(TOPNkeggPlot)
      pathWayID=unique(TCMNPobj$enrich$KEGG$Enrichs$ID)
      TOPNkeggPlot=min(TOPNkeggPlot,length(pathWayID))
      pathWayIDtoPlot=pathWayID[1:TOPNkeggPlot]
      NumPath=length(pathWayIDtoPlot)
      pathview(gene.data=TCMNPobj$selectGene$GeneID,pathway.id=pathWayIDtoPlot,species="hsa",out.suffix='PA')
      for (i in 1:NumPath){
        tempSheetName=pathWayIDtoPlot[i]
        sheetTemp=createSheet(wb, sheetName =tempSheetName)
        picfile=paste(tempSheetName,'.PA.png',sep='')
        addPicture(picfile,sheetTemp,scale=5)
        file.remove(picfile)
        file.remove(paste(tempSheetName,'.png',sep=''))
        file.remove(paste(tempSheetName,'.xml',sep=''))
      }
    }
  }
  ########
  if (PAenrich&(!is.null(TCMNPobj$enrich))){
    sheetPA=createSheet(wb, sheetName = "Reactome信号通路富集")
    xlsx.addTable(wb, sheetPA,TCMNPobj$enrich$PA$Enrichs)
    sheetPAPIC=createSheet(wb, sheetName = "Reactome信号通路富集PIC")
    plotFunction<-function(){
      p=Mdotplot(TCMNPobj$enrich$PA$OBject,showCategory=showCategory,font.size=PLOTfontsize)
      print(p)
    }
    xlsx.addPlot(wb, sheetPAPIC, plotFunction=plotFunction,width=PLOTwidth,height=PLOTheight)
    ##
  }
  ########
  ##
  if ((KEGGenrich)&(!is.null(diseaseGeneID))){
    sheetKEGGDisease=createSheet(wb, sheetName = "复方与疾病基因共关联KEGG通路")
    txt1=paste('复方与疾病基因共关联KEGG信号通路相似度：',TCMNPobj$Disease$CommonDiseaseKEGGSim,sep='')
    txt2=paste('复方与疾病基因共关联前',TopNIDCurve,'个KEGG信号通路曲线下面积AUC：',TCMNPobj$Disease$CommonKEGGauc,sep='')
    xlsx.addParagraph(wb, sheetKEGGDisease, txt1, fontSize=12, isItalic=TRUE,
                      fontColor="darkred", backGroundColor="gray",startRow=1,
                      colSpan=10, rowSpan=1)
    xlsx.addParagraph(wb, sheetKEGGDisease, txt2, fontSize=12, isItalic=TRUE,
                      fontColor="darkred", backGroundColor="gray",startRow=2,
                      colSpan=10, rowSpan=1)
    plotFunction<-function(){
      print(TCMNPobj$Disease$PlotCommonKEGG)
    }
    xlsx.addPlot(wb, sheetKEGGDisease, plotFunction=plotFunction,startRow=5)
  }
  ###########
  if ((PAenrich)&(!is.null(diseaseGeneID))){
    sheetPADisease=createSheet(wb, sheetName = "复方与疾病基因共关联Reactome通路")
    txt1=paste('复方与疾病基因共关联Reactome信号通路相似度：',TCMNPobj$Disease$CommonDiseasePASim,sep='')
    txt2=paste('复方与疾病基因共关联前',TopNIDCurve,'个Reactome信号通路曲线下面积AUC：',TCMNPobj$Disease$CommonPAauc,sep='')
    xlsx.addParagraph(wb, sheetPADisease, txt1, fontSize=12, isItalic=TRUE,
                      fontColor="darkred", backGroundColor="gray",startRow=1,
                      colSpan=10, rowSpan=1)
    xlsx.addParagraph(wb, sheetPADisease, txt2, fontSize=12, isItalic=TRUE,
                      fontColor="darkred", backGroundColor="gray",startRow=2,
                      colSpan=10, rowSpan=1)
    plotFunction<-function(){
      print(TCMNPobj$Disease$PlotCommonPA)
    }
    xlsx.addPlot(wb, sheetPADisease, plotFunction=plotFunction,startRow=5)
  }
  ###########
  if (DOenrich&(!is.null(TCMNPobj$enrich))){
    sheetDO=createSheet(wb, sheetName = "疾病本体富集")
    xlsx.addTable(wb, sheetDO,TCMNPobj$enrich$DO$Enrichs)
    sheetDOPIC=createSheet(wb, sheetName = "疾病本体富集PIC")
    plotFunction<-function(){
      p=Mdotplot(TCMNPobj$enrich$DO$OBject,showCategory=showCategory,font.size=PLOTfontsize)
      print(p)
    }
    xlsx.addPlot(wb, sheetDOPIC, plotFunction=plotFunction,width=PLOTwidth,height=PLOTheight)
  }
  ###
  if (MESHenrich&(!is.null(TCMNPobj$enrich))){
    sheetMESH=createSheet(wb, sheetName = "主题词富集")
    xlsx.addTable(wb, sheetMESH,TCMNPobj$enrich$MESH$Enrichs)
    sheetMESHPIC=createSheet(wb, sheetName = "主题词富集PIC")
    plotFunction<-function(){
      p=Mbarplot(TCMNPobj$enrich$MESH$OBject,showCategory=showCategory,font.size=PLOTfontsize)
      print(p)
    }
    xlsx.addPlot(wb, sheetMESHPIC, plotFunction=plotFunction,width=PLOTwidth,height=PLOTheight)
  }
  #############################
  if (!dir.exists('Report'))dir.create('Report')
  if (is.null(Savefilename))Savefilename='指定靶标群的机制分析'
  saveWorkbook(wb, paste('Report/',Savefilename,'.xlsx',sep=''))
}

ReportTCMNP = function(CustomherbTargetDat=NULL,herbList,herbListName='TCM',QEDset=0.2,geneScore=800,targetDatabase=c('HIT','TCMID','STITCH','TCMSP','CUSTOM'),diseaseGeneID=NULL,geneSelectPv=0.01,GOMFenrich=T,GOBPenrich=T,GOCCenrich=T,KEGGenrich=T,PAenrich=T,DOenrich=T,MESHenrich=T,meshcategory='C',qvalueCutoff=0.05,TOPNkeggPlot=1,showCategory=5,font.size=12,width=1024,height=950,TopNIDCurve=10,Savefilename){
  require(r2excel)
  require(xlsx)
  require(data.table)
  require(pathview)
  if (!is.null(CustomherbTargetDat)){
    CustomherbTargetDat=fread(CustomherbTargetDat,encoding='UTF-8')
  }
  TCMNPobj=TCMNP(CustomherbTargetDat=CustomherbTargetDat,herbList=herbList,herbListName=herbListName,QEDset=QEDset,geneScore=geneScore,targetDatabase=targetDatabase,diseaseGeneID=diseaseGeneID,geneSelectPv=geneSelectPv,GOMFenrich=GOMFenrich,GOBPenrich=GOBPenrich,GOCCenrich=GOCCenrich,KEGGenrich=KEGGenrich,PAenrich=PAenrich,DOenrich=DOenrich,MESHenrich=MESHenrich,meshcategory=meshcategory,qvalueCutoff=qvalueCutoff,TopNIDCurve=TopNIDCurve)
  wb <- createWorkbook(type="xlsx")
  sheet1 <- createSheet(wb, sheetName = "复方靶标筛选结果")
  sheet2 <- createSheet(wb, sheetName = "复方活性化合物筛选结果")
  sheet3 <- createSheet(wb, sheetName = "药物-成分-靶标检索结果")
  ##################
  xlsx.addHeader(wb, sheet1, " 复方靶标筛选结果", level=2, underline=1)
  xlsx.addLineBreak(sheet1, 2)
  paragraph1=paste('输入的药物列表：',toString(herbList),'||检索到的药物列表：',toString(TCMNPobj$actHerb),'||检索的靶标数据库：',toString(targetDatabase),'||设定的基因靶标关联阈值：',geneScore,'||设定的基因靶标显著性阈值：',geneSelectPv)
  xlsx.addParagraph(wb, sheet1, paragraph1, fontSize=12, isItalic=TRUE,
                    fontColor="darkred", backGroundColor="gray",
                    colSpan=10, rowSpan=3)
  xlsx.addTable(wb, sheet1,TCMNPobj$selectGene,startRow=8)
  #################
  xlsx.addHeader(wb, sheet2, " 复方活性化合物筛选结果", level=2, underline=1)
  xlsx.addLineBreak(sheet2, 2)
  paragraph2=paste('输入的药物列表：',toString(herbList),'||检索到的药物列表：',toString(TCMNPobj$actHerb),'||检索的靶标数据库：',toString(targetDatabase),'||设定的化合物类药性指数阈值：',QEDset,'||设定的基因靶标显著性阈值：',geneSelectPv)
  xlsx.addParagraph(wb, sheet2, paragraph2, fontSize=12, isItalic=TRUE,
                    fontColor="darkred", backGroundColor="gray",
                    colSpan=10, rowSpan=3)
  xlsx.addTable(wb, sheet2,TCMNPobj$selectChem,startRow=8)
  #################
  xlsx.addHeader(wb, sheet3, " 药物-成分-靶标检索结果", level=2, underline=1)
  xlsx.addLineBreak(sheet3, 2)
  paragraph3=paste('输入的药物列表：',toString(herbList),'||检索到的药物列表：',toString(TCMNPobj$actHerb),'||检索的靶标数据库：',toString(targetDatabase),'||设定的基因靶标关联阈值：',geneScore,'||设定的化合物类药性指数阈值：',QEDset,'||设定的基因靶标显著性阈值：',geneSelectPv)
  xlsx.addParagraph(wb, sheet3, paragraph3, fontSize=12, isItalic=TRUE,
                    fontColor="darkred", backGroundColor="gray",
                    colSpan=10, rowSpan=3)
  xlsx.addTable(wb, sheet3,TCMNPobj$HerbTarget,startRow=8)
  ##################
  if (GOMFenrich&(!is.null(TCMNPobj$enrich))){
    sheetGOMF=createSheet(wb, sheetName = "基因本体分子功能富集")
    xlsx.addTable(wb, sheetGOMF,TCMNPobj$enrich$GOmf$Enrichs)
    sheetGOMFPIC=createSheet(wb, sheetName = "基因本体分子功能富集PIC")
    plotFunction<-function(){
      p=Mdotplot(TCMNPobj$enrich$GOmf$OBject,showCategory=showCategory,font.size=font.size)
      print(p)
      }
    xlsx.addPlot(wb, sheetGOMFPIC, plotFunction=plotFunction,width=width,height=height)
  }
  ###
  if (GOBPenrich&(!is.null(TCMNPobj$enrich))){
    sheetGOBP=createSheet(wb, sheetName = "基因本体生物过程富集")
    xlsx.addTable(wb, sheetGOBP,TCMNPobj$enrich$GObp$Enrichs)
    sheetGOBPPIC=createSheet(wb, sheetName = "基因本体生物过程富集PIC")
    plotFunction<-function(){
      p=Mdotplot(TCMNPobj$enrich$GObp$OBject,showCategory=showCategory,font.size=font.size)
      print(p)
      }
    xlsx.addPlot(wb, sheetGOBPPIC, plotFunction=plotFunction,width=width,height=height)
  }
  ###
  if (GOCCenrich&(!is.null(TCMNPobj$enrich))){
    sheetGOCC=createSheet(wb, sheetName = "基因本体细胞组件富集")
    xlsx.addTable(wb, sheetGOCC,TCMNPobj$enrich$GOcc$Enrichs)
    sheetGOCCPIC=createSheet(wb, sheetName = "基因本体细胞组件富集PIC")
    plotFunction=function(){
      p=Mdotplot(TCMNPobj$enrich$GOcc$OBject,showCategory=showCategory,font.size=font.size)
      print(p)
      }
    xlsx.addPlot(wb, sheetGOCCPIC, plotFunction=plotFunction,width=width,height=height)
  }
  ####
  if ((GOMFenrich|GOBPenrich|GOCCenrich)&(!is.null(diseaseGeneID))){
    sheetGODisease=createSheet(wb, sheetName = "复方与疾病基因共关联GO")
    txt1=paste('复方与疾病基因共关联GO相似度：',TCMNPobj$Disease$CommonDiseaseGOSim,sep='')
    txt2=paste('复方与疾病基因共关联前',TopNIDCurve,'个GO本体曲线下面积AUC：',TCMNPobj$Disease$CommonGOauc,sep='')
    xlsx.addParagraph(wb, sheetGODisease, txt1, fontSize=12, isItalic=TRUE,
                      fontColor="darkred", backGroundColor="gray",startRow=1,
                      colSpan=10, rowSpan=1)
    xlsx.addParagraph(wb, sheetGODisease, txt2, fontSize=12, isItalic=TRUE,
                      fontColor="darkred", backGroundColor="gray",startRow=2,
                      colSpan=10, rowSpan=1)
    plotFunction<-function(){
      print(TCMNPobj$Disease$PlotCommonGO)
      }
    xlsx.addPlot(wb, sheetGODisease, plotFunction=plotFunction,startRow=5)
  }
  ####
  ###
  if (KEGGenrich&(!is.null(TCMNPobj$enrich))){
    sheetKEGG=createSheet(wb, sheetName = "KEGG信号通路富集")
    xlsx.addTable(wb, sheetKEGG,TCMNPobj$enrich$KEGG$Enrichs)
    sheetKEGGPIC=createSheet(wb, sheetName = "KEGG信号通路富集PIC")
    plotFunction<-function(){
      p=Mdotplot(TCMNPobj$enrich$KEGG$OBject,showCategory=showCategory,font.size=font.size)
      print(p)
      }
    xlsx.addPlot(wb, sheetKEGGPIC, plotFunction=plotFunction,width=width,height=height)
    ##
    if (TOPNkeggPlot!=0){
      TOPNkeggPlot=as.numeric(TOPNkeggPlot)
      pathWayID=unique(TCMNPobj$enrich$KEGG$Enrichs$ID)
      TOPNkeggPlot=min(TOPNkeggPlot,length(pathWayID))
      pathWayIDtoPlot=pathWayID[1:TOPNkeggPlot]
      NumPath=length(pathWayIDtoPlot)
      pathview(gene.data=TCMNPobj$selectGene$GeneID,pathway.id=pathWayIDtoPlot,species="hsa",out.suffix='PA')
      for (i in 1:NumPath){
        tempSheetName=pathWayIDtoPlot[i]
        sheetTemp=createSheet(wb, sheetName =tempSheetName)
        picfile=paste(tempSheetName,'.PA.png',sep='')
        addPicture(picfile,sheetTemp,scale=5)
        file.remove(picfile)
        file.remove(paste(tempSheetName,'.png',sep=''))
        file.remove(paste(tempSheetName,'.xml',sep=''))
      }
    }
  }
  ########
  if (PAenrich&(!is.null(TCMNPobj$enrich))){
    sheetPA=createSheet(wb, sheetName = "Reactome信号通路富集")
    xlsx.addTable(wb, sheetPA,TCMNPobj$enrich$PA$Enrichs)
    sheetPAPIC=createSheet(wb, sheetName = "Reactome信号通路富集PIC")
    plotFunction<-function(){
      p=Mdotplot(TCMNPobj$enrich$PA$OBject,showCategory=showCategory,font.size=font.size)
      print(p)
    }
    xlsx.addPlot(wb, sheetPAPIC, plotFunction=plotFunction,width=width,height=height)
    ##
  }
  ########
  ##
  if ((KEGGenrich)&(!is.null(diseaseGeneID))){
    sheetKEGGDisease=createSheet(wb, sheetName = "复方与疾病基因共关联KEGG通路")
    txt1=paste('复方与疾病基因共关联KEGG信号通路相似度：',TCMNPobj$Disease$CommonDiseaseKEGGSim,sep='')
    txt2=paste('复方与疾病基因共关联前',TopNIDCurve,'个KEGG信号通路曲线下面积AUC：',TCMNPobj$Disease$CommonKEGGauc,sep='')
    xlsx.addParagraph(wb, sheetKEGGDisease, txt1, fontSize=12, isItalic=TRUE,
                      fontColor="darkred", backGroundColor="gray",startRow=1,
                      colSpan=10, rowSpan=1)
    xlsx.addParagraph(wb, sheetKEGGDisease, txt2, fontSize=12, isItalic=TRUE,
                      fontColor="darkred", backGroundColor="gray",startRow=2,
                      colSpan=10, rowSpan=1)
    plotFunction<-function(){
      print(TCMNPobj$Disease$PlotCommonKEGG)
      }
    xlsx.addPlot(wb, sheetKEGGDisease, plotFunction=plotFunction,startRow=5)
  }
  ###########
  if ((PAenrich)&(!is.null(diseaseGeneID))){
    sheetPADisease=createSheet(wb, sheetName = "复方与疾病基因共关联Reactome通路")
    txt1=paste('复方与疾病基因共关联Reactome信号通路相似度：',TCMNPobj$Disease$CommonDiseasePASim,sep='')
    txt2=paste('复方与疾病基因共关联前',TopNIDCurve,'个Reactome信号通路曲线下面积AUC：',TCMNPobj$Disease$CommonPAauc,sep='')
    xlsx.addParagraph(wb, sheetPADisease, txt1, fontSize=12, isItalic=TRUE,
                      fontColor="darkred", backGroundColor="gray",startRow=1,
                      colSpan=10, rowSpan=1)
    xlsx.addParagraph(wb, sheetPADisease, txt2, fontSize=12, isItalic=TRUE,
                      fontColor="darkred", backGroundColor="gray",startRow=2,
                      colSpan=10, rowSpan=1)
    plotFunction<-function(){
      print(TCMNPobj$Disease$PlotCommonPA)
    }
    xlsx.addPlot(wb, sheetPADisease, plotFunction=plotFunction,startRow=5)
  }
  ###########
  if (DOenrich&(!is.null(TCMNPobj$enrich))){
    sheetDO=createSheet(wb, sheetName = "疾病本体富集")
    xlsx.addTable(wb, sheetDO,TCMNPobj$enrich$DO$Enrichs)
    sheetDOPIC=createSheet(wb, sheetName = "疾病本体富集PIC")
    plotFunction<-function(){
      p=Mdotplot(TCMNPobj$enrich$DO$OBject,showCategory=showCategory,font.size=font.size)
      print(p)
      }
    xlsx.addPlot(wb, sheetDOPIC, plotFunction=plotFunction,width=width,height=height)
  }
  ###
  if (MESHenrich&(!is.null(TCMNPobj$enrich))){
    sheetMESH=createSheet(wb, sheetName = "主题词富集")
    xlsx.addTable(wb, sheetMESH,TCMNPobj$enrich$MESH$Enrichs)
    sheetMESHPIC=createSheet(wb, sheetName = "主题词富集PIC")
    plotFunction<-function(){
      p=Mbarplot(TCMNPobj$enrich$MESH$OBject,showCategory=showCategory,font.size=font.size)
      print(p)
      }
    xlsx.addPlot(wb, sheetMESHPIC, plotFunction=plotFunction,width=width,height=height)
  }
  #############################
  if (!dir.exists('Report'))dir.create('Report')
  saveWorkbook(wb, paste('Report/',Savefilename,'.xlsx',sep=''))
}

RWanalysis = function(SeedGeneID,ScoreGeneID,scoreStand=1000,randomtimes=1000,EdgeWeight=FALSE,gamma=0.7,TopNCal=FALSE,TopN=1:10,Plot=FALSE,testM=1,ppiBinaryNet=ppiBinaryNet){
  #SeedGeneID-->seedGeneID
  #ScoreGeneID-->待评估的geneID
  #scoreStand-->score的扩大倍数
  #TopNCal-->是否要计算前N个gene，使用此项功能前ScoreGeneID必须按重要性降序排列
  #ppiBinaryNet
  #testM:1-->wilcox;testM:2-->perm
  require(data.table)
  require(plyr)
  require(stringr)
  require(igraph)
  require("org.Hs.eg.db")
  require(pracma)
  if (!exists('ppiBinaryNet')){
    load('db/ppiNetData.db')
	#load('db/ppiNetData.db')
  }
  ##
  if (sum(str_detect(SeedGeneID,',| |;'))>0){
    SeedGeneID=as.character(unlist(str_split(SeedGeneID,',| |;')))
    SeedGeneID=str_replace_all(SeedGeneID,'"','')
    SeedGeneID=str_replace_all(SeedGeneID,"'",'')
    SeedGeneID=SeedGeneID[SeedGeneID!='']
    SeedGeneID=unique(str_trim(SeedGeneID,side='both'))
  }
  SeedGeneID=as.character(unique(SeedGeneID))
  SeedGeneID=SeedGeneID[SeedGeneID%in%V(ppiBinaryNet)$name]
  ##
  if (sum(str_detect(ScoreGeneID,',| |;'))>0){
    ScoreGeneID=as.character(unlist(str_split(ScoreGeneID,',| |;')))
    ScoreGeneID=str_replace_all(ScoreGeneID,'"','')
    ScoreGeneID=str_replace_all(ScoreGeneID,"'",'')
    ScoreGeneID=ScoreGeneID[ScoreGeneID!='']
    ScoreGeneID=unique(str_trim(ScoreGeneID,side='both'))
  }
  ScoreGeneID=as.character(unique(ScoreGeneID))
  ScoreGeneID=ScoreGeneID[ScoreGeneID%in%V(ppiBinaryNet)$name]
  ##
  if (!TopNCal){
    RWR_Seed=randomWalk(ppiBinaryNet,SeedGeneID,EdgeWeight=EdgeWeight,gamma=gamma)*scoreStand
    TopNN=1:length(ScoreGeneID)
    RWR_Score=numeric()
    for (i in TopNN){
      tempN=TopNN[i]
      RWR_Score[i]=sum(RWR_Seed[ScoreGeneID[1:tempN]],na.rm=T)/tempN
      }
    #######
    TopNNauc=trapz(TopNN,RWR_Score)
    GeneScore=RWR_Seed[ScoreGeneID]
    TotalScore=sum(RWR_Seed[ScoreGeneID])*length(unique(SeedGeneID))
    if (Plot){
      RWRscore=data.frame(TopN=TopNN,geneScore=RWR_Score)
      #RWRscoreMelt=melt(RWRscore,id='TopN')
      colnames(RWRscore)=c('TopN','Score')
      ggplot(RWRscore)+aes(x=TopN,y=Score)+geom_point()+geom_line()+xlab('TopN Genes')+ylab('RWR Score')
      ggsave('RWRscore.png')
    }

  }else{
    RWR_Seed=randomWalk(ppiBinaryNet,SeedGeneID,EdgeWeight=EdgeWeight,gamma=gamma)*scoreStand
    #TopNN=TopN
    RWR_Score=numeric()
    for (i in 1:length(TopN)){
      tempN=TopN[i]
      RWR_Score[i]=sum(RWR_Seed[ScoreGeneID[1:tempN]],na.rm=T)/tempN
    }
    #######
    TopNNauc=trapz(TopN,RWR_Score)
    GeneScore=RWR_Seed[ScoreGeneID]
    TotalScore=sum(RWR_Seed[ScoreGeneID])*length(unique(SeedGeneID))
    if (Plot){
      RWRscore=data.frame(TopN=TopN,geneScore=RWR_Score)
      #RWRscoreMelt=melt(RWRscore,id='TopN')
      colnames(RWRscore)=c('TopN','Score')
      ggplot(RWRscore)+aes(x=TopN,y=Score)+geom_point()+geom_line()+xlab('TopN Genes')+ylab('RWR Score')
      ggsave('RWRscore.png')
    }
  }
  ###################
  tempRandomScore=numeric()
  #tempRandomScore2=numeric()
  #RWR_ScoregeneID=randomWalk(ppiBinaryNet,ScoreGeneID,EdgeWeight=EdgeWeight,gamma=gamma)*scoreStand
  #TotalScore_ScoreID=sum(RWR_ScoregeneID[SeedGeneID])*length(unique(ScoreGeneID))
  seedNum=123456
  for (j in 1:randomtimes){
    set.seed(seedNum+j)
    tempRandomSeed1=sample(V(ppiBinaryNet)$name,length(ScoreGeneID))
    #tempRandomSeed2=sample(V(ppiBinaryNet)$name,length(SeedGeneID))
    tempRandomScore[j]=sum(RWR_Seed[tempRandomSeed1])*length(unique(SeedGeneID))
    #tempRandomScore2[i]=sum(RWR_ScoregeneID[tempRandomSeed2])*length(unique(ScoreGeneID))
  }
  #tempRandomScore2[is.na(tempRandomScore2)]=0
  tempRandomScore[is.na(tempRandomScore)]=0
  if (testM==1){
  Pvalue1=wilcox.test(tempRandomScore,mu=TotalScore,alternative='less')$p.value
  }else{
  Pvalue1=sum(tempRandomScore>TotalScore)/randomtimes  
  }
  #Pvalue1=sum(tempRandomScore>TotalScore)/randomtimes  
  #Pvalue2=sum(tempRandomScore2>TotalScore_ScoreID)/randomtimes
  #Pvalue2=wilcox.test(tempRandomScore2,mu=TotalScore,alternative='less')$p.value
  #Pvalue=max(Pvalue1,Pvalue2)
  Pvalue=Pvalue1
  rm('ppiBinaryNet')
  gc()  
  return(list(SeedGeneID=SeedGeneID,NumSeedGeneIDPPI=length(SeedGeneID),ScoreGeneID=ScoreGeneID,NumScoreGeneIDPPI=length(ScoreGeneID),GeneScore=GeneScore,TopNauc=TopNNauc,TotalScore=TotalScore,Pvalue=Pvalue,RandomScoreMedian=paste(mean(tempRandomScore),' [',quantile(tempRandomScore,0.025),',',quantile(tempRandomScore,0.975),']',sep='')))
}

selectFun = function(protein){
  require(UniProt.ws)
  if (!exists('up')){
    load('db/UniProtwsOB.RData')
  }
  columns='ENTREZ_GENE'
  kt='ENSEMBL_PROTEIN'
  if (class(try({b=select(up,protein,columns,kt)},silent = T))!='try-error'){
    a=b
    return(a$ENTREZ_GENE[1])
  }else{
    a='NA'
    return(a)
  }
  #a=select(up,protein,columns,kt)

}

sim = function(x,y){
  ##x,y--->two vectors not continurous!
  ##output-->list(result=result,sim_score=sim_score)
  common_item=intersect(x,y)
  sim_score=length(common_item)/min(length(x),length(y))
  sim_score=round(sim_score,digits = 3)
  if (length(common_item)==0){
    result_txt='non-overlapping'
    sim_score=0
    result=result_txt
  }else{
    result_txt=paste(common_item,collapse = ',')
    result=paste('S=',sim_score,'->',result_txt,sep='')
  }
  return(list(result=result,sim_score=sim_score))
}

SmiToSDF = function(SMI,type='C'){
  #SMI='SMI.csv'-->colnames=c('id','smiles') or SMItext-->Named character vector
  #type='F'-->file for SMI or 'C'-->character for SMItext
  #outPut:'Molelular.sdf'
  library(data.table)
  library(plyr)
  library(stringr)
  library(clipr)
  library(igraph)
  library(readr)
  library(ggplot2)
  library(ChemmineR)
  if (type=='F'){
    SMI1=fread(SMI)
    SMI2=SMI1[smiles!='',]
    smiles=SMI2$smiles
    names(smiles)=SMI2$id
    SMISdfSets=ChemmineR::smiles2sdf(smiles)
    SMISdfSets=SMISdfSets[validSDF(SMISdfSets)]
    #sum(validSDF(SMISdfSets))==length(SMISdfSets)
    datablock(SMISdfSets)<-data.frame(id=SMISdfSets@ID)
    write.SDF(SMISdfSets,'Molelular.sdf',cid=T)
  }else if (type=='C'){
    SMISdfSets=ChemmineR::smiles2sdf(SMI)
    SMISdfSets=SMISdfSets[validSDF(SMISdfSets)]
    #sum(validSDF(SMISdfSets))==length(SMISdfSets)
    datablock(SMISdfSets)<-data.frame(id=SMISdfSets@ID)
    write.SDF(SMISdfSets,'Molelular.sdf',cid=T)
  }else{
    print('ERROR:type not found!')
  }
}

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

TCMNP = function(CustomherbTargetDat=NULL,herbList,herbListName='TCM',QEDset=0.2,geneScore=400,targetDatabase=c('HIT','TCMID','STITCH','TCMSP','CUSTOM'),diseaseGeneID=NULL,geneSelectPv=0.05,GOMFenrich=T,GOBPenrich=T,GOCCenrich=T,KEGGenrich=T,PAenrich=T,DOenrich=T,MESHenrich=T,meshcategory='C',qvalueCutoff=0.05,TopNIDCurve=10,AllHerbListData=AllHerbListData,TarDB='db/Tar.db'){
  #diseaseGeneID:data.frame('class','geneID') or character
  ##CustomherbTargetDat:colname-->c(herb,cid,chemical_name,geneID,score,inchikey,smiles,database,QED_DES)
  require(data.table)
  require(stringr)
  require(RSQLite)
  require(plyr)
  require(ggplot2)
  require(pracma)
  ##############################################
  if (!is.null(diseaseGeneID)){
    if (class(diseaseGeneID)%in%c('data.table','data.frame')){
      diseaseGeneID$geneID=as.character(diseaseGeneID$geneID)
      DGeneID=diseaseGeneID$geneID
      names(DGeneID)=diseaseGeneID$class
      DiseaseName=unique(names(DGeneID))
      DiseaseNum=length(DiseaseName)
    }else{
      if (is.null(names(diseaseGeneID))){
        diseaseGeneID=as.character(unlist(str_split(diseaseGeneID,',| |;')))
        diseaseGeneID=str_replace_all(diseaseGeneID,'"','')
        diseaseGeneID=str_replace_all(diseaseGeneID,"'",'')
        diseaseGeneID=diseaseGeneID[diseaseGeneID!='']
        diseaseGeneID=unique(str_trim(diseaseGeneID,side='both'))
        DGeneID=diseaseGeneID
        names(DGeneID)=rep('Disease',length(diseaseGeneID))
        DiseaseName='Disease'
        DiseaseNum=1
      }else{
        diseaseGeneID=as.character(unlist(str_split(diseaseGeneID,',| |;')))
        diseaseGeneID=str_replace_all(diseaseGeneID,'"','')
        diseaseGeneID=str_replace_all(diseaseGeneID,"'",'')
        diseaseGeneID=diseaseGeneID[diseaseGeneID!='']
        diseaseGeneID=unique(str_trim(diseaseGeneID,side='both'))
        DGeneID=diseaseGeneID
        DiseaseName=unique(names(DGeneID))
        DiseaseNum=length(DiseaseName)
      }
    }
  }
  #############################################################################################################################自定义输入靶标
  if (!is.null(CustomherbTargetDat)){
    standColname=c('herb','cid','chemical_name','geneID','score','inchikey','smiles','database','QED_DES')
    colNameCustom=colnames(CustomherbTargetDat)
    if (sum(colNameCustom%in%standColname)==9){
      subD=CustomherbTargetDat
    }else{
      OutcolNameCustom=standColname[!standColname%in%colNameCustom]
      OutSideCustomData=data.frame(matrix(NA,ncol=length(OutcolNameCustom),nrow=nrow(CustomherbTargetDat)))
      colnames(OutSideCustomData)=OutcolNameCustom
      if ('score'%in%OutcolNameCustom)OutSideCustomData$score=9999
      if ('database'%in%OutcolNameCustom)OutSideCustomData$database='Custom'
      if ('QED_DES'%in%OutcolNameCustom)OutSideCustomData$QED_DES=1
      TempsubD=cbind(CustomherbTargetDat,OutSideCustomData)
      subD=TempsubD[,c(standColname,setdiff(colNameCustom,standColname)),with=F]
    }
    #####################################
    actHerb=unique(subD$herb)
    subD=subD[QED_DES>=QEDset&score>=geneScore,]
    subD$geneID=as.character(subD$geneID)
    tempR=geneScoreCal(unique(subD[,.(herb,id=cid,geneID=geneID)]),Pcutoff = geneSelectPv)
    selectGeneID=tempR$selectGeneId###selectGeneID
    if (length(selectGeneID)!=0){
      selectGeneName=geneTrans(selectGeneID,type='ID')
      tempR2=chemScoreCal(unique(subD[,.(id=cid,geneID=geneID)]),geneScore = tempR$geneScore,selectGeneId=selectGeneID)
      selectCID=tempR2$selectChemId#####selectCID
      subDD=unique(subD[,.(cid,chemical_name)])
      setkey(subDD,cid)
      selectCName=subDD[selectCID,chemical_name,mult="first"]
      subD2=subD[cid%in%selectCID&geneID%in%selectGeneID,]###subD2
      if (is.null(diseaseGeneID)){
        Order=paste('geneCluster=list(',herbListName,'=selectGeneID)',sep='')
      }else{
        if (DiseaseNum==1){
          Order=paste('geneCluster=list(',herbListName,'=selectGeneID,',DiseaseName,'=DGeneID[names(DGeneID)==DiseaseName])',sep='')
        }
        ##
        if (DiseaseNum==2){
          Order=paste('geneCluster=list(',herbListName,'=selectGeneID,',DiseaseName[1],'=DGeneID[names(DGeneID)==DiseaseName[1]],',DiseaseName[2],'=DGeneID[names(DGeneID)==DiseaseName[2]])',sep='')
        }
        if (DiseaseNum==3){
          Order=paste('geneCluster=list(',herbListName,'=selectGeneID,',DiseaseName[1],'=DGeneID[names(DGeneID)==DiseaseName[1]],',DiseaseName[2],'=DGeneID[names(DGeneID)==DiseaseName[2]],',DiseaseName[3],'=DGeneID[names(DGeneID)==DiseaseName[3]])',sep='')
        }
      }
      eval(parse(text=Order))
      ###################################
      if (KEGGenrich){
        if (class(try({EnrichKEGG=EnrichMentAnalysis(geneCluster,type='KEGG',qvalueCutoff=qvalueCutoff)},silent=T))=='try-error')EnrichKEGG=NULL
      }else{
        EnrichKEGG=NULL
      }
      #################################
      if (PAenrich){
        if (class(try({EnrichPA=EnrichMentAnalysis(geneCluster,type='ReactomePA',qvalueCutoff=qvalueCutoff)},silent=T))=='try-error')EnrichPA=NULL
      }else{
        EnrichPA=NULL
      }
      ############
      if (GOBPenrich){
        if (class(try({enrichGObp=EnrichMentAnalysis(geneCluster,type='GO',ontGO='BP',qvalueCutoff=qvalueCutoff)},silent=T))=='try-error')enrichGObp=NULL
      }else{
        enrichGObp=NULL
      }
      #############
      if (GOMFenrich){
        if (class(try({enrichGOmf=EnrichMentAnalysis(geneCluster,type='GO',ontGO='MF',qvalueCutoff=qvalueCutoff)},silent=T))=='try-error')enrichGOmf=NULL
      }else{
        enrichGOmf=NULL
      }
      ############
      if (GOCCenrich){
        if (class(try({enrichGOcc=EnrichMentAnalysis(geneCluster,type='GO',ontGO='CC',qvalueCutoff=qvalueCutoff)},silent=T))=='try-error')enrichGOcc=NULL
      }else{
        enrichGOcc=NULL
      }
      ############
      if (DOenrich){
        data(EG2DOLite)
        data(DOLiteTerm)
        data(DOLite2EG)
        if (class(try({enrichmentDO=EnrichMentAnalysis(geneCluster,type='DO',qvalueCutoff=qvalueCutoff)},silent=T))=='try-error')enrichmentDO=NULL
      }else{
        enrichmentDO=NULL
      }
      ############
      if (MESHenrich){
        if (class(try({enrichMESH=EnrichMentAnalysis(geneCluster,type='MESH',meshcategory=meshcategory,qvalueCutoff=qvalueCutoff)},silent=T))=='try-error')enrichMESH=NULL
      }else{
        enrichMESH=NULL
      }
      ##############################################################################
      if (!is.null(diseaseGeneID)&!is.null(EnrichKEGG)&(herbListName%in%EnrichKEGG$Enrichs$Cluster)){
        tempD=EnrichKEGG$Enrichs
        tempD$Cluster=as.character(tempD$Cluster)
        tempD$ID=as.character(tempD$ID)
        item=unique(tempD$Cluster)
        diseaseItem=item[!item%in%herbListName]
        CommonDiseaseKEGGSim=numeric()
        herbListNameItem=tempD$ID[tempD$Cluster==herbListName]
        if (length(diseaseItem)>0){
          length(CommonDiseaseKEGGSim)=length(diseaseItem)
          names(CommonDiseaseKEGGSim)=diseaseItem
          for (i in 1:length(diseaseItem)){
            diseaseItemtemp=tempD$ID[tempD$Cluster==diseaseItem[i]]
            CommonDiseaseKEGGSim[i]=length(intersect(herbListNameItem,diseaseItemtemp))/min(length(herbListNameItem),length(diseaseItemtemp))
          }
          diseaseToCurve=diseaseItem[1]
          diseasetNameItem=tempD$ID[tempD$Cluster==diseaseToCurve]
          TopNIDCurve=min(TopNIDCurve,length(herbListNameItem))
          common_KeggDisaeae=numeric()
          for (k in 1:TopNIDCurve){
            tempC=herbListNameItem[1:k]
            common_KeggDisaeae[k]=length(intersect(tempC,diseasetNameItem))/k
          }
          common_Keggdata=data.frame(ID=1:TopNIDCurve,Disease=common_KeggDisaeae)
          common_Keggdata=as.data.table(common_Keggdata)
          common_Keggdata=melt(common_Keggdata,id='ID',variable='Disease',value='Percentage')
          common_Keggdata=as.data.frame(common_Keggdata)
          CommonKEGGauc=trapz(1:TopNIDCurve,common_KeggDisaeae)
          PlotCommonKEGG=ggplot(common_Keggdata)+aes(x=ID,y=Percentage,group=Disease,colour=Disease,shape=Disease,fill=Disease)+geom_line()+geom_point(size=2)+xlab(paste('Top N enriched-Kegg pathway in ',herbListName))+ylab('Percentage of common enriched-Kegg pathway')+theme(legend.text=element_text(size=10),legend.title=element_text(size=10))+theme(axis.title.x=element_text(size=8),axis.title.y=element_text(size=8))+theme(axis.text.x=element_text(size=6),axis.text.y=element_text(size=6))
        }
        else{
          CommonDiseaseKEGGSim=NULL
          CommonKEGGauc=NULL
          PlotCommonKEGG=NULL
        }
        ######
      }else{
        CommonDiseaseKEGGSim=NULL
        CommonKEGGauc=NULL
        PlotCommonKEGG=NULL
      }
      ###############################################################################
      ###########ReactomePA#############################
      ##############################################################################
      if (!is.null(diseaseGeneID)&!is.null(EnrichPA)&(herbListName%in%EnrichPA$Enrichs$Cluster)){
        tempD=EnrichPA$Enrichs
        tempD$Cluster=as.character(tempD$Cluster)
        tempD$ID=as.character(tempD$ID)
        item=unique(tempD$Cluster)
        diseaseItem=item[!item%in%herbListName]
        CommonDiseasePASim=numeric()
        herbListNameItem=tempD$ID[tempD$Cluster==herbListName]
        if (length(diseaseItem)>0){
          length(CommonDiseasePASim)=length(diseaseItem)
          names(CommonDiseasePASim)=diseaseItem
          for (i in 1:length(diseaseItem)){
            diseaseItemtemp=tempD$ID[tempD$Cluster==diseaseItem[i]]
            CommonDiseasePASim[i]=length(intersect(herbListNameItem,diseaseItemtemp))/min(length(herbListNameItem),length(diseaseItemtemp))
          }
          diseaseToCurve=diseaseItem[1]
          diseasetNameItem=tempD$ID[tempD$Cluster==diseaseToCurve]
          TopNIDCurve=min(TopNIDCurve,length(herbListNameItem))
          common_PADisaeae=numeric()
          for (k in 1:TopNIDCurve){
            tempC=herbListNameItem[1:k]
            common_PADisaeae[k]=length(intersect(tempC,diseasetNameItem))/k
          }
          common_PAdata=data.frame(ID=1:TopNIDCurve,Disease=common_PADisaeae)
          common_PAdata=as.data.table(common_PAdata)
          common_PAdata=melt(common_PAdata,id='ID',variable='Disease',value='Percentage')
          common_PAdata=as.data.frame(common_PAdata)
          CommonPAauc=trapz(1:TopNIDCurve,common_PADisaeae)
          PlotCommonPA=ggplot(common_PAdata)+aes(x=ID,y=Percentage,group=Disease,colour=Disease,shape=Disease,fill=Disease)+geom_line()+geom_point(size=2)+xlab(paste('Top N enriched-Reactome pathway in ',herbListName))+ylab('Percentage of common enriched-Reactome pathway')+theme(legend.text=element_text(size=10),legend.title=element_text(size=10))+theme(axis.title.x=element_text(size=8),axis.title.y=element_text(size=8))+theme(axis.text.x=element_text(size=6),axis.text.y=element_text(size=6))
        }
        else{
          CommonDiseasePASim=NULL
          CommonPAauc=NULL
          PlotCommonPA=NULL
        }
        ######
      }else{
        CommonDiseasePASim=NULL
        CommonPAauc=NULL
        PlotCommonPA=NULL
      }
      ################PA END
      ################################################################################
      BlankGO=(!is.null(enrichGObp))|(!is.null(enrichGOcc))|(!is.null(enrichGOmf))
      BlankherbListName=(herbListName%in%enrichGObp$Enrichs$Cluster)|(herbListName%in%enrichGOcc$Enrichs$Cluster)|(herbListName%in%enrichGOmf$Enrichs$Cluster)
      if (!is.null(diseaseGeneID)&BlankGO&BlankherbListName){
        tempD=rbind(enrichGObp$Enrichs,enrichGOcc$Enrichs,enrichGOmf$Enrichs)
        tempD=arrange(tempD,qvalue)
        tempD$Cluster=as.character(tempD$Cluster)
        tempD$ID=as.character(tempD$ID)
        item=unique(tempD$Cluster)
        diseaseItem=item[!item%in%herbListName]
        CommonDiseaseGOSim=numeric()
        herbListNameItem=unique(tempD$ID[tempD$Cluster==herbListName])
        if (length(diseaseItem)>0){
          length(CommonDiseaseGOSim)=length(diseaseItem)
          names(CommonDiseaseGOSim)=diseaseItem
          for (i in 1:length(diseaseItem)){
            diseaseItemtemp=unique(tempD$ID[tempD$Cluster==diseaseItem[i]])
            CommonDiseaseGOSim[i]=length(intersect(herbListNameItem,diseaseItemtemp))/min(length(herbListNameItem),length(diseaseItemtemp))
          }
          diseaseToCurve=diseaseItem[1]
          diseasetNameItem=unique(tempD$ID[tempD$Cluster==diseaseToCurve])
          TopNIDCurve=min(TopNIDCurve,length(herbListNameItem))
          common_GODisaeae=numeric()
          for (k in 1:TopNIDCurve){
            tempC=herbListNameItem[1:k]
            common_GODisaeae[k]=length(intersect(tempC,diseasetNameItem))/k
          }
          common_GOdata=data.frame(ID=1:TopNIDCurve,Disease=common_GODisaeae)
          common_GOdata=as.data.table(common_GOdata)
          common_GOdata=melt(common_GOdata,id='ID',variable='Disease',value='Percentage')
          common_GOdata=as.data.frame(common_GOdata)
          CommonGOauc=trapz(1:TopNIDCurve,common_GODisaeae)
          PlotCommonGO=ggplot(common_GOdata)+aes(x=ID,y=Percentage,group=Disease,colour=Disease,shape=Disease,fill=Disease)+geom_line()+geom_point(size=2)+xlab(paste('Top N enriched-GO in ',herbListName))+ylab('Percentage of common enriched-GO')+theme(legend.text=element_text(size=10),legend.title=element_text(size=10))+theme(axis.title.x=element_text(size=8),axis.title.y=element_text(size=8))+theme(axis.text.x=element_text(size=6),axis.text.y=element_text(size=6))
        }else{
          CommonDiseaseGOSim=NULL
          CommonGOauc=NULL
          PlotCommonGO=NULL
        }
        ##############
      }else{
        CommonDiseaseGOSim=NULL
        CommonGOauc=NULL
        PlotCommonGO=NULL
      }
      
      #result=list(selectGene=data.table(GeneID=selectGeneID,GeneName=selectGeneName),selectChem=data.table(CID=selectCID,ChemName=selectCName),HerbTarget=subD2,actHerb=actHerb,enrich=list(KEGG=EnrichKEGG,PA=EnrichPA,GObp=enrichGObp,GOmf=enrichGOmf,GOcc=enrichGOcc,DO=enrichmentDO,MESH=enrichMESH),Disease=list(CommonDiseaseKEGGSim=CommonDiseaseKEGGSim,CommonKEGGauc=CommonKEGGauc,PlotCommonKEGG=PlotCommonKEGG,CommonDiseasePASim=CommonDiseasePASim,CommonPAauc=CommonPAauc,PlotCommonPA=PlotCommonPA,CommonDiseaseGOSim=CommonDiseaseGOSim,CommonGOauc=CommonGOauc,PlotCommonGO=PlotCommonGO),geneScore=tempR,chemScore=tempR2)
	  result=list(selectGene=data.table(GeneID=selectGeneID,GeneName=selectGeneName,geneScore=tempR$geneScore),selectChem=data.table(CID=selectCID,ChemName=selectCName,chemScore=tempR2$chemScore[tempR2$chemScore>0]),HerbTarget=subD2,actHerb=actHerb,enrich=list(KEGG=EnrichKEGG,PA=EnrichPA,GObp=enrichGObp,GOmf=enrichGOmf,GOcc=enrichGOcc,DO=enrichmentDO,MESH=enrichMESH),Disease=list(CommonDiseaseKEGGSim=CommonDiseaseKEGGSim,CommonKEGGauc=CommonKEGGauc,PlotCommonKEGG=PlotCommonKEGG,CommonDiseasePASim=CommonDiseasePASim,CommonPAauc=CommonPAauc,PlotCommonPA=PlotCommonPA,CommonDiseaseGOSim=CommonDiseaseGOSim,CommonGOauc=CommonGOauc,PlotCommonGO=PlotCommonGO),geneScore=tempR,chemScore=tempR2)
    }else{
      result=list(selectGene=data.table(GeneID=character(),GeneName=character()),selectChem=data.table(CID=character(),ChemName=character()),HerbTarget=data.table(),actHerb=actHerb,enrich=NULL,Disease=NULL)
    }
    ###############
    
    #####################################
  }else{###非自定义输入靶标
    #########################################################################################################
    if (!exists('AllHerbListData')){
      #ALL=fread('db/AllHerbListData.csv')
      AllHerbListData=fread('db/AllHerbListData.csv',encoding='UTF-8')
      #load('db/AllHerbListData.RData')
    }
    herbListALL=AllHerbListData$herb
    herbList2=unlist(str_split(herbList,',| |;'))
    actHerb=herbList2[herbList2%in%herbListALL]
    actHerbID=AllHerbListData[AllHerbListData$herb%in%actHerb,herbID]
    if (length(actHerbID)!=0){
      paID2=paste("'",actHerbID,"'",sep='',collapse = ',')
      paID3=paste("(",paID2,")",sep='')
      herbIDOrder=paste("herbID IN ",paID3,sep='')
      OrderherbID=paste("SELECT * FROM STITCH_HIT_TCMID_TCMSP_Tar where",herbIDOrder,sep=' ')
      tmp <- dbConnect(SQLite(), TarDB)
      #tmp <- dbConnect(SQLite(), 'program/library/base/R/Tar.db')
      res <- dbSendQuery(tmp, OrderherbID)
      subD0 <- fetch(res, n =-1)
      dbDisconnect(tmp)
      subD=join(subD0,AllHerbListData,by='herbID')
      subD=as.data.table(subD)
      subD=subD[,!'herbID',with=F]
      subD=subD[,c('herb',setdiff(colnames(subD),'herb')),with=F]
      subD=subD[QED_DES>=QEDset&score>=geneScore&database%in%targetDatabase,]
      rm(subD0)
      gc()
      subD$geneID=as.character(subD$geneID)
      tempR=geneScoreCal(unique(subD[,.(herb,id=cid,geneID=geneID)]),Pcutoff = geneSelectPv)
      selectGeneID=tempR$selectGeneId###selectGeneID
      if (length(selectGeneID)!=0){
        selectGeneName=geneTrans(selectGeneID,type='ID')
        tempR2=chemScoreCal(unique(subD[,.(id=cid,geneID=geneID)]),geneScore = tempR$geneScore,selectGeneId=selectGeneID)
        selectCID=tempR2$selectChemId#####selectCID
        subDD=unique(subD[,.(cid,chemical_name)])
        setkey(subDD,cid)
        selectCName=subDD[selectCID,chemical_name,mult="first"]
        subD2=subD[cid%in%selectCID&geneID%in%selectGeneID,]###subD2
        if (is.null(diseaseGeneID)){
          Order=paste('geneCluster=list(',herbListName,'=selectGeneID)',sep='')
        }else{
          if (DiseaseNum==1){
            Order=paste('geneCluster=list(',herbListName,'=selectGeneID,',DiseaseName,'=DGeneID[names(DGeneID)==DiseaseName])',sep='')
          }
          ##
          if (DiseaseNum==2){
            Order=paste('geneCluster=list(',herbListName,'=selectGeneID,',DiseaseName[1],'=DGeneID[names(DGeneID)==DiseaseName[1]],',DiseaseName[2],'=DGeneID[names(DGeneID)==DiseaseName[2]])',sep='')
          }
          if (DiseaseNum==3){
            Order=paste('geneCluster=list(',herbListName,'=selectGeneID,',DiseaseName[1],'=DGeneID[names(DGeneID)==DiseaseName[1]],',DiseaseName[2],'=DGeneID[names(DGeneID)==DiseaseName[2]],',DiseaseName[3],'=DGeneID[names(DGeneID)==DiseaseName[3]])',sep='')
          }
        }
        eval(parse(text=Order))
        ###################################
        if (KEGGenrich){
          if (class(try({EnrichKEGG=EnrichMentAnalysis(geneCluster,type='KEGG',qvalueCutoff=qvalueCutoff)},silent=T))=='try-error')EnrichKEGG=NULL
        }else{
          EnrichKEGG=NULL
        }
        ##############
        if (PAenrich){
          if (class(try({EnrichPA=EnrichMentAnalysis(geneCluster,type='ReactomePA',qvalueCutoff=qvalueCutoff)},silent=T))=='try-error')EnrichPA=NULL
        }else{
          EnrichPA=NULL
        }
        ############
        if (GOBPenrich){
          if (class(try({enrichGObp=EnrichMentAnalysis(geneCluster,type='GO',ontGO='BP',qvalueCutoff=qvalueCutoff)},silent=T))=='try-error')enrichGObp=NULL
        }else{
          enrichGObp=NULL
        }
        #############
        if (GOMFenrich){
          if (class(try({enrichGOmf=EnrichMentAnalysis(geneCluster,type='GO',ontGO='MF',qvalueCutoff=qvalueCutoff)},silent=T))=='try-error')enrichGOmf=NULL
        }else{
          enrichGOmf=NULL
        }
        ############
        if (GOCCenrich){
          if (class(try({enrichGOcc=EnrichMentAnalysis(geneCluster,type='GO',ontGO='CC',qvalueCutoff=qvalueCutoff)},silent=T))=='try-error')enrichGOcc=NULL
        }else{
          enrichGOcc=NULL
        }
        ############
        if (DOenrich){
          data(EG2DOLite)
          data(DOLiteTerm)
          data(DOLite2EG)
          if (class(try({enrichmentDO=EnrichMentAnalysis(geneCluster,type='DO',qvalueCutoff=qvalueCutoff)},silent=T))=='try-error')enrichmentDO=NULL
        }else{
          enrichmentDO=NULL
        }
        ############
        if (MESHenrich){
          if (class(try({enrichMESH=EnrichMentAnalysis(geneCluster,type='MESH',meshcategory=meshcategory,qvalueCutoff=qvalueCutoff)},silent=T))=='try-error')enrichMESH=NULL
        }else{
          enrichMESH=NULL
        }
        ##############################################################################
        if (!is.null(diseaseGeneID)&!is.null(EnrichKEGG)&(herbListName%in%EnrichKEGG$Enrichs$Cluster)){
          tempD=EnrichKEGG$Enrichs
          tempD$Cluster=as.character(tempD$Cluster)
          tempD$ID=as.character(tempD$ID)
          item=unique(tempD$Cluster)
          diseaseItem=item[!item%in%herbListName]
          CommonDiseaseKEGGSim=numeric()
          herbListNameItem=tempD$ID[tempD$Cluster==herbListName]
          if (length(diseaseItem)>0){
            length(CommonDiseaseKEGGSim)=length(diseaseItem)
            names(CommonDiseaseKEGGSim)=diseaseItem
            for (i in 1:length(diseaseItem)){
              diseaseItemtemp=tempD$ID[tempD$Cluster==diseaseItem[i]]
              CommonDiseaseKEGGSim[i]=length(intersect(herbListNameItem,diseaseItemtemp))/min(length(herbListNameItem),length(diseaseItemtemp))
            }
            diseaseToCurve=diseaseItem[1]
            diseasetNameItem=tempD$ID[tempD$Cluster==diseaseToCurve]
            TopNIDCurve=min(TopNIDCurve,length(herbListNameItem))
            common_KeggDisaeae=numeric()
            for (k in 1:TopNIDCurve){
              tempC=herbListNameItem[1:k]
              common_KeggDisaeae[k]=length(intersect(tempC,diseasetNameItem))/k
            }
            common_Keggdata=data.frame(ID=1:TopNIDCurve,Disease=common_KeggDisaeae)
            common_Keggdata=as.data.table(common_Keggdata)
            common_Keggdata=melt(common_Keggdata,id='ID',variable='Disease',value='Percentage')
            common_Keggdata=as.data.frame(common_Keggdata)
            CommonKEGGauc=trapz(1:TopNIDCurve,common_KeggDisaeae)
            PlotCommonKEGG=ggplot(common_Keggdata)+aes(x=ID,y=Percentage,group=Disease,colour=Disease,shape=Disease,fill=Disease)+geom_line()+geom_point(size=2)+xlab(paste('Top N enriched-Kegg pathway in ',herbListName))+ylab('Percentage of common enriched-Kegg pathway')+theme(legend.text=element_text(size=10),legend.title=element_text(size=10))+theme(axis.title.x=element_text(size=8),axis.title.y=element_text(size=8))+theme(axis.text.x=element_text(size=6),axis.text.y=element_text(size=6))
          }
          else{
            CommonDiseaseKEGGSim=NULL
            CommonKEGGauc=NULL
            PlotCommonKEGG=NULL
          }
          ######
        }else{
          CommonDiseaseKEGGSim=NULL
          CommonKEGGauc=NULL
          PlotCommonKEGG=NULL
        }
        ###############################################################################
        ###########ReactomePA#############################
        ##############################################################################
        if (!is.null(diseaseGeneID)&!is.null(EnrichPA)&(herbListName%in%EnrichPA$Enrichs$Cluster)){
          tempD=EnrichPA$Enrichs
          tempD$Cluster=as.character(tempD$Cluster)
          tempD$ID=as.character(tempD$ID)
          item=unique(tempD$Cluster)
          diseaseItem=item[!item%in%herbListName]
          CommonDiseasePASim=numeric()
          herbListNameItem=tempD$ID[tempD$Cluster==herbListName]
          if (length(diseaseItem)>0){
            length(CommonDiseasePASim)=length(diseaseItem)
            names(CommonDiseasePASim)=diseaseItem
            for (i in 1:length(diseaseItem)){
              diseaseItemtemp=tempD$ID[tempD$Cluster==diseaseItem[i]]
              CommonDiseasePASim[i]=length(intersect(herbListNameItem,diseaseItemtemp))/min(length(herbListNameItem),length(diseaseItemtemp))
            }
            diseaseToCurve=diseaseItem[1]
            diseasetNameItem=tempD$ID[tempD$Cluster==diseaseToCurve]
            TopNIDCurve=min(TopNIDCurve,length(herbListNameItem))
            common_PADisaeae=numeric()
            for (k in 1:TopNIDCurve){
              tempC=herbListNameItem[1:k]
              common_PADisaeae[k]=length(intersect(tempC,diseasetNameItem))/k
            }
            common_PAdata=data.frame(ID=1:TopNIDCurve,Disease=common_PADisaeae)
            common_PAdata=as.data.table(common_PAdata)
            common_PAdata=melt(common_PAdata,id='ID',variable='Disease',value='Percentage')
            common_PAdata=as.data.frame(common_PAdata)
            CommonPAauc=trapz(1:TopNIDCurve,common_PADisaeae)
            PlotCommonPA=ggplot(common_PAdata)+aes(x=ID,y=Percentage,group=Disease,colour=Disease,shape=Disease,fill=Disease)+geom_line()+geom_point(size=2)+xlab(paste('Top N enriched-Reactome pathway in ',herbListName))+ylab('Percentage of common enriched-Reactome pathway')+theme(legend.text=element_text(size=10),legend.title=element_text(size=10))+theme(axis.title.x=element_text(size=8),axis.title.y=element_text(size=8))+theme(axis.text.x=element_text(size=6),axis.text.y=element_text(size=6))
          }
          else{
            CommonDiseasePASim=NULL
            CommonPAauc=NULL
            PlotCommonPA=NULL
          }
          ######
        }else{
          CommonDiseasePASim=NULL
          CommonPAauc=NULL
          PlotCommonPA=NULL
        }
        ################PA END
        ################################################################################
        ################
        BlankGO=(!is.null(enrichGObp))|(!is.null(enrichGOcc))|(!is.null(enrichGOmf))
        BlankherbListName=(herbListName%in%enrichGObp$Enrichs$Cluster)|(herbListName%in%enrichGOcc$Enrichs$Cluster)|(herbListName%in%enrichGOmf$Enrichs$Cluster)
        if (!is.null(diseaseGeneID)&BlankGO&BlankherbListName){
          tempD=rbind(enrichGObp$Enrichs,enrichGOcc$Enrichs,enrichGOmf$Enrichs)
          tempD=arrange(tempD,qvalue)
          tempD$Cluster=as.character(tempD$Cluster)
          tempD$ID=as.character(tempD$ID)
          item=unique(tempD$Cluster)
          diseaseItem=item[!item%in%herbListName]
          CommonDiseaseGOSim=numeric()
          herbListNameItem=unique(tempD$ID[tempD$Cluster==herbListName])
          if (length(diseaseItem)>0){
            length(CommonDiseaseGOSim)=length(diseaseItem)
            names(CommonDiseaseGOSim)=diseaseItem
            for (i in 1:length(diseaseItem)){
              diseaseItemtemp=unique(tempD$ID[tempD$Cluster==diseaseItem[i]])
              CommonDiseaseGOSim[i]=length(intersect(herbListNameItem,diseaseItemtemp))/min(length(herbListNameItem),length(diseaseItemtemp))
            }
            diseaseToCurve=diseaseItem[1]
            diseasetNameItem=unique(tempD$ID[tempD$Cluster==diseaseToCurve])
            TopNIDCurve=min(TopNIDCurve,length(herbListNameItem))
            common_GODisaeae=numeric()
            for (k in 1:TopNIDCurve){
              tempC=herbListNameItem[1:k]
              common_GODisaeae[k]=length(intersect(tempC,diseasetNameItem))/k
            }
            common_GOdata=data.frame(ID=1:TopNIDCurve,Disease=common_GODisaeae)
            common_GOdata=as.data.table(common_GOdata)
            common_GOdata=melt(common_GOdata,id='ID',variable='Disease',value='Percentage')
            common_GOdata=as.data.frame(common_GOdata)
            CommonGOauc=trapz(1:TopNIDCurve,common_GODisaeae)
            PlotCommonGO=ggplot(common_GOdata)+aes(x=ID,y=Percentage,group=Disease,colour=Disease,shape=Disease,fill=Disease)+geom_line()+geom_point(size=2)+xlab(paste('Top N enriched-GO in ',herbListName))+ylab('Percentage of common enriched-GO')+theme(legend.text=element_text(size=10),legend.title=element_text(size=10))+theme(axis.title.x=element_text(size=8),axis.title.y=element_text(size=8))+theme(axis.text.x=element_text(size=6),axis.text.y=element_text(size=6))
          }else{
            CommonDiseaseGOSim=NULL
            CommonGOauc=NULL
            PlotCommonGO=NULL
          }
          ##############
        }else{
          CommonDiseaseGOSim=NULL
          CommonGOauc=NULL
          PlotCommonGO=NULL
        }
        
        result=list(selectGene=data.table(GeneID=selectGeneID,GeneName=selectGeneName,geneScore=tempR$geneScore),selectChem=data.table(CID=selectCID,ChemName=selectCName,chemScore=tempR2$chemScore[tempR2$chemScore>0]),HerbTarget=subD2,actHerb=actHerb,enrich=list(KEGG=EnrichKEGG,PA=EnrichPA,GObp=enrichGObp,GOmf=enrichGOmf,GOcc=enrichGOcc,DO=enrichmentDO,MESH=enrichMESH),Disease=list(CommonDiseaseKEGGSim=CommonDiseaseKEGGSim,CommonKEGGauc=CommonKEGGauc,PlotCommonKEGG=PlotCommonKEGG,CommonDiseasePASim=CommonDiseasePASim,CommonPAauc=CommonPAauc,PlotCommonPA=PlotCommonPA,CommonDiseaseGOSim=CommonDiseaseGOSim,CommonGOauc=CommonGOauc,PlotCommonGO=PlotCommonGO),geneScore=tempR,chemScore=tempR2)
      }else{
        result=list(selectGene=data.table(GeneID=character(),GeneName=character()),selectChem=data.table(CID=character(),ChemName=character()),HerbTarget=data.table(),actHerb=actHerb,enrich=NULL,Disease=NULL)
      }
      ###############
    }else{
      result=list(selectGene=data.table(GeneID=character(),GeneName=character()),selectChem=data.table(CID=character(),ChemName=character()),HerbTarget=data.table(),actHerb=actHerb,enrich=NULL,Disease=NULL)
    }
    #######################################
  }
  return(result)
}

tempcatch = function(x){
  library(stringr)
  if (class(try(x$a$.attrs['href'],silent=T))=='try-error'){
    temp=x$a['href']
  }else{
    temp=x$a$.attrs['href']
  }
  return(str_extract(temp,pattern='[0-9]+'))
}

theme_dose = function(font.size=14){
  theme_bw() %+replace%
    theme(axis.text.x = element_text(colour = "black",
                                     size = font.size, vjust =1 ),
          axis.title.x = element_text(colour="black",
                                      size = font.size),
          axis.text.y = element_text(colour = "black",
                                     size = font.size, hjust =1 ),
          axis.title.y = element_text(colour="black",
                                      size = font.size, angle=90)
    )
}

ToDataFrame = function(x, itemSep = ",", setStart = "{", setEnd = "}", linebreak = NULL,...) {
  require(arules)
  n_of_itemsets <- length(x)
  if(n_of_itemsets == 0) return(invisible(NULL))
  ## Nothing to inspect here ...

  l <- labels(x, itemSep, setStart, setEnd)
  if(is.null(linebreak))
    linebreak <- any(nchar(l) > options("width")$width*2/3)

  if(!linebreak) {
    out <- data.frame(items = l)
    if(nrow(quality(x)) > 0) out <- cbind(out, quality(x))
    #print(out, right = FALSE)

  } else {


    ## number of rows + fix empty itemsets
    items <- unlist(lapply(as(items(x), "list"),
                           FUN = function(x) if(length(x) == 0) "" else x))
    n_of_items <- size(items(x))
    n_of_items[n_of_items == 0] <- 1

    ## calculate begin and end positions
    entry_beg_pos <- cumsum(c(1, n_of_items[-n_of_itemsets]))
    entry_end_pos <- entry_beg_pos+n_of_items-1

    ## prepare output
    n_of_rows <- sum(n_of_items)
    quality <- quality(x)
    ## Output.
    out <- matrix("", nrow = n_of_rows+1, ncol = 2 + NCOL(quality))

    ## Column 1: itemset nr.
    tmp <- rep.int("", n_of_rows + 1)
    tmp[entry_beg_pos+1] <- c(1:n_of_itemsets)
    out[,1] <- format(tmp)

    ## Column 2: items in the item sets, one per line.
    pre <- rep.int(" ", n_of_rows)
    pre[entry_beg_pos] <- rep.int(setStart, length(entry_beg_pos))
    post <- rep.int(itemSep, n_of_rows)
    post[entry_end_pos] <- rep.int(setEnd, length(entry_end_pos))
    out[, 2] <- format(c("items",
                         paste(pre, unlist(items), post, sep = "")))


    ## Remaining columns: quality measures.
    for(i in seq(length = NCOL(quality))) {
      tmp <- rep.int("", n_of_rows + 1)
      tmp[1] <- names(quality)[i]
      tmp[entry_end_pos + 1] <- format(quality[[i]])
      out[, i + 2] <- format(tmp, justify = "right")
    }

    ## Output.
    cat(t(out), sep = c(rep.int(" ", NCOL(out) - 1), "\n"))
  }
  return(out)
}

toVnaNet = function(net,filename='net.vna',NodeAtr=F,NodeData){
  ##net:igraphObject,V(net)$name,E(net)$weight
  ##filename:*.vna
  ##NodeAtr:whether to import Node Attributes DataFrame
  ##NodeData:Node Attributes DataFrame
  require(igraph)
  require(readr)
  require(data.table)
  a=as_data_frame(net,  what  =  c("both"))
  if (NodeAtr){
    Node=NodeData
  }else{
    if (ncol(a$vertices)==0){
      a$vertices=data.frame(Node=as.character(1:nrow(a$vertices)))
    }
    Node=a$vertices
  }
  #####################
  if (ncol(a$edges)==2){
    a$edges$weight=1
  }
  Edge=as.data.table(a$edges)
  Edge=unique(Edge)
  Edge=as.data.frame(Edge)
  Edge$from=as.character(Edge$from)
  Edge$to=as.character(Edge$to)
  write_lines('*Node data',filename)
  write_csv(Node,filename,append=T,col_names=T)
  write_lines('*Tie data',filename,append=T)
  write_csv(Edge,filename,append=T,col_names=T)
}

TransInchikey = function(molecularFile,molecularType='SMI'){
  #molecularType='SMI':#molecularFile:SMI.csv-->colnames(molecularFile)=c('id','smiles')
  #molecularType='SDF':#molecularFile:sdf-->datablocktag(sdf,'id')含id列
  #molecularType='SMItxt':SMI character text
  library(data.table)
  library(plyr)
  library(stringr)
  library(clipr)
  library(igraph)
  library(readr)
  library(ggplot2)
  library(ChemmineR)
  if (molecularType=='SMI'){
    SMI=fread(molecularFile)
    SMI2=SMI[smiles!='',]
    smiles=SMI2$smiles
    names(smiles)=SMI2$id
    SMISdfSets=ChemmineR::smiles2sdf(smiles)
    SMISdfSets=SMISdfSets[validSDF(SMISdfSets)]
    #sum(validSDF(SMISdfSets))==length(SMISdfSets)
    datablock(SMISdfSets)<-data.frame(id=SMISdfSets@ID)
    CID=datablocktag(SMISdfSets,'id')
    Nmol=length(CID)
    inchikey=character()
    length(inchikey)=Nmol
    names(inchikey)=CID
    for (i in 1:Nmol){
      TempSdf=write.SDF(SMISdfSets[i],'TempSdf.sdf',cid=T)
      TempName=CID[i]
      #Outputfilename=paste(TempName,'mol2',sep='.')
      Dos=paste('obabel -isdf TempSdf.sdf -oinchikey','-Oinchikey.txt',sep=' ')
      shell(Dos,intern=T)
      if (file.exists('inchikey.txt')){
        key=read_lines('inchikey.txt')
        inchikey[TempName]=key
      }else{
        inchikey[TempName]=NA
      }
    }
  }else if (molecularType=='SDF'){
    SMISdfSets=read.SDFset(molecularFile)
    Nmol=length(SMISdfSets)
    CID=datablocktag(SMISdfSets,'id')
    inchikey=character()
    length(inchikey)=Nmol
    names(inchikey)=CID
    for (i in 1:Nmol){
      TempSdf=write.SDF(SMISdfSets[i],'TempSdf.sdf',cid=T)
      TempName=CID[i]
      #Outputfilename=paste(TempName,'mol2',sep='.')
      Dos=paste('obabel -isdf TempSdf.sdf -oinchikey','-Oinchikey.txt',sep=' ')
      shell(Dos,intern=T)
      if (file.exists('inchikey.txt')){
        key=read_lines('inchikey.txt')
        inchikey[TempName]=key
      }else{
        inchikey[TempName]=NA
      }
    }
  }else if (molecularType=='SMItxt'){
    SMISdfSets=ChemmineR::smiles2sdf(molecularFile)
    #SMISdfSets=SMISdfSets[validSDF(SMISdfSets)]
    #sum(validSDF(SMISdfSets))==length(SMISdfSets)
    datablock(SMISdfSets)<-data.frame(id=paste('temp',1:length(SMISdfSets),sep=''))
    CID=datablocktag(SMISdfSets,'id')
    Nmol=length(CID)
    inchikey=character()
    length(inchikey)=Nmol
    names(inchikey)=CID
    for (i in 1:Nmol){
      TempSdf=write.SDF(SMISdfSets[i],'TempSdf.sdf',cid=T)
      TempName=CID[i]
      #Outputfilename=paste(TempName,'mol2',sep='.')
      Dos=paste('obabel -isdf TempSdf.sdf -oinchikey','-Oinchikey.txt',sep=' ')
      shell(Dos,intern=T)
      if (file.exists('inchikey.txt')){
        key=read_lines('inchikey.txt')
        inchikey[TempName]=key
      }else{
        inchikey[TempName]=NA
      }
    }
  }else{
    print('ERROR:molecularType not found!')
    inchikey=character()
  }
return(as.data.frame(inchikey,stringsAsFactors = FALSE))
}

TransMol = function(sdfFile){
  library(data.table)
  library(plyr)
  library(stringr)
  library(igraph)
  library(ggplot2)
  library(ChemmineR)
  library(PaDEL)
  SDF=read.SDFset(sdfFile)
  Nmol=length(SDF)
  CID=datablocktag(SDF,'id')
  for (i in 1:Nmol){
    TempSdf=write.SDF(SDF[i],'TempSdf.sdf',cid=T)
    TempName=CID[i]
    Outputfilename=paste(TempName,'mol2',sep='.')
    Dos=paste('obabel -isdf TempSdf.sdf -oMOL2','-O',Outputfilename,'--gen3D',sep=' ')
    shell(Dos,intern=T)
  }
}

upDateDataTargetFromStitch = function(dataSource='db/HITDatabase_herb_Chem_TargetALL.csv',STITCH_CTdb='db/STITCH5_CTdb.db',STITCH_CKEYdb='db/STITCH5_CKEYdb.db',outputfile=NULL){
  #dataSource:HITDatabase_herb_Chem_TargetALL.csv,TCMIDDatabase_herb_Chem_TargetALL.csv,TCMSPDatabase_herb_Chem_TargetALL.csv【colnames:"herb","cid","chemical_name","geneID","QED_DES","inchikey","smiles","database","score"】
  #chem_inchikey:inchikey
  #STITCH_CTdb:'STITCH5_CTdb.db',colnames=c('cid','ENSP','score','geneID','inchikey')
  #STITCH_CKEYdb:'STITCH_CKEYdb.db',colnames=c('inchikey','source_cid','stereo_chemical_id','flat_chemical_id')
  #StitchSearchout：data.frame(cid,ENSP,geneID,score,inchikey)
  #result:outputfile同STITCH_HIT_TCMID_TCMSP_Tar.csv【colnames:"herb","chemical_name","geneID","score","inchikey","smiles","database","QED_DES","cid"】
  require(RSQLite)
  require(stringr)
  require(data.table)
  require(foreign)
  require(plyr)
  options(stringsAsFactors = F)
  judgeTarget=function(x,databaseNameSet){
    if (length(unique(x))>1){
      result='STITCH'
    }else{
      result=databaseNameSet
    }
    return(result)
  }
  COLname=c("herb","cid","chemical_name","geneID","score","QED_DES","inchikey","smiles","database")
  dataSource=fread(dataSource,encoding='UTF-8')
  dataSource=as.data.frame(dataSource)
  dataSource=dataSource[,colnames(dataSource)%in%COLname]
  NotINCOLname=COLname[!COLname%in%colnames(dataSource)]
  if (length(NotINCOLname)>0){
    for (i in NotINCOLname){
      if (i=='score'){
        dataSource[,i]=9999
      }else if (i=='database'){
        dataSource[,i]='CUSTOM'
      }else{
        dataSource[,i]=NA
      }
    }
  }
  DatsourceName=unique(dataSource$database)[1]
  chemKey=unique(dataSource$inchikey)
  SearchSTITCH=findStitchTargetAll(chem_inchikey=chemKey,scoreSet=0,STITCH_CTdb=STITCH_CTdb,STITCH_CKEYdb=STITCH_CKEYdb)
  SearchSTITCH$database='STITCH'
  SearchSTITCH=as.data.frame(SearchSTITCH)
  SearchSTITCHsub=SearchSTITCH[,c('inchikey','geneID','score','database')]
  dataSourcesub=dataSource[,c('inchikey','geneID','score','database')]
  dataSourceCombine=rbind(dataSourcesub,SearchSTITCHsub)
  dataSourceCombine=as.data.table(dataSourceCombine)
  dataSourceCombine=dataSourceCombine[,.(score=min(score),database=judgeTarget(score,databaseNameSet=DatsourceName)),by=c('inchikey','geneID')]
  dataSourceTarget=dataSource[,c("herb","cid","chemical_name","QED_DES","inchikey","smiles")]
  dataSourceTarget=as.data.table(dataSourceTarget)
  dataSourceTarget=unique(dataSourceTarget)
  dataSourceTarget=join(dataSourceTarget,dataSourceCombine[,c('inchikey','geneID','score','database')],by='inchikey')
  dataSourceTarget=dataSourceTarget[!is.na(dataSourceTarget$geneID),]
  dataSourceTarget=unique(dataSourceTarget)
  if (is.null(outputfile)){
    UpTime=as.character(Sys.time())
    UpTime=str_replace_all(UpTime,' ','_')
    UpTime=str_replace_all(UpTime,':','_')
    fileName=paste(UpTime,'_upDated_',DatsourceName,'.csv',sep='')
  }else{
    fileName=outputfile
  }
  write.csv(dataSourceTarget,fileName,row.names = F)
}

UpdateHerbName = function(TargetFile='db/STITCH_HIT_TCMID_TCMSP_Tar.csv',HerbNameData){
  require(data.table)
  require(RSQLite)
  require(plyr)
  require(lubridate)
  require(stringr)
  options(stringsAsFactors = F)
  #HerbNameData-->colnames('oldName','newName');newName中标0的进行删除操作！
  #TargetFile-->由UpDateTargetData得到，需要进一步规整药名herb列
  #outFile-->新的STITCH_HIT_TCMID_TCMSP_Tar.db-->D:\PA2.1Plus\program\library\base\R\Tar.db；###新的AllHerbListData-->within(pa.db)##新的STITCH_HIT_TCMID_TCMSP_Tar.csv-->For updata later in the future
  HerbName=fread(HerbNameData,encoding = 'UTF-8')
  HerbName=HerbName[!is.na(HerbName$newName),]
  HerbName=HerbName[HerbName$newName!='',]
  excludName=unique(HerbName$oldName[HerbName$newName=='0'])
  HerbName=HerbName[HerbName$newName!='0',]
  oldName=HerbName$oldName
  newName=HerbName$newName
  replaceNameVector=newName
  names(replaceNameVector)=oldName
  TargetFile=fread(TargetFile,encoding='UTF-8')
  TargetFile=TargetFile[!TargetFile$herb%in%excludName,]
  TargetFile$herb=HerbName$newName[chmatch(TargetFile$herb,HerbName$oldName)]
  #TargetFile$herb=str_replace_all(TargetFile$herb,pattern =replaceNameVector )
  UpTime=as.character(Sys.time())
  UpTime=str_replace_all(UpTime,' ','_')
  UpTime=str_replace_all(UpTime,':','_')
  ALL=TargetFile
  ALL=unique(ALL,by=c('herb','inchikey','geneID','database'))
  AllHerbList=unique(ALL$herb)
  HerbID=paste('h',1:length(AllHerbList),sep='')
  AllHerbListData=data.table(herbID=HerbID,herb=AllHerbList)
  save(list='AllHerbListData',file=paste(UpTime,'AllHerbListData.RData',sep='_'))###AllHerbListData-->within(pa.db)
  ########
  NumCompound=length(unique(ALL$inchikey))
  cid=paste('c',1:NumCompound,sep='')
  CompoundListData=data.table(cid2=cid,inchikey=unique(ALL$inchikey))
  ALL2=plyr::join(ALL,AllHerbListData,by='herb')
  ALL3=plyr::join(ALL2,CompoundListData,by='inchikey')
  ALL3=ALL3[,!'cid',with=F]
  ALL3=plyr::rename(ALL3,c('cid2'='cid'))
  ALLtoSave=ALL3[,!'herbID',with=F]
  write.csv(ALLtoSave,paste(UpTime,'STITCH_HIT_TCMID_TCMSP_Tar.csv',sep='_'),row.names = F)##STITCH_HIT_TCMID_TCMSP_Tar.csv-->For updata later in the future
  ALL4=ALL3[,!'herb',with=F]
  ColNameALL3=c('herbID','cid','chemical_name',"geneID","score","inchikey","smiles","database",'QED_DES')
  if (!setequal(ColNameALL3,colnames(ALL4))){
    print('ERROR:colnames')
    stop()
  }
  ALL4=ALL4[,ColNameALL3,with=F]
  write.csv(ALL4,paste(UpTime,'STITCH_HIT_TCMID_TCMSP_Tar_db.csv',sep='_'),row.names = F)
  ########
  tmp <- dbConnect(SQLite(), paste(UpTime,'STITCH_HIT_TCMID_TCMSP_Tar.db',sep='_'))##STITCH_HIT_TCMID_TCMSP_Tar.db-->D:\PA2.1Plus\program\library\base\R\Tar.db
  dbWriteTable(tmp,'STITCH_HIT_TCMID_TCMSP_Tar',ALL4,append=T)
  dbSendQuery(tmp,'create index index_herbID on STITCH_HIT_TCMID_TCMSP_Tar (herbID)')
  dbDisconnect(tmp)
}

updateKEGGdb = function(OLDKEGG_data='E:/DATABASE/MYpro/BanXiaXieXin/result/ToolBox/db/KEGG_data.csv'){
  ##KEGG_data.db-->KEGG_data.csv[colnames:PathwayID,geneID,Annotatiion]
  ##OLDKEGG_data='E:/DATABASE/MYpro/BanXiaXieXin/result/ToolBox/db/KEGG_data.csv'
  ##KEGG_data.csv-->KEGG_data.db
  library(KEGGgraph)
  library(KEGGREST)
  library(stringr)
  library(data.table)
  options(stringsAsFactors = F)
  OldKEGG_data=fread(OLDKEGG_data)
  res <- keggList("pathway", "hsa")
  KEGG_data=data.frame()
  for (i in 1:length(res)){
    print(i)
    tempName=names(res[i])
    tempName=str_replace_all(tempName,'path:','')
    b=keggGet(tempName, "kgml")
    bnode=parseKGML(b)
    bnode2=nodes(bnode)
    bnode2=bnode2[sapply(bnode2,function(x){x@type})=='gene']
    bnode3=unique(unlist(lapply(bnode2, function(x){x@name})))
    geneID_Temp=translateKEGGID2GeneID(bnode3, organism="hsa")
    KEGG_dataTemp=data.frame(PathwayID=tempName,geneID=geneID_Temp,Annotatiion=res[i])
    KEGG_data=rbind(KEGG_data,KEGG_dataTemp)
  }
  UpTime=as.character(Sys.time())
  UpTime=str_replace_all(UpTime,' ','_')
  UpTime=str_replace_all(UpTime,':','_')
  print(paste0('update paths:',toString(setdiff(unique(KEGG_data$PathwayID),unique(OldKEGG_data$PathwayID)))))
  print(paste0('update gens:',toString(setdiff(unique(KEGG_data$geneID),unique(OldKEGG_data$geneID)))))
  print(paste0('update links:',nrow(KEGG_data)-nrow(OldKEGG_data)))
  write.csv(KEGG_data,file=paste(UpTime,'KEGG_data.csv',sep='_'),row.names=F)
}

updateKEGGsqlite = function(hsa00001json='E:/DATABASE/MYpro/BanXiaXieXin/result/ToolBox/db/KEGG.db/hsa00001_190502.json',KEGGsqlite="E:/DATABASE/MYpro/BanXiaXieXin/result/ToolBox/db/KEGG.db/KEGG.sqlite"){
  library("RSQLite")
  library(data.table)
  library(jsonlite)
  library(stringr)
  ##KEGG.sqlite-->KEGG.db/extdata/KEGG.sqlite
  ##hsa00001json-->https://www.genome.jp/kegg-bin/get_htext?hsa00001+3101
  ##output:list(PathGeneData,PathGeneDataSub,KEGG_data):PathGeneDataSub:dataFrame(pathway2gene-->KEGG.sqlite)。PathGeneData:dataFrame(HSApathdata,colnames:pathway_id,gene_or_orf_id)。KEGG_data-->KEGG_data.csv用于updatePathSimM2
  options(stringsAsFactors = F)
  a=read_json(path=hsa00001json, simplifyVector = TRUE)
  PathGeneData=data.frame()
  class1name=a[[2]]$name
  for (i in 1:length(class1name)){
    class2name=a[[2]]$children[[i]]$name
    for (j in 1:length(class2name)){
      pathName=a[[2]]$children[[i]]$children[[j]]$name
      for (k in 1:length(pathName)){
        pathName2=pathName[k]
        if (str_detect(pathName2,'PATH:hsa')){
          pathid0=str_extract(pathName2,'\\[.*\\]')
          pathid=str_replace_all(pathid0,'\\[PATH:|\\]','')
          tempid=str_replace(pathid,'hsa','')
          pathName2ok=str_replace_all(pathName2,'\\[.*\\]','')
          pathName2ok=str_trim(str_replace_all(pathName2ok,tempid,''),'both')
          ########
          tempa=a[[2]]$children[[i]]$children[[j]]$children[[k]]$name
          tempa2=sapply(tempa,function(x){y=unlist(str_split(x,';'))[1];y=unlist(str_split(y,' '))[1];return(y)})
          names(tempa2)=NULL
          tempPathGene=data.frame(gene_or_orf_id=tempa2)
          tempPathGene$pathway_id=pathid
          tempPathGene$pathway_name=pathName2ok
          tempPathGene$class1_name=class1name[i]
          tempPathGene$class2_name=class2name[j]
          PathGeneData=rbind(PathGeneData,tempPathGene)
        }
      }
    }
  }
  ###############
  con <- dbConnect(RSQLite::SQLite(), dbname=KEGGsqlite)
  #dbListTables(con)
  pathway2name0=dbGetQuery(con,"select * from pathway2name")
  PathGeneData_pathway2name=PathGeneData[,c('pathway_id','pathway_name')]
  colnames(PathGeneData_pathway2name)=c('path_id','path_name')
  PathGeneData_pathway2name=as.data.table(PathGeneData_pathway2name)
  PathGeneData_pathway2name=unique(PathGeneData_pathway2name)
  PathGeneData_pathway2name$path_id=str_replace_all(PathGeneData_pathway2name$path_id,'hsa','')
  pathway2name0=as.data.table(pathway2name0)
  pathway2name=rbind(pathway2name0,PathGeneData_pathway2name)
  pathway2name=unique(pathway2name,by='path_id')
  pathway2name=as.data.frame(pathway2name)
  dbWriteTable(con,name='pathway2name',value=pathway2name,overwrite=T)
  #####
  pathway2gene0=dbGetQuery(con,"select * from pathway2gene")
  pathway2gene_excludHSA=pathway2gene0[!str_detect(pathway2gene0$pathway_id,'hsa'),]
  PathGeneDataSub=rbind(PathGeneData[,colnames(pathway2gene0)],pathway2gene_excludHSA)
  dbWriteTable(con,name='pathway2gene',value=PathGeneDataSub,overwrite=T)
  dbDisconnect(con)
  KEGG_data=PathGeneData[,c('pathway_id','gene_or_orf_id','pathway_name')]
  colnames(KEGG_data)=c('PathwayID','geneID','Annotatiion')
  write.csv(KEGG_data,'KEGG_data.csv',row.names = F)
  return(list(PathGeneData=PathGeneData,PathGeneDataSub=PathGeneDataSub,KEGG_data=KEGG_data))
}

updatePPiNet = function(HIPPIEfiletxt='E:/DATABASE/HIPPIE-PPI database-Human Integrated Protein-Protein Interaction rEference/hippie_current_v2.2.txt'){
##ppiNetData.db-->包含ppiNet,ppiBinaryNet
  ##output-->ppiNetData.Rdata-->ppiNetData.db
  library(data.table)
  library(igraph)
  library(stringr)
  #ppi_data=fread('E://DATABASE//HIPPIE-PPI database-Human Integrated Protein-Protein Interaction rEference//hippie_current.txt')
  #HIPPIEfiletxt='E:/DATABASE/HIPPIE-PPI database-Human Integrated Protein-Protein Interaction rEference/hippie_current_v2.2.txt'
  ppi_data=fread(HIPPIEfiletxt)
  colnames(ppi_data)=c('geneName1','geneId1','geneName2','geneId2','score','annotation')
  ppi_data_sub=ppi_data[score>0,]
  ppi_data_sub=ppi_data_sub[geneId1!=geneId2,]
  ppi_data_sub=as.data.frame(ppi_data_sub)
  ppi_data_sub_forNet=subset(ppi_data_sub,select=c('geneId1','geneId2','score'))
  ppiNet=graph_from_data_frame(ppi_data_sub_forNet,directed=F)
  ppiBinaryNet=graph_from_data_frame(ppi_data_sub_forNet[,1:2],directed=F)
  rm(ppi_data,ppi_data_sub)
  UpTime=as.character(Sys.time())
  UpTime=str_replace_all(UpTime,' ','_')
  UpTime=str_replace_all(UpTime,':','_')
  save(ppiNet,ppiBinaryNet,file=paste(UpTime,'ppiNetData.Rdata',sep='_'))
}

updatePPIpath = function(ppiNetdata='ppiNetData.Rdata'){
  ##PPIpath.db-->ppiBinaryNet,ppiNet,ppiNetMatrix,ppiNetMatrix2,ppiNetMatrix3
  ##PPIpath.Rdata-->PPIpath.db
  ##ppiNetdata='2019-03-25_10_38_18_ppiNetData.Rdata'
  library(linkcomm)
  library(ggplot2)
  library(expm)
  library(igraph)
  library(stringr)
  load(ppiNetdata)
  ppiNetMatrix=as_adjacency_matrix(ppiBinaryNet)
  ppiNetMatrix=as.matrix(ppiNetMatrix)
  ppiNetMatrix2=ppiNetMatrix%^%2
  ppiNetMatrix3=ppiNetMatrix%^%3
  diag(ppiNetMatrix)=0
  diag(ppiNetMatrix2)=0
  diag(ppiNetMatrix3)=0
  UpTime=as.character(Sys.time())
  UpTime=str_replace_all(UpTime,' ','_')
  UpTime=str_replace_all(UpTime,':','_')
  save(ppiBinaryNet,ppiNet,ppiNetMatrix,ppiNetMatrix2,ppiNetMatrix3,file=paste(UpTime,'PPIpath.Rdata',sep='_'))
}

UpDateTargetData = function(oldTarData='db/STITCH_HIT_TCMID_TCMSP_Tar.csv',newTarData=NULL){
  require(data.table)
  require(RSQLite)
  require(plyr)
  require(lubridate)
  require(stringr)
  #oldTarData:"herb","chemical_name","geneID","score","inchikey","smiles","database","QED_DES","cid"
  ##herbID与cid将重新编号！
  ##Outfile:STITCH_HIT_TCMID_TCMSP_Tar.db-->D:\PA2.1Plus\program\library\base\R\Tar.db；###AllHerbListData-->within(pa.db)##STITCH_HIT_TCMID_TCMSP_Tar.csv-->For updata later in the future
  oldD=fread(oldTarData)
  if (is.null(newTarData)){
    newD=data.table()
  }else{
    newD=fread(newTarData)
	if (sum(colnames(oldD)%in%colnames(newD))==ncol(oldD)){
	newD=newD[,colnames(oldD),with=F]
	}else{
	print('ERROR:colnames for new Data')
	stop()
	}
}
  UpTime=as.character(Sys.time())
  UpTime=str_replace_all(UpTime,' ','_')
  UpTime=str_replace_all(UpTime,':','_')
  ALL=rbind(oldD,newD)
  ALL=unique(ALL,by=c('herb','inchikey','geneID','database'))
  AllHerbList=unique(ALL$herb)
  HerbID=paste('h',1:length(AllHerbList),sep='')
  AllHerbListData=data.table(herbID=HerbID,herb=AllHerbList)
  save(list='AllHerbListData',file=paste(UpTime,'AllHerbListData.RData',sep='_'))###AllHerbListData-->within(pa.db)
  ########
  NumCompound=length(unique(ALL$inchikey))
  cid=paste('c',1:NumCompound,sep='')
  CompoundListData=data.table(cid2=cid,inchikey=unique(ALL$inchikey))
  ALL2=plyr::join(ALL,AllHerbListData,by='herb')
  ALL3=plyr::join(ALL2,CompoundListData,by='inchikey')
  ALL3=ALL3[,!'cid',with=F]
  ALL3=plyr::rename(ALL3,c('cid2'='cid'))
  ALLtoSave=ALL3[,!'herbID',with=F]
  write.csv(ALLtoSave,paste(UpTime,'STITCH_HIT_TCMID_TCMSP_Tar.csv',sep='_'),row.names = F)##STITCH_HIT_TCMID_TCMSP_Tar.csv-->For updata later in the future
  ALL4=ALL3[,!'herb',with=F]
  ColNameALL3=c('herbID','cid','chemical_name',"geneID","score","inchikey","smiles","database",'QED_DES')
  if (!setequal(ColNameALL3,colnames(ALL4))){
    print('ERROR:colnames')
    stop()
  }
  ALL4=ALL4[,ColNameALL3,with=F]
  write.csv(ALL4,paste(UpTime,'STITCH_HIT_TCMID_TCMSP_Tar_db.csv',sep='_'),row.names = F)
  ########
  tmp <- dbConnect(SQLite(), paste(UpTime,'STITCH_HIT_TCMID_TCMSP_Tar.db',sep='_'))##STITCH_HIT_TCMID_TCMSP_Tar.db-->D:\PA2.1Plus\program\library\base\R\Tar.db
  dbWriteTable(tmp,'STITCH_HIT_TCMID_TCMSP_Tar',ALL4,append=T)
  dbSendQuery(tmp,'create index index_herbID on STITCH_HIT_TCMID_TCMSP_Tar (herbID)')
  dbSendQuery(tmp,'create index index_inchikey on STITCH_HIT_TCMID_TCMSP_Tar (inchikey)')
  dbDisconnect(tmp)
}

WNS = function(pathID1,cid1Genetemp,pathID2,cid2Genetemp,keggList,WWIDis,ppiDis){
  require(igraph)
  require(data.table)
  Npath1=length(pathID1)
  Npath2=length(pathID2)
  pathDis=matrix(0,nrow=Npath1,ncol=Npath2)
  rownames(pathDis)=pathID1
  colnames(pathDis)=pathID2
  for (i in 1:Npath1){
    for (j in 1:Npath2){
      path1=pathID1[i]
      path2=pathID2[j]
      if (path1!=path2){
        pathDis[path1,path2]=WWIDis[path1,path2]
      }else{
        tempGene=keggList[[path1]]
        cid1geneIN=cid1Genetemp[cid1Genetemp%in%tempGene]
        cid1geneIN=cid1geneIN[cid1geneIN%in%rownames(ppiDis)]
        cid2geneIN=cid2Genetemp[cid2Genetemp%in%tempGene]
        cid2geneIN=cid2geneIN[cid2geneIN%in%rownames(ppiDis)]
        if(length(cid1geneIN)==0|length(cid2geneIN)==0){
          pathDis[path1,path2]=Inf
        }else{
          ######################################
          ####oldVersion######
          tempPPIdis=ppiDis[cid1geneIN,cid2geneIN]
          tempPPIdis=sum(tempPPIdis)
          pathDis[path1,path2]=exp(-tempPPIdis/(length(cid1geneIN)*length(cid2geneIN)))
          ##############################
          ###newVersion####
          #tempPPIdis=ppiDis[cid1geneIN,cid2geneIN]
          #pathDis[path1,path2]=mean(exp(-tempPPIdis/(length(cid1geneIN)*length(cid2geneIN))))##newVersion
        }
      }
    }
  }
  ######################################
  ####oldVersion######
  SumPathDis=sum(pathDis)##
  Score=exp(-SumPathDis/(length(pathID1)*length(pathID2)))
  ##############################
  ###newVersion####
  #Score=mean(exp(-pathDis/(length(pathID1)*length(pathID2))))##newVersion
  return(Score)
}

WNSscore = function(chem_target,ppiBinaryNet=ppiBinaryNet,pathSimM2=pathSimM2,keggPath=keggPath){
  #chem_target:data.frame-->colnames(chem_target)=c('id','geneID')
  library(data.table)
  library(plyr)
  library(stringr)
  library(clipr)
  library(igraph)
  library(ggplot2)
  library(ChemmineR)
  library(pracma)
  chem_target$geneID=as.character(chem_target$geneID)
  chem_target=as.data.table(chem_target)
  if (!exists('ppiBinaryNet')){
    #load('db/ppiNetData.RData')
    load('db/ppiNetData.db')
  }
  ppiDis=distances(ppiBinaryNet)######ppiDis
  if (!exists('pathSimM2')){
    #load('db/pathSimM2.RData')
    load('db/pathSimM2.db')
  }
  WWI=graph_from_adjacency_matrix(pathSimM2,mode='undirected',diag = F)
  degreePathID=igraph::degree(WWI)
  WWI2=WWI-names(degreePathID)[degreePathID==0]
  DrugId=unique(chem_target$id)
  chem_target$geneID=as.character(chem_target$geneID)
  ##ppiDis,WWIDis
  WWIDis=distances(WWI2)
  if (!exists('keggPath')){
    #keggPath=fread('db/KEGG_data.csv')
    keggPath=fread('db/KEGG_data.db',encoding='UTF-8')
  }
  keggPath$geneID=as.character(keggPath$geneID)
  keggPath=keggPath[PathwayID%in%V(WWI2)$name,]
  PathID=unique(keggPath$PathwayID)
  kegGeneID=unique(keggPath$geneID)
  keggList=list()
  length(keggList)=length(PathID)
  names(keggList)=PathID
  for (i in PathID){
    keggList[[i]]=keggPath[PathwayID==i,geneID]
  }
  #####
  WNSMatrix=matrix(0,nrow=length(DrugId),ncol=length(DrugId))
  rownames(WNSMatrix)=DrugId
  colnames(WNSMatrix)=DrugId
  for (i in 1:length(DrugId)){
    j=i+1
    while(j<=length(DrugId)){
      cid1=DrugId[i]
      cid1Genetemp=chem_target[id==cid1,geneID]
      cid1Genetemp=cid1Genetemp[cid1Genetemp%in%rownames(ppiDis)]
      keggid1=keggPath[geneID%in%cid1Genetemp,unique(PathwayID)]
      cid2=DrugId[j]
      cid2Genetemp=chem_target[id==cid2,geneID]
      cid2Genetemp=cid2Genetemp[cid2Genetemp%in%rownames(ppiDis)]
      keggid2=keggPath[geneID%in%cid2Genetemp,unique(PathwayID)]
      if (length(keggid1)==0|length(keggid2)==0){
        WNSMatrix[cid1,cid2]=0
      }else{
        WNSMatrix[cid1,cid2]=WNS(pathID1=keggid1,cid1Genetemp,pathID2=keggid2,cid2Genetemp,keggList,WWIDis,ppiDis)
      }
      j=j+1
    }
  }
  rm(list=c('ppiBinaryNet','WWIDis','ppiDis'))
  WNSMatrix=asSymmetric(WNSMatrix)
  return(WNSMatrix)
}
