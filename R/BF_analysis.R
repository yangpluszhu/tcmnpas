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

