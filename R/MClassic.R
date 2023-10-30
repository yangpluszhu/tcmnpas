MClassic = function(Cname,classicF_Matrix,TopN){
  require(data.table)
  require(stringr)
  Nyao=length(Cname)
  classicFlength=rowSums(classicF_Matrix)
  subclassicF_Matrix=classicF_Matrix[,Cname,drop=F]
  similarity=rowSums(subclassicF_Matrix)/classicFlength
  names(similarity)=rownames(classicF_Matrix)
  similarity=similarity[similarity>0]
  similarity=sort(similarity,decreasing = T)
  similarity=round(similarity,digits = 2)
  if (length(similarity)>=TopN){
    similarity2=similarity[1:TopN]
    subM=classicF_Matrix[names(similarity2),Cname,drop=F]
    Comember=apply(subM,1,function(x){paste(Cname[which(x==1)],collapse = ',')})
    temp=paste(names(similarity2),'S=',similarity2,'[',Comember,']',sep='')
    Sim=paste(temp,collapse = ';')
    MaxSim=max(similarity2)
    MeanSim=mean(similarity2)
  }else if (length(similarity)>0){
    similarity2=similarity
    subM=classicF_Matrix[names(similarity2),Cname,drop=F]
    Comember=apply(subM,1,function(x){paste(Cname[which(x==1)],collapse = ',')})
    temp=paste(names(similarity2),'S=',similarity2,'[',Comember,']',sep='')
    Sim=paste(temp,collapse = ';')
    MaxSim=max(similarity2)
    MeanSim=mean(similarity2)
  }else{
    Sim=NA
    MaxSim=0
    MeanSim=0
  }
  return(list(Sim=Sim,MaxSim=MaxSim,MeanSim=MeanSim))
}

