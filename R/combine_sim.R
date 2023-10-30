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

