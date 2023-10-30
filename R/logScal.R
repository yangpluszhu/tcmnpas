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

