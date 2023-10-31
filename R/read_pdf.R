read_pdf = function(pdffile,sep=FALSE){ 
  require(stringr) 
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
 
