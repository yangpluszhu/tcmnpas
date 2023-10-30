pdftotxt = function(pdffile){
  output=str_replace(pdffile,'\\.pdf','.txt')
  dos=paste('pdftotext',pdffile,'-enc UTF-8')
  shell(dos)
  #Sys.sleep(1)
  return(output)
}

