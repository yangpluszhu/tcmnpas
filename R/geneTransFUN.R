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
 
