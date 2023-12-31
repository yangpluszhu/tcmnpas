getChemFromHerb = function(HerbList,QEDset=0.2,Database=c('HIT','TCMID','TCMSP','CUSTOM'),herb_ChemDataAll=herb_ChemDataAll,Lipinski_Viosset=0,Veber_Viosset=0){
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
  HerbListData=herb_ChemDataAll[herb_ChemDataAll$herb%in%HerbList&herb_ChemDataAll$QED_DES>=QEDset&herb_ChemDataAll$database%in%Database&herb_ChemDataAll$Lipinski_Vios<=Lipinski_Viosset&herb_ChemDataAll$Veber_Vios<=Veber_Viosset,]
  return(HerbListData)
}
 
