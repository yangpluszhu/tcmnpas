loadALLFunctions = function(){
  setwd('E:/DATABASE/MYpro/BanXiaXieXin/result/ToolBox/shinyFunctions')
  codes=list.files('functions/',pattern='.R')
  codes=codes[!codes%in%c('TempCode.R','loadALLFunctions.R')]
  for (i in codes){
    source(paste('functions/',i,sep=''),encoding = 'UTF-8')
  }
}

