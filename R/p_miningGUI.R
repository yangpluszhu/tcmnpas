p_miningGUI = function(请输入基本方分析的支持度阈值=0.03,请输入基本方的最小含药数=8,请输入基本方的期望含药数=10,基本方是否包含西药='否'){
  require(RSQLite)
  require(lubridate)
  require(stringr)
  require(data.table)
  require(foreign)
  require(fgui)
  #detach(package:DescTools)
  #prescription <- guiGetSafe("prescription")
  tS0.9=请输入基本方分析的支持度阈值
  tS0.9=as.numeric(tS0.9)
  min_yao=请输入基本方的最小含药数
  min_yao=as.numeric(min_yao)
  threshold_leastYao=请输入基本方的期望含药数
  threshold_leastYao=as.numeric(threshold_leastYao)
  if (基本方是否包含西药=='是'){
    onlyHerb=FALSE
  }else if (基本方是否包含西药=='否'){
    onlyHerb=TRUE
  }
  ####
  ##
  p_mining_Result=p_mining(tS0.9=tS0.9,min_yao=min_yao,threshold_leastYao=threshold_leastYao,onlyHerb=onlyHerb)
  txt='恭喜您！核心药物及核心处方挖掘分析完毕！请在Repor文件夹中Tab4_1核心药物分析结果.csv;Tab4_2核心处方分析结果.csv检视结果！'
  #print(txt)
  outputfile=file('OutPut_message.txt','wb')
  writeChar(txt,outputfile)
  close(outputfile)
  shell('OutPut_message.txt',wait=FALSE)
}

