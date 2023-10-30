import_press = function(){
  nameList=readLines('program/tempList.txt')
  nameList=as.character(nameList)
  setListElements("再选择要查询的药物", nameList)
}

