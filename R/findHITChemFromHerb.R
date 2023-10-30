findHITChemFromHerb = function(HerbName,type='CH'){
  #HerbName
  #type:type of HerbName-->chinese(CH) or PingYin(PY)
  require(data.table)
  require(stringr)
  if (!exists('HitHerbChem')){
    load('db/HITdatabase.RData')
  }
  if (type=='CH'){
    data=HitHerbChem[herbCH%like%HerbName,.(name,herbCH)]
  }else if (type=='PY'){
    data=HitHerbChem[herbCH==HerbName,.(name,herbCH)]
  }
  return(data)
}

