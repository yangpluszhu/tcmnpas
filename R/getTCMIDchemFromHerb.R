getTCMIDchemFromHerb = function(herbName){
  #herbName:herbChineseName
  #TCMIDdatabase:TCMID_data:colnames:c("chemical_name","pub_id","herb_name","herb_CHname")
  require(data.table)
  require(stringr)
  if (!exists('TCMID_data')){
    load('db/TCMIDdatabase.RData')
  }
  data=TCMID_data[herb_CHname%like%herbName,.(chemical_name,pub_id,herb_CHname)]
  if (nrow(data)!=0){
    data$InputName=herbName
  }
  return(data)
}

