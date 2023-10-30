make_net = function(database){
  library(data.table)
  library(plyr)
  library(stringr)
  library(igraph)
  ##database-->data.table:colnames-->c("drug1","drug2","id","influence","actionId","measures","level","literatureLevel", "mechanism","literature","assess")
  #####
  database=as.data.table(database)
  #setwd('E:\\DATABASE\\数据库\\151027\\prepared_data')
  #database=fread('drug_interact_151212.csv',header=T,stringsAsFactors=F)
  #setkey(database,drug1,drug2)
  #database=unique(database)
  net=graph_from_edgelist(as.matrix(database[,.(drug1,drug2)]),directed=F)
  E(net)$id=database[,actionId]
  #net=simplify(net,remove.multiple =TRUE)
  #Drug_data=fread('drugs_prepared.csv',head=T,stringsAsFactors=F)
  return(net)
}

