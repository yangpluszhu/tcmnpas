FiltedPstat = function(){
  require(data.table)
  prescription2=as.data.table(prescription)
  NumYao=length(unique(prescription2$drug))
  NumF=length(unique(prescription2$id))
  NumDoctor=length(unique(prescription2$doctor))
  Numpatient=length(unique(prescription2$patient))
  result=data.frame(处方数=NumF,药物数=NumYao,患者数=Numpatient,开方医生部门数=NumDoctor)
  return(result)
}

