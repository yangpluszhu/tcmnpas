P_BasicStat = function(sex=TRUE,age=TRUE,doctor=TRUE,drug=TRUE,diagnosis=TRUE,TCMdiagnosis=TRUE,cost=TRUE,type=TRUE,Manufactor=TRUE){
  require(data.table)
  require(plyr)
  require(stringr)
  require(ggplot2)
  require(DescTools)
  prescription=as.data.table(prescription)
  prescription$yaoId=as.character(prescription$yaoId)
  prescription$id=as.character(prescription$id)
  prescription$admission_number=as.character(prescription$admission_number)
  prescription$patient=as.character(prescription$patient)
  prescription$doctorCode=as.character(prescription$doctorCode)
  prescription_Merge=merge(prescription,mine_data,by.x='yaoId',by.y='idm')
  #############
  prescriptionUnique=prescription_Merge[,.(doctor=unique(doctor),cost=sum(cost)/length(unique(diagnosis))),by='id']##doctor,cost
  prescriptionBYpatient=prescription_Merge[,.(sex=unique(sex),age=unique(age)),by='patient']##sex,age
  prescriptionBYdrug=unique(prescription_Merge,by=c('id','drug'))##drug,type,Manufactor
  prescriptionBYDiag=unique(prescription_Merge,by=c('id','diagnosis'))##diagnosis
  prescriptionBYTCMDiag=unique(prescription_Merge,by=c('id','TCMdiagnosis'))##diagnosis
  name1=c('doctor','cost')
  name2=c('sex','age')
  name3=c('drug','type','Manufactor')
  name4=c('diagnosis')
  name5=c('TCMdiagnosis')
  D1=c(doctor,cost)
  D2=c(sex,age)
  D3=c(drug,type,Manufactor)
  D4=diagnosis
  D5=TCMdiagnosis
  #wrd=GetNewWrd()
  NumPrescription=nrow(unique(prescription,by=c('id','doctor','doctorCode')))
  if (sum(D1)>0){
    #wrd=GetNewWrd()
    #Desc(prescriptionUnique[,name1[D1],with=F],plotit=T,wrd=wrd,maxrows=0.999,ord='desc',digits=3)
	Desc(prescriptionUnique[,name1[D1],with=F],plotit=T,maxrows=0.999,ord='desc',digits=3)
    if (doctor){
      doctorResult=prescriptionUnique[,.(.N),by='doctor']
      doctorResult=doctorResult[,.(doctor,Fre=N,Per=N/NumPrescription)][order(-Fre)]
      doctorResult=doctorResult[,.(doctor,Fre,Per,CumPer=cumsum(Per))]
      write.csv(doctorResult,'Report/Tab13_1开方医生部门汇总.csv',row.names=F)
    }
    if (cost){
      costResult=prescriptionUnique[,.(MinCost=min(cost),Cost25=quantile(cost)[2],Cost50=quantile(cost)[3],Cost75=quantile(cost)[4],MaxCost=max(cost))]
      write.csv(costResult,'Report/Tab13_2处方金额分布结果.csv',row.names = F)
    }
  }
  ##
  if (sum(D2)>0){
    #wrd=GetNewWrd()
    #Desc(prescriptionBYpatient[,name2[D2],with=F],plotit=T,wrd=wrd,maxrows=0.999,ord='desc',digits=3)
	Desc(prescriptionBYpatient[,name2[D2],with=F],plotit=T,maxrows=0.999,ord='desc',digits=3)
    if (sex){
      sexResult=prescriptionBYpatient[,.(.N),by='sex']
      sexResult=sexResult[,.(sex,Fre=N,Per=N/(nrow(prescriptionBYpatient)))][order(-Fre)]
      sexResult=sexResult[,.(sex,Fre,Per,CumPer=cumsum(Per))]
      write.csv(sexResult,'Report/Tab13_3患者性别汇总.csv',row.names=F)
    }
    if (age){
      ageResult=prescriptionBYpatient[,.(MinAge=min(age),Age25=quantile(age)[2],Age50=quantile(age)[3],Age75=quantile(age)[4],MaxAge=max(age))]
      write.csv(ageResult,'Report/Tab13_4患者年龄汇总.csv',row.names = F)
    }
  }
  ##
  if (sum(D3)>0){
    #wrd=GetNewWrd()
    Desc(prescriptionBYdrug[,name3[c(F,type,F)],with=F],plotit=T,maxrows=0.999,ord='desc',digits=3)
    if (drug){
      drugResult=prescriptionBYdrug[,.(.N),by='drug']
      drugResult=drugResult[,.(drug,Fre=N,Per=N/NumPrescription)][order(-Fre)]
      drugResult=drugResult[,.(drug,Fre,Per)]
	  #Desc(prescriptionBYdrug[,name3[D3],with=F],plotit=T,wrd=wrd,maxrows=0.999,ord='desc',digits=3)
      write.csv(drugResult,'Report/Tab13_5药物频次汇总.csv',row.names=F)
    }
    if (type){
      typeResult=prescriptionBYdrug[,.(.N),by='type']
      typeResult=typeResult[,.(type,Fre=N)][order(-Fre)]
      typeResult=typeResult[,.(type,Fre,Per=Fre/sum(Fre))]
      typeResult=typeResult[,.(type,Fre,Per,CumPer=cumsum(Per))]
      write.csv(typeResult,'Report/Tab13_6药物分类汇总.csv',row.names=F)
    }
    if (Manufactor){
      ManufactorResult=prescriptionBYdrug[,.(.N),by='Manufactor']
      ManufactorResult=ManufactorResult[,.(Manufactor,Fre=N)][order(-Fre)]
      ManufactorResult=ManufactorResult[,.(Manufactor,Fre,Per=Fre/sum(Fre))]
      ManufactorResult=ManufactorResult[,.(Manufactor,Fre,Per=Fre/sum(Fre),CumPer=cumsum(Per))]
      write.csv(ManufactorResult,'Report/Tab13_7药物生产厂家汇总.csv',row.names=F)
    }
  }
  ##
  if (sum(D4)>0){
    #wrd=GetNewWrd()
    #Desc(as.character(prescriptionBYDiag[,name4[D4],with=F]),plotit=T,wrd=wrd,maxrows=0.999,ord='desc',digits=3)
    diagnosisResult=prescriptionBYDiag[,.(.N),by='diagnosis']
    diagnosisResult=diagnosisResult[,.(diagnosis,Fre=N,Per=N/NumPrescription)][order(-Fre)]
    diagnosisResult=diagnosisResult[,.(diagnosis,Fre,Per)]
    write.csv(diagnosisResult,'Report/Tab13_8诊断疾病频次汇总.csv',row.names=F)
  }
  ##
  if (sum(D5)>0){
    #wrd=GetNewWrd()
    #Desc(as.character(prescriptionBYTCMDiag[,name5[D4],with=F]),plotit=T,wrd=wrd,maxrows=0.999,ord='desc',digits=3)
    TCMdiagnosisResult=prescriptionBYTCMDiag[,.(.N),by='TCMdiagnosis']
    TCMdiagnosisResult=TCMdiagnosisResult[,.(TCMdiagnosis,Fre=N,Per=N/NumPrescription)][order(-Fre)]
    TCMdiagnosisResult=TCMdiagnosisResult[,.(TCMdiagnosis,Fre,Per)]
    write.csv(TCMdiagnosisResult,'Report/Tab13_9中医诊断频次汇总.csv',row.names=F)
  }
  detach(package:DescTools)
}

