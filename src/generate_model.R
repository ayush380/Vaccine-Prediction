#generating Model for Random Forest
#different mhcs will have different rfs
classifier<-function()
{
  library('caret')
  set.seed(1)
  data<-read.csv('./aminopos.csv')
  data<-data[,c(3,8:17)]
  mhcw<-unique(data$mhc)
  for(i in 3:11)
  {
    data[,i]<-as.factor(data[,i])
  }
  for(i in mhcw)
  {
    data_mhc<-subset(data,mhc==i)
    tryCatch({
    data_index<-createDataPartition(y=as.factor(data_mhc$binder),p=0.80,list=FALSE)
    data_train<-data_mhc[data_index,]
    data_test<-data_mhc[-data_index,]
    clas<-train(y=as.factor(data_train$binder),x=data_train[,3:11],data=data_train,method='rf')
    fname<-paste(gsub("[:*]","",i),"-",9,sep="")
    saveRDS(clas,paste('./classifiers/',fname,'.rds',sep=""))
    },error=function(e){print(e)})
  }

}