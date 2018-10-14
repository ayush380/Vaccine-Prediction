tester<-function(mhc2)
{
  library(caret)
  set.seed(1)
  data<-read.csv('./aminopos.csv')
  data<-data[,c(3,8:17)]
  data_mhc<-subset(data,mhc==mhc2)
  for(i in 3:11)
  {
    data_mhc[,i]<-as.factor(data_mhc[,i])
  }
  acc<-0
  oacc<-0
  count<-0
  print(1)
    tryCatch({
      index<-createDataPartition(y=as.factor(data_mhc$binder),p=0.8,list=FALSE)
    data_train<-data_mhc[index,]
    data_test<-data_mhc[-index,]
    fname<-paste(gsub("[:*]","",mhc2),"-",9,sep="")
    clas<-readRDS(paste("./classifiers/",fname,".rds",sep=""))
    print(paste("predicting for",mhc2))
    p<-predict(clas,data_test)
    c<-confusionMatrix(p,as.factor(data_test$binder))
    acc<-(c$overall[1]*100)
    ret<-acc
   
    return(ret)
    
    #oacc<-oacc+acc
    #count<-count+1
    },error=function(e){print(e)})
  #print(paste("overall accuracy:",oacc/count))
  
}