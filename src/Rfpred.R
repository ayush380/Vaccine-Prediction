#Random forest Pridiction

rfpred<-function(mhc2="HLA-A*01:01",seq)
{
  library(randomForest)
  data<-read.csv("./aminopos.csv")
  data<-data[,c(3,8:17)]
  
  data[1,1]<-mhc2
  fname<-paste(gsub("[:*]","",mhc2),"-",9,sep="")
  tryCatch({
    clas<-readRDS(paste("./classifiers/",fname,".rds",sep=""))
  },error=function(e){
    print(e)
    return(0)
  })
  #data<-data.frame(pos1=0,pos2=0,pos3=0,pos4=0,pos5=0,pos6=0,pos7=0,pos8=0,pos9=0)
  
  
  protein_array<-c('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y')
  #factor(protein_array,levels=c('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T', 'V','W','Y'))
    for(j in 1:9)
    {
      y<-substr(seq,j,j)
      data[1,j+2]<-y
      #data[1,j]<-as.character(factor((data[1,j]),levels=c('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T', 'V','W','Y')))   
    }
  
  
  for(i in 3:11)
  {
    data[,i]<-as.factor(data[,i])
  }
  data<-data[1,]
  #print(data)
  res<-predict(clas,data)
  #print(res)
  return(res)
}