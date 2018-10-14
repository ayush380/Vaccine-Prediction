create_matrix<-function(data_train)
  {
  #data<-read.csv(file="bdata.20130222.mhci.txt",sep="\t")
  #data<-subset(data,peptide_length==9)
  #data$mhc<-as.character(data$mhc)
  #binders_data<-subset(data,meas<=500)
  #binders_data['binder']=1
  #non_binders_data<-subset(data,meas>500)
  #non_binders_data['binder']=0
  #data2<-rbind(binders_data,non_binders_data)
  #data_index<-createDataPartition(y=data2$binder,p=0.80,list=FALSE)
  #data_train<-data2[data_index,]
  #data_test<-data2[-data_index,]
  mhcs<-unique(data_train$mhc)
  protein_array<-c('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y')
  #test_dataframe<-data.frame()
  for (i in mhcs)
  {
    matrixname<-" "
    #print(i)
    
    data_mhc<-subset(data_train,mhc==i)
    #print(data_mhc)
    dimension<-dim(data_mhc)
    #print(dimension)
    mat<-matrix(rep(0,180),nrow = 20,ncol = 9)
    limit<-as.numeric(dimension[1])
    for(j in 1:limit)
    {
      #print(j)
      if(data_mhc[j,7]==1)
      {
        #print("sequ is")
        sequ<-as.character(data_mhc[j,4])
        #print(sequ)
        for(pos in 1:9)
        {
          posi<-0
         # print("x is")
          x<-substr(sequ,pos,pos)
         # print(x)
          for(k in 1:20)
          {
            if(protein_array[k]==x)
            {
              posi<-k
            }
          }
          #print(posi)
          mat[posi,pos]<-mat[posi,pos]+1
        }      
      }
     

    }
    mat<-(mat/dimension[1])
    #print(mat)
    #path<-"C:/Users/Ayush Srivastava/Desktop/Minor/matrices/"
    
    #mhcw <- gsub( "HLA-([AB])-([0-9])", "HLA-\\1\\2", as.character(i) )
    matrixname <- paste(gsub("[:*]","",i),"-",9,sep="")
  #  matrixname
   # getwd()
    
    matrixname<-paste(matrixname,".txt",sep="")
    #setwd("C:/Users/Ayush Srivastava/Desktop/Minor/matrices")
    write.table(mat, file=paste("./matrices/",matrixname,sep=""), row.names=FALSE, col.names=FALSE)
   
  }
  
}