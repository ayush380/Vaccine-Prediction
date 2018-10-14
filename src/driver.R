driver<-function(his="HLA-A*01:01",hcv)
{
  library(caret)
  library(e1071)
  l<-9
  x <- as.character(hcv)
  start <- seq_len( nchar(x)-l+1 )
  end <- start+l-1
  m<-substring(x,start,end)
  source('./src/predict_matrix.R')
  source('./src/Rfpred.R')
  res_mat<-vector()
  res_rf<-vector()
  count<-0
  for(i in m)
  {
    count<-count+1
    print(count)
    print(paste("testing for",i))
    mat_res<-predi(mhc=his,seq=i)
    rf_res<-rfpred(mhc=his,seq=i)
    if(mat_res==1)
    {
    res_mat<-c(res_mat,i)
    }
    if(rf_res==1)
    {
     res_rf<-c(res_rf,i)
    }
  }
  print("WORKING")
  #send res_mat to UI
  source('./src/pred.R')
  res_smm<-binders(x=hcv,mhc=his)
  #send res_smm to UI
  #send res_rf to UI
  print("from smm")
  v<-as.character(res_smm$peptide)
  print("from matrix")
  print(res_mat)
  print('from random forest')
  print(res_rf)
  combined<-c(v,res_rf,res_mat)
  combined<-unique(combined)
  #setwd('../')
  return(combined)
}
Acc<-function(his,hcv)
{
  
  source('./src/Tester.R')
  source('./src/tested2.R')
  source('./src/tester3.R')
  
  test_mat<-0
  test_smm<-0
  test_rf<-0
  #test_mat<-main()
  test_smm<-smm_test(his)
  test_rf<-tester(his)
  
  acc<-data.frame(RF=test_rf ,SMM=test_smm,prob=73.113)
  return(acc)
}