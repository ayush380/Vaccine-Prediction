subset_human<-function(data)
{
  x<-colnames(data)
  k <- 0;
  library(stringr)
  reg<-"HLA+.";
  for(i in 6:29)
  {
    
    if(!str_detect(x[i],"HLA+."))
    {
      data <- data[,-(i-k)];
      k <- k + 1;
    }
  }
  return(data);
  
}


typechanger<-function(data)
{
  tempdata<-data;
  tempdata$HLA.A.0201 <- as.character(tempdata$HLA.A.0201)
  for(i in 1:180)
  {
    
    if(tempdata[i,6] == "<69444.444444")
    {
      tempdata[i,6] = "69444.444444"
    }
  }
  tempdata$HLA.A.0201 <- as.numeric(tempdata$HLA.A.0201)
  
  tempdata$HLA.A.3201 <- as.character(tempdata$HLA.A.3201)
  for(i in 1:180)
  {
    
    if(tempdata[i,9] == "<72949.62")
    {
      tempdata[i,9] = "72949.62"
    }
  }
  tempdata$HLA.A.3201 <- as.numeric(tempdata$HLA.A.3201)
  
  
  tempdata$HLA.B.5801 <- as.character(tempdata$HLA.B.5801)
  for(i in 1:180)
  {
    
    if(tempdata[i,11] == "<78125")
    {
      tempdata[i,11] = "78125"
    }
  }
  tempdata$HLA.B.5801 <- as.numeric(tempdata$HLA.B.5801)
  
  
  tempdata$HLA.B.5802 <- as.character(tempdata$HLA.B.5802)
  for(i in 1:180)
  {
    
    if(tempdata[i,12] == "<78125")
    {
      tempdata[i,12] = "78125"
    }
  }
  tempdata$HLA.B.5802 <- as.numeric(tempdata$HLA.B.5802)
  
  tempdata$HLA.B.0702 <- as.character(tempdata$HLA.B.0702)
  for(i in 1:180)
  {
    
    if(tempdata[i,13] == "<77464.788732")
    {
      tempdata[i,13] = "77464.788732"
    }
  }
  tempdata$HLA.B.0702 <- as.numeric(tempdata$HLA.B.0702)
  
  
  tempdata$HLA.B.3501 <- as.character(tempdata$HLA.B.3501)
  for(i in 1:180)
  {
    
    if(tempdata[i,14] == "<77419.354839")
    {
      tempdata[i,14] = "77419.354839"
    }
  }
  tempdata$HLA.B.3501 <- as.numeric(tempdata$HLA.B.3501)
  
  tempdata$HLA.B.5301 <- as.character(tempdata$HLA.B.5301)
  for(i in 1:180)
  {
    
    if(tempdata[i,15] == "<78151.260504")
    {
      tempdata[i,15] = "78151.260504"
    }
  }
  tempdata$HLA.B.5301 <- as.numeric(tempdata$HLA.B.5301)
  
  
  tempdata$HLA.B.5401 <- as.character(tempdata$HLA.B.5401)
  for(i in 1:180)
  {
    
    if(tempdata[i,16] == "<78125")
    {
      tempdata[i,16] = "78125"
    }
  }
  tempdata$HLA.B.5401 <- as.numeric(tempdata$HLA.B.5401)
  
  
  tempdata$HLA.B.5101 <- as.character(tempdata$HLA.B.5101)
  for(i in 1:180)
  {
    
    if(tempdata[i,17] == "<77464.788732")
    {
      tempdata[i,17] = "77464.788732"
    }
  }
  tempdata$HLA.B.5101 <- as.numeric(tempdata$HLA.B.5101)
  
  
  tempdata$HLA.A.6802 <- as.character(tempdata$HLA.A.6802)
  for(i in 1:180)
  {
    
    if(tempdata[i,19] == "<76923.076923")
    {
      tempdata[i,19] = "76923.076923"
    }
  }
  tempdata$HLA.A.6802 <- as.numeric(tempdata$HLA.A.6802)
  
  tempdata$HLA.B.2705 <- as.character(tempdata$HLA.B.2705)
  for(i in 1:180)
  {
    
    if(tempdata[i,20] == "<77906.976744")
    {
      tempdata[i,20] = "77906.976744"
    }
  }
  tempdata$HLA.B.2705 <- as.numeric(tempdata$HLA.B.2705)
  
  data<-tempdata;
  return(data);
  
}


energy_formula<-function(data,HLA,pos,value)
{
  tempdata<-data
  total<-0
  x<-as.integer(pos/20)+1;
  for(i in 1:180)
  {
    if(tempdata[i,5]==x)
    {
      total<- total + log10(tempdata[i,HLA]);
    }
  }
  total<-total/20;
  return(log10(value)-total)
}


evaluate_affinity<-function(data)
{
  tempdata<-data;
  energy<-matrix(rep(0,(180*15)),180,15)
  for(i in 6:20)
  {
    for(j in 1:180)
    {
      energy[j,i-5]<-energy_formula(data,i,j,tempdata[j,i]);
    }
  }
  return(energy);
}



summation<-function(i,j,energymat,data,arr)
{
  total<-0;
  for(k in 1:15)
  {
    x<-0;
    for(l in 1:9)
    {
      total<- total + ((energymat[i+(x*20),k]-arr[i%%20])*(energymat[j+(x*20),k]-arr[j%%20]));
      x<-x+1;
    }
  }
  print(total/(15*9));
}


covariance<-function(data,energymat)
{
  arr<-c(rep(0,20));
  cov<-matrix(rep(0,400),20,20)
  for(i in 1:180)
  {
    total<-0;
    for(j in 6:20)
    {
      total<-total + energymat[i,j-5];
    }
    if(i%%20==0)
    {
      arr[20]<-arr[20] + total;
    }
    else
    {
      arr[(i%%20)]<-arr[i%%20] + total;
    }
  }
  
  
  for(i in 1:20)
  {
    arr[i]=arr[i]/135;
  }
  
  
  for(i in 1:20)
  {
    for(j in 1:20 )
    {
      to<-0;
      for(k in 1:15)
      {
        
        for(l in 1:9)
        {
          
          if(i%%20==0 & j%%20!=0)
          {
            pro1<-(energymat[i+((l-1)*20),k]-arr[20]);
            
            pro2<-(energymat[j+((l-1)*20),k]-arr[j%%20]);
            
            to<- to + (pro1*pro2);
          }
          
          else if(j%%20==0 & i%%20!=0)
          {
            pro1<-(energymat[i+((l-1)*20),k]-arr[i%%20]);
            
            pro2<-(energymat[j+((l-1)*20),k]-arr[20]);
            
            to<- to + (pro1*pro2);
          }
          
          else if(i%%20==0 & j%%20==0)
          {
            pro1<-(energymat[i+((l-1)*20),k]-arr[20]);
            
            pro2<-(energymat[j+((l-1)*20),k]-arr[20]);
            
            to<- to + (pro1*pro2);
          }
          
          else
          {
            pro1<-(energymat[i+((l-1)*20),k]-arr[i%%20]);
            
            pro2<-(energymat[j+((l-1)*20),k]-arr[j%%20]);
            
            to<- to + (pro1*pro2);
          }
        }
      }
      cov[i,j]<-to/135;
    }
  }
  
  return(cov);
}
create_smm<-function(data){
  
}

main<-function()
{
  data<-read.csv("../dataset.txt",sep="\t");
  
  data1<-subset_human(data);
  
  data2<-typechanger(data1);
  
  energymat<-evaluate_affinity(data2);
  
  covmat<-covariance(data2,energymat);
  
  return_list<-list("subse"=data1,"typechanger"=data2,"energymatrix"=energymat,"covariance"=covmat)
  return(return_list);
}
