#Generate Dateset for applying random forest 
create_Dataset<-function(data)
{
  dimension<-dim(data)
  limit<-dimension[1]
  x<-rep(0,limit)
  data$pos1<-x
  data$pos2<-x
  data$pos3<-x
  data$pos4<-x
  data$pos5<-x
  data$pos6<-x
  data$pos7<-x
  data$pos8<-x
  data$pos9<-x
  
  protein_array<-c('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y')
  for(i in 1:limit)
  {
      sequ<-data[i,4]
      for(j in 1:9)
      {
        y<-substr(sequ,j,j)
        data[i,7+j]<-y
      }
    
  }
  fname<-path('./aminopos.csv')
  write.csv(data,fname)
}