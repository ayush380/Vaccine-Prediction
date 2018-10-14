#prediction using probability method 
predi<-function(mhc="HLA-A*01:01",seq)
{
  protein_array<-c('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y')
  
  matrixname <- paste(gsub("[:*]","",mhc),"-",9,sep="")
  filename<-paste("./matrices/",matrixname,".txt",sep="")
  mat<-as.matrix(read.table(file=filename,nrows = 20))
  #print(mat)
  res<-0
  for(pos in 1:9)
  {
    x<-substr(seq,pos,pos)
    posi<-0
    for(k in 1:20)
    {
      if(protein_array[k]==x)
      {
        posi<-k
      }
    }
    res<-res+mat[posi,pos]
  }
  if(res>0.5)
    return(1)
  else
    return(0)
}