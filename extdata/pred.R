#calculate binders using SMM

smmMatrix <- function( mhc="HLA-A-02:01", l=9 ){
  mhc <- gsub( "HLA-([AB])-([0-9])", "HLA-\\1\\2", as.character(mhc) )
  matrixname <- paste(gsub("[:*]","",mhc),"-",l,sep="")
  matrixfile<-paste(matrixname,".txt",sep="")
  
  if( matrixfile == "" ){
    stop( paste( "SMM matrix for MHC",mhc,"and peptide length",l,"not found!") )
  } else {
    M <- as.matrix( read.table(matrixfile, skip=1, nrows=20, row.names=1) )
    colnames(M) <- NULL
    Mc <- as.numeric( read.table(matrixfile, skip=21) )
    return( list(M=M, c=Mc) )
  }
}

binders <- function( x, mhc="HLA-A-02:01", l=9, ic50.threshold=500,
                     quantile.threshold=NULL,
                     include.peptide=TRUE, method="smm" ){
  
  x <- as.character(x)
  start <- seq_len( nchar(x)-l+1 )
  end <- start+l-1
  ic50 <- smm( substring(x,start,end), mhc, l )
  if( is.null( quantile.threshold ) || !is.finite( quantile.threshold ) ||
      quantile.threshold < 0 || quantile.threshold > 1 ){
    i <- which( ic50 < ic50.threshold )
  } else {
    i <- which( ic50 <= quantile( ic50, quantile.threshold ) )
  }
  start <- start[i]
  end <- end[i]
  if( include.peptide ){
    return(data.frame( peptide=substring(x,start,end) ))
  } else {
    return(data.frame( start=start, end=end, ic50=ic50[i]))
  }
}


smm <- function( x=c("SLYNTVATL","SYFPEITHI"), mhc="HLA-A-02:01",
                 output.IC50=TRUE ){
  if( length(mhc)>1 && ( length(mhc) != length(x) ) ){
    stop( "If 'mhc' is a vector, it must have the same length as 'peptides'!" )
  }
  if( length(mhc) == 1 && length(x) > 1 ){
    mhc <- rep.int( mhc, length(x) )
  }
  pred <- function( i ){
    l <-  nchar(x[i])
    M <- smmMatrix( mhc[i], l )
    M$c + sum( sapply( seq_len(l), function(j) M$M[substring(x[i],j,j),j] ) )
  }
  v <- sapply( seq_along(x), pred )
  if( output.IC50 ){
    return( 10^v )
  } else {
    return( v )
  }
}


anchorPositions <- function( mhc="HLA-A-02:01", l=9, k=2 ){
  M <- smmMatrix( mhc, l )$M
  return(sort(tail(order(apply( M, 2, var )),k)))
}


plotBindingMotif <- function( mhc="HLA-A-02:01",
                              l=9, motif.matrix=NULL, width=.5,
                              space=1, col=rainbow(20), main=paste(mhc,", ",l,"-mers",sep=""),
                              ... ){
  if( is.null( motif.matrix ) ){
    M <- 10^(-smmMatrix( mhc, l )$M)
  } else {
    M <- motif.matrix
    l <- ncol(M)
  }
  M <- scale(M,center=FALSE,scale=colSums(M))
  colEntropies <- apply(M, 2, function(x){
    x <- x*log(x)
    x[is.nan(x)] <- 0
    log(20) + sum(x)
  })
  M <- scale(M,center=FALSE,scale=1/colEntropies)
  barplot(M,col=col,width=width,space=space,main=main,...)
  for( i in seq_len( l ) ){
    ypos <- cumsum(M[,i])-M[,i]/2
    showl <- M[,i]>1.5*strheight("M")
    if( sum(showl) > 0 ){
      text( i*width*(1+space)-width/2, ypos[showl], rownames(M)[showl] )
    }
  }
  axis(1,at=(1:l)*(width*(1+space))-width/2,labels=1:l)
}
