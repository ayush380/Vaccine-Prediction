source('pred.R')
bind <- function( x, mhc="HLA-A-02:01", l=9, ic50.threshold=500,
                     quantile.threshold=NULL,
                     include.peptide=TRUE, method="smm" ){
  
  x <- as.character(x)

  ic50 <- smm( x, mhc, l )
  if( is.null( quantile.threshold ) || !is.finite( quantile.threshold ) ||
      quantile.threshold < 0 || quantile.threshold > 1 ){
    i <- which( ic50 < ic50.threshold )
  } else {
    i <- which( ic50 <= quantile( ic50, quantile.threshold ) )
  }
 return(i)
}