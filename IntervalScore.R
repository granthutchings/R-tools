interval_score = function(y,alpha=.05,y.samp=NULL,y.conf.int=NULL,terms=F){
  if(!is.null(y.samp)){
    n = dim(y.samp)[3]
    n.y = dim(y.samp)[1]
    y.conf.int = apply(y.samp,c(1,3),quantile,c(.025,.975))
  } else{
    n = dim(y.conf.int)[3]
    n.y = dim(y.conf.int)[2]
  }
  # make sure the field data shape matches the samples shape
  if(!all(dim(y)==c(n.y,n)))
    dim(y) = c(n.y,n)
  
  width = array(dim=c(n.y,n))
  penalty = array(dim=c(n.y,n))
  i_score = array(dim=c(n.y,n))
  
  for(k in 1:n){ # loop over pred locations
    # i_score = (u-l)+2/alpha(l-x)I(x<l)+2/alpha(x-u)I(x>u)
    width[,k] = (y.conf.int[2,,k]-y.conf.int[1,,k])
    penalty[,k] = (2/alpha)*(y.conf.int[1,,k] - y[,k])*as.integer(y[,k]<y.conf.int[1,,k]) +
      (2/alpha)*(y[,k] - y.conf.int[2,,k])*as.integer(y[,k]>y.conf.int[2,,k])
    i_score[,k] = width[,k] + penalty[,k]
  }
  if(terms){
    return(list(i_score=i_score,width=width,penalty=penalty))
  } else{
    return(i_score)
  }
}
