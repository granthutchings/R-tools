unit_xform = function(X,X.min=NULL,X.range=NULL){
  X = as.matrix(X)
  if(is.null(X.min)){
    X.min = apply(X,2,min)
  }
  if(is.null(X.range)){
    Xmax = apply(X,2,max)
    X.range = Xmax - X.min
  }
  if(!isTRUE(all.equal(X.range,rep(0,length(X.range))))){
    Xtrans = t( (t(X) - X.min) / X.range )
  } else{
    # just return X if dummyx
    Xtrans = X
  }
  return(list(orig=X,
              trans=Xtrans,
              min=X.min,
              range=X.range))
  # return(Xtrans)
}