lhs_empirco = function(data,nsample){
  # % s=lhs_empirco(data,nsample)
  # % perform lhs on multivariate empirical distribution
  # % with correlation
  # % Input:
  #   %   data    : data matrix (ndata,nvar)
  # %   nsample : no. of samples
  # % Output:
  #   %   s       : random sample (nsample,nvar)
  # %   Budiman (2003)
  m = nrow(data)
  nvar = ncol(data)
  corr = cor(data)
  rc = rank_corr(corr,nsample) # % induce correlation
  s = matrix(nrow=nsample,ncol=nvar)
  for(j in 1:nvar){
    r=rc[,j]
    # % draw random no.
    u=runif(nsample)
    # % calc. percentile
    p=(r-u)/nsample
    # % inverse from empirical distribution
    s[,j]=quantile(data[,j],p)
  }
  return(s)
}


rank_corr = function(corr,nsample){
  # % rc=rank_corr(corr,nsample)
  # % induce rank correlation
  # % method of Iman & Conover
  # % Iman, R. L., and W. J. Conover. 1982. A Distribution-free Approach to
  # % Inducing Rank Correlation Among Input Variables.
  # % Communications in Statistics B 11:311-334.
  # % Input:
  #   %   corr    : correlation matrix of the variables (nvar,nvar)
  # %   nsample : no. of samples
  # % Output:
  #   %   rc       : rank (nsample,nvar)
  # %   Budiman (2004)
  nvar = ncol(corr)
  # % induce data with correlation
  # R = qnorm(lhs::maximinLHS(nsample,nvar))
  R = qnorm(tgp::lhs(nsample,matrix(rep(c(0,1),nvar),ncol=2,byrow=T)))
  T = cor(R)
  P = t(chol(corr))
  Q = t(chol(T + 1e-8*diag(nvar)))
  # % use modified cholesky for corr matrix that is not quite positive definite
  # %[L,D,E]=mchol(corr)  
  # %P=L*sqrt(D)
  # %[L,D,E]=mchol(T)  
  # %Q=L*sqrt(D)
  S = P %*% forwardsolve(Q,diag(nvar))
  RB= R %*% t(S)
  rc = array(dim=dim(RB))
  for(j in 1:nvar){    
    # % rank RB
    i = order(RB[,j])                                                            
    rc[i,j] = 1:nsample
    # RB[,j] = order(RB[,j])
  }
  return(rc)
}
