# y - observed data array of size (n x d_y) where n is the number of observations and d_y is the length of the functional (multivariate) response
# z - array of draws from predictive distribution (n_samp x n x n_y) where n_samp is the number of samples from the predictive distribution
# n_cores - option for parallel computing with mclapply() over n. Depending on the dataset, the overhead associated with this might make it slower than sapply()
energy_score = function(y,z,n_cores=0){
  # The energy score is composed of two terms, a goodness of fit term which is fast to compute
  G = es_goodness_fit(y,z)
  # and a UQ term which is very intensive to compute because it requires a pairwise norm over samples from the predictive distribution
  UQ = es_uncertainty(z,n_cores)
  # Energy score is the difference of these terms ( G - UQ )
  return(list(G=G,UQ=UQ))
}

energy_score_2sample = function(y,z){
  # The energy score is composed of two terms, a goodness of fit term which is fast to compute
  G = es_goodness_fit_2sample(y,z)
  # and a UQ term which is very intensive to compute because it requires a pairwise norm over samples from the predictive distribution
  UQ = es_uncertainty_2sample(z)
  # Energy score is the difference of these terms ( G - UQ )
  ES = G-UQ
  return(list(ES=ES,G=G,UQ=UQ))
}

# this term is fast to compute - but may still benefit from parallel computing like in es_uncertainty
es_goodness_fit = function(y,z){
  d = dim(z)
  n_samp = d[1]
  n = d[2]
  
  colMeans(sapply(1:n,function(i) sapply(1:n_samp, function(k) norm(z[k,i,]-y[i,],type='2'))))
}

# this term is fast to compute - but may still benefit from parallel computing like in es_uncertainty
es_goodness_fit_2sample = function(y,z){
  n_samp = nrow(z)
  n = nrow(y)
  
  2*mean(sapply(1:n,function(i) sapply(1:n_samp, function(k) norm(z[k,]-y[i,],type='2'))))
}
es_uncertainty_2sample = function(z,n_cores){
  d = dim(z)
  n_samp = d[1]
  n = d[2]

  mean(as.matrix(distances::distances(z)))
}

# this term can be slow to compute because of the potentially large pairwise distance matrix
es_uncertainty = function(z,n_cores){
  d = dim(z)
  n_samp = d[1]
  n = d[2]
  
  if(n_cores>0){
    unlist(parallel::mclapply(1:n,function(i) (1/(2*n_samp^2))*sum(as.matrix(distances::distances(z[,i,]))),mc.cores=n_cores))
  } else{
    sapply(1:n,function(i) (1/(2*n_samp^2))*sum(as.matrix(distances::distances(z[,i,]))))
  }
}

# TEST
# set.seed(11)
# n = 4
# d_y = 10
# n_samp = 2
# y_mean = c(0,1,2,3)
# z_sd = c(.01,.01,.1,.2)
# y = matrix(c(rnorm(d_y,y_mean[1],.1),
#              rnorm(d_y,y_mean[2],.1),
#              rnorm(d_y,y_mean[3],.1),
#              rnorm(d_y,y_mean[4],.1)),nrow=n,byrow=T)
# z = array(dim=c(n_samp,n,d_y))
# for(i in 1:n){
#   for(k in 1:d_y){
#     if(i!=2){
#       z[,i,k] = rnorm(n_samp,y[i,k],z_sd[i]) # unbiased z
#     } else{
#       z[,i,k] = rnorm(n_samp,y[i,k]+.1,z_sd[i]) # biased z
#     }
#   }
# }
# matplot(t(y),type='l',col='black',ylim=c(-.5,3.5),lty=1)
# for(i in 1:n){
#   matplot(t(z[,i,]),type='l',col='darkorange',add=T,lty=2:6)
# }
# 
# ES = energy_score(y,z)
# ES$G - ES$UQ
# # 4 test cases ordered from y axis ~ 0-3 on the plot
# # 1. unbiased low noise samples of z - low ES
# # 2. biased low noise samples of z - high ES
# # 3. unbiased with medium noise - mid ES
# # 4. unbiased with high noise - high ES similar to biased low noise case
# ES$G
# # Goodness of fit is about the same for the median noise unbiased and low noise bias
# ES$UQ
# UQ term is a function of the noise in the predictive samples. Note that because the UQ term is subtracted from G, the higher noise samples are favorable for this term,
# although it doesn't necessarily lead to the lowest ES because the goodness of fit term increases.

# Note: I've always found this score to be odd, and a bit counter intuitive with the uncertainty term. It's a negative oriented score, and its a subtraction of the two terms, not an addition,
# so ES favors models with a larger UQ score for fixed Goodness of fit. The UQ term, i think, is a mean of pairwise distances, so it seems that it wants a kind of space
# filling property in the predictive samples, favoring a large pairwise distance.
