y = rgamma(1000,5,1)
library(MASS)
true = fitdistr(y, "gamma", start=list(shape=1))$estimate

y = rpois(10000,10)
true = fitdistr(y, "Poisson")$estimate

# learn shape
log_post = function(param,data){
  # sum(dgamma(y,param,1,log = T))
  sum(dpois(y,param,log=T))
}

prop.sd.log = .01
nsamp = 10000
par = numeric(nsamp)
cur = 10
accept = numeric(nsamp)
for(i in 1:nsamp){
  
  prop = exp(rnorm(1,log(cur),prop.sd.log))
  
  acc = log_post(prop) - log_post(cur)
  if(log(runif(1))<acc){
    par[i] = prop 
    cur = prop
    accept[i] = 1
  } else{
    par[i] = cur
  }
}
mean(accept)
hist(par)
mean(par) # BIASED
true

log_post_jac = function(param){
  # sum(dgamma(y,param,1,log = T)) - log(param)
  sum(dpois(y,param,log = T)) - log(param)
}

prop.sd.log = .01
nsamp = 10000
par = numeric(nsamp)
cur = 10
accept = numeric(nsamp)
for(i in 1:nsamp){
  
  prop = exp(rnorm(1,log(cur),prop.sd.log))
  
  acc = log_post_jac(prop) - log_post_jac(cur)
  if(log(runif(1))<acc){
    par[i] = prop 
    cur = prop
    accept[i] = 1
  } else{
    par[i] = cur
  }
}
# hist(par)
mean(par) # BIASED
