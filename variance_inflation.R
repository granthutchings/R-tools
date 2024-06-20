variance_inflation = function(x, delta = .5, n_dens_est = 1000, bounds=c(-Inf,Inf)){
  d = density(x, from = bounds[1], to = bounds[2], n=n_dens_est)
  d$y = d$y^delta # power up density
  d$bw
  cdf = cumsum(d$y)/sum(d$y)
  u = runif(length(x))
  samples = approx(cdf, d$x, xout = u)$y
  return(samples)
}

xx4 = variance_inflation(xx,bounds=c(0,1))
hist(xx4, breaks=seq(0, 1, by=0.01), xlim=c(0, 1), freq=FALSE)
