package.names = read.csv('~/Desktop/R tools/cran_packages.txt',header = F)

for(i in 1:length(package.names$V1)){
  install.packages(package.names$V1[i])
}

# duqling for interesting computer experiments functions
devtools::install_github("knrumsey/duqling")

# MH adaptive for mcmc
devtools::install_url("https://cran.r-project.org/src/contrib/Archive/MHadaptive/MHadaptive_1.1-8.tar.gz")

# quack for ecp
devtools::install_github("knrumsey/quack")
