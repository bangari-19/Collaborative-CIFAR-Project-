# import modules 
source("sim_modules.R")

# load in contemp and lagged matrices 
data <- read.csv("/Users/reneehlozek/Code/CIFAR_Network/GIMME/220429.FinalMatrices/1.csv")
size = 6
start <- rep(0, size)
steps = 100
start_index = size+1
stop = length(data)
matContemp <- data[, start_index:stop]
matLagged <- data[, 1:size]

# covariances
covContemp <- matrix(rnorm(size**2), nrow=size, ncol=size)
covLagged <- matrix(rnorm(size**2), nrow=size, ncol=size)
ampContemp = 0.1
ampLagged = 0.1
ampMeasure = 1.0
measureCov <- ampMeasure * diag(size)

# masks
maskContemp =  make_mask(matContemp, contemp=T)
maskLagged = make_mask(matLagged, contemp=F)

# clipping 
# if (require(TSstudio)){

# } else { install.packages('https://cran.r-project.org/src/contrib/TSstudio_0.1.7.tar.gz', repos = NULL, type="source",deps=TRUE)}
# library(TSstudio)
# generate timeseries
library(ggplot2)
samples = generate_timeseries(start, steps, ampContemp, matContemp, covContemp, ampLagged, matLagged, covLagged, measureCov, save=F)
data=as.data.frame(samples)

plot(samples[,1], type='l')