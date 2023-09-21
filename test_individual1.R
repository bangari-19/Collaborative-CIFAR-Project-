# import modules 
source("~/Desktop/CIFAR/R code/sim_modules.R")

# load in contemp and lagged matrices 
data <- read.csv("~/Desktop/CIFAR/220429.FinalMatrices/1.csv")
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
clip_indices <- c(0,1)
clip_mins <- c(0.5, 0.7)
clip_maxs <- c(0.8,1.3)

# clipping outliers
clip_sigma = 2

for (i in 0:50) {
  ts_clipped = clip_outliers(ts, clip_sigma, ampMeasure, debug=T, start, ampContemp, matContemp, covContemp, ampLagged, matLagged, covLagged, covMeasure)
  ts_clipped
}

ts = generate_timeseries(start, steps, ampContemp, matContemp, covContemp, ampLagged, matLagged, covLagged, covMeasure, debug=F)
ts
plot(ts, xlab = "time [days]", ylab = "value", main = 'original timeseries')

ts_clipped = clip_outliers(ts, clip_sigma, ampMeasure, debug=T, start, ampContemp, matContemp, covContemp, ampLagged, matLagged, covLagged, covMeasure)
ts_clipped
plot(ts_clipped, xlab = "time [days]", ylab = "value", main = 'clipped timeseries')

# clipping 
if (require(TSstudio)){

} else { install.packages('https://cran.r-project.org/src/contrib/TSstudio_0.1.7.tar.gz', repos = NULL, type="source",deps=TRUE)}
library(TSstudio)
# generate timeseries

library(ggplot2)
#samples = generate_timeseries(start, steps, ampContemp, matContemp, covContemp, ampLagged, matLagged, covLagged, measureCov, save=F)
#data=as.data.frame(samples)

#plot(samples[,1], type='l')

df <- data.frame(
  days = 1:(steps-1),
  original = ts[, 1],  # only param 1
  clipped = ts_clipped[, 1]
)

ggplot(df, aes(x = days)) +
  geom_line(aes(y = original, color = "Original"), linetype = "solid") +
  geom_line(aes(y = clipped, color = "Clipped"), linetype = "dashed") +
  xlab("time [days]") +
  ylab("value") +
  ggtitle("original vs. clipped timeseries") +
  scale_color_manual(values = c(Original = "blue", Clipped = "red")) +
  scale_linetype_manual(values = c(Original = "solid", Clipped = "dashed")) +
  labs(color = "Timeseries")
