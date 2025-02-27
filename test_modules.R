source("sim_modules.R")
# --------------------------TEST FOR ONE INDIVIDUAL------------------------------- #
data <- read.csv("~/Desktop/CIFAR/220429.FinalMatrices/1.csv")
size = 6
start <- rep(0, size)
steps = 100
start_index = size+1
stop = length(data)
matContemp <- data[, start_index:stop]
matLagged <- data[, 1:size]
covContemp <- matrix(rnorm(size**2), nrow=size, ncol=size)
covLagged <- matrix(rnorm(size**2), nrow=size, ncol=size)
ampContemp = 0.01
ampLagged = 0.01
ampMeasure = 1.0
measureCov <- ampMeasure * diag(size)
maskContemp =  make_mask(matContemp, contemp=T)
maskLagged = make_mask(matLagged, contemp=F)
clip_indices <- c(0,1)
clip_mins <- c(0.5, 0.7)
clip_maxs <- c(0.8,1.3)
clip_sigma = 3  # smaller clip_sigma --> more clipping

ts = generate_timeseries(start, steps, ampContemp, matContemp, covContemp, ampLagged, matLagged, covLagged, covMeasure, debug=F)
ts_clipped = clip_outliers(ts, clip_sigma, ampMeasure, debug=F, start, ampContemp, matContemp, covContemp, ampLagged, matLagged, covLagged, covMeasure)

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
  ggtitle("Original vs. Clipped Timeseries") +
  scale_color_manual(values = c(Original = "blue", Clipped = "red")) +
  scale_linetype_manual(values = c(Original = "solid", Clipped = "dashed")) +
  labs(color = "Timeseries") +
  theme_minimal()
