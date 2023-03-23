# install and load required packages
install.packages('plot.matrix')
install.packages('matlib')
install.packages('MASS')
install.packages('stats')

library('plot.matrix')
library('matlib')
library('MASS')
library('stats')


plot_matrix <- function(mat, save=FALSE, name) {
  plot(mat)
  if (save) {
    fname = paste0(name, '.png')
    png(file = fname)
    plot(mat)
    dev.off()
  }
}


save_matrix <- function(mat, name) {
  fname = paste0(name, '.csv')
  write.csv(mat, file = fname)
}


make_mask <- function(mat, contemp) {
  size = dim(mat)[1]  # returns number of rows 
  mask <- matrix(1, size, size)  
  if (contemp) {
    diag = rep(0, size)  # replicate element 0 of vector 'size' times
    diag(mask) <- diag  # fill diagonal with zeros
  } else {
    indices = which(mat == 0, arr.ind = T)  # arr.ind outputs indices in row, col format
    mask[indices] = 0  # indices where matrix = 0, make mask = 0 
  }
  return(mask)
}


# Draw a coefficient matrix from the mean matrix and the covariance
coeff_draw_from_cov <- function(amplitude, mat, cov, mask) {
  size <- dim(mat)[1]  
  random_mat <- amplitude * matrix(rnorm(size**2), nrow = size)
  sample_step <- cov %*% random_mat  # generate noise
  step_mat <- mat + sample_step  # add noise 
  noisy_mat <- step_mat %*% mask  # mask out terms 
  return(noisy_mat)
}


# Generate timeseries
generate_timeseries <- function(start, len, contempamp, contempmat, contempcov, lagamp, lagmat, lagcov, measurecov, save=F){
  size = length(start)  # 'start' is a vector 
  samples = matrix(0, nrow=len, ncol=size)  # np.zeros((len, size))
  samples[1, ] = start  # day 1 data, index starts at 1 in R
  lagmask = make_mask(lagmat, contemp=FALSE)
  contempmask = make_mask(contempmat, contemp=TRUE)
  print(lagmask)
  print(contempmask)
  
  noisy_lagmat = coeff_draw_from_cov(lagamp, lagmat, lagcov, lagmask)
  noisy_contempmat = coeff_draw_from_cov(contempamp, contempmat, contempcov, contempmask)
  print(noisy_lagmat)
  print(noisy_contempmat)
  
  if (save) {
    save_matrix(noisy_lagmat, 'noisy_lagmat')
    plot_matrix(noisy_lagmat, TRUE, 'noisy_lagmat')
    save_matrix(noisy_contempmat, 'noisy_contempmat')
    plot_matrix(noisy_contempmat, TRUE, 'noisy_contempmat')
  }
  
  days <- c(2:len)  # days = (2, 3, ..., len)
  for (i in days) {
    Ylag_lagmat = samples[i-1, ] %*% noisy_lagmat 
    inverse = inv(diag(size) - noisy_contempmat)  # diag(size) gives identity 
    print(Ylag_lagmat)
    print(inverse)
    samples[i, ] = (Ylag_lagmat %*% inverse) 
    print(samples[i, ])
    measurement_noise <- mvrnorm(n = 1, mu = matrix(0, size, size), Sigma = measurecov)  # generates random samples from a multivariate normal distribution
    print(measurement_noise)
    samples[i, ] = samples[i, ] + measurement_noise
    print(samples[i, ])
  }
  #samples = samples 
  samples = start[-1, ]  # get all data except for day 1 data  
  
  return (samples)
}

# clip_timeseries
clip_timeseries <- function(timeseries, indices_to_clip, min_vec, max_vec){
  # copy timeseries matrix 
  ts = matrix(timeseries, nrow=nrow(timeseries), ncol=ncol(timeseries))
  for (i in indices_to_clip) {
    min = min_vec[i]
    max = max_vec[i]
    min_inds = which(timeseries[, i] < min, arr.ind = T)
    ts[min_inds, i] = min
    max_inds = which(timeseries[, i] > max, arr.ind = T)
    ts[max_inds, i] = max
  }
  
  return (ts)
}

# clip_outliers

# test for individual 1  
data <- read.csv("~/Desktop/Astro/CIFAR Project/220429.FinalMatrices/1.csv")
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

# generate timeseries
samples = generate_timeseries(start, steps, ampContemp, matContemp, covContemp, ampLagged, matLagged, covLagged, measureCov, save=F)

  
