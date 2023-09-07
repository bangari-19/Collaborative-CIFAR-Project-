# SIM_MODULES

plot_matrix <- function(mat, save=FALSE, name=NULL) {
  if(require('plot.matrix')){}
  else {install.packages('https://cran.r-project.org/src/contrib/plot.matrix_1.6.2.tar.gz',repos = NULL, type="source",deps=TRUE)}

  library(plot.matrix)
  plot_matrix(mat)
  if (save) {
    plot(mat)  # create plot
    png(file = name)  # save plot as png
    dev.off()  # close plot 
  }
}

save_matrix <- function(mat, name) {
  write.csv(mat, file = name)
}


make_mask <- function(mat, contemp) {
  "If lagged then mask all zero entries 
  If contemp then mask diagonal with zeros in addition"
  size = nrow(mat)  # number of rows
  mask <- matrix(1, size, size)
  if (contemp) {
    diag(mask) <- 0  # fill diagonal of mask with zeros
  }
  
  inds = which(mat == 0, arr.ind = T)  # find indices where mat == 0 in (row, col) format
  mask[inds] = 0  # make mask 0 at those indices 
  
  return(mask)
}


coeff_draw_from_cov <- function(amplitude, mat, cov, mask) {
  size <- nrow(cov)
  sample_step <- cov %*% (amplitude * matrix(rnorm(size * size), nrow = size, ncol = size))
  stepmat <- mat + sample_step
  noisy_mat <- stepmat * mask
  return(noisy_mat)
}


generate_timeseries <- function(start, len, contempamp, contempmat, contempcov, lagamp, lagmat, lagcov, measurecov, save=FALSE) {
  size <- length(start)  # length of 'start' vector
  samples <- matrix(0, nrow = len, ncol = size)
  samples[1,] <- start  # day 1 data (index start at 1 in R)
  lagmask <- make_mask(lagmat, contemp=FALSE)
  contempmask <- make_mask(contempmat, contemp=TRUE)
  
  noisy_lagmat <- coeff_draw_from_cov(lagamp, lagmat, lagcov, lagmask)
  noisy_contempmat <- coeff_draw_from_cov(contempamp, contempmat, contempcov, contempmask)
  
  if (save) {
    save_matrix(noisy_lagmat, 'noisy_lagmat.csv')
    plot_matrix(noisy_lagmat, TRUE, 'noisy_lagmat.png')
    save_matrix(noisy_contempmat, 'noisy_contempmat.csv')
    plot_matrix(noisy_contempmat, TRUE, 'noisy_contempmat.png')
  }
  
if(require('mvtnorm')){}
else {install.packages('https://cran.r-project.org/src/contrib/mvtnorm_1.2-3.tar.gz', repos = NULL, type="source",deps=TRUE)}
library(mvtnorm)
# print(noisy_contempmat)
# print(class(solve(diag(size) - noisy_contempmat)))
# print(class(start))
# print(class(samples))
# print(class(as.matrix(noisy_contempmat)))
# print(class(noisy_lagmat))
  for (i in 2:len) {
    samples[i,] <- samples[i-1,] %*% as.matrix(noisy_lagmat) %*% solve(diag(size) - as.matrix(noisy_contempmat)) + start
    samples[i,] <- samples[i,] + rmvnorm(1, rep(0, size), measurecov)  
    # rmvnorm generates one random sample from the multivar norm dist with mean 0 and covariance 'measurecov'
  }
  
  samples <- samples[-1,] - start
  return(samples)
}


clip_timeseries <- function(timeseries, indices_to_clip, min_vec, max_vec) {
  
  ts <- timeseries  # copy of timeseries
  
  for (i in indices_to_clip) {
    min_val <- min_vec[i]
    max_val <- max_vec[i]
    
    min_inds <- which(timeseries[, i] < min_val)  # indices where the ith column is less than min_val
    ts[min_inds, i] <- min_val  # at those indices insert min_val
    
    max_inds <- which(timeseries[, i] > max_val)
    ts[max_inds, i] <- max_val
  }
  
  return(ts)
}


# #TEST FOR INDIVIDUAL 1 
# data <- read.csv("/Users/reneehlozek/Code/CIFAR_Network/GIMME/220429.FinalMatrices/1.csv")
# len = 10  # number of days for timeseries 
# start <- rep(0, 6)  # starting value for day 1
# size = 6
# start_ind = size + 1
# stop_ind = length(data)
# matContemp <- data[, start_ind : stop_ind]
# matLagged <- data[, 1 : size]
# matContemp <- as.matrix(matContemp)  # convert from data format to matrix format
# matLagged <- as.matrix(matLagged)
# matContemp
# matLagged

# plot_matrix(matLagged)
# plot_matrix(matContemp)

# maskContemp = make_mask(matContemp, TRUE)
# maskContemp
# maskLagged = make_mask(matLagged, FALSE)
# maskLagged

# covContemp = matrix(rnorm(size * size), nrow=size, ncol=size)
# covContemp
# covLagged = matrix(rnorm(size * size), nrow=size, ncol=size)
# covLagged

# ampContemp = 0.1
# ampLagged = 0.1
# ampMeasure = 1
# covMeasure = ampMeasure * diag(size)
# covMeasure

# ts = generate_timeseries(start, len, ampContemp, matContemp, covContemp, ampLagged, matLagged, covLagged, covMeasure)
# ts  # variables on column and data for each day on rows 

# clip_indices <- c(0,1)
# clip_mins <- c(0.5, 0.7)
# clip_maxs <- c(0.8,1.3)

# ts_clip = clip_timeseries(ts, clip_indices, clip_mins, clip_maxs)
# ts_clip
