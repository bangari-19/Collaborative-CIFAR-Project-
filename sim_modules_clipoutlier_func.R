# SIM_MODULES

plot_matrix <- function(mat, save=FALSE, name=NULL) {
  library('plot.matrix')
  plot(mat)
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
  
  library('mvtnorm')
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


clip_outliers <- function(timeseries, sigma, measure_amp, debug, start, contempamp, contempmat, contempcov, lagamp, lagmat, lagcov, measurecov) {

  ts <- timeseries
  
  # Iterate through all columns of timeseries
  for (i in 1:ncol(timeseries)) {
    
    min_val <- -sigma * sqrt(measure_amp)  
    max_val <- sigma * sqrt(measure_amp)  
    min_inds <- which(timeseries[, i] < min_val)  # indices where the column values are less than min_val
    
    if (debug && length(min_inds) > 0) {
      cat("bad vals min", timeseries[min_inds, i], "inds", min_inds, "\n")
    }
    
    countmax <- 2
    
    # Loop over the min outlier indices
    for (minind in min_inds) {
      count <- 0
      tmp <- generate_timeseries(ts[minind - 1, ], 3, contempamp, contempmat, contempcov, lagamp, lagmat, lagcov, measurecov, save = FALSE)[1]
      
      if (debug) {
        cat("param", i, tmp[i], "tmp[i]", "min =", min_val, "count", count, "minind", minind, "\n")
      }
      
      # Generate new timeseries for that outlier index until outlier within bound or count > countmax
      while (ts[minind, i] < min_val && count < countmax) {
        tmp <- generate_timeseries(ts[minind - 1, ], 3, contempamp, contempmat, contempcov, lagamp, lagmat, lagcov, measurecov, save = FALSE)[1]
        
        if (debug) {
          cat("iteration", count, "for param", i, "and ind", minind, "in min", "\n")
          print(ts[minind - 1:minind + 2, i])
          cat("param", i, "before", "\n")
        }
        
        # Update the value in ts
        ts[minind, i] <- tmp[i]
        count <- count + 1
        
        if (debug) {
          print(ts[minind - 1:minind + 2, i])
          cat("after", count, countmax, "\n")
          cat("-----\n")
        }
        
        if (count == (countmax - 1)) {
          ts[minind, i] <- min_val * (1 - 0.01 * runif(1))
          if (debug) {
            cat(ts[minind, i], "failing due to count", "\n")
          }
        }
      }
      
      if (debug) {
        cat("while loop complete", "\n")
        cat(ts[minind, i], "new value for param", i, "\n")
        cat("=======\n")
      }
    }
    
    
    max_inds <- which(timeseries[, i] > max_val)
    
    if (debug && length(max_inds) > 0) {
      cat("bad vals max", timeseries[max_inds, i], "inds", max_inds, "\n")
    }
    
    for (maxind in max_inds) {
      count <- 0
      tmp <- generate_timeseries(ts[maxind - 1, ], 3, contempamp, contempmat, contempcov, lagamp, lagmat, lagcov, measurecov, save = FALSE)[1]
      
      if (debug) {
        cat("param", i, tmp[i], "tmp[i]", "max =", max_val, "count =", count, "maxind", maxind, "\n")
      }
      
      while (ts[maxind, i] > max_val && count < countmax) {
        
        if (debug) {
          cat("iteration", count, "for param", i, "and ind", maxind, "in max", "\n")
          print(ts[maxind - 1:maxind + 2, i])
          cat("ts before", "\n")
          cat(tmp[i], "tmp[i]", "max =", max_val, "count =", count, "\n")
        }
        
        ts[maxind, i] <- tmp[i]
        count <- count + 1
        
        if (debug) {
          print(ts[maxind - 1:maxind + 2, i])
          cat("after", count, countmax, "\n")
          cat("-----\n")
        }
        
        if (count == (countmax - 1)) {
          ts[maxind, i] <- max_val * (1 - 0.01 * runif(1))  
          if (debug) {
            cat(ts[maxind, i], "failing due to count", "\n")
          }
        }
      }
      
      if (debug) {
        cat("while loop complete", "\n")
        cat(ts[maxind, i], "new value for param", i, "\n")
        cat("=======\n")
      }
    }
    
    if (debug) {
      if (length(min_inds) == 0 && length(max_inds) == 0) {
        cat("no clipping needed for param", i, "\n") 
      } else {
        cat("mean for param", i, "before clip:", mean(timeseries[, i]), ", mean after clip:", mean(ts[, i]), "\n")
      }
    }
  }
  
  return(ts)
}
































  
  
#TEST FOR INDIVIDUAL 1 
data <- read.csv("~/Desktop/CIFAR/220429.FinalMatrices/1.csv")
clip_sigma = 1
debug = FALSE
len = 10  # number of days for timeseries 
start <- rep(0, 6)  # starting value for day 1
size = 6
start_ind = size + 1
stop_ind = length(data)
matContemp <- data[, start_ind : stop_ind]
matLagged <- data[, 1 : size]
matContemp <- as.matrix(matContemp)  # convert from data format to matrix format
matLagged <- as.matrix(matLagged)
matContemp
matLagged

plot_matrix(matLagged)
plot_matrix(matContemp)

maskContemp = make_mask(matContemp, TRUE)
maskContemp
maskLagged = make_mask(matLagged, FALSE)
maskLagged

covContemp = matrix(rnorm(size * size), nrow=size, ncol=size)
covContemp
covLagged = matrix(rnorm(size * size), nrow=size, ncol=size)
covLagged

ampContemp = 0.1
ampLagged = 0.1
ampMeasure = 1
covMeasure = ampMeasure * diag(size)
covMeasure

ts = generate_timeseries(start, len, ampContemp, matContemp, covContemp, ampLagged, matLagged, covLagged, covMeasure)
ts  # variables on column and data for each day on rows 

clip_indices <- c(0,1)
clip_mins <- c(0.5, 0.7)
clip_maxs <- c(0.8,1.3)

ts_clip = clip_timeseries(ts, clip_indices, clip_mins, clip_maxs)
ts_clip = clip_outliers(ts, clip_sigma, ampMeasure, debug, start, ampContemp, matContemp, covContemp, ampLagged, matLagged, covLagged, covMeasure)
ts_clip

plot(ts[ , 2])
plot(ts_clip[ , 2])
