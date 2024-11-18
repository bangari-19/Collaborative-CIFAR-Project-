syncrng.box.muller <- function(mu, sigma, n, seed=0, rng=NULL) {
  "Generate random samples from a normal distribution using the Box-Muller transform.
  
  Args:
     mu: The mean of the normal distribution.
     sigma: The standard deviation of the normal distribution.
     n: The number of random samples to generate.
     seed: An optional seed value for the random number generator (default is 0).
     rng: An optional instance of SyncRNG. If NULL, a new instance will be created using the seed.
  
   Returns:
     A vector of n random samples from the specified normal distribution."
  
  rng <- if (is.null(rng)) universal_rng else rng
  
  two.pi <- 2 * pi
  ngen <- ceiling(n / 2)
  out <- replicate(2 * ngen, 0.0)
  
  for (i in 1:ngen) {
    u1 <- 0.0
    u2 <- 0.0
    
    while (u1 == 0) { u1 <- rng$rand(); }
    while (u2 == 0) { u2 <- rng$rand(); }
    
    mag <- sigma * sqrt(-2.0 * log(u1))
    z0 <- mag * cos(two.pi * u2) + mu
    z1 <- mag * sin(two.pi * u2) + mu
    
    out[2*i - 1] = z0;
    out[2*i] = z1;
  }
  return(out[1:n]);
}


generate_random_matrix <- function(size, mu, sigma, rng=NULL) {
  "Generate a matrix with elements sampled from a normal distribution.

  Args:
    size: The dimension of the square matrix (number of rows and columns).
    mu: The mean of the normal distribution.
    sigma: The standard deviation of the normal distribution.
    seed: A seed value for the random number generator.

  Returns:
    A size x size matrix with elements sampled from the specified normal distribution."
  
  rng <- if (is.null(rng)) universal_rng else rng
  random_numbers <- syncrng.box.muller(mu, sigma, size * size, rng = rng)
  matrix(random_numbers, nrow = size, ncol = size)
}


generate_mvn_samples <- function(mean_vector, cov_matrix, rng=NULL) {
  "Generate a random vector from a multivariate normal distribution.

  Args:
    size: The dimension of the multivariate normal distribution (length of the mean vector).
    mean_vector: The mean vector of the multivariate normal distribution.
    cov_matrix: The covariance matrix of the multivariate normal distribution.
    seed: A seed value for the random number generator to ensure reproducibility.

  Returns:
    A random vector sampled from the specified multivariate normal distribution."
  
  # Cholesky decomposition of the covariance matrix
  # Transpose to get lower triangular matrix
  rng <- if (is.null(rng)) universal_rng else rng
  size <- length(mean_vector)
  L <- t(chol(cov_matrix))
  Z <- syncrng.box.muller(0, 1, size, rng = rng)
  mvn_samples <- L %*% Z + mean_vector
  return(mvn_samples)
}
 

plot_matrix <- function(mat, save=FALSE, name=NULL) {
  # Plots a matrix using plot.matrix library and optionally saves the plot to a file.
  # 
  # Args:
  #   mat: Matrix to be plotted.
  #   save: Logical indicating whether to save the plot.
  #   name: Filename for saving the plot, required if save is TRUE.
  
  if(require('plot.matrix')){}
  else {install.packages('https://cran.r-project.org/src/contrib/plot.matrix_1.6.2.tar.gz',repos = NULL, type="source",deps=TRUE)}

  library(plot.matrix)
  plot(mat)
  
  if (save) {
    png(file = name)  
    dev.off()   
  }
}


save_matrix <- function(mat, name) {
  "Saves a matrix to a CSV file.

  Args:
    mat: Matrix to be saved.
    name: Filename for the saved CSV."
  
  write.csv(mat, file = name)
}


make_mask <- function(mat, contemp) {
  "Generates a mask for a matrix. Zeros in the original matrix and optionally the diagonal are masked.

  Args:
    mat: Matrix for which to generate the mask.
    contemp: Logical indicating whether to also mask the diagonal.

  Returns:
    A mask matrix with the same dimensions as 'mat'."
  
  size = nrow(mat)  
  mask <- matrix(1, size, size)
  if (contemp) {
    diag(mask) <- 0  
  }
  
  inds = which(mat == 0, arr.ind = T)  # find indices where mat == 0 in (row, col) format
  mask[inds] = 0  # make mask 0 at those indices 
  
  return(mask)
}


coeff_draw_from_cov <- function(amplitude, mat, cov, mask) {
  "Samples a new matrix using a given covariance structure, amplitude, and mask.

  Args:
    amplitude: Scaling factor for the covariance-based sample.
    mat: Base matrix to which the sampled changes are applied.
    cov: Covariance matrix used for sampling.
    mask: Mask matrix to apply to the sampled matrix.

  Returns:
    A matrix representing the noisy version of 'mat' after applying the sampled changes and mask."
  
  size <- nrow(cov)
  sample_step <- cov %*% (amplitude * matrix(rnorm(size * size), nrow = size, ncol = size))
  stepmat <- mat + sample_step
  noisy_mat <- stepmat * mask
  return(noisy_mat)
}


generate_timeseries <- function(start, len, contempamp, contempmat, contempcov, lagamp, lagmat, lagcov, measurecov, save=FALSE, debug=FALSE) {
  "Generates a time series matrix based on initial values, amplitude factors, and covariance structures.
  The function simulates the evolution of a set of variables over time, considering both contemporary
  and lagged influences, and can optionally save intermediate results and plots for analysis.

  Args:
    start: Numeric vector of initial values for the time series variables.
    len: Integer specifying the length of the time series (number of time points).
    contempamp: Numeric, amplitude factor for contemporary effects.
    contempmat: Matrix/Data frame, specifies the contemporary relationships between variables.
    contempcov: Matrix, covariance structure for contemporary effects.
    lagamp: Numeric, amplitude factor for lagged effects.
    lagmat: Matrix/Data frame, specifies the lagged relationships between variables.
    lagcov: Matrix, covariance structure for lagged effects.
    measurecov: Matrix, measurement error covariance structure.
    save: Logical, if TRUE, saves intermediate matrices and plots to files.
    debug: Logical, if TRUE, prints debug information during the simulation.

  Returns:
    A matrix representing the time series of variables over the specified time period."
  
  size <- length(start)  
  samples <- matrix(0, nrow = len, ncol = size)
  samples[1,] <- start  # index start at 1 in R
  lagmask <- make_mask(lagmat, contemp=FALSE)
  contempmask <- make_mask(contempmat, contemp=TRUE)
  
  # Convert type from data.frame to matrix 
  noisy_lagmat <- as.matrix(coeff_draw_from_cov(lagamp, lagmat, lagcov, lagmask))
  noisy_contempmat <- as.matrix(coeff_draw_from_cov(contempamp, contempmat, contempcov, contempmask))
  
  if (save) {
    save_matrix(noisy_lagmat, 'noisy_lagmat.csv')
    plot_matrix(noisy_lagmat, TRUE, 'noisy_lagmat.png')
    save_matrix(noisy_contempmat, 'noisy_contempmat.csv')
    plot_matrix(noisy_contempmat, TRUE, 'noisy_contempmat.png')
  }
  
  if(require('mvtnorm')){}
  else {install.packages('https://cran.r-project.org/src/contrib/mvtnorm_1.2-3.tar.gz', repos = NULL, type="source",deps=TRUE)}
  library(mvtnorm)

    for (i in 2:len) {
      samples[i,] <- samples[i-1,] %*% as.matrix(noisy_lagmat) %*% as.matrix(solve(diag(size) - as.matrix(noisy_contempmat))) + start
      samples[i,] <- samples[i,] + rmvnorm(1, rep(0, size), measurecov)  
      # rmvnorm generates one random sample from the multivariate normal distribution with mean 0 and covariance 'measurecov'
      
      if (debug) {
        cat(sprintf("Day: %d\n", i))
        cat("sample[i, ]:\n"); print(samples[i, ])
        cat(sprintf("length(samples[i-1, ]): %d\n", length(samples[i-1,])))
        cat(sprintf("dim of noisy_lagmat: %s\n", toString(dim(noisy_lagmat))))
        cat(sprintf("dimension of (identity - noisy_contempmat): %s\n", toString(dim(diag(size) - noisy_contempmat))))
        cat("noisy_contempmat:\n"); print(noisy_contempmat)
        cat(sprintf("class of (identity - noisy_contempmat): %s\n", class(solve(diag(size) - noisy_contempmat))))
        cat(sprintf("class of start: %s\n", class(start)))
        cat(sprintf("class of samples: %s\n", class(samples)))
        cat(sprintf("class of noisy_contempmat: %s\n", class(as.matrix(noisy_contempmat))))
        cat(sprintf("class of noisy_lagmat: %s\n", class(noisy_lagmat)))
      }
    }
    
    samples <- samples[-1,] - start  # Remove first row of samples and subtract 'start' from all rows
    return(samples)
}


clip_timeseries <- function(timeseries, indices_to_clip, min_vec, max_vec) {
  "Clips values in the time series that are outside the bounds defined by min_vec and max_vec.

  Args:
    timeseries: Matrix representing the time series data.
    indices_to_clip: Indices of columns in the time series to check for clipping.
    min_vec: Vector of minimum allowable values for each column.
    max_vec: Vector of maximum allowable values for each column.

  Returns:
    A matrix with the same dimensions as 'timeseries', with outliers clipped."
  
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
  "Clips outliers in a time series matrix by replacing outlying values based on a specified sigma threshold.
  Outliers are identified in each column of the time series and replaced with new data generated considering
  the sigma value and the underlying time series model parameters.

  Args:
    timeseries: Matrix representing the time series data.
    sigma: Numeric threshold for identifying outliers based on standard deviation.
    measure_amp: Numeric, amplitude factor for measurement errors.
    debug: Logical, if TRUE, enables printing of debug information.
    start: Numeric vector of initial values for the time series.
    contempamp: Numeric, amplitude factor for contemporary effects.
    contempmat: Matrix/Data frame specifying contemporary relationships between variables.
    contempcov: Matrix, covariance structure for contemporary effects.
    lagamp: Numeric, amplitude factor for lagged effects.
    lagmat: Matrix/Data frame specifying the lagged relationships between variables.
    lagcov: Matrix, covariance structure for lagged effects.
    measurecov: Matrix, covariance structure for measurement errors.

  Returns:
    A matrix representing the clipped time series, with outliers adjusted based on the specified criteria."
  
  ts <- timeseries
  # Iterate through all columns of timeseries
  for (i in 1:ncol(timeseries)) {
    min_val <- -sigma * sqrt(measure_amp)  
    max_val <- sigma * sqrt(measure_amp)  
    min_inds <- which(timeseries[, i] < min_val)  # indices where the column values are less than min_val
    
    if (debug && length(min_inds) > 0) {
      cat("bad vals min", timeseries[min_inds, i], "inds", min_inds, "\n")
    }
    
    countmax <- 2  # Number of replacement attempts for outliers
    # Loop over the min outlier indices
    for (minind in min_inds) {
      # Check if minind is within bounds
      if (minind >= 2 && minind <= (nrow(ts) - 2)) {
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
      } else {
        if (debug) {
          cat("minind out of bounds:", minind, "\n")
        }
      }
    }
    
    max_inds <- which(timeseries[, i] > max_val)
    
    if (debug && length(max_inds) > 0) {
      cat("bad vals max", timeseries[max_inds, i], "inds", max_inds, "\n")
    }
    
    for (maxind in max_inds) {
      # Check if maxind is within bounds
      if (maxind >= 2 && maxind <= (nrow(ts) - 2)) {
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
      } else {
        if (debug) {
          cat("maxind out of bounds:", maxind, "\n")
        }
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


