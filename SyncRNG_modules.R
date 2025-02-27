# Using SyncRNG to replicate 
library('SyncRNG')
library(MASS)
library('mvtnorm')
library(ggplot2)

universal_seed <- 123
universal_rng <- SyncRNG(seed=universal_seed)

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
  rng <- if (is.null(rng)) universal_rng else rng
  size <- length(mean_vector)
  L <- t(chol(cov_matrix))  # Returns upper triangular matrix by default. Transpose to match Python. 
  Z <- syncrng.box.muller(0, 1, size, rng = rng)
  mvn_samples <- as.vector(L %*% Z + mean_vector)
  return(mvn_samples)
}

# mean_vector <- c(1, 1)
# cov_matrix <- diag(c(2, 2))  # Zero on off-diagonals --> no correlation between x and y
# print(generate_mvn_samples(mean_vector=mean_vector, cov_matrix=cov_matrix))
# 
# print(generate_random_matrix(2, 0, 1))

################################## EXAMPLES ####################################
###################### CHECK IF OUTPUT MATCHES PYTHON ##########################

# mu <- 0
# sigma <- 1
# size <- 3
# mean_vector <- c(1, 2, 3)
# cov_matrix <- matrix(c(1, 0.5, 0.3,
#                        0.5, 1, 0.4,
#                        0.3, 0.4, 1), nrow=3, ncol=3)
# 
# # Generate random numbers
# random_samples <- syncrng.box.muller(mu, sigma, 10)
# print("Random Samples:")
# print(random_samples)
# 
# # Generate a random matrix
# random_matrix <- generate_random_matrix(size, mu, sigma)
# print("Random Matrix:")
# print(random_matrix)
# 
# # Generate multivariate normal samples
# mvn_sample <- generate_mvn_samples(mean_vector, cov_matrix)
# print("Multivariate Normal Sample:")
# print(mvn_sample)

#################### COMPARE R WITH SYNCRNG FOR BASIC PLOTS ####################
# n_samples <- 1000
# s <- SyncRNG(seed=42)
# sync_uniforms <- replicate(n_samples, s$rand())
# base_uniforms <- runif(n_samples)
# 
# par(mfrow = c(2, 1))
# hist(sync_uniforms, breaks = 30, main = "SyncRNG Uniform Distribution", xlab = "Value")
# hist(base_uniforms, breaks = 30, main = "Base R Uniform Distribution", xlab = "Value")
# 
# # Compare Normal Random Numbers
# sync_normals <- syncrng.box.muller(0, 1, n_samples, seed = 123)
# base_normals <- rnorm(n_samples)
# 
# par(mfrow = c(2, 1))
# hist(sync_normals, breaks = 30, main = "SyncRNG Normal Distribution", xlab = "Value")
# hist(base_normals, breaks = 30, main = "Base R Normal Distribution", xlab = "Value")


#################### COMPARE MVN SAMPLES FROM SYNCRNG AND RMVNORM ##############
# Plot rmvnorm vs syncrng samples to compare multivariate normal distribution
# mean_vector <- c(1, 1)
# cov_matrix <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)  # Zero on off-diagonals --> no correlation between x and y
# num_samples <- 1000
# 
# samples_rmvnorm <- rmvnorm(num_samples, mean_vector, cov_matrix)
# samples_syncrng <- matrix(0, nrow = num_samples, ncol = size)
# 
# for (i in 1:num_samples) {
#   samples_syncrng[i, ] <- generate_mvn_samples(mean_vector, cov_matrix)
# }
# 
# samples_rmvnorm_df <- data.frame(x = samples_rmvnorm[,1], y = samples_rmvnorm[,2])
# samples_syncrng_df <- data.frame(x = samples_syncrng[,1], y = samples_syncrng[,2])
# 
# samples_rmvnorm_df$source <- "rmvnorm"
# samples_syncrng_df$source <- "syncrng"
# combined_df <- rbind(samples_rmvnorm_df, samples_syncrng_df)
# 
# ggplot(combined_df, aes(x = x, y = y, color = source)) +
#   geom_point(alpha = 0.5) +
#   labs(title = "Comparison of Multivariate Normal Distribution Samples",
#        x = "X-axis",
#        y = "Y-axis") +
#   theme_minimal() +
#   scale_color_manual(values = c("rmvnorm" = "blue", "syncrng" = "red"))


#################### PLOT MVN SAMPLES FROM ONLY SYNCRNG ########################
# mean_vector <- c(0, 0)
# cov_matrix <- matrix(c(1, 0.9, 0.9, 1), nrow = 2)
# n_samples <- 100
# samples <- matrix(0, nrow = 2, ncol = n_samples)
# 
# for (i in 1:n_samples) {
#   samples[, i] <- generate_mvn_samples(mean_vector, cov_matrix)
# }
# 
# plot(samples[1, ], samples[2, ], main = "Bivariate normal with variance 1, covariance 0.9", asp = 1)
