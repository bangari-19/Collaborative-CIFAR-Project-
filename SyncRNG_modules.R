# Using SyncRNG to replicate 
library('SyncRNG')
library(MASS)
library('mvtnorm')
library(ggplot2)

# Print random integers 
# s <- SyncRNG(seed=42)
# for (i in 1:10)
#   cat(s$randi(), '\n')


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

  if (is.null(rng)) {
    rng <- SyncRNG(seed=seed)
  }
  
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


generate_random_matrix <- function(size, mu, sigma, seed) {
  "Generate a matrix with elements sampled from a normal distribution.

  Args:
    size: The dimension of the square matrix (number of rows and columns).
    mu: The mean of the normal distribution.
    sigma: The standard deviation of the normal distribution.
    seed: A seed value for the random number generator.

  Returns:
    A size x size matrix with elements sampled from the specified normal distribution."
  
  random_numbers <- syncrng.box.muller(mu, sigma, size * size, seed)
  matrix(random_numbers, nrow = size, ncol = size)
}


generate_mvn_samples <- function(size, mean_vector, cov_matrix, seed) {
  "Generate a random vector from a multivariate normal distribution.

  Args:
    size: The dimension of the multivariate normal distribution (length of the mean vector).
    mean_vector: The mean vector of the multivariate normal distribution.
    cov_matrix: The covariance matrix of the multivariate normal distribution.
    seed: A seed value for the random number generator to ensure reproducibility.

  Returns:
    A random vector sampled from the specified multivariate normal distribution."
  
  L <- chol(cov_matrix)  # Cholesky decomposition of the covariance matrix
  random_numbers <- syncrng.box.muller(0, 1, size, seed)  # Generate standard normal random numbers
  random_matrix <- matrix(random_numbers, nrow = 1, ncol = size)  # Convert to matrix form
  mvn_samples <- t(t(random_matrix %*% L) + mean_vector)  # Transform to multivariate normal
  
  return(mvn_samples)
}


################################## EXAMPLES ###################################
# Histogram from syncrng
mu = 3
sigma = 1
n = 10000
seed = 42
random_samples <- syncrng.box.muller(mu, sigma, n, seed)
hist(random_samples, breaks=100)

# Sample a random matrix
random_matrix = generate_random_matrix(size=6, mu, sigma, seed)
print(random_matrix)

# Sample from mvn distribution 
size = 6
mean_vector = rep(0, 6)
cov_matrix = diag(rep(0.1, 6))
random_matrix_mvn = generate_mvn_samples(size=6, rep(0, 6), cov_matrix, seed)
print(random_matrix_mvn)

# Plot rmvnorm vs syncrng samples to compare multivariate normal distribution
size <- 2  # For 2D plot
mean_vector <- c(0, 0)
cov_matrix <- diag(c(1, 1))  # Zero on off-diagonals --> no correlation between x and y 
num_samples <- 10000

samples_rmvnorm <- rmvnorm(num_samples, mean_vector, cov_matrix)
samples_syncrng <- matrix(0, nrow = num_samples, ncol = size)

seed <- 123
for (i in 1:num_samples) {
  samples_syncrng[i, ] <- generate_mvn_samples(size, mean_vector, cov_matrix, seed + i)
}

samples_rmvnorm_df <- data.frame(x = samples_rmvnorm[,1], y = samples_rmvnorm[,2])
samples_syncrng_df <- data.frame(x = samples_syncrng[,1], y = samples_syncrng[,2])

samples_rmvnorm_df$source <- "rmvnorm"
samples_syncrng_df$source <- "syncrng"
combined_df <- rbind(samples_rmvnorm_df, samples_syncrng_df)

ggplot(combined_df, aes(x = x, y = y, color = source)) +
  geom_point(alpha = 0.5) +
  labs(title = "Comparison of Multivariate Normal Distribution Samples",
       x = "X-axis",
       y = "Y-axis") +
  theme_minimal() +
  scale_color_manual(values = c("rmvnorm" = "blue", "syncrng" = "red"))
