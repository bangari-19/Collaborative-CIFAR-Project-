rm(list=ls())

# Prompt user for input parameters for matrix generation
nvar_input <- readline(prompt = "Enter the number of variables (e.g., 6): ")
AR_input <- readline(prompt = "Enter the autoregressive coefficient (e.g., 0.5): ")
dens_input <- readline(prompt = "Enter the density of connections (e.g., 0.2): ")
p_group_input <- readline(prompt = "Enter the proportion of group paths (e.g., 0.5): ")
p_con_input <- readline(prompt = "Enter the proportion of contemporaneous paths (e.g., 0.5): ")
con_b_input <- readline(prompt = "Enter the contemporaneous path coefficient (e.g., 0.3 or -0.3): ")
lag_b_input <- readline(prompt = "Enter the lagged path coefficient (e.g., 0.3 or -0.3): ")

# Convert to correct data types
nvar <- as.integer(nvar_input)
AR <- as.numeric(AR_input)
dens <- as.numeric(dens_input)
p.group <- as.numeric(p_group_input)
p.con <- as.numeric(p_con_input)
con.b <- as.numeric(con_b_input)
lag.b <- as.numeric(lag_b_input)

# Display user input
cat("You have entered the following parameters:\n")
cat("Number of Variables:", nvar, "\n")
cat("Autoregressive Coefficient:", AR, "\n")
cat("Density of Connections:", dens, "\n")
cat("Proportion of Group Paths:", p.group, "\n")
cat("Proportion of Contemporaneous Paths:", p.con, "\n")
cat("Contemporaneous Path Coefficient:", con.b, "\n")
cat("Lagged Path Coefficient:", lag.b, "\n")

# Generate matrices
mat.generate.asw <- function(p.con, nvar, AR, dens, p.group, con.b, lag.b) {
  A <- matrix(0, ncol = nvar, nrow = nvar)
  Phi <- matrix(0, ncol = nvar, nrow = nvar)
  
  indices <- which(Phi == 0, arr.ind = TRUE)
  row.col <- c(17, 24, 31, 28)
  diag.set <- c(8, 22, 36)
  
  n.p.1 <- row.col[1:round(p.con*length(row.col))]
  n.p.2 <- row.col[(length(n.p.1)+1):length(row.col)]
  
  grp.con <- n.p.1
  grp.lag <- n.p.2
  grp.diag <- diag.set
  
  Phi[indices[grp.lag,]] <- lag.b
  A[indices[grp.con,]] <- con.b
  Phi[indices[grp.diag,]] <- AR
  
  all <- cbind(Phi, A)
  ind.pres <- which(all != 0, arr.ind = TRUE)
  
  all.lvl <- matrix(NA, ncol = ncol(all), nrow = nrow(all))
  all.lvl[ind.pres] <- "grp"
  
  res <- list(sub1 = all, lvl1 = all.lvl)
  return(res)
}

# Generate time series
ts.generate.asw <- function (mat, lvl, t, dens, p.group, con.b, lag.b, p.con) {
  repeat {
    repeat {
      v <- ncol(mat) / 2
      Phi <- mat[, 1:v]
      A <- mat[, (v+1):(v*2)]
      indices.A <- which(A == 0, arr.ind = TRUE)
      indices.Phi <- which(Phi == 0, arr.ind = TRUE)
      indices.A <- indices.A[which(indices.A[,1] != indices.A[,2]),]
      pos.all <- (v*v-v) * 2
      cnt.all <- dens * p.group * pos.all
      rand <- sample(c(round(cnt.all/2), (round(cnt.all)-round(cnt.all/2))), size=2, replace=FALSE)
      row.col.A <- sample(1:nrow(indices.A), rand[1], replace = FALSE)
      row.col.Phi <- sample(1:nrow(indices.Phi), rand[2], replace = FALSE)
      Phi[indices.Phi[row.col.Phi,]] <- lag.b
      A[indices.A[row.col.A,]] <- con.b
      break
    }
    
    st <- (t + 50)
    noise <- matrix(rnorm(v * st, 0, 1), v)
    I <- diag(v)
    time <- matrix(0, nrow = v, ncol = (st+1))
    time1 <- matrix(0, nrow = v, ncol = st)
    
    for (i in 1:st) {
      time1[,i] <- solve(I - A) %*% (Phi %*% time[,i] + noise[,i])
      time[,i+1] <- time1[,i]
    }
    
    time1 <- time1[, (51:(50+t))]
    series <- t(time1)
    paths <- cbind(Phi, A)
    if (abs(max(series, na.rm = TRUE)) < 20 & abs(min(series, na.rm = TRUE)) > .01 & abs(min(series, na.rm = TRUE)) < 20) break
  }
  
  lvl[is.na(lvl) & paths != 0] <- "ind"
  
  list <- list("series" = series, "paths" = paths, "levels" = lvl)
  return(list)
}

# Generate matrices using user-defined parameters
res <- mat.generate.asw(
  p.con = p.con, 
  nvar = nvar, 
  AR = AR,
  p.group = p.group, 
  dens = dens,
  con.b = con.b, 
  lag.b = lag.b
)

# Generate time series
# Can also make t user-defined
out <- ts.generate.asw(
  mat = res$sub1,
  lvl = res$lvl1,
  t   = 100,
  p.group = p.group, 
  dens = dens,
  con.b = con.b, 
  lag.b = lag.b,
  p.con = p.con
)

# Display generated matrices and time series
cat("\nGenerated Matrix (sub1):\n")
print(res$sub1)
cat("\nConnection Levels (lvl1):\n")
print(res$lvl1)
cat("\nGenerated Time Series:\n")
print(out$series)

# Load necessary library
library(ggplot2)

# Convert time series to a data frame
ts_data <- as.data.frame(out$series)
ts_data$Time <- 1:nrow(ts_data)

# Create the plot using a format similar to your preferred style
ggplot(ts_data, aes(x = Time)) +
  geom_line(aes(y = ts_data[, 1], color = "Var 1"), linetype = "solid") +
  geom_line(aes(y = ts_data[, 2], color = "Var 2"), linetype = "solid") +
  geom_line(aes(y = ts_data[, 3], color = "Var 3"), linetype = "solid") +
  geom_line(aes(y = ts_data[, 4], color = "Var 4"), linetype = "solid") +
  geom_line(aes(y = ts_data[, 5], color = "Var 5"), linetype = "solid") +
  geom_line(aes(y = ts_data[, 6], color = "Var 6"), linetype = "solid") +
  xlab("Time") +
  ylab("Value") +
  ggtitle("Generated Time Series for All Parameters") +
  scale_color_manual(values = c("Var 1" = "blue", "Var 2" = "red", "Var 3" = "green", 
                                "Var 4" = "purple", "Var 5" = "orange", "Var 6" = "brown")) +
  scale_linetype_manual(values = c("Var 1" = "solid", "Var 2" = "dashed", "Var 3" = "solid", 
                                   "Var 4" = "dashed", "Var 5" = "solid", "Var 6" = "solid")) +
  labs(color = "Variables") +
  theme_minimal()
