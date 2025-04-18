# Script to generate and store simulated pre-clipped and clipped timeseries data for all individuals 

# Set working directory
setwd("/Users/reneehlozek/Code/Collaborative-CIFAR-Project-/")
library('SyncRNG')
seed=123456
s <- SyncRNG(seed=seed)
universal_rng <- SyncRNG(seed=seed)
# Load scripts
source("input.R")
source("sim_modules.R")


# Initialize input variables
size <- size
pathin <- pathin
pathout <- pathout
cols <- 2 * size
steps <- steps 
print(steps)
num_participants <- num_participants
num_iterations <- num_iterations
save <- save
covContempName <- covContempName
covLaggedName <- covLaggedName
measurecovName <- measurecovName
ampContemp <- ampContemp
ampLagged <- ampLagged
ampMeasure <- ampMeasure
maskZero <- maskZero
clipSamples <- clipSamples
clipIndices <- clipIndices
clipMins <- clipMins
clipMaxs <- clipMaxs
clipOutliers <- clipOutliers
clipSigma <- clipSigma

print(paste("Sim Directory is", pathin))
print(paste("reading from", num_participants, "participants", steps, "timesteps", ampContemp, ampLagged, ampMeasure, maskZero))
printpath <- file.path(pathout, "sims", paste0("numParticipants_", num_participants))
print(paste("Main dir is", printpath))

#dev.off()

# Make directory for sims
# Example: /Users/mac/Desktop/CIFAR/220429.FinalMatrices_sims/sims/numParticipants_15/steps_10/Contemp_randn_amp_0.1/Lagged_randn_amp_0.1/Measure_diag_amp_1/Mask_TRUE/clipsigma_3
savepath <- file.path(pathout, "sims")
savepath <- file.path(savepath, paste0("numParticipants_", num_participants))
savepath <- file.path(savepath, paste0("steps_", steps))
savepath <- file.path(savepath, paste0("Contemp_", covContempName, "_amp_", ampContemp))
savepath <- file.path(savepath, paste0("Lagged_", covLaggedName, "_amp_", ampLagged))
savepath <- file.path(savepath, paste0("Measure_", measurecovName, "_amp_", ampMeasure))
savepath <- file.path(savepath, paste0("Mask_", maskZero))
savepath <- file.path(savepath, paste0("clipsigma_", clipSigma))
cat("Save path is", savepath, "\n")


# Make directories if they do not exist
if (!dir.exists(savepath)) {
  dir.create(savepath, recursive = TRUE)
  # if (debug) { cat("making", savepath, "\n") }
  }

# Check and create clip and preclip directories if clip_outliers is TRUE
if (clipOutliers) {
  # if (debug) {cat("Making clip/preclip dirs\n")}
  savepathprecliporig <- file.path(savepath, "preclip/")
  savepathcliporig <- file.path(savepath, "clip/")
  dir.create(savepathprecliporig, showWarnings = FALSE, recursive = TRUE)
  dir.create(savepathcliporig, showWarnings = FALSE, recursive = TRUE)
}

# Loop over iterations
for (j in 1:num_iterations) {
  # Set unique seed for each iteration to control randomness 
  seed = seed+j
  s <- SyncRNG(seed=seed)
  universal_rng <- SyncRNG(seed=seed)
  
  # if (debug) {
  #   cat("making rep dirs\n")
  # }
  
  savepathclip <- file.path(savepathcliporig, paste0("rep_", j, "/"))
  savepathpreclip <- file.path(savepathprecliporig, paste0("rep_", j, "/"))
  dir.create(savepathclip, showWarnings = FALSE, recursive = TRUE)
  dir.create(savepathpreclip, showWarnings = FALSE, recursive = TRUE)
  
  # Loop over participants
  for (number in 1:num_participants) {
    # if (debug) {
    #   cat("Writing sims for participant", number, "and sim iteration", j, "\n")
    # }
    
    # printpath(pathin)
    # print(huh)
    # Load in input contemp and lagged coeff matrix data
    file <- file.path(pathin, paste0("ind_",number, ".csv"))

    data <- read.csv(file, header = TRUE, sep = ",")
    matContemp <- as.matrix(data[, (size + 1):cols])
    matLagged <- as.matrix(data[, 1:size])
    
    # Generate random covariance matrices
    if (covContempName == "randn") {
      covContemp <- matrix(rnorm(size * size), nrow = size, ncol = size)
    } else {
      covContemp <- matrix(1, nrow = size, ncol = size)
    }
    
    if (covLaggedName == 'randn') {
      covLagged <- matrix(rnorm(size * size), nrow = size, ncol = size)
    } else {
      covLagged <- matrix(1, nrow = size, ncol = size)
    }
    
    # Create measureCov matrix 
    if (measurecovName == 'diag') {
      measureCov <- ampMeasure * diag(size)
    } else {
      measureCov <- matrix(0, nrow = size)
    }
    
    # Apply mask
    maskContemp <- ifelse(maskZero, make_mask(matContemp, contemp = TRUE), matrix(1, nrow = size, ncol = size))
    maskLagged <- ifelse(maskZero, make_mask(matLagged, contemp = FALSE), matrix(1, nrow = size, ncol = size))
    
    # Generate timeseries data
    samples <- generate_timeseries(start, steps, ampContemp, matContemp, covContemp,ampLagged, matLagged, covLagged, measureCov, save)
    # print(nrow(samples))
    # print(ncol(samples))
    
    max_clip_threshold <- sqrt(diag(measureCov)[1]) * 3
    
    # if (debug) {
    #   cat(max_clip_threshold, 'max', "\n")
    #   for (p in 1:ncol(samples)) {
    #     inds_pre <- which(samples[, p] > max_clip_threshold)
    #     cat("pre clip:", inds_pre, samples[inds_pre, p], p, "\n")
    #   }
    #   cat("---- TEST pre clip -----\n")
    # }
    
    # Clip time series data if specified
    if (clipSamples) {
      # Checking the clipping
      samples_clip <- clip_timeseries(samples, clip_indices, clip_mins, clip_maxs)
      samples <- samples_clip
      # if (debug) {
      #   cat('saved clip\n')
      # }
    }
    
    if (clipOutliers) {
      # if (debug) {
      #   cat('Clipping samples if needed\n')
      # }
      
      savefile <- file.path(savepathpreclip, paste0('ind_', number, '_preclip.txt'))
      write.table(samples, file = savefile, sep = ',', row.names = FALSE, col.names = FALSE)
      
      samples_clip <- clip_outliers(
        samples, clip_sigma, ampMeasure, debug, start,
        ampContemp, matContemp, covContemp, ampLagged, matLagged,
        covLagged, measureCov
      )
      
      samples <- samples_clip
      
      #  if (debug) {
      #   for (p in 1:ncol(samples)) {
      #     inds_post <- which(samples[, p] > max_clip_threshold)
      #     cat("post clip:", inds_post, samples[inds_post, p], p, "\n")
      #   }
      #   cat("---- TEST post clip -----\n")
      #   cat('saved clip\n')
      # }
    
      savefile <- file.path(savepathclip, paste0('ind_', number, '.txt'))
    }
    
    write.table(samples, file = savefile, sep = ',', row.names = FALSE, col.names = FALSE)
  }
}



