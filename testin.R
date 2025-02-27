# Initialize input variables
size <- 6
pathin <- '/Users/reneehlozek/Dropbox/CIFAR_Sims/230222_syncrng_replaceNF/'
pathout <- '/Users/reneehlozek/Dropbox/CIFAR_Sims/230222_syncrng_replaceNF_sims/'
cols <- 2 * size
steps <- 200 
num_participants <- 150
num_iterations <- 10
save <- FALSE
covContempName <- 'randn'
covLaggedName <- 'randn'
measurecovName <- 'diag'
ampContemp <- 0.01
ampLagged <- 0.01
ampMeasure <- 1.0
maskZero <- TRUE
clipSamples <- FALSE
clipIndices <- c(0)
clipMins <- c(0.5,0.7)
clipMaxs <- c(0.8, 1.3)
clipOutliers <- TRUE
clipSigma <- 20