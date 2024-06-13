library(ggplot2)

# Function to construct file paths based on parameters
# Example:
# "/Users/mac/Desktop/CIFAR/220429.FinalMatrices_sims/sims/numParticipants_10/steps_100/Contemp_randn_amp_0.01/Lagged_randn_amp_0.01/Measure_diag_amp_1/Mask_TRUE/clipsigma_3/preclip/rep_1/ind_5_preclip.txt"

construct_filepath <- function(base_path, num_participants, steps, contemp_type, contemp_amp, lagged_type, lagged_amp, measure_type, measure_amp, mask, clipsigma, participant_number, iteration, preclipped = FALSE) {
  # Build the directory path with the correct separators between configurations and values
  directory <- paste(paste("numParticipants", num_participants, sep='_'),
                     paste("steps", steps, sep='_'),
                     paste(contemp_type, "amp", contemp_amp, sep="_"),
                     paste(lagged_type, "amp", lagged_amp, sep="_"),
                     paste(measure_type, "amp", measure_amp, sep="_"),
                     paste("Mask", ifelse(mask, "TRUE", "FALSE"), sep="_"),
                     paste("clipsigma", clipsigma, sep="_"), 
                     sep="/")
  
  # Choose the correct subfolder based on whether the data is preclipped or clipped
  subfolder <- ifelse(preclipped, "preclip", "clip")
  
  # Construct the full file path using the base path and all directory components
  folder <- file.path(base_path, "sims", directory, subfolder, paste0("rep_", iteration))
  
  # Define the filename depending on whether it is preclipped or clipped data
  filename_suffix <- ifelse(preclipped, "_preclip", "")
  filename <- paste0("ind_", participant_number, filename_suffix, ".txt")
  
  # Combine the folder path and filename into the full file path
  file.path(folder, filename)
}

# Function to plot timeseries data for a specified variable for multiple individuals
plot_timeseries <- function(base_path, num_participants, steps, contemp_type, contemp_amp, lagged_type, lagged_amp, measure_type, measure_amp, mask, clipsigma, column, iterations, participant_numbers, show_preclip=TRUE, show_clip=TRUE) {
  # Initialize an empty data frame for plotting
  df <- data.frame()
  
  for (iteration in iterations) {
    for (participant_number in participant_numbers) {
      # Construct file paths for preclipped and clipped data
      preclip_path <- construct_filepath(base_path, num_participants, steps, contemp_type, contemp_amp, lagged_type, lagged_amp, measure_type, measure_amp, mask, clipsigma, participant_number, iteration, TRUE)
      clip_path <- construct_filepath(base_path, num_participants, steps, contemp_type, contemp_amp, lagged_type, lagged_amp, measure_type, measure_amp, mask, clipsigma, participant_number, iteration, FALSE)
      print(preclip_path)
      # Read the pre-clipped and clipped data
      if (show_preclip && file.exists(preclip_path)) {
        preclip_data <- read.table(preclip_path, sep=",", header=FALSE)
        if (nrow(preclip_data) > 0) {
          preclip_df <- data.frame(
            days = 1:nrow(preclip_data),
            value = preclip_data[[column]],
            type = "Preclip",
            participant = paste("Participant", participant_number),
            iteration = iteration
          )
          df <- rbind(df, preclip_df)
        }
      }
      
      if (show_clip && file.exists(clip_path)) {
        clip_data <- read.table(clip_path, sep=",", header=FALSE)
        if (nrow(clip_data) > 0) {
          clip_df <- data.frame(
            days = 1:nrow(clip_data),
            value = clip_data[[column]],
            type = "Clip",
            participant = paste("Participant", participant_number),
            iteration = iteration
          )
          df <- rbind(df, clip_df)
        }
      }
    }
  }
  
  if (nrow(df) == 0) {
    stop("No data found for the specified parameters.")
  }
  
  # Plotting
  p <- ggplot(df, aes(x = days, y = value, color = as.factor(participant), linetype = type, group = interaction(participant, type))) +
    geom_line() +
    xlab("Time [days]") +
    ylab("Value") +
    ggtitle(paste("Timeseries for Variable (", column,") for Individuals (", length(participant_numbers), ")")) +
    theme_minimal() +
    scale_color_viridis_d() +  # Using viridis color scale for better readability
    scale_linetype_manual(values = c("Preclip" = "solid", "Clip" = "dotted")) +
    labs(color = "Participant", linetype = "Timeseries")
  
  print(p)
  
}

# Example usage
# Provide the base path and other parameters
base_path <- "/Users/mac/Desktop/CIFAR/220429.FinalMatrices_sims"  # Update this to your base path
num_participants <- 10
steps <- 100
contemp_type <- "Contemp_randn"
contemp_amp <- 0.01
lagged_type <- "Lagged_randn"
lagged_amp <- 0.01
measure_type <- "Measure_diag"
measure_amp <- 1
mask <- TRUE
clipsigma <- 3
column <- 1  # Specify the column index you are interested in (1-based index for R)
iterations <- 1:5  # Specify the iterations you want to plot
participant_numbers <- 1:3  # Specify the participant numbers

# Plot the timeseries for the specified column
plot_timeseries(base_path, num_participants, steps, contemp_type, contemp_amp, lagged_type, lagged_amp, measure_type, measure_amp, mask, clipsigma, column, iterations, participant_numbers, show_preclip=TRUE, show_clip=TRUE)

