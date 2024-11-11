# Load necessary libraries
library(tidyverse)

setwd("C:/Users/magol001/OneDrive - Wageningen University & Research/Git_AMD/GammarusPulex_DEB-GUTS-T/results/quick_check/Hn2150")

# Define the exposure concentration groups
exposure_concentrations <- rep(seq(0, 1.2, by = 0.2), each = 10)

# Define model versions and filenames
model_versions <- c("SD_IMI", "SDT_IMI", "IT_IMI", "ITT_IMI")
model_versions <- c("SD_FPF", "SDT_FPF", "IT_FPF", "ITT_FPF")
filenames <- paste0(model_versions, ".txt")

# Initialize an empty list to store mean data
mean_data_list <- list()

# Loop through each file to read data, calculate means, and store in the list
for (i in seq_along(filenames)) {
  # Read the data file
  data <- read.table(filenames[i], header = TRUE, col.names = c("embryos", "juveniles", "adults", "juv.adults"))
  
  # Add a column for exposure concentration
  data$exposure_concentration <- exposure_concentrations
  
  # Calculate mean values for each exposure concentration
  mean_data <- data %>%
    group_by(exposure_concentration) %>%
    summarize(across(c("embryos", "juveniles", "adults", "juv.adults"), mean)) %>%
    mutate(model_version = model_versions[i])
  
  # Append to list
  mean_data_list[[i]] <- mean_data
}

# Combine all mean data into a single data frame
mean_data_all <- bind_rows(mean_data_list)

# Plot the mean values for each model version and concentration
ggplot(mean_data_all, aes(x = exposure_concentration)) +
  # Separate plots for each column of data
  geom_line(aes(y = juv.adults, color = model_version), size = 1) +
  labs(title = "Mean Values by Exposure Concentration and Model Version",
       x = "Exposure Concentration",
       y = "Mean Values",
       color = "Model Version") +
  theme_minimal() +
  #facet_wrap(~ model_version, scales = "free_y") +
  theme(legend.position = "right")

