# Load necessary libraries
library(readr)
library(ggplot2)
library(gghalves)
library(gridExtra)  # For displaying plots

# Read the CSV file with skipping the header row
data <- read.csv("/home/joeh/rdl_helicon_topspeeds_and_increase_decrease.csv", header = TRUE)

# Filter out rows with missing values
data <- data[complete.cases(data), ]

# Filter data to include only values under 200
data <- data[data$top_speed_difference < 200, ]

# Check data types and unique values
str(data$top_speed_difference)
unique(data$sleep_condition)

# Unique genotypes in the dataset
genotypes <- unique(data$genotype)

# Function to perform t-test for each genotype
perform_t_test <- function(genotype) {
  # Subset data for the current genotype
  genotype_data <- data[data$genotype == genotype, ]
  
  # Extract top_speed_difference for sd and rested sleep conditions
  sd_speeds <- genotype_data$top_speed_difference[genotype_data$sleep_condition == "sd"]
  rested_speeds <- genotype_data$top_speed_difference[genotype_data$sleep_condition == "rested"]
  
  # Check if there are enough observations for t-test
  if (length(sd_speeds) < 2 || length(rested_speeds) < 2) {
    cat("Insufficient observations for genotype", genotype, "\n")
    return()  # Skip t-test if there are not enough observations
  }
  
  # Perform t-test
  t_test_result <- t.test(sd_speeds, rested_speeds)
  
  # Print results
  cat("Genotype:", genotype, "\n")
  cat("t-test p-value:", t_test_result$p.value, "\n")
  cat("\n")
}

# Perform t-test for each genotype
for (genotype in genotypes) {
  perform_t_test(genotype)
}

# Function to create raincloud plots for each genotype
plot_raincloud <- function(genotype) {
  # Subset data for the current genotype
  genotype_data <- data[data$genotype == genotype, ]
  
  # Create the raincloud plot only if there are observations
  if (nrow(genotype_data) > 0) {
    # Create the raincloud plot
    p <- ggplot(genotype_data, aes(x = genotype, y = top_speed_difference, color = sleep_condition)) +
      geom_point(position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.1), size = 2, alpha = 0.7) +
      labs(title = paste( genotype),
           x = "Genotype", y = "Top Speed Difference") +
      theme_minimal() +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))  # Adjust margins
    
    return(p)
  } else {
    return(NULL)
  }
}

# Create a list to store plots
plots <- list()

# Create raincloud plot for each genotype and store in the list
for (genotype in genotypes) {
  plot <- plot_raincloud(genotype)
  if (!is.null(plot)) {
    plots[[genotype]] <- plot
  }
}

# Display plots
grid.arrange(grobs = plots, ncol = 1)
