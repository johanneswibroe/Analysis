# Load necessary libraries
library(readr)
library(ggplot2)
library(dplyr)
library(lme4)
library(emmeans)
library(gridExtra)

# Read the CSV file
data <- read.csv("/home/joeh/Helicon_RDL_FEMALE/female_helicon_rdl_2800MS.csv", header = TRUE)

# Filter out rows with missing values
data <- data[complete.cases(data), ]

# Filter to include the relevant genotypes
data <- data[data$genotype %in% c("24B11_RDL", "24B11_G4", "RDL_RNAi"), ]

# Reorder the genotype factor levels to make "24B11_RDL" last
data$genotype <- factor(data$genotype, levels = c("24B11_G4", "RDL_RNAi", "24B11_RDL"))

# Function to flag outliers based on z-score
flag_outliers <- function(x, cutoff = 2) {
  z_scores <- abs(scale(x))
  z_scores < cutoff
}

# Flag outliers for each genotype and sleep condition combination
data <- data %>%
  group_by(genotype, sleep_condition) %>%
  mutate(is_not_outlier = flag_outliers(abs_difference)) %>%
  ungroup()

# Create a new column for cleaned data, replacing outliers with NA
data$abs_difference_clean <- ifelse(data$is_not_outlier, data$abs_difference, NA)

# Fit a linear model without random effects, using the cleaned data
model <- lm(abs_difference_clean ~ sleep_condition * genotype, data = data)

# Summary of the model
summary(model)

# Post-hoc analysis with p-value adjustment
emmeans_results <- emmeans(model, specs = pairwise ~ sleep_condition | genotype)
emmeans_results_df <- as.data.frame(emmeans_results$contrasts)
emmeans_results_df$p.value.original <- emmeans_results_df$p.value
emmeans_results_df$p.value.adjusted <- p.adjust(emmeans_results_df$p.value, method = "BH")

# Display the results with original and adjusted p-values
print(emmeans_results_df)

# Plot the interaction using cleaned data
interaction_plot <- ggplot(data, aes(x = sleep_condition, y = abs_difference_clean, color = genotype, group = genotype)) +
  stat_summary(fun = mean, geom = "line", position = position_dodge(width = 0.2), size = 1, na.rm = TRUE) +
  stat_summary(fun = mean, geom = "point", position = position_dodge(width = 0.2), size = 3, na.rm = TRUE) +
  labs(title = "Interaction Plot of Sleep Condition and Genotype on Absolute Difference (Outliers Removed)",
       x = "Sleep Condition",
       y = "Absolute Difference",
       color = "Genotype") +
  theme_minimal()

# Display the interaction plot
print(interaction_plot)

# Function to create pairwise comparison plot for each genotype
plot_pairwise_comparison <- function(genotype) {
  # Subset data for the current genotype
  genotype_data <- data[data$genotype == genotype, ]
  
  # Create the boxplot
  p_boxplot <- ggplot(genotype_data, aes(x = sleep_condition, y = abs_difference_clean, fill = sleep_condition)) +
    geom_boxplot(alpha = 0.7) +
    labs(title = paste("Genotype:", genotype, "(Outliers Removed)"),
         x = "Sleep Condition", y = "Absolute Difference (rested - SD)") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_fill_manual(values = c("rested" = "blue", "SD" = "red"))
  
  return(p_boxplot)
}

# List of genotypes to plot
genotypes <- unique(data$genotype)

# Create a list to store plots
plots <- list()

# Create pairwise comparison plots for each genotype and store in the list
for (genotype in genotypes) {
  plot <- plot_pairwise_comparison(genotype)
  plots[[genotype]] <- plot
}

# Display boxplots
grid.arrange(grobs = plots, ncol = 1)

# Perform post-hoc comparisons for specific genotypes with adjusted p-values
posthoc_comparisons <- emmeans(model, specs = pairwise ~ sleep_condition | genotype)
posthoc_comparisons_df <- as.data.frame(posthoc_comparisons$contrasts)
posthoc_comparisons_df$p.value.original <- posthoc_comparisons_df$p.value
posthoc_comparisons_df$p.value.adjusted <- p.adjust(posthoc_comparisons_df$p.value, method = "BH")

# Display post-hoc comparison results with original and adjusted p-values
print(posthoc_comparisons_df)