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

# Function to check normality
check_normality <- function(data) {
  data %>%
    group_by(genotype, sleep_condition) %>%
    summarise(
      shapiro_p = shapiro.test(abs_difference_clean)$p.value,
      .groups = "drop"
    )
}

# Perform normality checks
normality_results <- check_normality(data)
print("Normality tests:")
print(normality_results)

# Determine if we should use parametric or non-parametric tests
use_parametric <- all(normality_results$shapiro_p > 0.05)

if (use_parametric) {
  print("Using parametric tests.")
  
  # Fit a linear model
  model <- lm(abs_difference_clean ~ sleep_condition * genotype, data = data)
  print(summary(model))
  
  # Post-hoc analysis
  emmeans_results <- emmeans(model, specs = pairwise ~ sleep_condition | genotype)
  posthoc_results <- as.data.frame(emmeans_results$contrasts)
  
  # Test for differences in rested condition
  rested_data <- data[data$sleep_condition == "rested", ]
  rested_model <- aov(abs_difference_clean ~ genotype, data = rested_data)
  rested_results <- TukeyHSD(rested_model)
  print("Differences in rested condition:")
  print(rested_results)
} else {
  print("Using non-parametric tests.")
  
  # Perform Wilcoxon tests for each genotype
  posthoc_results <- data %>%
    group_by(genotype) %>%
    summarise(
      p_value = wilcox.test(abs_difference_clean ~ sleep_condition)$p.value,
      .groups = "drop"
    )
  
  # Test for differences in rested condition
  rested_data <- data[data$sleep_condition == "rested", ]
  rested_results <- kruskal.test(abs_difference_clean ~ genotype, data = rested_data)
  print("Differences in rested condition:")
  print(rested_results)
  
  # If Kruskal-Wallis test is significant, perform pairwise Wilcoxon tests
  if (rested_results$p.value < 0.05) {
    rested_pairwise <- pairwise.wilcox.test(rested_data$abs_difference_clean, rested_data$genotype, p.adjust.method = "BH")
    print("Pairwise comparisons in rested condition:")
    print(rested_pairwise)
  }
}

# Adjust p-values for posthoc results
posthoc_results$p.value.adjusted <- p.adjust(posthoc_results$p_value, method = "BH")

# Display the results with original and adjusted p-values
print("Post-hoc results:")
print(posthoc_results)

# Calculate summary statistics
summary_stats <- data %>%
  group_by(genotype, sleep_condition) %>%
  summarise(
    mean = mean(abs_difference_clean, na.rm = TRUE),
    sd = sd(abs_difference_clean, na.rm = TRUE),
    median = median(abs_difference_clean, na.rm = TRUE),
    n = sum(!is.na(abs_difference_clean)),
    .groups = "drop"
  )
print("Summary statistics:")
print(summary_stats)

# Function to create pairwise comparison plot for each genotype
plot_pairwise_comparison <- function(genotype) {
  # Subset data for the current genotype
  genotype_data <- data[data$genotype == genotype, ]
  
  # Get p-value for this genotype
  p_value <- posthoc_results$p.value.adjusted[posthoc_results$genotype == genotype]
  sig <- if(p_value < 0.05) "*" else ""
  
  # Create the boxplot
  p_boxplot <- ggplot(genotype_data, aes(x = sleep_condition, y = abs_difference_clean, fill = sleep_condition)) +
    geom_boxplot(alpha = 0.7) +
    labs(title = paste("Genotype:", genotype),
         x = "Sleep Condition", 
         y = "Difference Pre-Post Stimulus (mm/s)") +
    annotate("text", x = 1.5, y = max(genotype_data$abs_difference_clean, na.rm = TRUE), 
             label = paste("p =", round(p_value, 3), sig), size = 5) +
    theme_minimal(base_size = 16) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 18),
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 16)
    ) +
    scale_fill_manual(values = c("rested" = "blue", "SD" = "red"))
  
  return(p_boxplot)
}

# Create and display plots
plots <- list()
for (genotype in levels(data$genotype)) {
  plot <- plot_pairwise_comparison(genotype)
  plots[[genotype]] <- plot
}
grid.arrange(grobs = plots, ncol = 3)
