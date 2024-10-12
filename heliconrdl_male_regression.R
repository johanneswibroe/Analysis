# Load necessary libraries
library(readr)
library(ggplot2)
library(gghalves)
library(gridExtra)
library(dplyr)
library(lme4)
library(emmeans)

# Read the CSV file
data <- read.csv("/home/joeh/rdl_helicon_topspeeds_and_increase_decrease.csv", header = TRUE)

# Filter out rows with missing values
data <- data[complete.cases(data), ]

# Convert sleep_condition and genotype to factors
data$sleep_condition <- factor(data$sleep_condition)
data$genotype <- factor(data$genotype)

# Function to remove outliers
remove_outliers <- function(x) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = TRUE)
  H <- 1.5 * IQR(x, na.rm = TRUE)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

# Remove outliers for each group
data <- data %>%
  group_by(genotype, sleep_condition) %>%
  mutate(top_speed_difference = remove_outliers(top_speed_difference)) %>%
  ungroup()

# Remove NA values created by outlier removal
data <- data[complete.cases(data), ]

# Function to check normality and homogeneity of variance
check_assumptions <- function(data) {
  # Shapiro-Wilk test for each group
  normality_tests <- data %>%
    group_by(genotype, sleep_condition) %>%
    summarise(
      shapiro_p = shapiro.test(top_speed_difference)$p.value,
      .groups = "drop"
    )
  
  # Bartlett's test for homogeneity of variance
  bartlett_test <- bartlett.test(top_speed_difference ~ interaction(genotype, sleep_condition), data = data)
  
  list(normality = normality_tests, homogeneity = bartlett_test)
}

# Perform assumption checks
assumptions <- check_assumptions(data)
print("Normality tests:")
print(assumptions$normality)
print("Homogeneity of variance test:")
print(assumptions$homogeneity)

# Determine if we should use parametric or non-parametric tests
use_parametric <- all(assumptions$normality$shapiro_p > 0.05) && assumptions$homogeneity$p.value > 0.05

if (use_parametric) {
  # Fit a linear model
  model <- lm(top_speed_difference ~ sleep_condition * genotype, data = data)
  print(summary(model))
  
  # Post-hoc analysis
  emmeans_results <- emmeans(model, specs = pairwise ~ sleep_condition | genotype)
  emmeans_results_df <- as.data.frame(emmeans_results$contrasts)
} else {
  # Non-parametric analysis
  print("Using non-parametric tests due to violations of normality or homogeneity of variance.")
  
  # Kruskal-Wallis test
  kruskal_result <- kruskal.test(top_speed_difference ~ interaction(sleep_condition, genotype), data = data)
  print(kruskal_result)
  
  # Pairwise Wilcoxon tests
  pairwise_results <- pairwise.wilcox.test(data$top_speed_difference, 
                                           interaction(data$sleep_condition, data$genotype), 
                                           p.adjust.method = "none")
  
  # Convert pairwise results to a dataframe
  emmeans_results_df <- as.data.frame(pairwise_results$p.value)
  emmeans_results_df$comparison <- rownames(emmeans_results_df)
  emmeans_results_df <- tidyr::pivot_longer(emmeans_results_df, 
                                            cols = -comparison, 
                                            names_to = "group2", 
                                            values_to = "p.value")
  emmeans_results_df <- emmeans_results_df[!is.na(emmeans_results_df$p.value), ]
}

# Adjust p-values
emmeans_results_df$p.value.original <- emmeans_results_df$p.value
emmeans_results_df$p.value.adjusted <- p.adjust(emmeans_results_df$p.value, method = "fdr")

# Display the results with original and adjusted p-values
print(emmeans_results_df)

# Calculate significance levels
sig_levels <- emmeans_results_df %>%
  mutate(significance = case_when(
    p.value.adjusted < 0.001 ~ "***",
    p.value.adjusted < 0.01 ~ "**",
    p.value.adjusted < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Function to create raincloud plots for each genotype
plot_raincloud <- function(genotype) {
  # Subset data for the current genotype
  genotype_data <- data[data$genotype == genotype, ]
  
  # Get significance for this genotype
  sig <- sig_levels$significance[sig_levels$genotype == genotype]
  
  # Create the raincloud plot
  p <- ggplot(genotype_data, aes(x = genotype, y = top_speed_difference, fill = sleep_condition)) +
    geom_half_violin(side = "l", position = position_nudge(x = -.2), adjust = 1.5, trim = FALSE, scale = "width") +
    geom_point(aes(color = sleep_condition), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
    geom_boxplot(width = .1, outlier.shape = NA, alpha = 0.5) +
    labs(title = paste(genotype),
         x = "Genotype", y = "Top Speed Difference") +
    scale_fill_manual(values = c("sd" = "red", "rested" = "blue")) +
    scale_color_manual(values = c("sd" = "red", "rested" = "blue")) +
    annotate("text", x = 1, y = max(genotype_data$top_speed_difference), label = sig, size = 8) +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
  
  return(p)
}

# Create and display plots
plots <- list()
for (genotype in levels(data$genotype)) {
  plot <- plot_raincloud(genotype)
  plots[[genotype]] <- plot
}
grid.arrange(grobs = plots, ncol = 2)

# Print summary statistics
summary_stats <- data %>%
  group_by(genotype, sleep_condition) %>%
  summarise(
    mean = mean(top_speed_difference),
    sd = sd(top_speed_difference),
    median = median(top_speed_difference),
    n = n(),
    .groups = "drop"
  )
print(summary_stats)