# Load required libraries
library(readr)
library(dplyr)
library(ggplot2)
library(stats)

# Read the CSV file
data <- read_csv("/home/joeh/wildtypes/wildtype_with_top_speed3.csv")

# Print data structure
str(data)

# Ensure light_cond and genotype are factors
data$light_cond <- as.factor(data$light_cond)
data$genotype <- as.factor(data$genotype)

# Print unique values
print(unique(data$genotype))
print(unique(data$light_cond))

# Function to remove outliers
remove_outliers <- function(x) {
  mean_x <- mean(x, na.rm = TRUE)
  sd_x <- sd(x, na.rm = TRUE)
  x[abs(x - mean_x) > 2*sd_x] <- NA
  return(x)
}

# Remove outliers for each group
data <- data %>%
  group_by(genotype, light_cond) %>%
  mutate(
    top_speed_difference_clean = remove_outliers(top_speed_difference),
    pre_onset_top_speed_clean = remove_outliers(pre_onset_top_speed),
    post_onset_top_speed_clean = remove_outliers(post_onset_top_speed)
  ) %>%
  ungroup()

# Function to perform Wilcoxon rank-sum test
perform_test <- function(data, x, y) {
  test <- wilcox.test(data[[y]] ~ data[[x]], data = data, na.action = na.omit)
  return(list(p_value = test$p.value, statistic = test$statistic, test_type = "Wilcoxon"))
}

# Function to perform Wilcoxon signed-rank test (for paired data)
perform_paired_test <- function(data, x, y) {
  test <- wilcox.test(data[[x]], data[[y]], paired = TRUE, na.action = na.omit)
  return(list(p_value = test$p.value, statistic = test$statistic, test_type = "Wilcoxon Signed-Rank"))
}

# Analysis 1: Dark vs Light for each genotype (top_speed_difference)
genotypes <- unique(data$genotype)
test_results_1 <- list()
for (genotype in genotypes) {
  subset_data <- data[data$genotype == genotype, ]
  test_results_1[[genotype]] <- tryCatch({
    perform_test(subset_data, "light_cond", "top_speed_difference_clean")
  }, error = function(e) {
    list(p_value = NA, statistic = NA, test_type = "Error")
  })
}

# Apply FDR correction
p_values_1 <- sapply(test_results_1, function(x) x$p_value)
fdr_corrected_1 <- p.adjust(p_values_1, method = "fdr")
for (i in seq_along(test_results_1)) {
  test_results_1[[i]]$fdr_corrected_p <- fdr_corrected_1[i]
}

# Print results for Analysis 1
cat("Analysis 1: Dark vs Light for each genotype (top_speed_difference)\n")
for (genotype in names(test_results_1)) {
  cat("Genotype:", genotype, "\n")
  cat("P-value:", test_results_1[[genotype]]$p_value, "\n")
  cat("FDR-corrected P-value:", test_results_1[[genotype]]$fdr_corrected_p, "\n\n")
}

# Analysis 2: Genotype comparisons in Light condition
light_data <- data[data$light_cond == "light", ]
kruskal_light <- kruskal.test(top_speed_difference_clean ~ genotype, data = light_data)
cat("Analysis 2: Genotype comparisons in Light condition\n")
print(kruskal_light)
if (kruskal_light$p.value < 0.05) {
  pairwise_light <- pairwise.wilcox.test(light_data$top_speed_difference_clean, 
                                         light_data$genotype, 
                                         p.adjust.method = "fdr")
  print(pairwise_light)
}

# Analysis 3: Genotype comparisons in Dark condition
dark_data <- data[data$light_cond == "dark", ]
kruskal_dark <- kruskal.test(top_speed_difference_clean ~ genotype, data = dark_data)
cat("\nAnalysis 3: Genotype comparisons in Dark condition\n")
print(kruskal_dark)
if (kruskal_dark$p.value < 0.05) {
  pairwise_dark <- pairwise.wilcox.test(dark_data$top_speed_difference_clean, 
                                        dark_data$genotype, 
                                        p.adjust.method = "fdr")
  print(pairwise_dark)
}

# Analysis 4: Pre vs Post Onset Top Speed for each genotype in dark and light conditions
test_results_4 <- list()
for (genotype in genotypes) {
  for (condition in c("dark", "light")) {
    subset_data <- data[data$genotype == genotype & data$light_cond == condition, ]
    test_key <- paste(genotype, condition, sep = "_")
    test_results_4[[test_key]] <- tryCatch({
      perform_paired_test(subset_data, "pre_onset_top_speed_clean", "post_onset_top_speed_clean")
    }, error = function(e) {
      list(p_value = NA, statistic = NA, test_type = "Error")
    })
  }
}

# Apply FDR correction for Analysis 4
p_values_4 <- sapply(test_results_4, function(x) x$p_value)
fdr_corrected_4 <- p.adjust(p_values_4, method = "fdr")
for (i in seq_along(test_results_4)) {
  test_results_4[[i]]$fdr_corrected_p <- fdr_corrected_4[i]
}

# Print results for Analysis 4
cat("\nAnalysis 4: Pre vs Post Onset Top Speed for each genotype in dark and light conditions\n")
for (test_key in names(test_results_4)) {
  cat("Genotype and Condition:", test_key, "\n")
  cat("P-value:", test_results_4[[test_key]]$p_value, "\n")
  cat("FDR-corrected P-value:", test_results_4[[test_key]]$fdr_corrected_p, "\n\n")
}

# Plot 1: Dark vs Light for each genotype (top_speed_difference)
plot1 <- ggplot(data, aes(x = light_cond, y = top_speed_difference_clean, fill = light_cond)) +
  geom_boxplot() +
  facet_wrap(~ genotype) +
  theme_bw() +
  labs(title = "Dark vs Light for each Genotype",
       x = "Light Condition",
       y = "Top Speed Difference") +
  scale_fill_manual(values = c("dark" = "gray", "light" = "yellow")) +
  coord_cartesian(ylim = c(min(data$top_speed_difference_clean, na.rm = TRUE),
                           max(data$top_speed_difference_clean, na.rm = TRUE)))

# Display Plot 1
print(plot1)

# Save Plot 1
ggsave("dark_vs_light_by_genotype.png", plot1, width = 12, height = 6, dpi = 300)

# Plot 2: Genotype comparisons in Light condition
plot2 <- ggplot(light_data, aes(x = genotype, y = top_speed_difference_clean, fill = genotype)) +
  geom_boxplot() +
  theme_bw() +
  labs(title = "Genotype Comparisons in Light Condition",
       x = "Genotype",
       y = "Top Speed Difference") +
  coord_cartesian(ylim = c(min(data$top_speed_difference_clean, na.rm = TRUE),
                           max(data$top_speed_difference_clean, na.rm = TRUE)))

# Display Plot 2
print(plot2)

# Save Plot 2
ggsave("genotype_comparisons_light.png", plot2, width = 8, height = 6, dpi = 300)

# Plot 3: Genotype comparisons in Dark condition
plot3 <- ggplot(dark_data, aes(x = genotype, y = top_speed_difference_clean, fill = genotype)) +
  geom_boxplot() +
  theme_bw() +
  labs(title = "Genotype Comparisons in Dark Condition",
       x = "Genotype",
       y = "Top Speed Difference") +
  coord_cartesian(ylim = c(min(data$top_speed_difference_clean, na.rm = TRUE),
                           max(data$top_speed_difference_clean, na.rm = TRUE)))

# Display Plot 3
print(plot3)

# Save Plot 3
ggsave("genotype_comparisons_dark.png", plot3, width = 8, height = 6, dpi = 300)

# Calculate summary statistics
summary_stats <- data %>%
  group_by(genotype, light_cond) %>%
  summarise(
    mean_diff = mean(top_speed_difference_clean, na.rm = TRUE),
    sd_diff = sd(top_speed_difference_clean, na.rm = TRUE),
    median_diff = median(top_speed_difference_clean, na.rm = TRUE),
    n_diff = sum(!is.na(top_speed_difference_clean)),
    mean_pre = mean(pre_onset_top_speed_clean, na.rm = TRUE),
    sd_pre = sd(pre_onset_top_speed_clean, na.rm = TRUE),
    median_pre = median(pre_onset_top_speed_clean, na.rm = TRUE),
    n_pre = sum(!is.na(pre_onset_top_speed_clean)),
    mean_post = mean(post_onset_top_speed_clean, na.rm = TRUE),
    sd_post = sd(post_onset_top_speed_clean, na.rm = TRUE),
    median_post = median(post_onset_top_speed_clean, na.rm = TRUE),
    n_post = sum(!is.na(post_onset_top_speed_clean)),
    .groups = 'drop'
  )

# Print summary statistics
print(summary_stats)