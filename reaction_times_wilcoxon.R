# Read the CSV file
data <- read.csv("/home/joeh/Downloads/sleep_wake_gbd_reaction_times.csv")

# Perform Wilcoxon rank-sum tests
wake_vs_sd <- wilcox.test(
  data$Response.Time..ms.[data$cond == "WAKE"],
  data$Response.Time..ms.[data$cond == "SD"]
)

wake_vs_gbd <- wilcox.test(
  data$Response.Time..ms.[data$cond == "WAKE"],
  data$Response.Time..ms.[data$cond == "GBD"]
)

# Extract p-values
p_values <- c(wake_vs_sd$p.value, wake_vs_gbd$p.value)

# Apply FDR correction
adjusted_p_values <- p.adjust(p_values, method = "BH")

# Create a results data frame
results <- data.frame(
  Comparison = c("WAKE vs SD", "WAKE vs GBD"),
  W_statistic = c(wake_vs_sd$statistic, wake_vs_gbd$statistic),
  Raw_P_Value = p_values,
  Adjusted_P_Value = adjusted_p_values
)

# Print results
print("Wilcoxon Rank-Sum Test Results:")
print(results)

# Optional: Save results to a CSV file
write.csv(results, "sleep_wake_analysis_results_nonparametric.csv", row.names = FALSE)

# Additional summary statistics
print("\nSummary of Response Times by Condition:")
print(tapply(data$Response.Time..ms., data$cond, summary))

print("\nSample Sizes:")
print(table(data$cond))

# Calculate and print effect sizes (r = Z / sqrt(N))
n_wake <- sum(data$cond == "WAKE")
n_sd <- sum(data$cond == "SD")
n_gbd <- sum(data$cond == "GBD")

effect_size_wake_sd <- abs(qnorm(wake_vs_sd$p.value/2)) / sqrt(n_wake + n_sd)
effect_size_wake_gbd <- abs(qnorm(wake_vs_gbd$p.value/2)) / sqrt(n_wake + n_gbd)

print("\nEffect Sizes (r):")
print(paste("WAKE vs SD:", round(effect_size_wake_sd, 3)))
print(paste("WAKE vs GBD:", round(effect_size_wake_gbd, 3)))

# Calculate IQR for each condition
iqr_wake <- IQR(data$Response.Time..ms.[data$cond == "WAKE"])
iqr_sd <- IQR(data$Response.Time..ms.[data$cond == "SD"])
iqr_gbd <- IQR(data$Response.Time..ms.[data$cond == "GBD"])

print("\nInterquartile Ranges:")
print(paste("WAKE IQR:", round(iqr_wake, 2)))
print(paste("SD IQR:", round(iqr_sd, 2)))
print(paste("GBD IQR:", round(iqr_gbd, 2)))

# Generate summary for paper
print("\nSummary for paper:")
cat(sprintf("Response times were significantly different between WAKE and SD conditions (Wilcoxon rank-sum test: W = %d, p = %.3f, adjusted p = %.3f, r = %.3f), as well as between WAKE and GBD conditions (W = %d, p = %.3f, adjusted p = %.3f, r = %.3f). Median response times were shortest in the WAKE condition (Mdn = %.2f ms, IQR = %.2f ms, n = %d), followed b