# Load required libraries
library(dplyr)
library(ggplot2)

# Read the CSV file
data <- read.csv("/home/joeh/helicon_split_RDL/ellens_flies/helicon_rdl_overview_with_everything.csv")

# Filter out individual_slep_lat over 2000 for ALL analyses
data <- data %>% filter(is.na(individual_slep_lat) | individual_slep_lat <= 2000)

# Calculate the percentage of '1' reactions per genotype and sample size
reaction_summary <- data %>%
  group_by(genotype) %>%
  summarise(
    percentage_1 = mean(reaction == 1, na.rm = TRUE) * 100,
    n = n()  # Count number of observations for each genotype
  )

# Perform Fisher's exact tests
# Create contingency tables for each comparison

# 1. helicon_rdl vs helicon_W1118
helicon_rdl_data <- data %>% filter(genotype == "helicon_rdl")
helicon_W1118_data <- data %>% filter(genotype == "helicon_w1118")

# Count the number of 1s and 0s in each group
helicon_rdl_pos <- sum(helicon_rdl_data$reaction == 1, na.rm = TRUE)
helicon_rdl_neg <- sum(helicon_rdl_data$reaction == 0, na.rm = TRUE)
helicon_W1118_pos <- sum(helicon_W1118_data$reaction == 1, na.rm = TRUE)
helicon_W1118_neg <- sum(helicon_W1118_data$reaction == 0, na.rm = TRUE)

contingency_table1 <- matrix(
  c(
    helicon_rdl_pos, helicon_rdl_neg,
    helicon_W1118_pos, helicon_W1118_neg
  ),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(
    c("helicon_rdl", "helicon_w1118"),
    c("Positive", "Negative")
  )
)

# 2. helicon_rdl vs rdl_W1118
rdl_W1118_data <- data %>% filter(genotype == "rdl_w1118")

# Count the number of 1s and 0s
rdl_W1118_pos <- sum(rdl_W1118_data$reaction == 1, na.rm = TRUE)
rdl_W1118_neg <- sum(rdl_W1118_data$reaction == 0, na.rm = TRUE)

contingency_table2 <- matrix(
  c(
    helicon_rdl_pos, helicon_rdl_neg,
    rdl_W1118_pos, rdl_W1118_neg
  ),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(
    c("helicon_rdl", "rdl_W1118"),
    c("Positive", "Negative")
  )
)

# Print the contingency tables to verify they're correct
print("Contingency Table 1: helicon_rdl vs helicon_w1118")
print(contingency_table1)
print("Contingency Table 2: helicon_rdl vs rdl_w1118")
print(contingency_table2)

# Run Fisher's exact tests
fisher_test1 <- fisher.test(contingency_table1)
fisher_test2 <- fisher.test(contingency_table2)

# Print test results
cat("Fisher's Exact Test: helicon_rdl vs helicon_w1118\n")
cat("p-value:", fisher_test1$p.value, "\n")
cat("odds ratio:", fisher_test1$estimate, "\n\n")

cat("Fisher's Exact Test: helicon_rdl vs rdl_w1118\n")
cat("p-value:", fisher_test2$p.value, "\n")
cat("odds ratio:", fisher_test2$estimate, "\n\n")

# Create a data frame to store significance markers
sig_markers <- data.frame(
  comparison = c("helicon_rdl vs helicon_w1118", "helicon_rdl vs rdl_w1118"),
  p_value = c(fisher_test1$p.value, fisher_test2$p.value)
)

# Format p-values for the plot
sig_markers$label <- ifelse(
  sig_markers$p_value < 0.001, "***",
  ifelse(sig_markers$p_value < 0.01, "**",
         ifelse(sig_markers$p_value < 0.05, "*", "ns")
  )
)

# Define colors for each genotype
genotype_colors <- c(
  "helicon_rdl" = "#1b8a5a",
  "helicon_w1118" = "#fbb021",
  "rdl_w1118" = "#ee3e32"
)

# Plot the results using ggplot2 with color coding
p <- ggplot(reaction_summary, aes(x = genotype, y = percentage_1, fill = genotype)) +
  geom_bar(stat = "identity") +
  # Add sample size and percentage text above each bar
  geom_text(aes(label = paste0(round(percentage_1, 1), "% (n=", n, ")")), 
            vjust = -0.5, size = 5, fontface = "bold") +
  scale_fill_manual(values = genotype_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +  # Remove space at bottom, keep space at top for labels
  labs(
    x = "Genotype",
    y = "Percentage of Positive Reactions (1) (%)",
    title = "Percentage of Positive Reactions (1) per Genotype",
    subtitle = paste0(
      "Fisher's test: helicon_rdl vs helicon_w1118 (p = ", 
      format(fisher_test1$p.value, digits = 3), 
      "), helicon_rdl vs rdl_W1118 (p = ", 
      format(fisher_test2$p.value, digits = 3), ")"
    )
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14, face = "bold"),
    axis.text.y = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12, face = "bold")
  ) +  # Rotate x-labels for readability
  guides(fill = "none")  # Remove legend since colors match x-axis labels

# Add significance brackets
max_percentage <- max(reaction_summary$percentage_1)
bracket_height <- max_percentage * 1.1

# Add significance annotations
if (fisher_test1$p.value < 0.05) {
  p <- p + annotate("segment", 
                    x = which(reaction_summary$genotype == "helicon_rdl"), 
                    xend = which(reaction_summary$genotype == "helicon_w1118"),
                    y = bracket_height, 
                    yend = bracket_height) +
    annotate("text", 
             x = mean(c(which(reaction_summary$genotype == "helicon_rdl"), 
                        which(reaction_summary$genotype == "helicon_w1118"))), 
             y = bracket_height * 1.05, 
             label = sig_markers$label[1])
}

if (fisher_test2$p.value < 0.05) {
  p <- p + annotate("segment", 
                    x = which(reaction_summary$genotype == "helicon_rdl"), 
                    xend = which(reaction_summary$genotype == "rdl_w1118"),
                    y = bracket_height * 1.15, 
                    yend = bracket_height * 1.15) +
    annotate("text", 
             x = mean(c(which(reaction_summary$genotype == "helicon_rdl"), 
                        which(reaction_summary$genotype == "rdl_w1118"))), 
             y = bracket_height * 1.2, 
             label = sig_markers$label[2])
}

# Print the plot
print(p)

# Create scatter plot with filtered data
# Filter data for scatter plot (already filtered individual_slep_lat <= 2000)
scatter_data <- data %>%
  filter(thresh <= 99.9) %>%
  filter(!is.na(thresh) & !is.na(individual_slep_lat))

# Create scatter plot
scatter_plot <- ggplot(scatter_data, aes(x = individual_slep_lat, y = thresh, color = genotype)) +
  geom_point(alpha = 0.7, size = 4) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2, aes(fill = genotype)) +
  scale_color_manual(values = genotype_colors) +
  scale_fill_manual(values = genotype_colors) +
  labs(
    x = "Individual Sleep Latency (seconds)",
    y = "Volume Threshold (%)",
    title = "Scatter Plot: Volume Threshold vs Individual Sleep Latency by Genotype",
    subtitle = "Filtered: thresh ≤ 99.9, individual_slep_lat ≤ 2000",
    color = "Genotype"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12, face = "bold")
  )

# Print the scatter plot
print(scatter_plot)

# Create scatter plot with single correlation line for all genotypes
scatter_plot_combined <- ggplot(scatter_data, aes(x = individual_slep_lat, y = thresh)) +
  geom_point(aes(color = genotype), alpha = 0.7, size = 4) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2, color = "black", fill = "gray") +
  scale_color_manual(values = genotype_colors) +
  labs(
    x = "Individual Sleep Latency (seconds)",
    y = "Volume Threshold (%)",
    title = "Scatter Plot: Volume Threshold vs Individual Sleep Latency (Combined Correlation)",
    subtitle = "Single correlation line for all genotypes combined",
    color = "Genotype"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12, face = "bold")
  )

# Print the combined scatter plot
print(scatter_plot_combined)

# Create scatter plot correlating thresh with time_simp
time_scatter_data <- data %>%
  filter(thresh <= 99.9) %>%
  filter(!is.na(thresh) & !is.na(time_simp))

time_scatter_plot <- ggplot(time_scatter_data, aes(x = time_simp, y = thresh, color = genotype)) +
  geom_point(alpha = 0.7, size = 4) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2, aes(fill = genotype)) +
  scale_color_manual(values = genotype_colors) +
  scale_fill_manual(values = genotype_colors) +
  scale_x_continuous(
    breaks = 1:7,
    labels = c("1pm", "2pm", "3pm", "4pm", "5pm", "6pm", "7pm")
  ) +
  labs(
    x = "Time of Day",
    y = "Volume Threshold (%)",
    title = "Scatter Plot: Volume Threshold vs Time of Day by Genotype",
    subtitle = "Correlation between volume threshold and time of day",
    color = "Genotype"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12, face = "bold")
  )

# Print the time scatter plot
print(time_scatter_plot)

# Check the data distribution for helicon_w1118 in the time scatter plot
cat("\nData distribution for time scatter plot:\n")
time_scatter_summary <- time_scatter_data %>%
  group_by(genotype) %>%
  summarise(
    n = n(),
    time_range = paste(range(time_simp, na.rm = TRUE), collapse = " - "),
    .groups = 'drop'
  )
print(time_scatter_summary)

# Visualize how time_simp influences likelihood of reaction
# Create binned time_simp to better visualize the relationship
reaction_time_data <- data %>%
  filter(!is.na(time_simp) & !is.na(reaction)) %>%
  mutate(time_bin = cut(time_simp, breaks = 10, include.lowest = TRUE))

# Calculate reaction probability by time bin and genotype
reaction_by_time <- reaction_time_data %>%
  group_by(time_bin, genotype) %>%
  summarise(
    reaction_rate = mean(reaction == 1, na.rm = TRUE) * 100,
    n = n(),
    .groups = 'drop'
  ) %>%
  filter(n >= 5)  # Only show bins with at least 5 observations

# Create plot showing reaction likelihood vs time_simp
reaction_time_plot <- ggplot(reaction_by_time, aes(x = time_bin, y = reaction_rate, fill = genotype)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
  geom_text(aes(label = paste0(round(reaction_rate, 1), "% (n=", n, ")")), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5, size = 3, fontface = "bold") +
  scale_fill_manual(values = genotype_colors) +
  labs(
    x = "Time of Day",
    y = "Reaction Rate (%)",
    title = "Reaction Likelihood vs Time of Day by Genotype",
    subtitle = "Percentage of positive reactions (1) across time periods",
    fill = "Genotype"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14, face = "bold"),
    axis.text.y = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.position = "bottom"
  )

# Print the reaction time plot
print(reaction_time_plot)

# Alternative: Smooth line plot showing reaction probability vs time_simp
reaction_smooth_plot <- ggplot(data %>% filter(!is.na(time_simp) & !is.na(reaction)), 
                               aes(x = time_simp, y = as.numeric(reaction == 1) * 100, color = genotype)) +
  geom_smooth(method = "loess", se = TRUE, alpha = 0.2, aes(fill = genotype)) +
  scale_color_manual(values = genotype_colors) +
  scale_fill_manual(values = genotype_colors) +
  scale_x_continuous(
    breaks = 1:7,
    labels = c("1pm", "2pm", "3pm", "4pm", "5pm", "6pm", "7pm")
  ) +
  labs(
    x = "Time of Day",
    y = "Reaction Probability (%)",
    title = "Reaction Probability vs Time of Day (Smooth Trends)",
    subtitle = "Probability of positive reaction across time of day",
    color = "Genotype"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12, face = "bold")
  )

# Print the smooth reaction plot
print(reaction_smooth_plot)

# Create scatter plot for helicon_rdl only
helicon_only_data <- scatter_data %>%
  filter(genotype == "helicon_rdl")

helicon_scatter_plot <- ggplot(helicon_only_data, aes(x = individual_slep_lat, y = thresh)) +
  geom_point(alpha = 0.7, size = 4, color = "#1b8a5a") +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2, color = "#1b8a5a", fill = "#1b8a5a") +
  labs(
    x = "Individual Sleep Latency (seconds)",
    y = "Volume Threshold (%)",
    title = "Scatter Plot: Volume Threshold vs Individual Sleep Latency (helicon_rdl only)",
    subtitle = "Correlation for helicon_rdl genotype only"
  ) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14, face = "bold")
  )

# Print the helicon-only scatter plot
print(helicon_scatter_plot)

# LOESS plot with unclipped error bands
smoothed_data <- data %>%
  filter(!is.na(time_simp), !is.na(reaction)) %>%
  group_by(genotype) %>%
  do({
    fit <- loess(as.numeric(reaction == 1) ~ time_simp, data = ., span = 0.75)
    time_seq <- seq(min(.$time_simp), max(.$time_simp), length.out = 200)
    pred <- predict(fit, newdata = data.frame(time_simp = time_seq), se = TRUE)
    tibble(
      time_simp = time_seq,
      fit = pmin(pmax(pred$fit, 0), 1),             # clamp fitted line only
      ymin = pred$fit - 1.96 * pred$se.fit,         # error bands not clipped
      ymax = pred$fit + 1.96 * pred$se.fit,
      genotype = unique(.$genotype)
    )
  })

# Plot with extended y-axis to show full error bands
reaction_smooth_plot <- ggplot(smoothed_data, aes(x = time_simp, y = fit, color = genotype, fill = genotype)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.2, color = NA) +
  geom_line(size = 1.2) +
  scale_color_manual(values = genotype_colors) +
  scale_fill_manual(values = genotype_colors) +
  scale_x_continuous(
    breaks = 1:7,
    labels = c("1pm", "2pm", "3pm", "4pm", "5pm", "6pm", "7pm")
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0.1, 0.1)),   # allow headroom/footroom
    oob = scales::oob_keep                   # keep CI lines even outside limits
  ) +
  labs(
    x = "Time of Day",
    y = "Reaction Probability (%)",
    title = "Reaction Probability vs Time of Day (LOESS with Unclipped Error Bands)",
    subtitle = "Fitted line constrained to 0-100%, full confidence bands shown",
    color = "Genotype",
    fill = "Genotype"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12, face = "bold")
  )

# Print the plot
print(reaction_smooth_plot)