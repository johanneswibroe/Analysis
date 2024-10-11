# Load required libraries
library(ordinal)
library(dplyr)
library(ggplot2)
library(tidyr)
library(emmeans)
library(nortest)  # for Anderson-Darling normality test
library(effsize)  # for effect size calculation

# Import and prepare data
df <- read.csv("/home/joeh/Volume_threshold/still_awake_CS_vs_sleeping_CS/still_awake_cs_overview3.csv")
df$ReactionThreshold <- factor(df$ReactionThreshold, 
                               levels = c("2", "5", "10", "30", "100", "na"),
                               labels = c("50dB", "60dB", "70dB", "80dB", "90dB", "No response"),
                               ordered = TRUE)
df$SleepState <- factor(df$SleepState)
df$SoundType <- factor(df$SoundType)
df$ReactionThresholdNumeric <- as.numeric(as.character(factor(df$ReactionThreshold,
                                                              levels = c("50dB", "60dB", "70dB", "80dB", "90dB", "No response"),
                                                              labels = c(50, 60, 70, 80, 90, 100))))

# Check for normality
cat("\nNormality Tests:\n")
print(ad.test(df$ReactionThresholdNumeric))
for(sleep in levels(df$SleepState)) {
  for(sound in levels(df$SoundType)) {
    subset <- df$ReactionThresholdNumeric[df$SleepState == sleep & df$SoundType == sound]
    cat("\nAnderson-Darling test for", sleep, "and", sound, ":\n")
    print(ad.test(subset))
  }
}

# Check for homogeneity of variances using Fligner-Killeen test
cat("\nFligner-Killeen Test for Homogeneity of Variances:\n")
print(fligner.test(ReactionThresholdNumeric ~ interaction(SleepState, SoundType), data = df))

# Non-parametric tests
cat("\nKruskal-Wallis Test for overall differences:\n")
print(kruskal.test(ReactionThresholdNumeric ~ interaction(SleepState, SoundType), data = df))

cat("\nWilcoxon Rank Sum Tests:\n")

# Conduct pairwise Wilcoxon tests for each SoundType
pairwise_tests <- list()
sound_types <- levels(df$SoundType)

for (sound in sound_types) {
  test_result <- wilcox.test(ReactionThresholdNumeric ~ SleepState, data = df[df$SoundType == sound,])
  pairwise_tests[[sound]] <- test_result$p.value
}

# Combine the p-values to apply FDR correction only to the relevant comparisons
fdr_corrected_p_values <- p.adjust(unlist(pairwise_tests), method = "fdr")

# Combine results into a data frame for easier reading
results <- data.frame(
  Comparison = names(pairwise_tests),
  RawPValue = unlist(pairwise_tests),
  FDRCorrectedPValue = fdr_corrected_p_values
)
print(results)

# Ordinal logistic regression (this is appropriate for ordinal data regardless of normality)
model_main <- clm(ReactionThreshold ~ SleepState + SoundType, data = df)
cat("\nMain Effects Model Summary:\n")
print(summary(model_main))

model_interaction <- clm(ReactionThreshold ~ SleepState * SoundType, data = df)
emm <- emmeans(model_interaction, specs = ~ SleepState | SoundType)
within_sound <- contrast(emm, method = "pairwise", adjust = "fdr")
cat("\nPost-hoc Comparisons (FDR-corrected):\n")
print(within_sound)

# Visualization
plot_theme <- theme_minimal() +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18, face = "bold"),
    plot.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16, face = "bold")
  )

violin_plot <- ggplot(df, aes(x = interaction(SleepState, SoundType), y = ReactionThresholdNumeric, fill = SoundType)) +
  geom_violin(alpha = 0.7, width = 1.2) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.5, alpha = 0.5, binwidth = 1.5) +
  scale_y_continuous(breaks = c(50, 60, 70, 80, 90, 100), labels = c("50dB", "60dB", "70dB", "80dB", "90dB", "No response")) +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "Distribution of Reaction Thresholds",
       x = "Sleep State and Sound Type",
       y = "Reaction Threshold") +
  plot_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_cartesian(clip = "off")

print(violin_plot)
ggsave("reaction_threshold_violin_plot.png", plot = violin_plot, width = 14, height = 8)

# Additional data exploration
cat("\nCounts of each ReactionThreshold:\n")
print(table(df$ReactionThreshold))

cat("\nContingency table of SleepState and SoundType:\n")
print(table(df$SleepState, df$SoundType))

cat("\nProportion of 'No response' responses by SleepState and SoundType:\n")
df %>%
  group_by(SleepState, SoundType) %>%
  summarize(prop_no_response = mean(ReactionThreshold == "No response", na.rm = TRUE)) %>%
  print()

# Additional test: Compare noise vs courtship for sleeping flies
cat("\nWilcoxon Rank Sum Test for Noise vs Courtship in Sleeping Flies:\n")
sleeping_data <- df[df$SleepState == "Sleeping", ]
wilcox_result <- wilcox.test(ReactionThresholdNumeric ~ SoundType, data = sleeping_data)
fdr_corrected_p_value <- p.adjust(wilcox_result$p.value, method = "fdr")
cat("Raw p-value:", wilcox_result$p.value, "\n")
cat("FDR-corrected p-value:", fdr_corrected_p_value, "\n")

# Effect size calculation (Cliff's delta)
cliff_delta <- cliff.delta(sleeping_data$ReactionThresholdNumeric[sleeping_data$SoundType == "Noise"],
                           sleeping_data$ReactionThresholdNumeric[sleeping_data$SoundType == "Courtship"])
cat("\nEffect Size (Cliff's delta) for Noise vs Courtship in Sleeping Flies:\n")
print(cliff_delta)

# Visualization for the additional comparison
sleeping_plot <- ggplot(sleeping_data, aes(x = SoundType, y = ReactionThresholdNumeric, fill = SoundType)) +
  geom_violin(alpha = 0.7, width = 0.8) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.5, alpha = 0.5, binwidth = 1.5) +
  scale_y_continuous(breaks = c(50, 60, 70, 80, 90, 100), labels = c("50dB", "60dB", "70dB", "80dB", "90dB", "No response")) +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "Reaction Thresholds for Sleeping Flies: Noise vs Courtship",
       x = "Sound Type",
       y = "Reaction Threshold") +
  plot_theme

print(sleeping_plot)
ggsave("sleeping_flies_noise_vs_courtship_plot.png", plot = sleeping_plot, width = 10, height = 8)