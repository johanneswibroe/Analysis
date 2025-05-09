# Load necessary libraries
library(ggplot2)
library(ggdist)
library(dplyr)

# Load the data
data <- read.csv("/home/joeh/inc_vs_w1118/inc_vs_w1118_with_sleep_latency.csv")

# Clean the data: remove rows with missing or infinite values for `volume_threshold`
clean_data <- data %>%
  dplyr::filter(!is.na(volume_threshold) & is.finite(volume_threshold))

# Raincloud plot for volume threshold
ggplot(clean_data, aes(x = sleep_status, y = volume_threshold, fill = genotype)) +
  stat_halfeye(aes(color = genotype), adjust = .5, width = .6, .width = 0, justification = -.3, point_colour = NA) + 
  geom_boxplot(width = .25, outlier.shape = NA, alpha = 0.5) +
  geom_point(aes(color = genotype), position = position_jitter(width = .1), size = 1, alpha = 0.5) +
  scale_color_brewer(palette = "Set2") +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal() +
  labs(title = "Volume Threshold Distribution by Sleep Status and Genotype",
       y = "Volume Threshold")

# Reorder `sleep_status` factor levels explicitly
state_switch_prob <- state_switch_prob %>%
  dplyr::filter(!is.na(prob) & !is.na(se) & !is.na(sleep_status) & !is.na(genotype)) %>%
  mutate(sleep_status = factor(sleep_status, levels = c("awake_rested", "awake_sd", "sleeping")),
         genotype = factor(genotype, levels = c("w1118", "inc2")))

# Plot with corrected x-axis order
ggplot(state_switch_prob, aes(x = sleep_status, y = prob, fill = genotype)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = prob - se, ymax = prob + se), 
                position = position_dodge(width = 0.9), width = 0.2) +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal() +
  labs(title = "State-Switching Probability by Sleep Status and Genotype",
       y = "Probability of State Switch")

# Fisher's exact test for sleeping condition
sleeping_data <- data %>%
  dplyr::filter(sleep_status == "sleeping" & !is.na(state_switch)) %>%
  with(table(genotype, state_switch))
print("\nContingency table for sleeping condition:")
print(sleeping_data)
fisher_result_sleeping <- fisher.test(sleeping_data)
print("\nFisher's exact test result for sleeping condition:")
print(fisher_result_sleeping)

# Fisher's exact test for awake_rested condition
awake_rested_data <- data %>%
  dplyr::filter(sleep_status == "awake_rested" & !is.na(state_switch)) %>%
  with(table(genotype, state_switch))
print("\nContingency table for awake_rested condition:")
print(awake_rested_data)
fisher_result_awake_rested <- fisher.test(awake_rested_data)
print("\nFisher's exact test result for awake_rested condition:")
print(fisher_result_awake_rested)

# Fisher's exact test for awake_sd condition
awake_sd_data <- data %>%
  dplyr::filter(sleep_status == "awake_sd" & !is.na(state_switch)) %>%
  with(table(genotype, state_switch))
print("\nContingency table for awake_sd condition:")
print(awake_sd_data)
fisher_result_awake_sd <- fisher.test(awake_sd_data)
print("\nFisher's exact test result for awake_sd condition:")
print(fisher_result_awake_sd)

# NEW: Fisher's exact test for inc2 sleeping_gbd vs sleeping comparison
inc2_sleep_comparison <- data %>%
  dplyr::filter(genotype == "inc2", 
                sleep_status %in% c("sleeping_gbd", "sleeping") & !is.na(state_switch)) %>%
  with(table(sleep_status, state_switch))
print("\nContingency table for inc2 sleeping_gbd vs sleeping comparison:")
print(inc2_sleep_comparison)
fisher_result_inc2_sleep <- fisher.test(inc2_sleep_comparison)
print("\nFisher's exact test result for inc2 sleeping_gbd vs sleeping comparison:")
print(fisher_result_inc2_sleep)

# Print data structure check
print("\nData structure check (excluding NA):")
data_check <- data %>%
  dplyr::filter(!is.na(sleep_status) & !is.na(genotype) & !is.na(state_switch)) %>%
  group_by(sleep_status, genotype, state_switch) %>%
  summarize(count = n())
print(data_check)
