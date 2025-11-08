# R script to plot simulation results with no clamping (RSC = direct Poisson lambda)
# Install packages once: install.packages(c('tidyverse'))

library(tidyverse)

# Find the most recent simulation results file
parallel_files <- list.files(pattern = "parallel_sim_results_.*\\.csv$", path = ".", full.names = TRUE)

if(length(parallel_files) > 0) {
  file <- parallel_files[order(parallel_files, decreasing = TRUE)][1]
  cat("Using simulation file:", file, "\n")
} else {
  stop("No simulation results CSV found in this directory")
}

# Read data
df <- read_csv(file)

# Convert to appropriate types
df <- df %>% 
  mutate(
    Generation = as.integer(Generation),
    Rep = as.factor(Rep)
  )

cat("\n=== Data Summary ===\n")
cat("Total generations:", max(df$Generation), "\n")
cat("Total replicates:", length(unique(df$Rep)), "\n")
cat("Mean RSC range:", round(min(df$MeanRSC, na.rm=TRUE), 3), "to", round(max(df$MeanRSC, na.rm=TRUE), 3), "\n")
cat("Mean mates range:", round(min(df$MeanMates, na.rm=TRUE), 3), "to", round(max(df$MeanMates, na.rm=TRUE), 3), "\n")

# Create output directory for plots
if(!dir.exists("Plots")) {
  dir.create("Plots")
}

# Plot 1: RSC Evolution (Mean RSC across generations)
p1 <- ggplot(df, aes(x = Generation, y = MeanRSC, color = Rep)) +
  geom_line(linewidth = 1.2, alpha = 0.8) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.position = "right"
  ) +
  scale_y_continuous(
    name = "Mean RSC Phenotype (Lambda)",
    labels = scales::number_format(accuracy = 0.01)
  ) +
  labs(
    title = 'Evolution of Mean RSC (Risk of Sperm Competition)',
    subtitle = 'RSC directly serves as lambda (mean) for Poisson distribution\nNo clamping - RSC evolves freely',
    x = 'Generation',
    y = 'Mean RSC (Lambda)',
    color = 'Replicate'
  )

# Plot 2: Mean Mates Evolution
p2 <- ggplot(df, aes(x = Generation, y = MeanMates, color = Rep)) +
  geom_line(linewidth = 1.2, alpha = 0.8) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.position = "right"
  ) +
  scale_y_continuous(
    name = "Mean Mates per Female",
    labels = scales::number_format(accuracy = 0.01)
  ) +
  labs(
    title = 'Evolution of Mean Mates per Female',
    subtitle = 'Number of mates sampled from Poisson(RSC)\nNo clamping - can be 0 or any positive value',
    x = 'Generation',
    y = 'Mean Mates per Female',
    color = 'Replicate'
  )

# Plot 3: RSC vs Mates Relationship (should be approximately linear since RSC = lambda)
p3 <- ggplot(df, aes(x = MeanRSC, y = MeanMates)) +
  geom_point(alpha = 0.5, size = 1.5, color = "steelblue") +
  geom_smooth(method = "lm", se = TRUE, color = "red", linewidth = 1.2, linetype = "dashed") +
  geom_abline(intercept = 0, slope = 1, color = "black", linewidth = 1, linetype = "dotted") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  ) +
  scale_x_continuous(name = "Mean RSC (Lambda)", labels = scales::number_format(accuracy = 0.1)) +
  scale_y_continuous(name = "Mean Mates per Female", labels = scales::number_format(accuracy = 0.1)) +
  labs(
    title = "Relationship: RSC (Lambda) vs Mean Mates",
    subtitle = "Red line: linear fit | Black dotted line: y=x (theoretical)\nSince RSC = lambda, mates should approximate RSC",
    caption = "Points show each generation across all replicates"
  )

# Plot 4: RSC and Mates evolution side-by-side (faceted by replicate)
p4 <- df %>%
  select(Generation, Rep, MeanRSC, MeanMates) %>%
  pivot_longer(cols = c(MeanRSC, MeanMates), names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = ifelse(Metric == "MeanMates", "Mean Mates", "Mean RSC (Lambda)")) %>%
  ggplot(aes(x = Generation, y = Value, color = Metric)) +
  geom_line(alpha = 0.8, linewidth = 1) +
  facet_wrap(~Rep, scales = "free_y", ncol = 3) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "bottom",
    axis.text = element_text(size = 8)
  ) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  labs(
    title = "Evolution of RSC and Mates by Replicate",
    subtitle = "RSC (lambda) should track closely with Mean Mates",
    x = "Generation",
    y = "Value",
    color = "Metric"
  )

# Plot 5: Summary statistics across replicates
summary_data <- df %>%
  group_by(Generation) %>%
  summarise(
    MeanRSC = mean(MeanRSC, na.rm=TRUE),
    SE_RSC = sd(MeanRSC, na.rm=TRUE) / sqrt(n()),
    MeanMates = mean(MeanMates, na.rm=TRUE),
    SE_Mates = sd(MeanMates, na.rm=TRUE) / sqrt(n()),
    .groups = 'drop'
  )

p5 <- summary_data %>%
  pivot_longer(cols = c(MeanRSC, MeanMates), names_to = "Metric", values_to = "Value") %>%
  mutate(
    SE = ifelse(Metric == "MeanRSC", summary_data$SE_RSC, summary_data$SE_Mates),
    Metric = ifelse(Metric == "MeanMates", "Mean Mates", "Mean RSC (Lambda)")
  ) %>%
  ggplot(aes(x = Generation, y = Value, color = Metric, fill = Metric)) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = Value - SE, ymax = Value + SE), alpha = 0.2, color = NA) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.position = "right"
  ) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  labs(
    title = "Summary: RSC and Mates Evolution (Mean ± SE)",
    subtitle = "Average across all replicates with standard error",
    x = "Generation",
    y = "Value",
    color = "Metric",
    fill = "Metric"
  )

# Plot 6: Correlation over time
p6 <- df %>%
  group_by(Generation) %>%
  summarise(
    Correlation = cor(MeanRSC, MeanMates, use = "complete.obs"),
    .groups = 'drop'
  ) %>%
  ggplot(aes(x = Generation, y = Correlation)) +
  geom_line(color = "darkgreen", linewidth = 1.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 1) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  ) +
  scale_y_continuous(limits = c(0, 1.1), labels = scales::number_format(accuracy = 0.01)) +
  labs(
    title = "Correlation between RSC and Mates over Time",
    subtitle = "Red dashed line at y=1 (perfect correlation)\nShould be high since RSC = lambda for Poisson",
    x = "Generation",
    y = "Correlation (RSC, Mates)"
  )

# Save all plots
ggsave('Plots/rsc_evolution_no_clamp.png', p1, width = 10, height = 6, dpi = 300)
ggsave('Plots/mates_evolution_no_clamp.png', p2, width = 10, height = 6, dpi = 300)
ggsave('Plots/rsc_vs_mates_relationship.png', p3, width = 10, height = 7, dpi = 300)
ggsave('Plots/rsc_mates_by_replicate.png', p4, width = 14, height = 10, dpi = 300)
ggsave('Plots/rsc_mates_summary.png', p5, width = 10, height = 6, dpi = 300)
ggsave('Plots/rsc_mates_correlation.png', p6, width = 10, height = 6, dpi = 300)

cat("\n=== Statistical Analysis: RSC vs Mates ===\n")
cat("\nOverall relationship (all data):\n")
model_all <- lm(MeanMates ~ MeanRSC, data = df)
print(summary(model_all))

cat("\n=== Correlation Analysis ===\n")
overall_cor <- cor(df$MeanRSC, df$MeanMates, use = "complete.obs")
cat("Overall correlation (RSC, Mates):", round(overall_cor, 4), "\n")

cat("\n=== Final Generation Statistics ===\n")
final_gen <- max(df$Generation)
final_data <- df %>% filter(Generation == final_gen)
cat("Final generation:", final_gen, "\n")
cat("Mean RSC (final):", round(mean(final_data$MeanRSC, na.rm=TRUE), 4), "\n")
cat("Mean Mates (final):", round(mean(final_data$MeanMates, na.rm=TRUE), 4), "\n")
cat("Correlation (final gen):", round(cor(final_data$MeanRSC, final_data$MeanMates, use = "complete.obs"), 4), "\n")

cat("\n=== Plots Saved to Plots/ directory ===\n")
cat("1. rsc_evolution_no_clamp.png - RSC evolution across replicates\n")
cat("2. mates_evolution_no_clamp.png - Mates evolution across replicates\n")
cat("3. rsc_vs_mates_relationship.png - Scatter plot with linear fit\n")
cat("4. rsc_mates_by_replicate.png - Side-by-side by replicate\n")
cat("5. rsc_mates_summary.png - Summary with SE\n")
cat("6. rsc_mates_correlation.png - Correlation over time\n")

# Summary statistics table
cat("\n=== Summary Statistics by Replicate ===\n")
summary_table <- df %>%
  group_by(Rep) %>%
  summarize(
    mean_rsc = round(mean(MeanRSC, na.rm=TRUE), 3),
    final_rsc = round(last(MeanRSC), 3),
    rsc_range = paste(round(min(MeanRSC, na.rm=TRUE), 3), "to", round(max(MeanRSC, na.rm=TRUE), 3)),
    mean_mates = round(mean(MeanMates, na.rm=TRUE), 3),
    final_mates = round(last(MeanMates), 3),
    mates_range = paste(round(min(MeanMates, na.rm=TRUE), 3), "to", round(max(MeanMates, na.rm=TRUE), 3)),
    correlation = round(cor(MeanRSC, MeanMates, use = "complete.obs"), 3),
    .groups = 'drop'
  )
print(summary_table)

