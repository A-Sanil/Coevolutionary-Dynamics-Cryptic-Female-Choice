# R script to plot RSC vs Mates comparison over generations
# This verifies that RSC (lambda) correctly relates to mean mates per female

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

# Plot 1: RSC and Mates evolution overlaid (all replicates)
p1 <- df %>%
  select(Generation, Rep, MeanRSC, MeanMates) %>%
  pivot_longer(cols = c(MeanRSC, MeanMates), names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = ifelse(Metric == "MeanMates", "Mean Mates", "Mean RSC (Lambda)")) %>%
  ggplot(aes(x = Generation, y = Value, color = Metric, linetype = Metric)) +
  geom_line(alpha = 0.7, linewidth = 1) +
  facet_wrap(~Rep, scales = "free_y", ncol = 3) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "bottom",
    axis.text = element_text(size = 8)
  ) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  scale_color_manual(values = c("Mean RSC (Lambda)" = "darkgreen", "Mean Mates" = "darkblue")) +
  scale_linetype_manual(values = c("Mean RSC (Lambda)" = "solid", "Mean Mates" = "dashed")) +
  labs(
    title = "RSC (Lambda) vs Mean Mates Evolution by Replicate",
    subtitle = "RSC should track closely with Mean Mates since RSC = Poisson lambda\nSolid line: RSC | Dashed line: Mean Mates",
    x = "Generation",
    y = "Value",
    color = "Metric",
    linetype = "Metric"
  )

# Plot 2: Summary across all replicates (mean ± SE)
summary_data <- df %>%
  group_by(Generation) %>%
  summarise(
    MeanRSC = mean(MeanRSC, na.rm=TRUE),
    SE_RSC = sd(MeanRSC, na.rm=TRUE) / sqrt(n()),
    MeanMates = mean(MeanMates, na.rm=TRUE),
    SE_Mates = sd(MeanMates, na.rm=TRUE) / sqrt(n()),
    .groups = 'drop'
  )

p2 <- ggplot(summary_data, aes(x = Generation)) +
  geom_line(aes(y = MeanRSC, color = "Mean RSC (Lambda)"), linewidth = 1.2) +
  geom_ribbon(aes(ymin = MeanRSC - SE_RSC, ymax = MeanRSC + SE_RSC), 
              fill = "darkgreen", alpha = 0.2) +
  geom_line(aes(y = MeanMates, color = "Mean Mates"), linewidth = 1.2, linetype = "dashed") +
  geom_ribbon(aes(ymin = MeanMates - SE_Mates, ymax = MeanMates + SE_Mates), 
              fill = "darkblue", alpha = 0.2) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.position = "right"
  ) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  scale_color_manual(values = c("Mean RSC (Lambda)" = "darkgreen", "Mean Mates" = "darkblue")) +
  labs(
    title = "RSC (Lambda) vs Mean Mates: Summary Across All Replicates",
    subtitle = "Mean ± SE across replicates\nRSC (solid) should approximate Mean Mates (dashed) since RSC = Poisson lambda",
    x = "Generation",
    y = "Value",
    color = "Metric"
  )

# Plot 3: Scatter plot - RSC vs Mates (should be approximately y=x)
p3 <- ggplot(df, aes(x = MeanRSC, y = MeanMates)) +
  geom_point(alpha = 0.4, size = 1, color = "steelblue") +
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
    subtitle = "Red line: linear fit | Black dotted line: y=x (theoretical perfect match)\nSince RSC = lambda, mates should approximately equal RSC on average",
    caption = "Each point = one generation from one replicate"
  )

# Plot 4: Correlation over time
corr_data <- df %>%
  group_by(Generation) %>%
  summarise(
    Correlation = cor(MeanRSC, MeanMates, use = "complete.obs"),
    .groups = 'drop'
  )

p4 <- ggplot(corr_data, aes(x = Generation, y = Correlation)) +
  geom_line(color = "darkgreen", linewidth = 1.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 1) +
  geom_hline(yintercept = 0.95, linetype = "dotted", color = "orange", linewidth = 0.8) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  ) +
  scale_y_continuous(limits = c(0, 1.05), labels = scales::number_format(accuracy = 0.01)) +
  labs(
    title = "Correlation between RSC and Mates over Time",
    subtitle = "Red dashed line at y=1 (perfect correlation)\nOrange dotted line at y=0.95 (high correlation threshold)",
    x = "Generation",
    y = "Correlation (RSC, Mates)"
  )

# Plot 5: Difference (RSC - Mates) over time to see deviation
p5 <- df %>%
  mutate(Difference = MeanRSC - MeanMates) %>%
  group_by(Generation) %>%
  summarise(
    MeanDiff = mean(Difference, na.rm=TRUE),
    SE_Diff = sd(Difference, na.rm=TRUE) / sqrt(n()),
    .groups = 'drop'
  ) %>%
  ggplot(aes(x = Generation, y = MeanDiff)) +
  geom_line(color = "purple", linewidth = 1.2) +
  geom_ribbon(aes(ymin = MeanDiff - SE_Diff, ymax = MeanDiff + SE_Diff), 
              fill = "purple", alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  ) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  labs(
    title = "Difference: RSC - Mean Mates over Time",
    subtitle = "Mean difference across replicates (should be close to 0)\nPositive = RSC > Mates | Negative = RSC < Mates",
    x = "Generation",
    y = "Mean Difference (RSC - Mates)"
  )

# Save all plots
ggsave('Plots/rsc_vs_mates_by_replicate.png', p1, width = 14, height = 10, dpi = 300)
ggsave('Plots/rsc_vs_mates_summary.png', p2, width = 10, height = 6, dpi = 300)
ggsave('Plots/rsc_vs_mates_scatter.png', p3, width = 10, height = 7, dpi = 300)
ggsave('Plots/rsc_vs_mates_correlation.png', p4, width = 10, height = 6, dpi = 300)
ggsave('Plots/rsc_vs_mates_difference.png', p5, width = 10, height = 6, dpi = 300)

cat("\n=== Statistical Analysis ===\n")
cat("\nOverall correlation (all data):\n")
overall_cor <- cor(df$MeanRSC, df$MeanMates, use = "complete.obs")
cat("Correlation:", round(overall_cor, 4), "\n")

cat("\nLinear model: Mates ~ RSC\n")
model <- lm(MeanMates ~ MeanRSC, data = df)
print(summary(model))

cat("\nExpected: slope ≈ 1, intercept ≈ 0\n")
cat("Actual slope:", round(coef(model)[2], 4), "\n")
cat("Actual intercept:", round(coef(model)[1], 4), "\n")

cat("\n=== Final Generation Statistics ===\n")
final_gen <- max(df$Generation)
final_data <- df %>% filter(Generation == final_gen)
cat("Final generation:", final_gen, "\n")
cat("Mean RSC (final):", round(mean(final_data$MeanRSC, na.rm=TRUE), 4), "\n")
cat("Mean Mates (final):", round(mean(final_data$MeanMates, na.rm=TRUE), 4), "\n")
cat("Mean Difference (RSC - Mates):", round(mean(final_data$MeanRSC - final_data$MeanMates, na.rm=TRUE), 4), "\n")
cat("Correlation (final gen):", round(cor(final_data$MeanRSC, final_data$MeanMates, use = "complete.obs"), 4), "\n")

cat("\n=== Plots Saved to Plots/ directory ===\n")
cat("1. rsc_vs_mates_by_replicate.png - Side-by-side comparison by replicate\n")
cat("2. rsc_vs_mates_summary.png - Summary with SE across replicates\n")
cat("3. rsc_vs_mates_scatter.png - Scatter plot with y=x reference line\n")
cat("4. rsc_vs_mates_correlation.png - Correlation over time\n")
cat("5. rsc_vs_mates_difference.png - Difference (RSC - Mates) over time\n")

