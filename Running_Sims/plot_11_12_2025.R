# R script to plot all simulation data for 11/12/2025
# Comprehensive plotting including RSC vs mates, trait evolution, and other analyses
# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

# Get the most recent CSV file
csv_files <- list.files(pattern = "parallel_sim_results_.*\\.csv$", full.names = TRUE)
if (length(csv_files) == 0) {
  stop("No CSV files found! Make sure the simulation has run and created a CSV file.")
}

# Get the most recent file
latest_file <- csv_files[order(file.info(csv_files)$mtime, decreasing = TRUE)[1]]
cat("Reading data from:", latest_file, "\n")

# Read the data
df <- read.csv(latest_file)

# Create output directory for plots
# Using underscores for Windows compatibility (user requested "11/12/2025 plots")
plot_dir <- "11_12_2025_plots"
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir)
  cat("Created directory:", plot_dir, "\n")
}

# Convert Generation to numeric
df$Generation <- as.numeric(df$Generation)
df$Rep <- as.factor(df$Rep)

cat("\n=== Data Summary ===\n")
cat("Total generations:", max(df$Generation), "\n")
cat("Total replicates:", length(unique(df$Rep)), "\n")
cat("Mean RSC range:", round(min(df$MeanRSC, na.rm=TRUE), 3), "to", round(max(df$MeanRSC, na.rm=TRUE), 3), "\n")
cat("Mean mates range:", round(min(df$MeanMates, na.rm=TRUE), 3), "to", round(max(df$MeanMates, na.rm=TRUE), 3), "\n")

# =============================================================================
# PLOT 1: RSC vs Mates Scatter Plot
# =============================================================================
p1 <- ggplot(df, aes(x = MeanRSC, y = MeanMates)) +
  geom_point(alpha = 0.5, color = "steelblue", size = 1.5) +
  geom_smooth(method = "lm", se = TRUE, color = "red", linewidth = 1.2, linetype = "dashed") +
  geom_abline(intercept = 0, slope = 1, color = "black", linewidth = 1, linetype = "dotted") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  ) +
  scale_x_continuous(name = "Mean RSC (Lambda)", labels = number_format(accuracy = 0.1)) +
  scale_y_continuous(name = "Mean Mates per Female", labels = number_format(accuracy = 0.1)) +
  labs(
    title = "Relationship: RSC vs Mean Mates",
    subtitle = "Red line: linear fit | Black dotted line: y=x (theoretical perfect match)",
    caption = "Each point = one generation from one replicate"
  )

ggsave(file.path(plot_dir, "01_rsc_vs_mates_scatter.png"), p1, width = 10, height = 7, dpi = 300)
cat("Saved: 01_rsc_vs_mates_scatter.png\n")

# =============================================================================
# PLOT 2: RSC vs Mates by Replicate
# =============================================================================
p2 <- df %>%
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
  scale_y_continuous(labels = number_format(accuracy = 0.1)) +
  scale_color_manual(values = c("Mean RSC (Lambda)" = "darkgreen", "Mean Mates" = "darkblue")) +
  scale_linetype_manual(values = c("Mean RSC (Lambda)" = "solid", "Mean Mates" = "dashed")) +
  labs(
    title = "RSC vs Mean Mates Evolution by Replicate",
    subtitle = "Solid line: RSC | Dashed line: Mean Mates",
    x = "Generation",
    y = "Value",
    color = "Metric",
    linetype = "Metric"
  )

ggsave(file.path(plot_dir, "02_rsc_vs_mates_by_replicate.png"), p2, width = 14, height = 10, dpi = 300)
cat("Saved: 02_rsc_vs_mates_by_replicate.png\n")

# =============================================================================
# PLOT 3: RSC vs Mates Summary (Mean ± SE)
# =============================================================================
summary_data <- df %>%
  group_by(Generation) %>%
  summarise(
    MeanRSC = mean(MeanRSC, na.rm=TRUE),
    SE_RSC = sd(MeanRSC, na.rm=TRUE) / sqrt(n()),
    MeanMates = mean(MeanMates, na.rm=TRUE),
    SE_Mates = sd(MeanMates, na.rm=TRUE) / sqrt(n()),
    .groups = 'drop'
  )

p3 <- ggplot(summary_data, aes(x = Generation)) +
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
  scale_y_continuous(labels = number_format(accuracy = 0.01)) +
  scale_color_manual(values = c("Mean RSC (Lambda)" = "darkgreen", "Mean Mates" = "darkblue")) +
  labs(
    title = "RSC vs Mean Mates: Summary Across All Replicates",
    subtitle = "Mean ± SE across replicates",
    x = "Generation",
    y = "Value",
    color = "Metric"
  )

ggsave(file.path(plot_dir, "03_rsc_vs_mates_summary.png"), p3, width = 10, height = 6, dpi = 300)
cat("Saved: 03_rsc_vs_mates_summary.png\n")

# =============================================================================
# PLOT 4: Correlation between RSC and Mates over Time
# =============================================================================
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
  scale_y_continuous(limits = c(0, 1.05), labels = number_format(accuracy = 0.01)) +
  labs(
    title = "Correlation between RSC and Mates over Time",
    subtitle = "Red dashed line at y=1 (perfect correlation)",
    x = "Generation",
    y = "Correlation (RSC, Mates)"
  )

ggsave(file.path(plot_dir, "04_rsc_vs_mates_correlation.png"), p4, width = 10, height = 6, dpi = 300)
cat("Saved: 04_rsc_vs_mates_correlation.png\n")

# =============================================================================
# PLOT 5: Difference (RSC - Mates) over Time
# =============================================================================
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
  scale_y_continuous(labels = number_format(accuracy = 0.01)) +
  labs(
    title = "Difference: RSC - Mean Mates over Time",
    subtitle = "Mean difference across replicates (should be close to 0)",
    x = "Generation",
    y = "Mean Difference (RSC - Mates)"
  )

ggsave(file.path(plot_dir, "05_rsc_vs_mates_difference.png"), p5, width = 10, height = 6, dpi = 300)
cat("Saved: 05_rsc_vs_mates_difference.png\n")

# =============================================================================
# PLOT 6: Evolution of Mean Male and Female Traits
# =============================================================================
p6 <- df %>%
  pivot_longer(cols = c(MeanMale, MeanFemale), names_to = "Trait", values_to = "Value") %>%
  ggplot(aes(x = Generation, y = Value, color = Trait, group = interaction(Trait, Rep))) +
  geom_line(alpha = 0.6) +
  stat_smooth(aes(group = Trait), method = "loess", se = TRUE, linewidth = 1.2) +
  labs(
    title = "Evolution of Mean Male and Female Traits",
    subtitle = "Lines show individual replicates, smoothed lines show trends",
    x = "Generation",
    y = "Mean Trait Value",
    color = "Trait"
  ) +
  scale_color_manual(values = c("MeanMale" = "blue", "MeanFemale" = "red")) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(file.path(plot_dir, "06_traits_evolution.png"), p6, width = 12, height = 8, dpi = 300)
cat("Saved: 06_traits_evolution.png\n")

# =============================================================================
# PLOT 7: Evolution of Mean RSC
# =============================================================================
p7 <- df %>%
  ggplot(aes(x = Generation, y = MeanRSC, group = Rep)) +
  geom_line(alpha = 0.6, color = "purple") +
  stat_smooth(method = "loess", se = TRUE, color = "darkviolet", linewidth = 1.2) +
  labs(
    title = "Evolution of Mean RSC (Risk of Sperm Competition)",
    subtitle = "RSC phenotype evolution across generations",
    x = "Generation",
    y = "Mean RSC Phenotype"
  ) +
  theme_minimal()

ggsave(file.path(plot_dir, "07_rsc_evolution.png"), p7, width = 12, height = 8, dpi = 300)
cat("Saved: 07_rsc_evolution.png\n")

# =============================================================================
# PLOT 8: Evolution of Mean Mates per Female
# =============================================================================
p8 <- df %>%
  ggplot(aes(x = Generation, y = MeanMates, group = Rep)) +
  geom_line(alpha = 0.6, color = "darkgreen") +
  stat_smooth(method = "loess", se = TRUE, color = "forestgreen", linewidth = 1.2) +
  labs(
    title = "Evolution of Mean Mates per Female",
    subtitle = "Average number of mates per female across generations",
    x = "Generation",
    y = "Mean Mates per Female"
  ) +
  theme_minimal()

ggsave(file.path(plot_dir, "08_mates_evolution.png"), p8, width = 12, height = 8, dpi = 300)
cat("Saved: 08_mates_evolution.png\n")

# =============================================================================
# PLOT 9: Trait Correlation Evolution
# =============================================================================
p9 <- df %>%
  ggplot(aes(x = Generation, y = cor, group = Rep)) +
  geom_line(alpha = 0.6, color = "orange") +
  stat_smooth(method = "loess", se = TRUE, color = "darkorange", linewidth = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  labs(
    title = "Evolution of Male-Female Trait Correlation",
    subtitle = "Correlation between male and female traits",
    x = "Generation",
    y = "Correlation (r)"
  ) +
  theme_minimal()

ggsave(file.path(plot_dir, "09_correlation_evolution.png"), p9, width = 12, height = 8, dpi = 300)
cat("Saved: 09_correlation_evolution.png\n")

# =============================================================================
# PLOT 10: Sperm Count Evolution
# =============================================================================
p10 <- df %>%
  ggplot(aes(x = Generation, y = MeanCount, group = Rep)) +
  geom_line(alpha = 0.6, color = "brown") +
  stat_smooth(method = "loess", se = TRUE, color = "saddlebrown", linewidth = 1.2) +
  labs(
    title = "Evolution of Mean Sperm Count",
    subtitle = "Average sperm number across generations",
    x = "Generation",
    y = "Mean Sperm Count"
  ) +
  theme_minimal()

ggsave(file.path(plot_dir, "10_sperm_evolution.png"), p10, width = 12, height = 8, dpi = 300)
cat("Saved: 10_sperm_evolution.png\n")

# =============================================================================
# PLOT 11: Multi-panel Summary
# =============================================================================
p11 <- df %>%
  select(Generation, Rep, MeanMale, MeanFemale, MeanRSC, MeanMates, cor, MeanCount) %>%
  pivot_longer(cols = c(MeanMale, MeanFemale, MeanRSC, MeanMates, cor, MeanCount), names_to = "Metric", values_to = "Value") %>%
  ggplot(aes(x = Generation, y = Value)) +
  geom_line(alpha = 0.3, aes(group = interaction(Metric, Rep))) +
  stat_smooth(method = "loess", se = TRUE, aes(color = Metric), linewidth = 1.2) +
  facet_wrap(~ Metric, scales = "free_y", ncol = 2) +
  labs(
    title = "Summary of All Traits and Metrics",
    x = "Generation",
    y = "Value"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

ggsave(file.path(plot_dir, "11_summary_all_metrics.png"), p11, width = 14, height = 10, dpi = 300)
cat("Saved: 11_summary_all_metrics.png\n")

# =============================================================================
# PLOT 12: Trait Distributions (Final Generation)
# =============================================================================
final_gen <- max(df$Generation)
p12 <- df %>%
  filter(Generation == final_gen) %>%
  select(MeanMale, MeanFemale, MeanRSC, MeanMates, MeanCount) %>%
  pivot_longer(everything(), names_to = "Trait", values_to = "Value") %>%
  ggplot(aes(x = Value, fill = Trait)) +
  geom_histogram(bins = 30, alpha = 0.7) +
  facet_wrap(~ Trait, scales = "free", ncol = 2) +
  labs(
    title = paste("Distribution of Traits at Final Generation", final_gen),
    x = "Trait Value",
    y = "Frequency"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

ggsave(file.path(plot_dir, "12_final_generation_distributions.png"), p12, width = 14, height = 10, dpi = 300)
cat("Saved: 12_final_generation_distributions.png\n")

# =============================================================================
# PLOT 13: Mates Evolution by Replicate (Overlay)
# =============================================================================
p13 <- ggplot(df, aes(x = Generation, y = MeanMates, color = Rep)) +
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
    labels = number_format(accuracy = 0.01),
    breaks = pretty_breaks(n = 8)
  ) +
  labs(
    title = "Evolution of Mean Mates per Female Across Generations",
    subtitle = "Each line represents one replicate",
    x = "Generation",
    y = "Mean Mates per Female",
    color = "Replicate"
  )

ggsave(file.path(plot_dir, "13_mates_evolution_overlay.png"), p13, width = 10, height = 6, dpi = 300)
cat("Saved: 13_mates_evolution_overlay.png\n")

# =============================================================================
# PLOT 14: Mates Evolution by Replicate (Faceted)
# =============================================================================
p14 <- ggplot(df, aes(x = Generation, y = MeanMates, color = Rep)) +
  geom_line(linewidth = 1.2, alpha = 0.9) +
  facet_wrap(~Rep, ncol = 3, scales = "free_y") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    strip.text = element_text(face = "bold", size = 11),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 12),
    legend.position = "none",
    panel.grid.minor = element_line(color = "grey95", linewidth = 0.3)
  ) +
  scale_y_continuous(
    name = "Mean Mates per Female",
    labels = number_format(accuracy = 0.05),
    breaks = pretty_breaks(n = 6)
  ) +
  labs(
    title = "Evolution of Mean Mates per Female by Replicate",
    subtitle = "Each panel shows one replicate",
    x = "Generation",
    y = "Mean Mates per Female"
  )

ggsave(file.path(plot_dir, "14_mates_evolution_by_rep.png"), p14, width = 14, height = 12, dpi = 300)
cat("Saved: 14_mates_evolution_by_rep.png\n")

# =============================================================================
# Statistical Analysis
# =============================================================================
cat("\n=== Statistical Analysis: RSC vs Mates ===\n")
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
final_data <- df %>% filter(Generation == final_gen)
cat("Final generation:", final_gen, "\n")
cat("Mean RSC (final):", round(mean(final_data$MeanRSC, na.rm=TRUE), 4), "\n")
cat("Mean Mates (final):", round(mean(final_data$MeanMates, na.rm=TRUE), 4), "\n")
cat("Mean Difference (RSC - Mates):", round(mean(final_data$MeanRSC - final_data$MeanMates, na.rm=TRUE), 4), "\n")
cat("Correlation (final gen):", round(cor(final_data$MeanRSC, final_data$MeanMates, use = "complete.obs"), 4), "\n")

cat("\n=== Summary Statistics Table ===\n")
summary_table <- df %>%
  group_by(Rep) %>%
  summarize(
    mean_mates = round(mean(MeanMates), 3),
    final_mates = round(last(MeanMates), 3),
    mates_range = paste(round(min(MeanMates), 3), "to", round(max(MeanMates), 3)),
    mean_rsc = round(mean(MeanRSC), 3),
    final_rsc = round(last(MeanRSC), 3),
    rsc_range = paste(round(min(MeanRSC), 3), "to", round(max(MeanRSC), 3)),
    .groups = 'drop'
  )
print(summary_table)

cat("\n=== All plots saved to:", plot_dir, "===\n")
cat("Total plots created: 14\n")
cat("\nPlots created:\n")
cat("1. RSC vs Mates Scatter\n")
cat("2. RSC vs Mates by Replicate\n")
cat("3. RSC vs Mates Summary\n")
cat("4. RSC vs Mates Correlation\n")
cat("5. RSC vs Mates Difference\n")
cat("6. Traits Evolution\n")
cat("7. RSC Evolution\n")
cat("8. Mates Evolution\n")
cat("9. Correlation Evolution\n")
cat("10. Sperm Evolution\n")
cat("11. Summary All Metrics\n")
cat("12. Final Generation Distributions\n")
cat("13. Mates Evolution Overlay\n")
cat("14. Mates Evolution by Replicate\n")

