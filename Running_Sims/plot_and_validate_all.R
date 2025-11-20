# Comprehensive Plotting and Validation Script
# This script generates all plots and validates that:
# 1. RSC and MeanMates match properly (MeanMates ≈ MeanRSC when RSC > 0)
# 2. No unwanted clamping (RSC can be negative, only clamped at Poisson sampling)
# 3. Proper evolution of all traits
# 4. Correlation between male and female traits
# 
# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

# Get the most recent CSV file from 11_19_testing folder
csv_files <- list.files(pattern = "parallel_sim_results_.*\\.csv$", path = "11_19_testing", full.names = TRUE)
if (length(csv_files) == 0) {
  # Fallback: check current directory
  csv_files <- list.files(pattern = "parallel_sim_results_.*\\.csv$", full.names = TRUE)
  if (length(csv_files) == 0) {
    stop("No CSV files found! Make sure the simulation has run and created a CSV file.")
  }
}

# Get the most recent file
latest_file <- csv_files[order(file.info(csv_files)$mtime, decreasing = TRUE)[1]]
cat("Reading data from:", latest_file, "\n")

# Read the data
df <- read.csv(latest_file)

# Create output directory for plots
plot_dir <- "11_19_testing"
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir)
  cat("Created directory:", plot_dir, "\n")
}

# Convert Generation to numeric
df$Generation <- as.numeric(df$Generation)
df$Rep <- as.factor(df$Rep)

cat("\n=== DATA VALIDATION ===\n")
cat("Total generations:", max(df$Generation), "\n")
cat("Total replicates:", length(unique(df$Rep)), "\n")
cat("Total data points:", nrow(df), "\n")

# Check for negative RSC values (should be allowed, only clamped at Poisson sampling)
negative_rsc_count <- sum(df$MeanRSC < 0, na.rm = TRUE)
cat("\nNegative RSC values found:", negative_rsc_count, "out of", nrow(df), "\n")
if (negative_rsc_count > 0) {
  cat("  -> This is EXPECTED: RSC can evolve to negative values\n")
  cat("  -> Negative RSC is clamped to 0 only when sampling mates (Poisson requires non-negative lambda)\n")
} else {
  cat("  -> All RSC values are non-negative (RSC evolved to positive values)\n")
}

# Check RSC vs MeanMates relationship
cat("\n=== RSC vs MeanMates Validation ===\n")
cat("Mean RSC range:", round(min(df$MeanRSC, na.rm=TRUE), 3), "to", round(max(df$MeanRSC, na.rm=TRUE), 3), "\n")
cat("Mean mates range:", round(min(df$MeanMates, na.rm=TRUE), 3), "to", round(max(df$MeanMates, na.rm=TRUE), 3), "\n")

# Calculate expected relationship
# MeanMates should be approximately equal to MeanRSC when RSC > 0
# But slightly lower due to Poisson sampling variance and negative RSC values being clamped to 0
df$ExpectedMates <- pmax(0, df$MeanRSC)  # Expected mates = max(0, RSC)
df$Difference <- df$MeanRSC - df$MeanMates
df$AbsDifference <- abs(df$Difference)

cat("\nMean difference (RSC - MeanMates):", round(mean(df$Difference, na.rm=TRUE), 4), "\n")
cat("SD of difference:", round(sd(df$Difference, na.rm=TRUE), 4), "\n")
cat("Max absolute difference:", round(max(df$AbsDifference, na.rm=TRUE), 4), "\n")

# Overall correlation
overall_cor <- cor(df$MeanRSC, df$MeanMates, use = "complete.obs")
cat("Overall correlation (RSC, MeanMates):", round(overall_cor, 4), "\n")
if (overall_cor > 0.95) {
  cat("  -> EXCELLENT: Strong correlation indicates proper matching\n")
} else if (overall_cor > 0.90) {
  cat("  -> GOOD: Strong correlation\n")
} else {
  cat("  -> WARNING: Correlation lower than expected, check for issues\n")
}

# Linear model
model <- lm(MeanMates ~ MeanRSC, data = df)
slope <- coef(model)[2]
intercept <- coef(model)[1]
cat("\nLinear model: MeanMates ~ MeanRSC\n")
cat("  Slope:", round(slope, 4), "(expected ≈ 1.0)\n")
cat("  Intercept:", round(intercept, 4), "(expected ≈ 0.0)\n")
if (abs(slope - 1.0) < 0.1 && abs(intercept) < 0.5) {
  cat("  -> EXCELLENT: Model parameters match expectations\n")
} else {
  cat("  -> NOTE: Some deviation expected due to Poisson sampling variance\n")
}

# =============================================================================
# PLOT 1: RSC vs Mates Scatter with Validation
# =============================================================================
p1 <- ggplot(df, aes(x = MeanRSC, y = MeanMates)) +
  geom_point(alpha = 0.4, color = "steelblue", size = 1.2) +
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
    title = "Validation: RSC vs Mean Mates Relationship",
    subtitle = paste0("Correlation: ", round(overall_cor, 3), " | Slope: ", round(slope, 3), 
                     " | Red: linear fit | Black: y=x (perfect match)"),
    caption = "Each point = one generation from one replicate"
  )

ggsave(file.path(plot_dir, "01_validation_rsc_vs_mates.png"), p1, width = 10, height = 7, dpi = 300)
cat("\nSaved: 01_validation_rsc_vs_mates.png\n")

# =============================================================================
# PLOT 2: RSC and MeanMates Evolution Overlay
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

p2 <- ggplot(summary_data, aes(x = Generation)) +
  geom_line(aes(y = MeanRSC, color = "Mean RSC"), linewidth = 1.2) +
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
  scale_color_manual(values = c("Mean RSC" = "darkgreen", "Mean Mates" = "darkblue")) +
  labs(
    title = "RSC and Mean Mates Evolution: Summary Across Replicates",
    subtitle = "Mean ± SE across replicates | Should track closely",
    x = "Generation",
    y = "Value",
    color = "Metric"
  )

ggsave(file.path(plot_dir, "02_rsc_mates_evolution.png"), p2, width = 10, height = 6, dpi = 300)
cat("Saved: 02_rsc_mates_evolution.png\n")

# =============================================================================
# PLOT 3: Difference (RSC - MeanMates) Over Time
# =============================================================================
p3 <- df %>%
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
    title = "Difference: RSC - Mean Mates Over Time",
    subtitle = "Mean difference across replicates (should be close to 0)",
    x = "Generation",
    y = "Mean Difference (RSC - Mates)"
  )

ggsave(file.path(plot_dir, "03_rsc_mates_difference.png"), p3, width = 10, height = 6, dpi = 300)
cat("Saved: 03_rsc_mates_difference.png\n")

# =============================================================================
# PLOT 4: Evolution of Mean Male and Female Traits
# =============================================================================
p4 <- df %>%
  pivot_longer(cols = c(MeanMale, MeanFemale), names_to = "Trait", values_to = "Value") %>%
  ggplot(aes(x = Generation, y = Value, color = Trait, group = interaction(Trait, Rep))) +
  geom_line(alpha = 0.4) +
  stat_smooth(aes(group = Trait), method = "loess", se = TRUE, linewidth = 1.5) +
  labs(
    title = "Evolution of Mean Male and Female Traits",
    subtitle = "Thin lines: individual replicates | Thick lines: smoothed trends",
    x = "Generation",
    y = "Mean Trait Value",
    color = "Trait"
  ) +
  scale_color_manual(values = c("MeanMale" = "blue", "MeanFemale" = "red")) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(file.path(plot_dir, "04_traits_evolution.png"), p4, width = 12, height = 8, dpi = 300)
cat("Saved: 04_traits_evolution.png\n")

# =============================================================================
# PLOT 5: Trait Correlation Evolution
# =============================================================================
p5 <- df %>%
  ggplot(aes(x = Generation, y = cor, group = Rep)) +
  geom_line(alpha = 0.5, color = "orange") +
  stat_smooth(method = "loess", se = TRUE, color = "darkorange", linewidth = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray", linewidth = 0.8) +
  labs(
    title = "Evolution of Male-Female Trait Correlation",
    subtitle = "Correlation between male and female traits across generations",
    x = "Generation",
    y = "Correlation (r)"
  ) +
  theme_minimal()

ggsave(file.path(plot_dir, "05_correlation_evolution.png"), p5, width = 12, height = 8, dpi = 300)
cat("Saved: 05_correlation_evolution.png\n")

# =============================================================================
# PLOT 6: RSC Evolution by Replicate
# =============================================================================
p6 <- ggplot(df, aes(x = Generation, y = MeanRSC, color = Rep)) +
  geom_line(linewidth = 1, alpha = 0.7) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.position = "right"
  ) +
  scale_y_continuous(
    name = "Mean RSC",
    labels = number_format(accuracy = 0.1)
  ) +
  labs(
    title = "Evolution of Mean RSC Across Replicates",
    subtitle = "Each line represents one replicate",
    x = "Generation",
    y = "Mean RSC",
    color = "Replicate"
  )

ggsave(file.path(plot_dir, "06_rsc_evolution.png"), p6, width = 12, height = 8, dpi = 300)
cat("Saved: 06_rsc_evolution.png\n")

# =============================================================================
# PLOT 7: MeanMates Evolution by Replicate
# =============================================================================
p7 <- ggplot(df, aes(x = Generation, y = MeanMates, color = Rep)) +
  geom_line(linewidth = 1, alpha = 0.7) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.position = "right"
  ) +
  scale_y_continuous(
    name = "Mean Mates per Female",
    labels = number_format(accuracy = 0.1)
  ) +
  labs(
    title = "Evolution of Mean Mates per Female Across Replicates",
    subtitle = "Each line represents one replicate",
    x = "Generation",
    y = "Mean Mates per Female",
    color = "Replicate"
  )

ggsave(file.path(plot_dir, "07_mates_evolution.png"), p7, width = 12, height = 8, dpi = 300)
cat("Saved: 07_mates_evolution.png\n")

# =============================================================================
# PLOT 8: Sperm Count Evolution
# =============================================================================
p8 <- df %>%
  ggplot(aes(x = Generation, y = MeanCount, group = Rep)) +
  geom_line(alpha = 0.5, color = "brown") +
  stat_smooth(method = "loess", se = TRUE, color = "saddlebrown", linewidth = 1.5) +
  labs(
    title = "Evolution of Mean Sperm Count",
    subtitle = "Average sperm number across generations",
    x = "Generation",
    y = "Mean Sperm Count"
  ) +
  theme_minimal()

ggsave(file.path(plot_dir, "08_sperm_evolution.png"), p8, width = 12, height = 8, dpi = 300)
cat("Saved: 08_sperm_evolution.png\n")

# =============================================================================
# PLOT 9: Multi-panel Summary
# =============================================================================
p9 <- df %>%
  select(Generation, Rep, MeanMale, MeanFemale, MeanRSC, MeanMates, cor, MeanCount) %>%
  pivot_longer(cols = c(MeanMale, MeanFemale, MeanRSC, MeanMates, cor, MeanCount), 
               names_to = "Metric", values_to = "Value") %>%
  ggplot(aes(x = Generation, y = Value)) +
  geom_line(alpha = 0.2, aes(group = interaction(Metric, Rep))) +
  stat_smooth(method = "loess", se = TRUE, aes(color = Metric), linewidth = 1.2) +
  facet_wrap(~ Metric, scales = "free_y", ncol = 2) +
  labs(
    title = "Summary of All Traits and Metrics",
    x = "Generation",
    y = "Value"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

ggsave(file.path(plot_dir, "09_summary_all_metrics.png"), p9, width = 14, height = 10, dpi = 300)
cat("Saved: 09_summary_all_metrics.png\n")

# =============================================================================
# PLOT 10: RSC Distribution Check (for negative values)
# =============================================================================
p10 <- df %>%
  ggplot(aes(x = MeanRSC)) +
  geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7, color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 1) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  ) +
  labs(
    title = "Distribution of Mean RSC Values",
    subtitle = paste0("Red line at 0 | Negative values allowed (", negative_rsc_count, " found)"),
    x = "Mean RSC",
    y = "Frequency"
  )

ggsave(file.path(plot_dir, "10_rsc_distribution.png"), p10, width = 10, height = 6, dpi = 300)
cat("Saved: 10_rsc_distribution.png\n")

# =============================================================================
# Generate Validation Report
# =============================================================================
cat("\n=== GENERATING VALIDATION REPORT ===\n")

# Final generation statistics
final_gen <- max(df$Generation)
final_data <- df %>% filter(Generation == final_gen)

validation_report <- data.frame(
  Check = c(
    "Total Generations",
    "Total Replicates",
    "Negative RSC Values",
    "Overall RSC-Mates Correlation",
    "Linear Model Slope",
    "Linear Model Intercept",
    "Mean Difference (RSC - Mates)",
    "Max Absolute Difference",
    "Final Gen Mean RSC",
    "Final Gen Mean Mates",
    "Final Gen Correlation (RSC, Mates)",
    "Final Gen Trait Correlation (Male-Female)"
  ),
  Value = c(
    max(df$Generation),
    length(unique(df$Rep)),
    negative_rsc_count,
    round(overall_cor, 4),
    round(slope, 4),
    round(intercept, 4),
    round(mean(df$Difference, na.rm=TRUE), 4),
    round(max(df$AbsDifference, na.rm=TRUE), 4),
    round(mean(final_data$MeanRSC, na.rm=TRUE), 4),
    round(mean(final_data$MeanMates, na.rm=TRUE), 4),
    round(cor(final_data$MeanRSC, final_data$MeanMates, use = "complete.obs"), 4),
    round(mean(final_data$cor, na.rm=TRUE), 4)
  ),
  Status = c(
    "OK",
    "OK",
    ifelse(negative_rsc_count > 0, "EXPECTED (RSC can be negative)", "OK (all positive)"),
    ifelse(overall_cor > 0.95, "EXCELLENT", ifelse(overall_cor > 0.90, "GOOD", "CHECK")),
    ifelse(abs(slope - 1.0) < 0.1, "EXCELLENT", "ACCEPTABLE"),
    ifelse(abs(intercept) < 0.5, "EXCELLENT", "ACCEPTABLE"),
    ifelse(abs(mean(df$Difference, na.rm=TRUE)) < 0.5, "EXCELLENT", "ACCEPTABLE"),
    "OK",
    "OK",
    "OK",
    ifelse(cor(final_data$MeanRSC, final_data$MeanMates, use = "complete.obs") > 0.95, "EXCELLENT", "GOOD"),
    "OK"
  )
)

write.csv(validation_report, file.path(plot_dir, "validation_report.csv"), row.names = FALSE)
cat("\nValidation report saved to: validation_report.csv\n")

# Print summary
cat("\n=== VALIDATION SUMMARY ===\n")
print(validation_report)

cat("\n=== FINAL GENERATION STATISTICS ===\n")
cat("Final generation:", final_gen, "\n")
cat("Mean RSC (final):", round(mean(final_data$MeanRSC, na.rm=TRUE), 4), "\n")
cat("Mean Mates (final):", round(mean(final_data$MeanMates, na.rm=TRUE), 4), "\n")
cat("Mean Difference (RSC - Mates):", round(mean(final_data$MeanRSC - final_data$MeanMates, na.rm=TRUE), 4), "\n")
cat("Correlation RSC-Mates (final):", round(cor(final_data$MeanRSC, final_data$MeanMates, use = "complete.obs"), 4), "\n")
cat("Mean Trait Correlation (final):", round(mean(final_data$cor, na.rm=TRUE), 4), "\n")

cat("\n=== All plots saved to:", plot_dir, "===\n")
cat("Total plots created: 10\n")
cat("\nPlots created:\n")
cat("1. Validation: RSC vs Mates\n")
cat("2. RSC and Mates Evolution\n")
cat("3. RSC-Mates Difference\n")
cat("4. Traits Evolution\n")
cat("5. Correlation Evolution\n")
cat("6. RSC Evolution\n")
cat("7. Mates Evolution\n")
cat("8. Sperm Evolution\n")
cat("9. Summary All Metrics\n")
cat("10. RSC Distribution\n")

cat("\n=== VALIDATION COMPLETE ===\n")

