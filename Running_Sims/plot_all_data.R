# R script to plot all simulation data
# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)

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
plot_dir <- "Plots_AllData"
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir)
}

# Convert Generation to numeric if needed
df$Generation <- as.numeric(df$Generation)

# Plot 1: Evolution of Mean Male and Female Traits
p1 <- df %>%
  pivot_longer(cols = c(MeanMale, MeanFemale), names_to = "Trait", values_to = "Value") %>%
  ggplot(aes(x = Generation, y = Value, color = Trait, group = interaction(Trait, Rep))) +
  geom_line(alpha = 0.6) +
  stat_smooth(aes(group = Trait), method = "loess", se = TRUE, size = 1.2) +
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

ggsave(file.path(plot_dir, "01_traits_evolution.png"), p1, width = 12, height = 8, dpi = 300)
print("Saved: 01_traits_evolution.png")

# Plot 2: Evolution of Mean RSC
p2 <- df %>%
  ggplot(aes(x = Generation, y = MeanRSC, group = Rep)) +
  geom_line(alpha = 0.6, color = "purple") +
  stat_smooth(method = "loess", se = TRUE, color = "darkviolet", size = 1.2) +
  labs(
    title = "Evolution of Mean RSC (Risk of Sperm Competition)",
    subtitle = "RSC phenotype evolution across generations",
    x = "Generation",
    y = "Mean RSC Phenotype"
  ) +
  theme_minimal()

ggsave(file.path(plot_dir, "02_rsc_evolution.png"), p2, width = 12, height = 8, dpi = 300)
print("Saved: 02_rsc_evolution.png")

# Plot 3: Evolution of Mean Mates per Female
p3 <- df %>%
  ggplot(aes(x = Generation, y = MeanMates, group = Rep)) +
  geom_line(alpha = 0.6, color = "darkgreen") +
  stat_smooth(method = "loess", se = TRUE, color = "forestgreen", size = 1.2) +
  labs(
    title = "Evolution of Mean Mates per Female",
    subtitle = "Average number of mates per female across generations",
    x = "Generation",
    y = "Mean Mates per Female"
  ) +
  theme_minimal()

ggsave(file.path(plot_dir, "03_mates_evolution.png"), p3, width = 12, height = 8, dpi = 300)
print("Saved: 03_mates_evolution.png")

# Plot 4: Trait Correlation Evolution
p4 <- df %>%
  ggplot(aes(x = Generation, y = cor, group = Rep)) +
  geom_line(alpha = 0.6, color = "orange") +
  stat_smooth(method = "loess", se = TRUE, color = "darkorange", size = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  labs(
    title = "Evolution of Male-Female Trait Correlation",
    subtitle = "Correlation between male and female traits",
    x = "Generation",
    y = "Correlation (r)"
  ) +
  theme_minimal()

ggsave(file.path(plot_dir, "04_correlation_evolution.png"), p4, width = 12, height = 8, dpi = 300)
print("Saved: 04_correlation_evolution.png")

# Plot 5: RSC vs Mates Relationship
p5 <- df %>%
  ggplot(aes(x = MeanRSC, y = MeanMates)) +
  geom_point(alpha = 0.5, color = "steelblue") +
  geom_smooth(method = "lm", se = TRUE, color = "darkblue") +
  labs(
    title = "Relationship between RSC and Mean Mates",
    subtitle = "How RSC phenotype relates to actual mating rate",
    x = "Mean RSC Phenotype",
    y = "Mean Mates per Female"
  ) +
  theme_minimal()

ggsave(file.path(plot_dir, "05_rsc_vs_mates.png"), p5, width = 12, height = 8, dpi = 300)
print("Saved: 05_rsc_vs_mates.png")

# Plot 6: Sperm Count Evolution
p6 <- df %>%
  ggplot(aes(x = Generation, y = MeanCount, group = Rep)) +
  geom_line(alpha = 0.6, color = "brown") +
  stat_smooth(method = "loess", se = TRUE, color = "saddlebrown", size = 1.2) +
  labs(
    title = "Evolution of Mean Sperm Count",
    subtitle = "Average sperm number across generations",
    x = "Generation",
    y = "Mean Sperm Count"
  ) +
  theme_minimal()

ggsave(file.path(plot_dir, "06_sperm_evolution.png"), p6, width = 12, height = 8, dpi = 300)
print("Saved: 06_sperm_evolution.png")

# Plot 7: Multi-panel Summary
p7 <- df %>%
  select(Generation, MeanMale, MeanFemale, MeanRSC, MeanMates, cor, MeanCount) %>%
  pivot_longer(cols = -Generation, names_to = "Metric", values_to = "Value") %>%
  ggplot(aes(x = Generation, y = Value)) +
  geom_line(alpha = 0.3, aes(group = interaction(Metric, Rep))) +
  stat_smooth(method = "loess", se = TRUE, aes(color = Metric)) +
  facet_wrap(~ Metric, scales = "free_y", ncol = 2) +
  labs(
    title = "Summary of All Traits and Metrics",
    x = "Generation",
    y = "Value"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

ggsave(file.path(plot_dir, "07_summary_all_metrics.png"), p7, width = 14, height = 10, dpi = 300)
print("Saved: 07_summary_all_metrics.png")

# Plot 8: Trait Distributions (Final Generation)
final_gen <- max(df$Generation)
p8 <- df %>%
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

ggsave(file.path(plot_dir, "08_final_generation_distributions.png"), p8, width = 14, height = 10, dpi = 300)
print("Saved: 08_final_generation_distributions.png")

cat("\n=== All plots saved to:", plot_dir, "===\n")
cat("Total plots created: 8\n")

