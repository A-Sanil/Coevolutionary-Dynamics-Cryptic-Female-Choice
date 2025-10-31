# R script for analyzing and plotting simulation results

# --- 1. Setup ---
# Install and load necessary packages
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("here")) install.packages("here")

library(tidyverse)
library(here)

# --- 2. Load Data ---
# Find the most recent results file, prioritizing parallel runs
results_dir <- here("Running_Sims")
parallel_files <- list.files(path = results_dir, pattern = "parallel_sim_results_.*\\.csv$", full.names = TRUE)
quick_files <- list.files(path = results_dir, pattern = "quick_sim_results_.*\\.csv$", full.names = TRUE)

if(length(parallel_files) > 0) {
  file_paths <- parallel_files
  cat("Found parallel simulation files. Using the latest one.\n")
} else if(length(quick_files) > 0) {
  file_paths <- quick_files
  cat("No parallel files found. Using latest quick simulation file.\n")
} else {
  stop("No simulation results CSV found in ", results_dir)
}

latest_csv <- file_paths[which.max(file.info(file_paths)$mtime)]

cat("Loading data from:", latest_csv, "\n")
sim_data <- read_csv(latest_csv, show_col_types = FALSE)

# Convert Generation and Rep to factors/integers for plotting
sim_data <- sim_data %>% mutate(Generation = as.integer(Generation), Rep = as.factor(Rep))

# --- 3. Data Aggregation ---
# Calculate the mean and standard error across replicates for each generation
summary_data <- sim_data %>%
  group_by(Generation) %>%
  summarise(
    MeanMale = mean(MeanMale),
    SE_Male = sd(MeanMale) / sqrt(n()),
    MeanFemale = mean(MeanFemale),
    SE_Female = sd(MeanFemale) / sqrt(n()),
    MeanRSC = mean(MeanRSC),
    SE_RSC = sd(MeanRSC) / sqrt(n()),
    MeanCor = mean(cor, na.rm = TRUE),
    SE_Cor = sd(cor, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  )

# --- 4. Generate Plots ---

# Plot 1: Evolution of Mean Male and Female Traits
plot_traits <- ggplot(summary_data, aes(x = Generation)) +
  geom_line(aes(y = MeanMale, color = "Male Trait")) +
  geom_ribbon(aes(ymin = MeanMale - SE_Male, ymax = MeanMale + SE_Male), fill = "blue", alpha = 0.2) +
  geom_line(aes(y = MeanFemale, color = "Female Trait")) +
  geom_ribbon(aes(ymin = MeanFemale - SE_Female, ymax = MeanFemale + SE_Female), fill = "red", alpha = 0.2) +
  labs(
    title = "Evolution of Mean Male and Female Traits",
    x = "Generation",
    y = "Mean Trait Value",
    color = "Trait"
  ) +
  scale_color_manual(values = c("Male Trait" = "blue", "Female Trait" = "red")) +
  theme_minimal()

print(plot_traits)

# Plot 2: Evolution of Mean RSC (Risk of Sperm Competition)
plot_rsc <- ggplot(summary_data, aes(x = Generation)) +
  geom_line(aes(y = MeanRSC), color = "darkgreen") +
  geom_ribbon(aes(ymin = MeanRSC - SE_RSC, ymax = MeanRSC + SE_RSC), fill = "darkgreen", alpha = 0.2) +
  labs(
    title = "Evolution of Mean RSC (Risk of Sperm Competition)",
    x = "Generation",
    y = "Mean RSC Phenotype",
    subtitle = "Average RSC phenotype across all replicates, with standard error."
  ) +
  theme_minimal()

print(plot_rsc)

# Plot 3: Evolution of Male-Female Trait Correlation
plot_correlation <- ggplot(summary_data, aes(x = Generation)) +
  geom_line(aes(y = MeanCor), color = "purple") +
  geom_ribbon(aes(ymin = MeanCor - SE_Cor, ymax = MeanCor + SE_Cor), fill = "purple", alpha = 0.2) +
  labs(
    title = "Evolution of Male-Female Trait Correlation",
    x = "Generation",
    y = "Mean Correlation (cor)"
  ) +
  theme_minimal()

print(plot_correlation)

# Plot 4: Correlation vs. Mean RSC
# We'll use the raw data from the final generation for a clearer picture
final_gen_data <- sim_data %>% filter(Generation == max(Generation))

plot_cor_vs_rsc <- ggplot(final_gen_data, aes(x = MeanRSC, y = cor)) +
  geom_point(alpha = 0.8, color = "sienna") +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
  labs(
    title = "Male-Female Correlation vs. RSC in Final Generation",
    subtitle = "Each point represents one replicate",
    x = "Mean RSC Phenotype",
    y = "Male-Female Trait Correlation (cor)"
  ) +
  theme_minimal()

print(plot_cor_vs_rsc)

# Plot 5: Mean Sperm Count (MeanCount) across generations, by replicate
plot_sperm_count <- ggplot(sim_data, aes(x = Generation, y = MeanCount, color = Rep)) +
  geom_line(alpha = 0.8) +
  theme_minimal() +
  labs(
    title = 'Mean Sperm Count Across Generations (by Replicate)',
    x = 'Generation',
    y = 'Mean Sperm Count',
    color = 'Replicate'
    )

print(plot_sperm_count)

# --- 5. Save Plots ---
ggsave(here(results_dir, "summary_traits_evolution.png"), plot_traits, width = 8, height = 6)
ggsave(here(results_dir, "summary_rsc_evolution.png"), plot_rsc, width = 8, height = 6)
ggsave(here(results_dir, "summary_correlation_evolution.png"), plot_correlation, width = 8, height = 6)
ggsave(here(results_dir, "summary_cor_vs_rsc_final_gen.png"), plot_cor_vs_rsc, width = 8, height = 6)
ggsave(here(results_dir, "replicates_sperm_count.png"), plot_sperm_count, width = 8, height = 6)

cat("All plots have been generated and saved to the '", results_dir, "' directory.\n")