# R script for analyzing and plotting simulation results
# Generates all standard graphs from simulation data (based on plot_results.R)
# Modified to work with Child_RunModel.jl output

# --- 1. Setup ---
# Install and load necessary packages
if (!require("tidyverse")) install.packages("tidyverse")

library(tidyverse)

# --- 2. Load Data ---
# Find the most recent results file in csv subfolder
parallel_files <- list.files("csv", pattern = "parallel_sim_results_.*\\.csv$", full.names = TRUE)
runmodel_files <- list.files("csv", pattern = "runmodel_results_.*\\.csv$", full.names = TRUE)
quick_files <- list.files("csv", pattern = "quick_sim_results_.*\\.csv$", full.names = TRUE)
child_files <- list.files("csv", pattern = "child.*results.*\\.csv$", full.names = TRUE)
all_results_files <- list.files("csv", pattern = ".*results.*\\.csv$", full.names = TRUE)

all_files <- c(parallel_files, runmodel_files, quick_files, child_files, all_results_files)

if(length(all_files) == 0) {
  stop("No simulation results CSV found in current directory")
}

latest_csv <- all_files[which.max(file.info(all_files)$mtime)]
cat("Loading data from:", latest_csv, "\n")
sim_data <- read_csv(latest_csv, show_col_types = FALSE)

# Convert Generation and Rep to factors/integers for plotting
sim_data <- sim_data %>% mutate(Generation = as.integer(Generation), Rep = as.factor(Rep))

# --- 3. Data Aggregation ---
# Calculate the mean and standard error across replicates for each generation
summary_data <- sim_data %>%
  group_by(Generation) %>%
  summarise(
    MeanMale = mean(MeanMale, na.rm = TRUE),
    SE_Male = sd(MeanMale, na.rm = TRUE) / sqrt(n()),
    MeanFemale = mean(MeanFemale, na.rm = TRUE),
    SE_Female = sd(MeanFemale, na.rm = TRUE) / sqrt(n()),
    MeanRSC = mean(MeanRSC, na.rm = TRUE),
    SE_RSC = sd(MeanRSC, na.rm = TRUE) / sqrt(n()),
    MeanCount = mean(MeanCount, na.rm = TRUE),
    SE_Count = sd(MeanCount, na.rm = TRUE) / sqrt(n()),
    MeanCor = mean(cor, na.rm = TRUE),
    SE_Cor = sd(cor, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  )

# If MeanMates and MeanChildren columns exist, add them
if ("MeanMates" %in% names(sim_data)) {
  summary_data <- summary_data %>%
    left_join(
      sim_data %>%
        group_by(Generation) %>%
        summarise(
          MeanMates = mean(MeanMates, na.rm = TRUE),
          SE_Mates = sd(MeanMates, na.rm = TRUE) / sqrt(n()),
          MeanChildren = mean(MeanChildren, na.rm = TRUE),
          SE_Children = sd(MeanChildren, na.rm = TRUE) / sqrt(n()),
          .groups = 'drop'
        ),
      by = "Generation"
    )
}

# --- 4. Generate Plots ---

# Plot 1: Evolution of Mean Male and Female Traits
plot_traits <- ggplot(summary_data, aes(x = Generation)) +
  geom_line(aes(y = MeanMale, color = "Male"), linewidth = 1.2) +
  geom_ribbon(aes(ymin = MeanMale - SE_Male, ymax = MeanMale + SE_Male), 
              fill = "blue", alpha = 0.2) +
  geom_line(aes(y = MeanFemale, color = "Female"), linewidth = 1.2) +
  geom_ribbon(aes(ymin = MeanFemale - SE_Female, ymax = MeanFemale + SE_Female), 
              fill = "red", alpha = 0.2) +
  scale_color_manual(values = c("Male" = "blue", "Female" = "red")) +
  theme_minimal() +
  labs(
    title = "Evolution of Mean Male and Female Traits",
    subtitle = "Mean ± SE across replicates",
    x = "Generation",
    y = "Trait Value",
    color = "Sex"
  )

ggsave("graphs/plot_traits_evolution.png", plot_traits, width = 10, height = 6, dpi = 300)
cat("Saved: graphs/plot_traits_evolution.png\n")

# Plot 2: Mean Sperm Count Evolution
plot_sperm <- ggplot(summary_data, aes(x = Generation, y = MeanCount)) +
  geom_line(color = "steelblue", linewidth = 1.2) +
  geom_ribbon(aes(ymin = MeanCount - SE_Count, ymax = MeanCount + SE_Count), 
              fill = "steelblue", alpha = 0.2) +
  geom_point(color = "darkblue", size = 2) +
  theme_minimal() +
  labs(
    title = "Evolution of Mean Sperm Count",
    subtitle = "Mean ± SE across replicates",
    x = "Generation",
    y = "Mean Sperm Count"
  )

ggsave("graphs/plot_sperm_evolution.png", plot_sperm, width = 10, height = 6, dpi = 300)
cat("Saved: graphs/plot_sperm_evolution.png\n")

# Plot 3: Mean RSC Evolution
plot_rsc <- ggplot(summary_data, aes(x = Generation, y = MeanRSC)) +
  geom_line(color = "red", linewidth = 1.2) +
  geom_ribbon(aes(ymin = MeanRSC - SE_RSC, ymax = MeanRSC + SE_RSC), 
              fill = "red", alpha = 0.2) +
  geom_point(color = "darkred", size = 2) +
  theme_minimal() +
  labs(
    title = "Evolution of Mean RSC (Risk of Sperm Competition)",
    subtitle = "Mean ± SE across replicates",
    x = "Generation",
    y = "Mean RSC"
  )

ggsave("graphs/plot_rsc_evolution.png", plot_rsc, width = 10, height = 6, dpi = 300)
cat("Saved: graphs/plot_rsc_evolution.png\n")

# Plot 4: Trait Correlation Evolution
plot_cor <- ggplot(summary_data, aes(x = Generation, y = MeanCor)) +
  geom_line(color = "darkgreen", linewidth = 1.2) +
  geom_ribbon(aes(ymin = MeanCor - SE_Cor, ymax = MeanCor + SE_Cor), 
              fill = "darkgreen", alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_point(color = "forestgreen", size = 2) +
  theme_minimal() +
  labs(
    title = "Evolution of Male-Female Trait Correlation",
    subtitle = "Mean ± SE across replicates",
    x = "Generation",
    y = "Correlation (r)"
  )

ggsave("graphs/plot_correlation_evolution.png", plot_cor, width = 10, height = 6, dpi = 300)
cat("Saved: graphs/plot_correlation_evolution.png\n")

# Plot 5: All Replicates Overlay - Sperm Count
plot_all_sperm <- ggplot(sim_data, aes(x = Generation, y = MeanCount, color = Rep)) +
  geom_line(alpha = 0.6, linewidth = 0.8) +
  theme_minimal() +
  labs(
    title = "Mean Sperm Count: All Replicates",
    subtitle = "Each line = one replicate",
    x = "Generation",
    y = "Mean Sperm Count",
    color = "Replicate"
  )

ggsave("graphs/plot_all_replicates_sperm.png", plot_all_sperm, width = 12, height = 8, dpi = 300)
cat("Saved: graphs/plot_all_replicates_sperm.png\n")

# Plot 6: All Replicates Overlay - RSC
plot_all_rsc <- ggplot(sim_data, aes(x = Generation, y = MeanRSC, color = Rep)) +
  geom_line(alpha = 0.6, linewidth = 0.8) +
  theme_minimal() +
  labs(
    title = "Mean RSC: All Replicates",
    subtitle = "Each line = one replicate",
    x = "Generation",
    y = "Mean RSC",
    color = "Replicate"
  )

ggsave("graphs/plot_all_replicates_rsc.png", plot_all_rsc, width = 12, height = 8, dpi = 300)
cat("Saved: graphs/plot_all_replicates_rsc.png\n")

# Plot 7: Sperm vs RSC Relationship
plot_sperm_rsc <- ggplot(sim_data, aes(x = MeanRSC, y = MeanCount)) +
  geom_point(alpha = 0.5, size = 1.5, color = "steelblue") +
  geom_smooth(method = "lm", se = TRUE, color = "red", linewidth = 1.2) +
  theme_minimal() +
  labs(
    title = "Mean Sperm Count vs RSC Relationship",
    subtitle = "All data points across all generations and replicates",
    x = "Mean RSC",
    y = "Mean Sperm Count"
  )

ggsave("graphs/plot_sperm_vs_rsc.png", plot_sperm_rsc, width = 10, height = 8, dpi = 300)
cat("Saved: graphs/plot_sperm_vs_rsc.png\n")

# Plot 8: Mean Mates Evolution (if available)
if ("MeanMates" %in% names(summary_data)) {
  plot_mates <- ggplot(summary_data, aes(x = Generation, y = MeanMates)) +
    geom_line(color = "purple", linewidth = 1.2) +
    geom_ribbon(aes(ymin = MeanMates - SE_Mates, ymax = MeanMates + SE_Mates), 
                fill = "purple", alpha = 0.2) +
    geom_point(color = "darkviolet", size = 2) +
    theme_minimal() +
    labs(
      title = "Evolution of Mean Mates per Female",
      subtitle = "Mean ± SE across replicates",
      x = "Generation",
      y = "Mean Mates per Female"
    )
  
  ggsave("graphs/plot_mates_evolution.png", plot_mates, width = 10, height = 6, dpi = 300)
  cat("Saved: graphs/plot_mates_evolution.png\n")
}

# Plot 9: Mean Children Evolution (if available)
if ("MeanChildren" %in% names(summary_data)) {
  plot_children <- ggplot(summary_data, aes(x = Generation, y = MeanChildren)) +
    geom_line(color = "orange", linewidth = 1.2) +
    geom_ribbon(aes(ymin = MeanChildren - SE_Children, ymax = MeanChildren + SE_Children), 
                fill = "orange", alpha = 0.2) +
    geom_point(color = "darkorange", size = 2) +
    theme_minimal() +
    labs(
      title = "Evolution of Mean Children per Female",
      subtitle = "Mean ± SE across replicates",
      x = "Generation",
      y = "Mean Children per Female"
    )
  
  ggsave("graphs/plot_children_evolution.png", plot_children, width = 10, height = 6, dpi = 300)
  cat("Saved: graphs/plot_children_evolution.png\n")
}

cat("\n=== ALL PLOTS GENERATED ===\n")
cat("Total plots created: ", ifelse("MeanMates" %in% names(summary_data), 9, 7), "\n")

