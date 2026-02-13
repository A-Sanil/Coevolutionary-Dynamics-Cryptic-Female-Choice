# R script for analyzing children per generation and mates-to-children relationship
# Reads simulation results and generates plots showing:
# 1. Children born per generation
# 2. Relationship between mates and children
# 3. Comparison across child function types (parabola, asymptote, linear)

# --- 1. Setup ---
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("ggplot2")) install.packages("ggplot2")

library(tidyverse)
library(ggplot2)

# --- 2. Load Data ---
# Find the most recent results file in csv subfolder
csv_files <- list.files("csv", pattern = ".*results.*\\.csv$", full.names = TRUE)

if(length(csv_files) == 0) {
  stop("No simulation results CSV found in current directory")
}

latest_csv <- csv_files[which.max(file.info(csv_files)$mtime)]
cat("Loading data from:", latest_csv, "\n")
sim_data <- read_csv(latest_csv, show_col_types = FALSE)

# Check if required columns exist
required_cols <- c("Generation", "MeanMates", "MeanChildren", "TotalChildren", "Rep")
missing <- setdiff(required_cols, names(sim_data))
if (length(missing) > 0) {
  stop("Missing required columns: ", paste(missing, collapse = ", "))
}

# Convert Generation and Rep to appropriate types
sim_data <- sim_data %>% 
  mutate(
    Generation = as.integer(Generation), 
    Rep = as.factor(Rep),
    MeanMates = as.numeric(MeanMates),
    MeanChildren = as.numeric(MeanChildren),
    TotalChildren = as.numeric(TotalChildren)
  )

# Check if there's a child_model column (if running multiple models)
if ("child_model" %in% names(sim_data)) {
  sim_data$child_model <- as.factor(sim_data$child_model)
} else {
  # If no child_model column, assume single model (extract from filename if possible)
  sim_data$child_model <- "unknown"
}

# --- 3. Data Aggregation ---
# Calculate mean and standard error across replicates for each generation
summary_data <- sim_data %>%
  group_by(Generation) %>%
  summarise(
    MeanChildren = mean(MeanChildren, na.rm = TRUE),
    SE_Children = sd(MeanChildren, na.rm = TRUE) / sqrt(n()),
    TotalChildren = mean(TotalChildren, na.rm = TRUE),
    SE_TotalChildren = sd(TotalChildren, na.rm = TRUE) / sqrt(n()),
    MeanMates = mean(MeanMates, na.rm = TRUE),
    SE_Mates = sd(MeanMates, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  )

# If multiple child models, aggregate by model and generation
if ("child_model" %in% names(sim_data) && length(unique(sim_data$child_model)) > 1) {
  summary_by_model <- sim_data %>%
    group_by(Generation, child_model) %>%
    summarise(
      MeanChildren = mean(MeanChildren, na.rm = TRUE),
      SE_Children = sd(MeanChildren, na.rm = TRUE) / sqrt(n()),
      TotalChildren = mean(TotalChildren, na.rm = TRUE),
      SE_TotalChildren = sd(TotalChildren, na.rm = TRUE) / sqrt(n()),
      MeanMates = mean(MeanMates, na.rm = TRUE),
      SE_Mates = sd(MeanMates, na.rm = TRUE) / sqrt(n()),
      .groups = 'drop'
    )
}

# --- 4. Generate Plots ---

# Plot 1: Mean Children per Female Over Generations
plot_children_gen <- ggplot(summary_data, aes(x = Generation, y = MeanChildren)) +
  geom_line(color = "steelblue", linewidth = 1.2) +
  geom_ribbon(aes(ymin = MeanChildren - SE_Children, ymax = MeanChildren + SE_Children), 
              fill = "steelblue", alpha = 0.2) +
  geom_point(color = "darkblue", size = 2) +
  theme_minimal() +
  labs(
    title = "Mean Children per Female Over Generations",
    subtitle = "Mean ± SE across replicates",
    x = "Generation",
    y = "Mean Children per Female"
  )

ggsave("graphs/plot_children_per_generation.png", plot_children_gen, width = 10, height = 6, dpi = 300)
cat("Saved: graphs/plot_children_per_generation.png\n")

# Plot 2: Total Children Born Per Generation
plot_total_children <- ggplot(summary_data, aes(x = Generation, y = TotalChildren)) +
  geom_line(color = "darkgreen", linewidth = 1.2) +
  geom_ribbon(aes(ymin = TotalChildren - SE_TotalChildren, ymax = TotalChildren + SE_TotalChildren), 
              fill = "darkgreen", alpha = 0.2) +
  geom_point(color = "forestgreen", size = 2) +
  theme_minimal() +
  labs(
    title = "Total Children Born Per Generation",
    subtitle = "Mean ± SE across replicates",
    x = "Generation",
    y = "Total Children Born"
  )

ggsave("graphs/plot_total_children_per_generation.png", plot_total_children, width = 10, height = 6, dpi = 300)
cat("Saved: graphs/plot_total_children_per_generation.png\n")

# Plot 3: Mates vs Children Relationship (scatter plot)
plot_mates_children <- ggplot(sim_data, aes(x = MeanMates, y = MeanChildren)) +
  geom_point(alpha = 0.5, size = 1.5, color = "steelblue") +
  geom_smooth(method = "lm", se = TRUE, color = "red", linewidth = 1.2) +
  theme_minimal() +
  labs(
    title = "Relationship Between Mates and Children",
    subtitle = "All data points across all generations and replicates",
    x = "Mean Number of Mates per Female",
    y = "Mean Children per Female"
  )

ggsave("graphs/plot_mates_vs_children.png", plot_mates_children, width = 10, height = 8, dpi = 300)
cat("Saved: graphs/plot_mates_vs_children.png\n")

# Plot 4: Mates and Children Over Time (dual axis)
plot_combined <- ggplot(summary_data, aes(x = Generation)) +
  geom_line(aes(y = MeanMates, color = "Mates"), linewidth = 1.2) +
  geom_ribbon(aes(ymin = MeanMates - SE_Mates, ymax = MeanMates + SE_Mates, fill = "Mates"), 
              alpha = 0.2) +
  geom_line(aes(y = MeanChildren * 2, color = "Children (scaled)"), linewidth = 1.2) +
  geom_ribbon(aes(ymin = (MeanChildren - SE_Children) * 2, 
                  ymax = (MeanChildren + SE_Children) * 2, fill = "Children (scaled)"), 
              alpha = 0.2) +
  scale_y_continuous(
    name = "Mean Mates per Female",
    sec.axis = sec_axis(~ . / 2, name = "Mean Children per Female")
  ) +
  scale_color_manual(values = c("Mates" = "blue", "Children (scaled)" = "red")) +
  scale_fill_manual(values = c("Mates" = "blue", "Children (scaled)" = "red")) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(
    title = "Mates and Children Over Generations",
    subtitle = "Mean ± SE across replicates",
    x = "Generation",
    color = "Variable",
    fill = "Variable"
  )

ggsave("graphs/plot_mates_and_children_over_time.png", plot_combined, width = 12, height = 8, dpi = 300)
cat("Saved: graphs/plot_mates_and_children_over_time.png\n")

# Plot 5: All Replicates Overlay - Children
plot_all_children <- ggplot(sim_data, aes(x = Generation, y = MeanChildren, color = Rep)) +
  geom_line(alpha = 0.6, linewidth = 0.8) +
  theme_minimal() +
  labs(
    title = "Mean Children per Female: All Replicates",
    subtitle = "Each line = one replicate",
    x = "Generation",
    y = "Mean Children per Female",
    color = "Replicate"
  )

ggsave("graphs/plot_all_replicates_children.png", plot_all_children, width = 12, height = 8, dpi = 300)
cat("Saved: graphs/plot_all_replicates_children.png\n")

# Plot 6: If multiple child models, compare them
if (exists("summary_by_model") && length(unique(sim_data$child_model)) > 1) {
  plot_models_comparison <- ggplot(summary_by_model, aes(x = Generation, y = MeanChildren, color = child_model)) +
    geom_line(linewidth = 1.2) +
    geom_ribbon(aes(ymin = MeanChildren - SE_Children, ymax = MeanChildren + SE_Children, fill = child_model), 
                alpha = 0.2) +
    theme_minimal() +
    labs(
      title = "Mean Children per Female: Comparison Across Child Models",
      subtitle = "Mean ± SE across replicates",
      x = "Generation",
      y = "Mean Children per Female",
      color = "Child Model",
      fill = "Child Model"
    )
  
  ggsave("graphs/plot_child_models_comparison.png", plot_models_comparison, width = 12, height = 8, dpi = 300)
  cat("Saved: graphs/plot_child_models_comparison.png\n")
}

# --- 5. Statistical Summary ---
cat("\n=== STATISTICAL SUMMARY ===\n")
cat("Mean children per female (overall):", mean(sim_data$MeanChildren, na.rm = TRUE), "\n")
cat("Mean mates per female (overall):", mean(sim_data$MeanMates, na.rm = TRUE), "\n")
cat("Mean total children per generation (overall):", mean(sim_data$TotalChildren, na.rm = TRUE), "\n")

# Correlation between mates and children
correlation <- cor(sim_data$MeanMates, sim_data$MeanChildren, use = "complete.obs")
cat("Correlation between mates and children:", correlation, "\n")

cat("\n=== ALL PLOTS GENERATED ===\n")

