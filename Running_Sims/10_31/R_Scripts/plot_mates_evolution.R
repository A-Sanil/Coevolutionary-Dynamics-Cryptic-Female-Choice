# R script to plot the evolution of mean mates per female across generations
library(tidyverse)

# Find the most recent simulation results file with mates data
parallel_files <- list.files(pattern = "parallel_sim_results_.*_with_mates\\.csv$", path = "CSV/", full.names = TRUE)

if(length(parallel_files) > 0) {
  file <- parallel_files[order(parallel_files, decreasing = TRUE)][1]
  cat("Using parallel simulation file with mates:", file, "\n")
} else {
  stop("No simulation results CSV with mates data found")
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
cat("Mean mates range:", round(min(df$MeanMates), 3), "to", round(max(df$MeanMates), 3), "\n")

# Plot 1a: Mean mates evolution across all replicates (overlay view)
p1a <- ggplot(df, aes(x = Generation, y = MeanMates, color = Rep)) +
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
    labels = scales::number_format(accuracy = 0.01),
    breaks = scales::pretty_breaks(n = 8)
  ) +
  labs(
    title = 'Evolution of Mean Mates per Female Across Generations',
    subtitle = 'Higher values indicate increased sperm competition risk',
    x = 'Generation',
    y = 'Mean Mates per Female',
    color = 'Replicate'
  )

# Plot 1b: Faceted by replicate with precise y-axis increments
p1b <- ggplot(df, aes(x = Generation, y = MeanMates, color = Rep)) +
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
    labels = scales::number_format(accuracy = 0.05),
    breaks = scales::pretty_breaks(n = 6)
  ) +
  labs(
    title = 'Evolution of Mean Mates per Female by Replicate',
    subtitle = 'Each panel shows one replicate with optimized y-axis scale',
    x = 'Generation',
    y = 'Mean Mates per Female'
  )

# Plot 2: Summary with mean and SE across replicates
summary_data <- df %>%
  group_by(Generation) %>%
  summarise(
    MeanMates = mean(MeanMates),
    SE_Mates = sd(MeanMates) / sqrt(n()),
    .groups = 'drop'
  )

p2 <- ggplot(summary_data, aes(x = Generation)) +
  geom_line(aes(y = MeanMates), color = "darkblue", linewidth = 1.2) +
  geom_ribbon(aes(ymin = MeanMates - SE_Mates, ymax = MeanMates + SE_Mates), 
              fill = "darkblue", alpha = 0.2) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  ) +
  scale_y_continuous(
    name = "Mean Mates per Female",
    labels = scales::number_format(accuracy = 0.01)
  ) +
  labs(
    title = 'Mean Mates per Female Across Generations',
    subtitle = 'Average across all replicates with standard error',
    x = 'Generation',
    y = 'Mean Mates per Female'
  )

# Plot 3: Mates vs RSC relationship
p3 <- ggplot(df, aes(x = MeanRSC, y = MeanMates, color = Rep)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 1.2, linetype = "dashed") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.position = "right"
  ) +
  scale_x_continuous(name = "Mean RSC Phenotype", labels = scales::number_format(accuracy = 0.1)) +
  scale_y_continuous(name = "Mean Mates per Female", labels = scales::number_format(accuracy = 0.01)) +
  labs(
    title = "Relationship between RSC and Mates per Female",
    subtitle = "Black line shows overall trend across all generations",
    color = "Replicate",
    caption = "Higher RSC should theoretically lead to more mates"
  )

# Plot 4: Evolution comparison - RSC vs Mates
p4 <- df %>%
  select(Generation, Rep, MeanRSC, MeanMates) %>%
  pivot_longer(cols = c(MeanRSC, MeanMates), names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = ifelse(Metric == "MeanMates", "Mean Mates", "Mean RSC")) %>%
  ggplot(aes(x = Generation, y = Value, color = Metric)) +
  geom_line(alpha = 0.7, linewidth = 1) +
  facet_wrap(~Rep, scales = "free_y", ncol = 3) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "bottom"
  ) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  labs(
    title = "Evolution of RSC and Mates by Replicate",
    x = "Generation",
    y = "Value",
    color = "Metric",
    subtitle = "Each panel shows one replicate"
  )

# Save all plots
ggsave('Images/mates_evolution_overlay.png', p1a, width = 10, height = 6, dpi = 300)
ggsave('Images/mates_evolution_by_rep.png', p1b, width = 14, height = 12, dpi = 300)
ggsave('Images/mates_evolution_summary.png', p2, width = 10, height = 6, dpi = 300)
ggsave('Images/mates_vs_rsc_scatter.png', p3, width = 10, height = 7, dpi = 300)
ggsave('Images/mates_rsc_comparison.png', p4, width = 14, height = 10, dpi = 300)

cat("\n=== Statistical Analysis: Mates vs RSC ===\n")
cat("\nOverall relationship (all data):\n")
model_all <- lm(MeanMates ~ MeanRSC, data = df)
print(summary(model_all))

cat("\n=== Plots Saved ===\n")
cat("1. mates_evolution_overlay.png - All replicates overlaid\n")
cat("2. mates_evolution_by_rep.png - Individual replicate panels with precise scales\n")
cat("3. mates_evolution_summary.png - Average with SE\n")
cat("4. mates_vs_rsc_scatter.png - Relationship between RSC and mates\n")
cat("5. mates_rsc_comparison.png - Side-by-side evolution comparison\n")

# Summary statistics table
cat("\n=== Summary Statistics Table ===\n")
summary_table <- df %>%
  group_by(Rep) %>%
  summarize(
    mean_mates = round(mean(MeanMates), 3),
    final_mates = round(last(MeanMates), 3),
    mates_range = paste(round(min(MeanMates), 3), "to", round(max(MeanMates), 3)),
    mean_rsc = round(mean(MeanRSC), 3),
    .groups = 'drop'
  )
print(summary_table)

