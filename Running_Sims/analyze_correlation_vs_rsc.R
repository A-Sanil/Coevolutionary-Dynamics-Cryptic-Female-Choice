# R script to analyze how correlation between male and female traits
# varies based on the evolving RSC phenotype
# This addresses the question: How does sperm competition risk affect trait correlation?

library(tidyverse)

# Find the most recent simulation results file
parallel_files <- list.files(pattern = "parallel_sim_results_.*\\.csv$", path = ".", full.names = TRUE)
quick_files <- list.files(pattern = "quick_sim_results_.*\\.csv$", path = ".", full.names = TRUE)

if(length(parallel_files) > 0) {
  file <- parallel_files[order(parallel_files, decreasing = TRUE)][1]
  cat("Using parallel simulation file:", file, "\n")
} else if(length(quick_files) > 0) {
  file <- quick_files[order(quick_files, decreasing = TRUE)][1]
  cat("Using quick simulation file:", file, "\n")
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
cat("Mean RSC range:", round(min(df$MeanRSC), 3), "to", round(max(df$MeanRSC), 3), "\n")

# 1. Correlation vs RSC scatter plot (all generations)
# This shows the relationship between RSC and male-female correlation
p1 <- ggplot(df, aes(x = MeanRSC, y = cor, color = Rep)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 1.2) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.position = "right"
  ) +
  scale_x_continuous(
    name = "Mean RSC Phenotype",
    labels = scales::number_format(accuracy = 0.1)
  ) +
  scale_y_continuous(
    name = "Male-Female Trait Correlation",
    labels = scales::number_format(accuracy = 0.01)
  ) +
  labs(
    title = "Relationship between RSC and Male-Female Trait Correlation",
    subtitle = "Each point = one generation. Black line shows overall trend.",
    color = "Replicate",
    caption = "Higher RSC = higher sperm competition risk"
  )

# 2. Correlation vs RSC faceted by Replicate
# Shows if the relationship is consistent across replicates
p2 <- ggplot(df, aes(x = MeanRSC, y = cor)) +
  geom_point(alpha = 0.7, size = 1) +
  geom_smooth(method = "lm", se = TRUE, color = "red", linewidth = 1) +
  facet_wrap(~Rep, scales = "free_x") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    strip.text = element_text(face = "bold")
  ) +
  labs(
    title = "Correlation vs RSC by Replicate",
    x = "Mean RSC Phenotype",
    y = "Male-Female Trait Correlation",
    subtitle = "Each panel shows one replicate. Red line shows replicate-specific trend."
  )

# 3. Correlation and RSC evolution over time (side by side)
p3 <- df %>%
  select(Generation, Rep, cor, MeanRSC) %>%
  pivot_longer(cols = c(cor, MeanRSC), names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = ifelse(Metric == "cor", "Correlation", "Mean RSC")) %>%
  ggplot(aes(x = Generation, y = Value, color = Rep)) +
  geom_line(alpha = 0.8, linewidth = 1) +
  facet_wrap(~Metric, scales = "free_y", ncol = 1) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    strip.text = element_text(face = "bold", size = 12)
  ) +
  labs(
    title = "Evolution of Correlation and RSC Over Time",
    x = "Generation",
    y = "Value",
    color = "Replicate",
    subtitle = "Top: Male-Female trait correlation. Bottom: Mean RSC phenotype."
  )

# 4. Statistical analysis: Correlation vs RSC relationship
cat("\n=== Statistical Analysis ===\n")
cat("\nOverall relationship (all data):\n")
model_all <- lm(cor ~ MeanRSC, data = df)
print(summary(model_all))

cat("\nRelationship by replicate:\n")
models_by_rep <- df %>%
  group_by(Rep) %>%
  do(model = lm(cor ~ MeanRSC, data = .)) %>%
  mutate(
    slope = coef(model)[2],
    intercept = coef(model)[1],
    r_squared = summary(model)$r.squared,
    p_value = summary(model)$coefficients[2,4]
  ) %>%
  select(-model)

print(models_by_rep)

# 5. Binned analysis: Average correlation within RSC bins
# This helps see if there's a non-linear relationship
df_binned <- df %>%
  mutate(RSC_bin = cut(MeanRSC, breaks = 10, include.lowest = TRUE)) %>%
  group_by(RSC_bin) %>%
  summarize(
    mean_RSC = mean(MeanRSC),
    mean_cor = mean(cor),
    sd_cor = sd(cor),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(
    se_cor = sd_cor / sqrt(n),
    ci_lower = mean_cor - 1.96 * se_cor,
    ci_upper = mean_cor + 1.96 * se_cor
  )

p4 <- ggplot(df_binned, aes(x = mean_RSC, y = mean_cor)) +
  geom_point(size = 3, color = "darkblue") +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.1, color = "darkblue") +
  geom_smooth(method = "lm", se = TRUE, color = "red", linewidth = 1.2) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  ) +
  scale_x_continuous(name = "Mean RSC Phenotype (binned)", labels = scales::number_format(accuracy = 0.1)) +
  scale_y_continuous(name = "Mean Correlation (within bin)", labels = scales::number_format(accuracy = 0.01)) +
  labs(
    title = "Binned Analysis: Average Correlation by RSC Levels",
    subtitle = "Error bars show 95% confidence intervals. Red line shows linear trend.",
    caption = paste("Data divided into", nrow(df_binned), "bins")
  )

# Save all plots
ggsave('correlation_vs_rsc_scatter.png', p1, width = 10, height = 7, dpi = 300)
ggsave('correlation_vs_rsc_by_rep.png', p2, width = 12, height = 8, dpi = 300)
ggsave('correlation_and_rsc_evolution.png', p3, width = 10, height = 8, dpi = 300)
ggsave('correlation_vs_rsc_binned.png', p4, width = 10, height = 7, dpi = 300)

cat("\n=== Plots Saved ===\n")
cat("1. correlation_vs_rsc_scatter.png - Overall relationship\n")
cat("2. correlation_vs_rsc_by_rep.png - Relationship by replicate\n")
cat("3. correlation_and_rsc_evolution.png - Evolution over time\n")
cat("4. correlation_vs_rsc_binned.png - Binned analysis\n")

# Summary statistics table
cat("\n=== Summary Statistics Table ===\n")
summary_table <- df %>%
  group_by(Rep) %>%
  summarize(
    mean_RSC = round(mean(MeanRSC), 3),
    mean_cor = round(mean(cor), 4),
    cor_range = paste(round(min(cor), 4), "to", round(max(cor), 4)),
    .groups = "drop"
  )
print(summary_table)

