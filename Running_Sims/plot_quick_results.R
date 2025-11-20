# R script to plot quick_sim_results CSV produced by run_quick.jl
# Install packages once: install.packages(c('tidyverse'))

library(tidyverse)

# adjust path if needed - check for both quick_sim_results and parallel_sim_results
# Prioritize parallel_sim_results if available, otherwise use most recent
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

df <- read_csv(file)

# Convert Generation and Rep to factors/integers
df <- df %>% mutate(Generation = as.integer(Generation), Rep = as.factor(Rep))

# Example plots:
# 1) MeanMale and MeanFemale across generations (faceted by Rep)
p1 <- df %>% pivot_longer(cols = c('MeanMale','MeanFemale'), names_to = 'Trait', values_to = 'Value') %>%
  ggplot(aes(x = Generation, y = Value, color = Trait)) +
  geom_line(alpha = 0.8) +
  facet_wrap(~Rep, ncol = 1) +
  theme_minimal() +
  labs(title = 'Mean Male and Female trait across generations', y = 'Trait value')

# 2) MeanCount (sperm count) across generations
p2 <- ggplot(df, aes(x = Generation, y = MeanCount, color = Rep)) +
  geom_line() +
  theme_minimal() +
  labs(title = 'Mean sperm count across generations', y = 'Mean sperm count')

# 3) MeanRSC across generations - with precise scaling and detailed labels
# MeanRSC represents the mean Risk of Sperm Competition phenotype in the population
# RSC phenotype = sum of alleles (additive model), can be negative (natural evolution)
# Higher RSC values indicate higher average mating rates (females mate with more males)
# This affects the Poisson lambda parameter for number of mates per female
# Note: RSC can be negative (evolves freely), only clamped to 0 when sampling mates (Poisson requires non-negative)
p3 <- ggplot(df, aes(x = Generation, y = MeanRSC, color = Rep)) +
  geom_line(linewidth = 1.2, alpha = 0.8) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.position = "right"
  ) +
  scale_y_continuous(
    name = "Mean RSC Phenotype Value",
    labels = scales::number_format(accuracy = 0.01),  # Show 2 decimal places
    breaks = scales::pretty_breaks(n = 10)  # More breaks for better precision
  ) +
  labs(
    title = 'Evolution of Mean RSC (Risk of Sperm Competition)',
    subtitle = 'RSC determines mating rates: higher values = more mates per female\nRSC can evolve to negative values (natural evolution)',
    x = 'Generation',
    y = 'Mean RSC Phenotype Value',
    color = 'Replicate'
  )

# 4) MeanMates across generations - should match RSC closely
p4 <- ggplot(df, aes(x = Generation, y = MeanMates, color = Rep)) +
  geom_line(linewidth = 1.2, alpha = 0.8) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.position = "right"
  ) +
  scale_y_continuous(
    name = "Mean Mates per Female",
    labels = scales::number_format(accuracy = 0.01),
    breaks = scales::pretty_breaks(n = 10)
  ) +
  labs(
    title = 'Evolution of Mean Mates per Female',
    subtitle = 'Mean number of mates per female (should match RSC closely)',
    x = 'Generation',
    y = 'Mean Mates per Female',
    color = 'Replicate'
  )

# 5) RSC vs MeanMates scatter plot - validation
p5 <- ggplot(df, aes(x = MeanRSC, y = MeanMates)) +
  geom_point(alpha = 0.5, color = "steelblue", size = 1.5) +
  geom_smooth(method = "lm", se = TRUE, color = "red", linewidth = 1.2, linetype = "dashed") +
  geom_abline(intercept = 0, slope = 1, color = "black", linewidth = 1, linetype = "dotted") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  ) +
  scale_x_continuous(name = "Mean RSC", labels = scales::number_format(accuracy = 0.1)) +
  scale_y_continuous(name = "Mean Mates per Female", labels = scales::number_format(accuracy = 0.1)) +
  labs(
    title = 'Validation: RSC vs Mean Mates Relationship',
    subtitle = 'Red line: linear fit | Black dotted line: y=x (perfect match)',
    caption = 'Each point = one generation from one replicate'
  )

# Add summary statistics text
rsc_range <- paste0("Range: ", round(min(df$MeanRSC), 3), " - ", round(max(df$MeanRSC), 3))
rsc_mean <- paste0("Overall Mean: ", round(mean(df$MeanRSC), 3))

# Save plots
ggsave('quick_mean_traits.png', p1, width = 8, height = 6, dpi = 300)
ggsave('quick_mean_count.png', p2, width = 8, height = 4, dpi = 300)
ggsave('quick_mean_rsc.png', p3, width = 10, height = 6, dpi = 300)  # Wider for better visibility
ggsave('quick_mean_mates.png', p4, width = 10, height = 6, dpi = 300)
ggsave('quick_rsc_vs_mates.png', p5, width = 10, height = 7, dpi = 300)

# Print summary statistics
cat('\n=== RSC Summary Statistics ===\n')
cat('Overall Mean RSC:', round(mean(df$MeanRSC), 4), '\n')
cat('Overall Min RSC:', round(min(df$MeanRSC), 4), '\n')
cat('Overall Max RSC:', round(max(df$MeanRSC), 4), '\n')
cat('\nFinal Generation RSC by Replicate:\n')
final_gen <- max(df$Generation)
final_rsc <- df %>% 
  filter(Generation == final_gen) %>% 
  select(Rep, MeanRSC) %>% 
  arrange(Rep)
print(final_rsc)

# Calculate and print validation statistics
cat('\n=== RSC vs MeanMates Validation ===\n')
overall_cor <- cor(df$MeanRSC, df$MeanMates, use = "complete.obs")
cat('Overall correlation (RSC, MeanMates):', round(overall_cor, 4), '\n')
if (overall_cor > 0.95) {
  cat('  -> EXCELLENT: Strong correlation indicates proper matching\n')
} else if (overall_cor > 0.90) {
  cat('  -> GOOD: Strong correlation\n')
} else {
  cat('  -> WARNING: Correlation lower than expected\n')
}

model <- lm(MeanMates ~ MeanRSC, data = df)
slope <- coef(model)[2]
intercept <- coef(model)[1]
cat('Linear model: MeanMates ~ MeanRSC\n')
cat('  Slope:', round(slope, 4), '(expected ≈ 1.0)\n')
cat('  Intercept:', round(intercept, 4), '(expected ≈ 0.0)\n')

mean_diff <- mean(df$MeanRSC - df$MeanMates, na.rm=TRUE)
cat('Mean difference (RSC - MeanMates):', round(mean_diff, 4), '\n')

cat('\nPlots saved: quick_mean_traits.png, quick_mean_count.png, quick_mean_rsc.png, quick_mean_mates.png, quick_rsc_vs_mates.png\n')
