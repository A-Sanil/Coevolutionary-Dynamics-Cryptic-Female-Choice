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

# 3) MeanRSC across generations
p3 <- ggplot(df, aes(x = Generation, y = MeanRSC, color = Rep)) +
  geom_line() +
  theme_minimal() +
  labs(title = 'Mean RSC across generations', y = 'Mean RSC')

# Save plots
ggsave('quick_mean_traits.png', p1, width = 8, height = 6)
ggsave('quick_mean_count.png', p2, width = 8, height = 4)
ggsave('quick_mean_rsc.png', p3, width = 8, height = 4)

cat('Plots saved: quick_mean_traits.png, quick_mean_count.png, quick_mean_rsc.png\n')
