library(tidyverse)

# Find the most recent results file
files <- list.files(pattern = "quick_test_results_.*\\.csv$")
latest_file <- files[which.max(file.info(files)$mtime)]

# Read the data
data <- read_csv(latest_file)

# Create plots
p1 <- ggplot(data, aes(x = Generation, y = MeanRSC, group = Rep, color = factor(Rep))) +
  geom_line() +
  labs(title = "Evolution of RSC Trait",
       y = "Mean RSC Value",
       color = "Replicate") +
  theme_minimal()

p2 <- ggplot(data, aes(x = Generation, y = MeanMale, group = Rep, color = factor(Rep))) +
  geom_line() +
  labs(title = "Evolution of Male Trait",
       y = "Mean Male Trait Value",
       color = "Replicate") +
  theme_minimal()

p3 <- ggplot(data, aes(x = Generation, y = MeanFemale, group = Rep, color = factor(Rep))) +
  geom_line() +
  labs(title = "Evolution of Female Trait",
       y = "Mean Female Trait Value",
       color = "Replicate") +
  theme_minimal()

# Save plots
ggsave("quick_test_rsc_evolution.png", p1, width = 8, height = 6)
ggsave("quick_test_male_evolution.png", p2, width = 8, height = 6)
ggsave("quick_test_female_evolution.png", p3, width = 8, height = 6)

print("Plots saved as PNG files")