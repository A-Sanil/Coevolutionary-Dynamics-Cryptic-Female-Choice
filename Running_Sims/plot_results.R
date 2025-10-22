# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Check if the file exists
if (!file.exists("test_simulation_results.csv")) {
  stop("Error: test_simulation_results.csv not found in the current directory")
}

# Read the data
simulation_data <- read.csv("test_simulation_results.csv")

# Print the structure of the data to verify it loaded correctly
print("Data structure:")
str(simulation_data)

# Create plots directory if it doesn't exist
dir.create("plots", showWarnings = FALSE)

# 1. Plot Mean RSC over generations
p1 <- ggplot(data, aes(x=Generation, y=MeanRSC, color=factor(Rep))) +
  geom_line() +
  theme_minimal() +
  labs(title="Evolution of RSC Trait",
       x="Generation",
       y="Mean RSC",
       color="Replicate") +
  theme(legend.position="right")
ggsave("plots/rsc_evolution.png", p1, width=10, height=6)

# 2. Plot Male and Female trait means
p2 <- ggplot(data, aes(x=Generation)) +
  geom_line(aes(y=MeanMale, color="Male")) +
  geom_line(aes(y=MeanFemale, color="Female")) +
  facet_wrap(~Rep) +
  theme_minimal() +
  labs(title="Male and Female Trait Evolution",
       x="Generation",
       y="Mean Trait Value",
       color="Sex") +
  theme(legend.position="bottom")
ggsave("plots/trait_evolution.png", p2, width=12, height=8)

# 3. Plot correlation between male and female traits
p3 <- ggplot(data, aes(x=Generation, y=cor, color=factor(Rep))) +
  geom_line() +
  theme_minimal() +
  labs(title="Male-Female Trait Correlation",
       x="Generation",
       y="Correlation",
       color="Replicate") +
  theme(legend.position="right")
ggsave("plots/trait_correlation.png", p3, width=10, height=6)

# 4. Create summary statistics
summary_stats <- data %>%
  group_by(Rep) %>%
  summarize(
    mean_rsc = mean(MeanRSC),
    final_rsc = last(MeanRSC),
    mean_correlation = mean(cor),
    mean_male_trait = mean(MeanMale),
    mean_female_trait = mean(MeanFemale)
  )

# Save summary statistics
write.csv(summary_stats, "plots/summary_statistics.csv", row.names=FALSE)

print("Plots have been saved in the 'plots' directory")