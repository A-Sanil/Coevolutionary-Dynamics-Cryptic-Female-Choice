# Analysis for child_testing.jl simulations
# Reads the newest child_testing_results_*.csv (or path passed via INPUT env var)
# Produces population and mean-trait trajectories for each offspring curve.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(glue)
  library(tidyr)
})

# Determine script directory even when sourced or run via Rscript
this_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(...) NA_character_)
script_dir <- if (!is.na(this_file)) dirname(this_file) else getwd()
input_path <- Sys.getenv("INPUT")

# Helper to pick the most recent results file in the script directory
if (identical(input_path, "")) {
  files <- list.files(script_dir, pattern = "^child_testing_results_.*\\.csv$", full.names = TRUE)
  if (length(files) == 0) {
    stop("No child_testing_results_*.csv file found. Set INPUT env var or run child_testing.jl first.")
  }
  input_path <- files[order(file.info(files)$mtime, decreasing = TRUE)][1]
}

message("Reading: ", input_path)
raw <- read_csv(input_path, show_col_types = FALSE)

# Clean up and ensure expected columns exist
required_cols <- c("model", "rep", "generation", "population", "mean_trait",
                   "mean_offspring_per_pair", "mean_expected_offspring")
missing <- setdiff(required_cols, names(raw))
if (length(missing) > 0) {
  stop("Missing expected columns: ", paste(missing, collapse = ", "))
}

df <- raw %>%
  mutate(model = factor(model, levels = c("parabola", "asymptote", "mixed")))

# Population summary across replicates
pop_summary <- df %>%
  group_by(model, generation) %>%
  summarise(
    mean_pop = mean(population, na.rm = TRUE),
    median_pop = median(population, na.rm = TRUE),
    lo = quantile(population, 0.1, na.rm = TRUE),
    hi = quantile(population, 0.9, na.rm = TRUE),
    .groups = "drop"
  )

# Mean trait summary
trait_summary <- df %>%
  group_by(model, generation) %>%
  summarise(
    mean_trait = mean(mean_trait, na.rm = TRUE),
    lo = quantile(mean_trait, 0.1, na.rm = TRUE),
    hi = quantile(mean_trait, 0.9, na.rm = TRUE),
    .groups = "drop"
  )

# Plot population trajectories
p1 <- ggplot(pop_summary, aes(generation, mean_pop, color = model, fill = model)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, color = NA) +
  geom_line(size = 1) +
  labs(title = "Population size over time", y = "Population", x = "Generation") +
  theme_minimal(base_size = 12)

# Plot trait trajectories
p2 <- ggplot(trait_summary, aes(generation, mean_trait, color = model, fill = model)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, color = NA) +
  geom_line(size = 1) +
  labs(title = "Mean trait over time", y = "Trait value", x = "Generation") +
  theme_minimal(base_size = 12)

out_pop <- file.path(script_dir, "child_testing_population.png")
out_trait <- file.path(script_dir, "child_testing_trait.png")

ggsave(out_pop, p1, width = 7, height = 4, dpi = 300)
ggsave(out_trait, p2, width = 7, height = 4, dpi = 300)

message("Wrote plots:\n - ", out_pop, "\n - ", out_trait)

# Console summary for quick inspection
summary_tbl <- df %>%
  group_by(model) %>%
  summarise(final_pop = population[generation == max(generation)][1],
            final_trait = mean_trait[generation == max(generation)][1],
            avg_offspring = mean(mean_offspring_per_pair, na.rm = TRUE),
            .groups = "drop")

print(summary_tbl)
