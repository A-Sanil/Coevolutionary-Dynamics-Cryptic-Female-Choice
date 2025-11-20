# Master script to generate all plots for 11/19 testing
# This script runs all plotting scripts and generates comprehensive plots

cat("========================================\n")
cat("Generating All Plots for 11/19 Testing\n")
cat("========================================\n\n")

# Ensure output directory exists
if (!dir.exists("11_19_testing")) {
  dir.create("11_19_testing")
  cat("Created directory: 11_19_testing\n")
}

# Check if CSV file exists
csv_file <- "11_19_testing/parallel_sim_results_2025-11-12_175453.csv"
if (!file.exists(csv_file)) {
  stop("CSV file not found: ", csv_file)
}

cat("Found CSV file:", csv_file, "\n\n")

# Source the plotting scripts
cat("1. Running plot_and_validate_all.R...\n")
tryCatch({
  source("plot_and_validate_all.R")
  cat("   ✓ plot_and_validate_all.R completed\n\n")
}, error = function(e) {
  cat("   ✗ Error in plot_and_validate_all.R:", e$message, "\n\n")
})

cat("2. Running plot_quick_results.R...\n")
tryCatch({
  source("plot_quick_results.R")
  cat("   ✓ plot_quick_results.R completed\n\n")
}, error = function(e) {
  cat("   ✗ Error in plot_quick_results.R:", e$message, "\n\n")
})

cat("3. Running plot_11_12_2025.R...\n")
tryCatch({
  source("plot_11_12_2025.R")
  cat("   ✓ plot_11_12_2025.R completed\n\n")
}, error = function(e) {
  cat("   ✗ Error in plot_11_12_2025.R:", e$message, "\n\n")
})

# List all generated plots
cat("========================================\n")
cat("Generated Plots Summary\n")
cat("========================================\n")
plot_files <- list.files("11_19_testing", pattern = "\\.png$", full.names = FALSE)
if (length(plot_files) > 0) {
  cat("Total plots generated:", length(plot_files), "\n\n")
  for (i in seq_along(plot_files)) {
    cat(sprintf("%2d. %s\n", i, plot_files[i]))
  }
} else {
  cat("No plots found in 11_19_testing directory\n")
}

cat("\n========================================\n")
cat("All plots saved to: 11_19_testing/\n")
cat("========================================\n")

