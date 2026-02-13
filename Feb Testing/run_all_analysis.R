# Combined R script to run all analyses
# This script runs both analyze_children.R and plot_standard_results.R

cat(paste(rep("=", 60), collapse=""), "\n")
cat("Running all analyses\n")
cat(paste(rep("=", 60), collapse=""), "\n\n")

# Source the analysis scripts
cat("1. Running children analysis...\n")
source("analyze_children.R")

cat("\n2. Running standard plots...\n")
source("plot_standard_results.R")

cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("All analyses completed!\n")
cat(paste(rep("=", 60), collapse=""), "\n")

