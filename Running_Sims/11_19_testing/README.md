# 11/19 Testing Results

This folder contains simulation results and generated plots from testing on November 19, 2025.

## Contents

- `parallel_sim_results_2025-11-12_175453.csv` - Simulation results data
- Generated plots (PNG files) - Created by running plotting scripts

## Generating Plots

To generate all plots, run one of the following from the `Running_Sims` directory:

### Option 1: Generate all plots at once
```r
source("generate_all_plots_11_19.R")
```

### Option 2: Run individual plotting scripts
```r
# Comprehensive validation plots
source("plot_and_validate_all.R")

# Quick summary plots
source("plot_quick_results.R")

# Detailed analysis plots
source("plot_11_12_2025.R")
```

All plots will be saved to this `11_19_testing/` folder.

## Data Description

The CSV file contains the following columns:
- `MeanMale`, `MeanFemale` - Mean trait values for males and females
- `SDMale`, `SDFemale` - Standard deviations
- `cor` - Correlation between male and female traits
- `MeanCount`, `SDCount` - Mean and SD of sperm count
- `MeanRSC` - Mean Risk of Sperm Competition phenotype
- `MeanMates` - Mean number of mates per female
- `Generation` - Generation number
- `Rep` - Replicate number
- Selection coefficients: `BMale`, `GMale`, `BFemale`, `GFemale`, `BSperm`, `GSperm`, `GMF`, `GMS`, `GFS`
- Other metrics: `is`, `int`, `a`

## Expected Plots

After running the plotting scripts, you should see:
- Validation plots (RSC vs MeanMates relationship)
- Trait evolution plots
- Correlation evolution
- Sperm count evolution
- Summary multi-panel plots
- And more...

