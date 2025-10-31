# October 31, 2025 Simulation Results

## Overview
This folder contains the results of 10 simulations run for 100 generations each, examining the coevolutionary dynamics of cryptic female choice with evolving RSC (Risk of Sperm Competition).

## Folders

### CSV/
Contains 4 simulation result files:
- `parallel_sim_results_2025-10-31_144852.csv` - **Latest run** (10 reps × 100 gens)
- `parallel_sim_results_2025-10-31_144020.csv`
- `parallel_sim_results_2025-10-31_142310.csv`
- `parallel_sim_results_2025-10-31_131737.csv`

All CSV files contain the following columns:
- MeanMale, MeanFemale, SDMale, SDFemale
- cor (correlation between male and female traits)
- MeanCount, SDCount (sperm count)
- is (opportunity for selection)
- int (intercept from selection model)
- BMale, GMale, BFemale, GFemale, BSperm, GSperm (beta and gamma coefficients)
- GMF, GMS, GFS (gamma coefficients for interactions)
- a (selection strength parameter)
- MeanRSC (mean RSC phenotype)
- Generation, Rep

### Images/
Contains 12 plots generated from the simulation results:

**Quick Plots (from plot_quick_results.R):**
- `quick_mean_traits.png` - Male and female trait evolution across replicates
- `quick_mean_count.png` - Sperm count evolution
- `quick_mean_rsc.png` - RSC evolution

**Summary Plots (from analyze_results.R):**
- `summary_traits_evolution.png` - Mean trait evolution with SE
- `summary_rsc_evolution.png` - Mean RSC evolution with SE
- `summary_correlation_evolution.png` - Mean correlation evolution with SE
- `summary_cor_vs_rsc_final_gen.png` - Correlation vs RSC in final generation
- `replicates_sperm_count.png` - Sperm count by replicate

**Correlation Analysis Plots (from analyze_correlation_vs_rsc.R):**
- `correlation_vs_rsc_scatter.png` - Overall correlation vs RSC relationship
- `correlation_vs_rsc_by_rep.png` - Correlation vs RSC by replicate
- `correlation_and_rsc_evolution.png` - Evolution of both metrics
- `correlation_vs_rsc_binned.png` - Binned analysis

### R_Scripts/
Contains the R scripts used to generate all plots:
- `analyze_correlation_vs_rsc.R` - Correlation vs RSC analysis
- `analyze_results.R` - Summary statistics and evolution plots
- `plot_quick_results.R` - Quick visualization plots

## Simulation Parameters
- **Generations:** 100
- **Replicates:** 10
- **Population Size (N):** 200
- **Mean (μ):** 1.25
- **Variance (σ²):** (4×5²/40)^0.5 ≈ 1.58
- **Selection strength (a):** 1
- **Tradeoff mode:** true
- **RSC:** Evolving (starting at 0.25 phenotype)

## Key Findings
- Mean RSC phenotype range: 1.16 - 5.20
- Overall mean RSC: 2.15
- No significant correlation between RSC and male-female trait correlation (p = 0.73)
- RSC shows substantial evolutionary variation across replicates

## Running the Analysis
To regenerate plots, run from the Running_Sims directory:
```r
source("plot_quick_results.R")
source("analyze_results.R")
source("analyze_correlation_vs_rsc.R")
```

## Files
- `run_10_sims_100_gens.jl` - Julia script used to run simulations

