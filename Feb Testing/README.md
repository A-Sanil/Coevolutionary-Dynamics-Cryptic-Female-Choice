# Feb Testing - Child Function Simulations

This folder contains the combined model that integrates child function options (parabola, asymptote, linear) into the main coevolutionary dynamics model.

## Files

- **Child_RunModel.jl**: Main simulation code that combines the original RunModel.jl with child function options
- **run_simulations.jl**: Script to run simulations for all three child function types
- **analyze_children.R**: R script to analyze children per generation and mates-to-children relationships
- **plot_standard_results.R**: R script to generate standard plots (traits, sperm, RSC, etc.)

## Usage

### Running Simulations

1. **Quick test run** (small parameters):
```bash
cd "Feb Testing"
GENERATIONS=50 REPLICATES=3 POP_SIZE=100 NO_AUTO_ADDPROCS=1 julia run_simulations.jl
```

2. **Full simulation run**:
```bash
cd "Feb Testing"
GENERATIONS=100 REPLICATES=5 POP_SIZE=200 NO_AUTO_ADDPROCS=1 julia run_simulations.jl
```

3. **Run with parallel processing** (remove NO_AUTO_ADDPROCS):
```bash
cd "Feb Testing"
GENERATIONS=100 REPLICATES=5 POP_SIZE=200 julia run_simulations.jl
```

### Child Function Options

The child function can be set in the simulation by passing `child_model` parameter:
- `:parabola` - Offspring peaks at intermediate compatibility
- `:asymptote` - Offspring rises then levels off
- `:linear` - Offspring increases linearly with compatibility

### Generating Plots

After running simulations, generate plots:

1. **Children analysis** (children per generation, mates-to-children relationship):
```bash
cd "Feb Testing"
Rscript analyze_children.R
```

2. **Standard plots** (traits, sperm, RSC evolution):
```bash
cd "Feb Testing"
Rscript plot_standard_results.R
```

## Output Files

- `child_results_<model>_<timestamp>.csv`: Simulation results for each child model
- `plot_*.png`: Generated plots showing various relationships and evolution

## Parameters

Default parameters (can be modified in run_simulations.jl):
- Generations: 100
- Replicates: 5
- Population size: 200
- mu: 1.25
- var: (4*5^2/40)^0.5
- a: 1.0
- rsc: 0.25
- tradeoff: true

## Notes

- The model maintains a fixed population size by sampling from offspring
- Children are determined by compatibility between male and female traits
- The number of mates is determined by the RSC (Risk of Sperm Competition) trait
- All three child function types are tested and compared

