# ============================================================
# Simple Runner for RunModel_KBuffer.jl
# ============================================================
# Edit parameters below and run: julia --project=@. run_kbuffer_simple.jl
#
# Or on HPC with multiple cores, allow distributed processing:
#   julia --project=@. run_kbuffer_simple.jl

# ========== EDITABLE PARAMETERS ==========

# Population Parameters
N = 100              # Initial population size (total = 2*N: N males + N females)
generations = 20     # Number of generations to simulate
replicates = 1       # Number of replicate runs

# Carrying Capacity
K = N^2              # Carrying capacity options:
                     #   N^2     = 10,000 (for N=100)
                     #   5*N     = 500 (for N=100)
                     #   10*N    = 1,000 (for N=100)
                     #   Custom: set any value

# Sex Ratio at Capacity
maintain_sex_ratio = true   
                     # true:  Maintain 50/50 male/female split when pop > K
                     # false: Random sampling (sex ratio can drift)

# Progress Display
show_gui = true      # true:  Show live progress bar with population counts
                     # false: Print updates every 10 generations only

# ========== Advanced Parameters (usually keep defaults) ==========
mu = 5.0             # Mean trait value
var = 1.0            # Trait variance
a = 1.0              # Selection parameter
rsc = 0.25           # RSC (reproductive skew) parameter
tradeoff = true      # Enable tradeoff

# =====================================================

# Load simulation code
include("RunModel_KBuffer.jl")

using DataFrames, CSV, Dates

println("="^60)
println("Starting K-buffer simulation...")
println("="^60)
println("  Initial population (N): $N ($(2*N) total: $N males + $N females)")
println("  Generations: $generations")
println("  Replicates: $replicates")
println("  Carrying capacity (K): $K")
println("  Maintain 50/50 sex ratio: $maintain_sex_ratio")
println("  Progress display: $(show_gui ? "Live progress bar" : "Periodic updates")")
println("="^60)
println()

# Run simulation
results = runsim_serial(replicates, N, mu, var, a, rsc, tradeoff, generations, -1, K, maintain_sex_ratio, show_gui)

# Create DataFrame with all columns
data = DataFrame(results, [
    :MeanMale, :MeanFemale, :SDMale, :SDFemale, :cor, 
    :MeanCount, :SDCount, :is, :int, 
    :BMale, :GMale, :BFemale, :GFemale, :BSperm, :GSperm, 
    :GMF, :GMS, :GFS, :a, :MeanRSC, :Generation, :MeanMates,
    :population, :males, :females, :offspring, :Rep
])

# Save results
runstamp = Dates.format(now(), "yyyy-mm-ddTHH-MM-SS")
outdir = joinpath("CSV data", "run_" * runstamp)
isdir(outdir) || mkpath(outdir)
outfile = joinpath(outdir, "kbuffer_results_" * runstamp * ".csv")
CSV.write(outfile, data)

println()
println("="^60)
println("✓ Simulation complete!")
println("="^60)
println("Saved CSV to:")
println("  $outfile")
println("="^60)
