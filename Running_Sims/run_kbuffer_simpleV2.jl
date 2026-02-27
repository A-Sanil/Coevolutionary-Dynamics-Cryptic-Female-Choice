# ============================================================
# Simple Runner for RunModel_KBufferV2.jl
# ============================================================
# Run with:
#   julia --project=@. run_kbuffer_simpleV2.jl

# ========== EDITABLE PARAMETERS ==========
N = 500
generations = 100
replicates = 1

K = 100000
maintain_sex_ratio = true
show_gui = true

mu = 5.0
var = 1.0
a = 1.0
rsc = 0.25
tradeoff = true

# V2 switches
# false = original fixed 2-offspring behavior
# true  = dynamic offspring count from RSC-based Poisson target
# offspring_scale tunes expected offspring intensity in dynamic mode

dynamic_offspring_mode = true
offspring_scale = 1.0

# =====================================================

include("RunModel_KBufferV2.jl")

using DataFrames, CSV, Dates

println("="^60)
println("Starting K-buffer V2 simulation...")
println("="^60)
println("  Initial population (N): $N ($(2*N) total: $N males + $N females)")
println("  Generations: $generations")
println("  Replicates: $replicates")
println("  Carrying capacity (K): $K")
println("  Maintain 50/50 sex ratio: $maintain_sex_ratio")
println("  Progress display: $(show_gui ? "Live progress bar" : "Periodic updates")")
println("  Dynamic offspring mode: $dynamic_offspring_mode")
println("  Offspring scale: $offspring_scale")
println("="^60)
println()

results = runsim_serial(
    replicates, N, mu, var, a, rsc, tradeoff, generations, -1,
    K, maintain_sex_ratio, show_gui, dynamic_offspring_mode, offspring_scale
)

data = DataFrame(results, [
    :MeanMale, :MeanFemale, :SDMale, :SDFemale, :cor,
    :MeanCount, :SDCount, :is, :int,
    :BMale, :GMale, :BFemale, :GFemale, :BSperm, :GSperm,
    :GMF, :GMS, :GFS, :a, :MeanRSC, :Generation, :MeanMates,
    :population, :males, :females, :offspring,
    :offspring_desired, :offspring_realized, :offspring_dropped,
    :Rep
])

runstamp = Dates.format(now(), "yyyy-mm-ddTHH-MM-SS")
outdir = joinpath("CSV data", "run_v2_" * runstamp)
isdir(outdir) || mkpath(outdir)
outfile = joinpath(outdir, "kbuffer_v2_results_" * runstamp * ".csv")
CSV.write(outfile, data)

println()
println("="^60)
println("✓ V2 simulation complete!")
println("="^60)
println("Saved CSV to:")
println("  $outfile")
println("="^60)
