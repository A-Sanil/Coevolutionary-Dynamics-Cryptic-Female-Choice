# ============================================================
# Simple Runner for RunModel_KBufferV2.jl
# ============================================================
# Run with:
#   julia --project=@. run_kbuffer_simpleV2.jl

# ========== EDITABLE PARAMETERS ==========
N = 5000  # 5000 males + 5000 females = 10000 initial total
generations = 100
replicates = 5

K = 10000
maintain_sex_ratio = true
show_gui = false  # Turned off for parallel runs so progress bars don't overlap!

mu = 5.0
var = 1.0
a = 1.0
rsc = 0.25
tradeoff = true

# V2 switches
# false = original fixed 2-offspring behavior
# true  = dynamic offspring count from selected offspring mode
# offspring_scale tunes expected offspring intensity in Poisson mode

dynamic_offspring_mode = true

# Toggle which evolutionary curves to run by commenting out (adding a #) to the ones you want to skip.
modes_to_run = [
    #:poisson,   # Direct scaling linked to RSC
    #:logistic,  # S-curve plateau
    #:expdecay,  # Exponential decay (penalizes more mates immediately)
    :gaussian   # Gaussian peak (penalizes having more or less than the optimal mates)
]

offspring_scale = 1.0
max_mates = 100      # Increased from 5 to effectively remove the ceiling cap
max_offspring = 100  # Increased from 5 to allow offspring numbers to climb naturally
offspring_r = 1.4
offspring_c = 3.0
offspring_mu = 3.0
offspring_sigma = 1.0

# =====================================================

include(joinpath(@__DIR__, "RunModel_KBufferV2.jl"))

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
println("  Modes to run: $modes_to_run")
println("  Offspring scale: $offspring_scale")
println("  Max mates: $max_mates")
println("  Max offspring: $max_offspring")
println("  Offspring gaussian params: mu=$offspring_mu, sigma=$offspring_sigma")
println("="^60)
println()

runstamp = Dates.format(now(), "yyyy-mm-ddTHH-MM-SS")
outdir = joinpath(@__DIR__, "CSV data March", "run_v2_batch_" * runstamp)
isdir(outdir) || mkpath(outdir)

# Start a timer for the entire batch
total_start_time = time()

for mode in modes_to_run
    println("\n" * "="^60)
    println("▶ RUNNING MODE: $mode")
    println("="^60)

    println("Executing simulation (tracking performance)...")
    timed_run = @timed runsim(
        replicates, N, mu, var, a, rsc, tradeoff, generations, -1,
        K, maintain_sex_ratio, show_gui, dynamic_offspring_mode,
        mode, offspring_scale,
        max_mates, max_offspring,
        offspring_r, offspring_c,
        offspring_mu, offspring_sigma
    )
    
    results = timed_run.value
    run_time_seconds = timed_run.time
    run_memory_mb = timed_run.bytes / (1024^2)
    println("✓ Finished in $(round(run_time_seconds, digits=2)) seconds using $(round(run_memory_mb, digits=2)) MB of memory.")

    data = DataFrame(results, [
        :MeanMale, :MeanFemale, :SDMale, :SDFemale, :cor,
        :MeanCount, :SDCount, :is, :int,
        :BMale, :GMale, :BFemale, :GFemale, :BSperm, :GSperm,
        :GMF, :GMS, :GFS, :a, :MeanRSC, :Generation, :MeanMates,
        :population, :males, :females, :offspring,
        :offspring_desired, :offspring_realized, :offspring_dropped,
        :Rep
    ])

    # add cumulative total offspring born per replicate
    if :Rep in names(data)
      data = combine(groupby(data, :Rep)) do sub
        sub.TotalOffspringBorn = cumsum(sub.offspring)
        sub
      end
    else
      data.TotalOffspringBorn = cumsum(data.offspring)
    end

    # Save the performance stats directly into the CSV as new columns
    data[!, :Runtime_Seconds] .= run_time_seconds
    data[!, :Memory_MB] .= run_memory_mb

    outfile = joinpath(outdir, "kbuffer_v2_results_$(mode)_" * runstamp * ".csv")
    CSV.write(outfile, data)

    println("Saved CSV to: $outfile")
end

total_duration = time() - total_start_time

println("\n" * "="^60)
println("✓ ALL EXPERIMENTAL MODES COMPLETE!")
println("✓ Total Batch Runtime: $(round(total_duration, digits=2)) seconds ($(round(total_duration/60, digits=2)) minutes)")
println("="^60)
