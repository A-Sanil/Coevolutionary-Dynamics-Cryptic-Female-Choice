# ============================================================
# Simple Runner for RunModel_KBufferV2.jl
# ============================================================
# Run with:
#   julia --project=@. run_kbuffer_simpleV2.jl

# ========== EDITABLE PARAMETERS ==========
N = 500
generations = 5
replicates = 16

K = 1000
maintain_sex_ratio = true
show_gui = get(ENV, "SHOW_GUI", "0") == "1"

# HPC-friendly defaults: enforce single-threaded BLAS and pick worker count
using LinearAlgebra
BLAS.set_num_threads(1)

# Number of Julia worker processes to use for distributed runs.
# Prefer SLURM_CPUS_PER_TASK-1 when available, otherwise N_WORKERS env or Sys.CPU_THREADS-1
if !haskey(ENV, "N_WORKERS")
    slurm_cpus = tryparse(Int, get(ENV, "SLURM_CPUS_PER_TASK", ""))
    default_workers = if slurm_cpus != nothing && slurm_cpus > 1
        max(1, 15)
    else
        max(1, 15)
    end
    ENV["N_WORKERS"] = string(default_workers)
end

# Output root override (environment variable BRC_OUTPUT_ROOT already supported inside included file)
ENV["BRC_OUTPUT_ROOT"] = get(ENV, "BRC_OUTPUT_ROOT", get(ENV, "OUTPUT_ROOT", "CSV data"))

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
# Offspring mode should be a Symbol: :poisson, :logistic, :expdecay, or :gaussian
# default to :poisson for the original behavior
offspring_mode = :poisson

# progress and checkpointing
# how often to print status (generations)
progress_interval = parse(Int, get(ENV, "PROGRESS_INTERVAL", "250"))
# how often to write partial CSV checkpoints (0 = disabled)
checkpoint_interval = parse(Int, get(ENV, "CHECKPOINT_INTERVAL", "0"))

# execution mode
# true  = use distributed replicate-level parallelism (runsim)
# false = use single-process serial execution (runsim_serial)
use_parallel = true

# optional output root for cluster runs
# set BRC_OUTPUT_ROOT to a path like /global/scratch/users/<you>/kbuffer_outputs
output_root = get(ENV, "BRC_OUTPUT_ROOT", "CSV data")

# Smoke test mode: force a tiny serial run for quick validation.
smoke_test = get(ENV, "SMOKE_TEST", "0") == "1"
if smoke_test
    generations = 1
    replicates = 1
    use_parallel = false
    show_gui = false
end

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
println("  Parallel mode: $use_parallel")
println("  Smoke test mode: $smoke_test")
println("="^60)
println()

run_show_gui = use_parallel ? false : show_gui
run_mode = (use_parallel && replicates > 1) ? "parallel" : "serial"

elapsed_seconds = @elapsed begin
    global results = if use_parallel && replicates > 1
        runsim(
            replicates, N, mu, var, a, rsc, tradeoff, generations, -1,
            K, maintain_sex_ratio, run_show_gui, dynamic_offspring_mode, offspring_mode, offspring_scale,
            progress_interval, checkpoint_interval
        )
    else
        runsim_serial(
            replicates, N, mu, var, a, rsc, tradeoff, generations, -1,
            K, maintain_sex_ratio, run_show_gui, dynamic_offspring_mode, offspring_mode, offspring_scale,
            progress_interval, checkpoint_interval
        )
    end
end

data = DataFrame(results, [
    :MeanMale, :MeanFemale, :SDMale, :SDFemale, :cor,
    :MeanCount, :SDCount, :is, :int,
    :BMale, :GMale, :BFemale, :GFemale, :BSperm, :GSperm,
    :GMF, :GMS, :GFS, :a, :MeanRSC, :Generation, :MeanMates,
    :population, :males, :females, :offspring,
    :offspring_desired, :offspring_realized, :offspring_dropped,
    :Rep
])

sanitize_csv_dataframe!(data)

runstamp = Dates.format(now(), "yyyy-mm-ddTHH-MM-SS")
outdir = joinpath(output_root, "run_v2_" * runstamp)
isdir(outdir) || mkpath(outdir)
outfile = joinpath(outdir, "kbuffer_v2_results_" * runstamp * ".csv")
CSV.write(outfile, data)

println()
println("="^60)
println("✓ V2 simulation complete!")
println("="^60)
println("Run mode: $run_mode")
println("Elapsed seconds: $(round(elapsed_seconds, digits=3))")
println("Output root: $output_root")
println("Saved CSV to:")
println("  $outfile")
println("="^60)
