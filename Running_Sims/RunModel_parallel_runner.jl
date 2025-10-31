# Parallel runner for RunModel.jl
# Starts local worker processes and loads the main RunModel.jl on all workers.
# Controls:
#  - NO_AUTO_ADDPROCS=1  -> skip adding workers (useful on laptops)
#  - RUNWORKERS=<n>     -> explicitly set number of workers to add

using Distributed
# decide number of workers: environment override or use (CPU_THREADS - 1)
const DEFAULT_WORKERS = haskey(ENV, "RUNWORKERS") ? parse(Int, ENV["RUNWORKERS"]) : max(1, Sys.CPU_THREADS - 1)
if !haskey(ENV, "NO_AUTO_ADDPROCS")
  addprocs(DEFAULT_WORKERS)
  @info "Added $DEFAULT_WORKERS worker processes"
else
  @info "NO_AUTO_ADDPROCS set: not adding worker processes"
end

# Ensure required packages are available on all processes and load the model file on each worker.
# Using an absolute path to the RunModel.jl so workers can include it from the same filesystem.
mainfile = joinpath(@__DIR__, "RunModel.jl")
@info "Including model file on all processes: $mainfile"
@everywhere using Random, Distributions, StatsBase, GLM, DataFrames, CSV, SharedArrays, Dates
@everywhere include(mainfile)

# make CSV/DataFrame helpers available on the main process as well
using CSV, DataFrames

# At this point the functions `sim`, `runsim`, `runsim_serial` etc. are defined on all processes.
# To run a parallel simulation from this runner, you can either:
#  - call `runsim(reps,N,mu,var,a,rsc,tradeoff,gens)` from the REPL after starting julia with this file
#  - or modify this script to call `runsim(...)` with parameters you choose and write outputs to CSV

# Example (commented): a small parallel test run (uncomment to execute automatically)
# if get(ENV, "RUN_AUTO_TEST", "0") == "1"
#   @info "Running automatic parallel test run..."
#   results = runsim(4, 200, 1.25, (4*5^2/40)^0.5, 12.5, 0.25, true, 100)
#   data = DataFrame(results, [:MeanMale,:MeanFemale,:SDMale,:SDFemale,:cor,:MeanCount,:SDCount,:is,:int,:BMale,:GMale,:BFemale,:GFemale,:BSperm,:GSperm,:GMF,:GMS,:GFS,:a,:MeanRSC,:Generation,:Rep])
#   CSV.write(joinpath(@__DIR__, "parallel_test_results.csv"), data)
#   @info "Auto test complete."
# end

# Small helpers for writing CSVs safely (timestamp backups) or appending
function safe_csv_write(df::DataFrame, fname::AbstractString)
  # make parent dir if needed
  parent = dirname(fname)
  if parent != "" && !isdir(parent)
    mkpath(parent)
  end
  if isfile(fname)
    bak = fname * "." * Dates.format(now(), "yyyy-mm-dd_HHMMSS") * ".bak"
    cp(fname, bak; force=true)
    @info "Backed up existing '$fname' -> '$bak'"
  end
  CSV.write(fname, df)
  @info "Wrote CSV: $fname"
  return fname
end

function append_csv(df::DataFrame, fname::AbstractString)
  # if file missing, write header; otherwise append rows without header
  if !isfile(fname)
    safe_csv_write(df, fname)
  else
    open(fname, "a") do io
      CSV.write(io, df; writeheader=false)
    end
    @info "Appended $(size(df,1)) rows to $fname"
  end
  return fname
end
