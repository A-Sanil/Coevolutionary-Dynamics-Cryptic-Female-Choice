# Parallel run for longer simulations
# Uses distributed computing across multiple cores for faster execution
# Best for multiple replicates with many generations

using Dates
# Don't set NO_AUTO_ADDPROCS - we want parallel processing!
# include the main model (assumes RunModel.jl is in the same directory)
include(joinpath(@__DIR__, "RunModel.jl"))

# Parallel simulation parameters
reps = 8  # Under 10 as requested, good number for parallel execution
N = 200
mu = 1.25
var = (4*5^2/40)^0.5
a = 1.0
rsc = 0.25
tradeoff = true
gens = 200  # More generations to see evolution - parallel makes this feasible

@info "Running PARALLEL sim: reps=$reps N=$N gens=$gens"
@info "Using $(nworkers()) worker processes for parallel execution"

# Use parallel runsim function (distributes replicates across workers)
results = runsim(reps, N, mu, var, a, rsc, tradeoff, gens)

using DataFrames, CSV
data = DataFrame(results, [:MeanMale,:MeanFemale,:SDMale,:SDFemale,:cor,:MeanCount,:SDCount,:is,:int,:BMale,:GMale,:BFemale,:GFemale,:BSperm,:GSperm,:GMF,:GMS,:GFS,:a,:MeanRSC,:Generation,:Rep])

timestamp = Dates.format(now(), "yyyy-mm-dd_HHMMSS")
outfile = joinpath(@__DIR__, "parallel_sim_results_$(timestamp).csv")
CSV.write(outfile, data)
@info "Wrote parallel results to: $outfile"

println("\n" * repeat("=", 60))
println("SIMULATION SUMMARY:")
println(repeat("=", 60))
println("Replicates: $reps")
println("Generations: $gens")
println("Final generation RSC values by replicate:")
for rep in 1:reps
    rep_data = data[(data.Generation .== gens) .& (data.Rep .== rep), :]
    if nrow(rep_data) > 0
        println("  Rep $rep: MeanRSC = $(round(rep_data.MeanRSC[1], digits=3))")
    end
end
println(repeat("=", 60))
println("\nDone. CSV: $outfile")

