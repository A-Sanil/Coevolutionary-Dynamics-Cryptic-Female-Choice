# Small runner to execute 3 simulations quickly and save results
# Place this file in Running_Sims and run with `julia Running_Sims/run_3_sims.jl` from repo root

# Include the main model file (it contains addprocs and @everywhere declarations)
include("RunModel.jl")

# Parameters for a quick smoke-run; adjust as needed
reps = 3           # number of replicate simulations
N = 200            # population size (smaller for speed)
mu = 1.25
var = (4*5^2/40)^0.5
a = 1.0
rsc = 0.25
tradeoff = true
gens = 50          # number of generations per replicate (small for a quick run)

aprintln = println
println("Starting runsim with reps=$reps, N=$N, gens=$gens")

# Run the simulations (this uses the runsim function defined in RunModel.jl)
results = runsim(reps, N, mu, var, a, rsc, tradeoff, gens)

# Save results to CSV
using DataFrames, CSV
# Column names mirror those used elsewhere in the project
colnames = [:MeanMale,:MeanFemale,:SDMale,:SDFemale,:cor,:MeanCount,:SDCount,:is,:int,:BMale,:GMale,:BFemale,:GFemale,:BSperm,:GSperm,:GMF,:GMS,:GFS,:a,:rsc,:Generation,:Rep]

data = DataFrame(results, colnames)
outfile = joinpath(@__DIR__, "run_3_sims_results.csv")
CSV.write(outfile, data)
println("Wrote results to: $outfile")
println("Done.")
