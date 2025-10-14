# Single-process runner: avoids spawning workers and SharedArray
# Run with: julia Running_Sims/run_3_sims_single.jl

# Tell RunModel.jl not to auto addprocs
ENV["NO_AUTO_ADDPROCS"] = "1"

# Include the main model file
include("RunModel.jl")

using CSV, DataFrames

# Parameters
reps = 3
N = 200
mu = 1.25
var = (4*5^2/40)^0.5
a = 1.0
rsc = 0.25
tradeoff = true
gens = 50

all_results = Float64[]
results_mat = Array{Float64,2}(undef, reps*gens, 22)
row = 1
for rep in 1:reps
  println("Running rep $rep / $reps ...")
  df = sim(N, mu, var, a, rsc, tradeoff, gens)
  # sim returns a gens x 21 matrix; add Rep column to make 22 columns
  results_mat[(1+(rep-1)*gens):(rep*gens), 1:21] .= df
  results_mat[(1+(rep-1)*gens):(rep*gens), 22] .= rep
end

colnames = [:MeanMale,:MeanFemale,:SDMale,:SDFemale,:cor,:MeanCount,:SDCount,:is,:int,:BMale,:GMale,:BFemale,:GFemale,:BSperm,:GSperm,:GMF,:GMS,:GFS,:a,:rsc,:Generation,:Rep]

data = DataFrame(results_mat, colnames)
outfile = joinpath(@__DIR__, "run_3_sims_single_results.csv")
CSV.write(outfile, data)
println("Wrote: $outfile")
println("Done")
