# Quick serial run for testing and plotting
# Runs a small, fast simulation (serial) and writes a CSV with timestamped filename.

# Set environment variable to skip adding workers for quick serial run
ENV["NO_AUTO_ADDPROCS"] = "1"

using Dates
# include the main model (assumes RunModel.jl is in the same directory)
include(joinpath(@__DIR__, "RunModel.jl"))

# Quick parameters (small for speed)
reps = 1
N = 200
mu = 1.25
var = (4*5^2/40)^0.5
a = 1.0
rsc = 0.25
tradeoff = true
gens = 3

@info "Running quick serial sim: reps=$reps N=$N gens=$gens"
results = runsim_serial(reps, N, mu, var, a, rsc, tradeoff, gens)

using DataFrames, CSV
data = DataFrame(results, [:MeanMale,:MeanFemale,:SDMale,:SDFemale,:cor,:MeanCount,:SDCount,:is,:int,:BMale,:GMale,:BFemale,:GFemale,:BSperm,:GSperm,:GMF,:GMS,:GFS,:a,:MeanRSC,:Generation,:Rep])

timestamp = Dates.format(now(), "yyyy-mm-dd_HHMMSS")
outfile = joinpath(@__DIR__, "quick_sim_results_$(timestamp).csv")
CSV.write(outfile, data)
@info "Wrote quick results to: $outfile"

println("\n" * repeat("=", 60))
println("RSC VALUES BY GENERATION:")
println(repeat("=", 60))
for gen in 1:gens
    gen_data = data[data.Generation .== gen, :]
    if nrow(gen_data) > 0
        println("Generation $gen: MeanRSC = $(gen_data.MeanRSC[1])")
    end
end
println(repeat("=", 60))
println("\nDone. CSV: $outfile")
