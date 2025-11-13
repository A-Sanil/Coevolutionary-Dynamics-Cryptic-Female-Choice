 # Minimal serial runner to test phenotype behavior
 include(joinpath(@__DIR__, "RunModel.jl"))
 using Dates, CSV, DataFrames

reps = 2
N = 100
mu = 1.25
var = (4*5^2/40)^0.5
a = 1.0
rsc = 0.25
tradeoff = true
gens = 50

println("Running minimal serial sim: reps=$reps N=$N gens=$gens")
results = runsim_serial(reps, N, mu, var, a, rsc, tradeoff, gens)

data = DataFrame(results, [:MeanMale,:MeanFemale,:SDMale,:SDFemale,:cor,:MeanCount,:SDCount,:is,:int,:BMale,:GMale,:BFemale,:GFemale,:BSperm,:GSperm,:GMF,:GMS,:GFS,:a,:MeanRSC,:Generation,:Rep])
outfile = joinpath(@__DIR__, "minimal_results.csv")
CSV.write(outfile, data)
println("Wrote: ", outfile)
