using Random, Distributions, StatsBase, GLM, DataFrames, CSV, Dates

include("RunModel.jl") # Include the main model code

# Set test parameters (small values for quick testing)
gens = 100        # Reduced generations for quick test
N = 200           # Small population size
mu = 1.25         # Keep original mu
var = (4*5^2/40)^0.5  # Keep original var
a = 1             # Simple a value
tradeoff = true   # Use tradeoff mode
rsc = 0.25        # RSC parameter (now evolving)

# Run a quick test with 3 replicates
println("Starting quick test simulation...")
results = runsim_serial(3, N, mu, var, a, rsc, tradeoff, gens)
data = DataFrame(results, [:MeanMale,:MeanFemale,:SDMale,:SDFemale,:cor,:MeanCount,:SDCount,:is,:int,:BMale,:GMale,:BFemale,:GFemale,:BSperm,:GSperm,:GMF,:GMS,:GFS,:a,:MeanRSC,:Generation,:Rep])

# Save results with timestamp
timestamp = Dates.format(now(), "yyyy-mm-dd_HHMMSS")
outfile = "quick_test_results_$(timestamp).csv"
CSV.write(outfile, data)
println("Test completed! Results saved to: ", outfile)

# Print summary statistics
println("\nSummary of final generation across replicates:")
final_gens = data[data.Generation .== gens, :]
println("Mean RSC value: ", mean(final_gens.MeanRSC))
println("Mean male trait: ", mean(final_gens.MeanMale))
println("Mean female trait: ", mean(final_gens.MeanFemale))