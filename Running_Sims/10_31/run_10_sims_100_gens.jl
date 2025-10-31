using Random, Distributions, StatsBase, GLM, DataFrames, CSV, Dates

include("RunModel.jl") # Include the main model code

# Set simulation parameters
gens = 100        # 100 generations
reps = 10         # 10 simulations
N = 200           # Population size
mu = 1.25         # Keep original mu
var = (4*5^2/40)^0.5  # Keep original var
a = 1             # Simple a value
tradeoff = true   # Use tradeoff mode
rsc = 0.25        # RSC parameter (now evolving, but kept for compatibility)

# Run 10 simulations
println("Starting 10 simulations with 100 generations each...")
println("Parameters: N=$N, gens=$gens, reps=$reps, mu=$mu, var=$var, a=$a, tradeoff=$tradeoff")
results = runsim_serial(reps, N, mu, var, a, rsc, tradeoff, gens)

# Convert to DataFrame with proper column names
data = DataFrame(results, [:MeanMale,:MeanFemale,:SDMale,:SDFemale,:cor,:MeanCount,:SDCount,:is,:int,:BMale,:GMale,:BFemale,:GFemale,:BSperm,:GSperm,:GMF,:GMS,:GFS,:a,:MeanRSC,:Generation,:Rep])

# Save results with timestamp
timestamp = Dates.format(now(), "yyyy-mm-dd_HHMMSS")
outfile = "parallel_sim_results_$(timestamp).csv"
CSV.write(outfile, data)
println("Simulations completed! Results saved to: ", outfile)

# Print summary statistics
println("\n=== Summary Statistics ===")
final_gens = data[data.Generation .== gens, :]
println("Mean RSC value (final gen): ", round(mean(final_gens.MeanRSC), digits=4))
println("Mean male trait (final gen): ", round(mean(final_gens.MeanMale), digits=4))
println("Mean female trait (final gen): ", round(mean(final_gens.MeanFemale), digits=4))
println("Mean correlation (final gen): ", round(mean(final_gens.cor), digits=4))
println("\nRSC range across all runs: ", round(minimum(data.MeanRSC), digits=4), " - ", round(maximum(data.MeanRSC), digits=4))
println("Total rows: ", nrow(data))

