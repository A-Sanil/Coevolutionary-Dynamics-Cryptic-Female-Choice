# Script to run simulations for all three child function types
# Usage: julia run_simulations.jl
# Or set environment variables:
#   GENERATIONS=100 REPLICATES=5 POP_SIZE=200 julia run_simulations.jl

include("Child_RunModel.jl")

# Get parameters from environment or use defaults
gens = parse(Int, get(ENV, "GENERATIONS", "100"))
reps = parse(Int, get(ENV, "REPLICATES", "5"))
pop_size = parse(Int, get(ENV, "POP_SIZE", "200"))

# Standard parameters
mu = 1.25
var = (4*5^2/40)^0.5
a = 1.0
rsc = 0.25
tradeoff = true
d = -1

# Child models to test
child_models = [:parabola, :asymptote, :linear]

println("=" ^ 60)
println("Running simulations for child function testing")
println("=" ^ 60)
println("Parameters:")
println("  Generations: $gens")
println("  Replicates: $reps")
println("  Population size: $pop_size")
println("  Child models: $(join(string.(child_models), ", "))")
println("=" ^ 60)

# Run simulations for each child model
for child_mod in child_models
    println("\n" * "=" ^ 60)
    println("Running simulations for child model: $child_mod")
    println("=" ^ 60)
    
    # Run simulations (serial version for simplicity)
    results = runsim_serial(reps, pop_size, mu, var, a, rsc, tradeoff, gens, d, child_mod)
    
    # Create DataFrame with proper column names
    data = DataFrame(results, [:MeanMale,:MeanFemale,:SDMale,:SDFemale,:cor,:MeanCount,:SDCount,:is,:int,:BMale,:GMale,:BFemale,:GFemale,:BSperm,:GSperm,:GMF,:GMS,:GFS,:a,:MeanRSC,:Generation,:MeanMates,:MeanChildren,:TotalChildren,:Rep])
    
    # Add child_model column
    data.child_model = fill(string(child_mod), nrow(data))
    
    # Save to CSV in csv subfolder
    timestamp = Dates.format(now(), "yyyy-mm-dd_HHMMSS")
    outfile = joinpath("csv", "child_results_$(child_mod)_$(timestamp).csv")
    CSV.write(outfile, data)
    println("Saved results to: $outfile")
end

println("\n" * "=" ^ 60)
println("All simulations completed!")
println("=" ^ 60)

