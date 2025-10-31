using Random, Distributions, StatsBase, GLM, DataFrames, CSV, Dates

# Very simple test run
function quick_test()
    println("Starting minimal test simulation...")
    
    # Parameters
    N = 100           # Small population for testing
    gens = 50         # Very few generations
    mu = 1.25
    var = (4*5^2/40)^0.5
    a = 1
    tradeoff = true
    rsc = 0.25
    
    # Run simulation with minimal parameters
    try
        println("Running simulation...")
        results = sim(N, mu, var, a, rsc, tradeoff, gens)
        data = DataFrame(results, [:MeanMale,:MeanFemale,:SDMale,:SDFemale,:cor,:MeanCount,:SDCount,:is,:int,:BMale,:GMale,:BFemale,:GFemale,:BSperm,:GSperm,:GMF,:GMS,:GFS,:a,:MeanRSC,:Generation])
        
        # Add replicate column
        data[!,:Rep] .= 1
        
        # Save with timestamp in test_results directory
        timestamp = Dates.format(now(), "yyyy-mm-dd_HHMMSS")
        outdir = "test_results"
        isdir(outdir) || mkdir(outdir)
        outfile = joinpath(outdir, "minimal_test_results_$(timestamp).csv")
        CSV.write(outfile, data)
        println("Test completed! Results saved to: ", outfile)
        
        # Print key statistics
        println("\nFinal generation summary:")
        final_gen = data[data.Generation .== gens, :]
        println("Mean RSC: ", mean(final_gen.MeanRSC))
        println("Mean male trait: ", mean(final_gen.MeanMale))
        println("Mean female trait: ", mean(final_gen.MeanFemale))
        
        return true
    catch e
        println("Error in simulation: ", e)
        return false
    end
end

# Import core model code
include("RunModel.jl")

# Run test
quick_test()