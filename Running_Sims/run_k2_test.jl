using Dates, DataFrames, CSV

include("Child_sim.jl")
include("refactored_sim.jl")

function run_test()
  Nf = 50
  Nm = 50
  generations = 30
  mu = 0.0
  var = 1.0
  a = 1.0
  rsc = 1.0
  tradeoff = false

  dfall = sim_k2(Nf, Nm, mu, var, a, rsc, tradeoff, generations; buffer_multiplier=2, scale_param=1.0, csvpath=nothing, kids_per_fertilization=2)

  # map dfall columns to expected CSV columns used by plot_results.R
  # dfall columns: [MeanMale, MeanFemale, sdMale, sdFemale, cor, Meansperm, Stdsperm, is, coef1..coef10, a, MeanRSC, gen, MeanMates]
  gen_col = Int.(dfall[:,21])
  out = DataFrame(Generation=gen_col,
                  Rep=fill(1, length(gen_col)),
                  MeanMale=dfall[:,1],
                  MeanFemale=dfall[:,2],
                  MeanRSC=dfall[:,20],
                  MeanCount=dfall[:,6],
                  cor=dfall[:,5])

  # create main CSV data folder and timestamped run subfolder
  timestamp = Dates.format(now(), "yyyy-mm-ddTHH-MM-SS")
  main_dir = joinpath(@__DIR__, "CSV data")
  isdir(main_dir) || mkpath(main_dir)
  outdir = joinpath(main_dir, "run_" * timestamp)
  isdir(outdir) || mkpath(outdir)
  csvname = joinpath(outdir, "quick_sim_results_" * timestamp * ".csv")
  CSV.write(csvname, out)
  println("Wrote CSV: ", csvname)
  return csvname, outdir
end

if abspath(PROGRAM_FILE) == @__FILE__
  csvname, outdir = run_test()
  # run R plotting script in this directory
  # run plotting from the output folder so the quick_sim CSV is visible
  cd(outdir) do
    run(`Rscript ../../plot_results.R`)
  end
  println("Results (CSV + plots) written to: ", outdir)
end
