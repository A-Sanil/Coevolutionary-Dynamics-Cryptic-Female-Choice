using CSV, DataFrames, Statistics, Dates

infile = "Running_Sims/CSV data april/run_5x1000_20gens_parallel_2026-04-30_120751.csv"
println("READING: ", infile)
df = CSV.read(infile, DataFrame)
cols = [:MeanMale,:MeanFemale,:MeanCount,:BSperm,:MeanMates]
g = groupby(df, :Generation)
gens = sort(unique(df.Generation))
out = DataFrame(Generation=gens)
for c in cols
    m = combine(g, c => mean)
    s = combine(g, c => std)
    out[!, Symbol(string(c)*"_mean")] = m[!,2]
    out[!, Symbol(string(c)*"_sd")] = s[!,2]
end
summaryfile = joinpath("Running_Sims","CSV data april","summary_run_5x1000_20gens_parallel_" * Dates.format(now(),"yyyy-mm-dd_HHMMSS") * ".csv")
CSV.write(summaryfile, out)
println("WROTE_SUMMARY: ", summaryfile)
