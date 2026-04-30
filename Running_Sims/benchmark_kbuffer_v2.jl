using Distributed
using Dates
using Statistics
using BenchmarkTools
using DataFrames
using CSV

include("RunModel_KBufferV2.jl")

# -----------------------
# Benchmark configuration
# (overridable through environment variables)
# -----------------------
ns_env = get(ENV, "BENCH_NS", "100,200,400")
Ns = parse.(Int, split(ns_env, ','))
gens = parse(Int, get(ENV, "BENCH_GENS", "20"))
reps = parse(Int, get(ENV, "BENCH_REPS", "2"))
K_factor = parse(Float64, get(ENV, "BENCH_K_FACTOR", "4"))
maintain_sex_ratio = true
show_gui = false
dynamic_offspring_mode = true
offspring_scale = 1.0

mu = 5.0
trait_var = 1.0
a = 1.0
rsc = 0.25
tradeoff = true

samples = parse(Int, get(ENV, "BENCH_SAMPLES", "3"))
evals = parse(Int, get(ENV, "BENCH_EVALS", "1"))
bench_label = get(ENV, "BENCH_LABEL", "default")

function estimate_loglog_slope(xs::Vector{Float64}, ys::Vector{Float64})
    lx = log.(xs)
    ly = log.(ys)
    vx = Statistics.var(lx)
    if vx == 0.0
        return NaN
    end
    return cov(lx, ly) / vx
end

function bench_serial(N)
    K = Int(K_factor * N)
    t = @belapsed runsim_serial($reps, $N, $mu, $trait_var, $a, $rsc, $tradeoff, $gens, -1,
                                $K, $maintain_sex_ratio, $show_gui, $dynamic_offspring_mode, $offspring_scale) samples=samples evals=evals
    return t
end

function bench_parallel(N)
    K = Int(K_factor * N)
    t = @belapsed runsim($reps, $N, $mu, $trait_var, $a, $rsc, $tradeoff, $gens, -1,
                         $K, $maintain_sex_ratio, $show_gui, $dynamic_offspring_mode, $offspring_scale) samples=samples evals=evals
    return t
end

println("="^70)
println("Benchmarking KBuffer V2 (serial vs parallel)")
println("="^70)
println("Config: Ns=$(Ns), gens=$gens, reps=$reps, dynamic_offspring_mode=$dynamic_offspring_mode")
println("Workers available: ", nworkers())
println("="^70)

serial_times = Float64[]
parallel_times = Float64[]

for N in Ns
    println("Running benchmarks for N=$N ...")
    push!(serial_times, bench_serial(N))
    push!(parallel_times, bench_parallel(N))
end

speedups = serial_times ./ parallel_times

slope_serial = estimate_loglog_slope(Float64.(Ns), serial_times)
slope_parallel = estimate_loglog_slope(Float64.(Ns), parallel_times)

bench_df = DataFrame(
    N = Ns,
    workers = fill(nworkers(), length(Ns)),
    serial_seconds = serial_times,
    parallel_seconds = parallel_times,
    speedup = speedups,
)

runstamp = Dates.format(now(), "yyyy-mm-ddTHH-MM-SS")
outdir = joinpath("CSV data", "bench_v2_" * bench_label * "_" * runstamp * "_w" * string(nworkers()))
isdir(outdir) || mkpath(outdir)

csv_file = joinpath(outdir, "benchmark_results.csv")
CSV.write(csv_file, bench_df)

report_file = joinpath(outdir, "benchmark_report.txt")
open(report_file, "w") do io
    println(io, "KBuffer V2 Benchmark Report")
    println(io, "Date: ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))
    println(io, "")
    println(io, "Configuration")
    println(io, "- Ns: ", Ns)
    println(io, "- generations: ", gens)
    println(io, "- replicates: ", reps)
    println(io, "- maintain_sex_ratio: ", maintain_sex_ratio)
    println(io, "- dynamic_offspring_mode: ", dynamic_offspring_mode)
    println(io, "- offspring_scale: ", offspring_scale)
    println(io, "- workers: ", nworkers())
    println(io, "")
    println(io, "Timings (seconds)")
    for i in eachindex(Ns)
        println(io, "- N=", Ns[i], ": serial=", round(serial_times[i], digits=4),
                ", parallel=", round(parallel_times[i], digits=4),
                ", speedup=", round(speedups[i], digits=3), "x")
    end
    println(io, "")
    println(io, "Empirical scaling (log-log slope)")
    println(io, "- serial slope: ", round(slope_serial, digits=3))
    println(io, "- parallel slope: ", round(slope_parallel, digits=3))
    println(io, "")
    println(io, "Interpretation")
    println(io, "- Runtime approximately scales as O(N^p), where p is the fitted slope above.")
    println(io, "- This is empirical complexity for this configuration (not a strict symbolic proof).")
end

println("\nBenchmark complete.")
println("CSV:    ", csv_file)
println("Report: ", report_file)
println("serial slope ~ ", round(slope_serial, digits=3), " ; parallel slope ~ ", round(slope_parallel, digits=3))
