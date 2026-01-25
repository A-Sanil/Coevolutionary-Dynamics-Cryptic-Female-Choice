# Simple offspring-curve experiments for cryptic female choice project
# Three reproductive response curves:
#   1) parabola: offspring peaks at an intermediate compatibility
#   2) asymptote: offspring rises then levels off
#   3) mixed: blend of parabola + asymptote
# Output: CSV with population/trait trajectories across generations and replicates.

using Random, Statistics, Distributions, DataFrames, CSV, Dates

# Compatibility between a male and female trait (higher when traits are similar)
compatibility(m::Float64, f::Float64) = 1 / (1 + abs(m - f))

# Empirical modifiers from literature
const MHC_PENALTY = 0.10      # ~10% lower sperm success for related pairs (guppy study) citeturn2search3
const REL_THRESHOLD = 0.10    # treat |trait difference| < 0.1 as related
const OVARIAN_BOOST = 2.0     # ovarian fluid roughly doubles motile life span in salmonids citeturn2search4
const OVARIAN_THRESHOLD = 0.50 # boost starts once compatibility exceeds this level
const GL_OPT = 0.50           # gestation-length optimum (normalized trait) from sow study where extremes lowered TNB citeturn0search5
const GL_WIDTH = 0.20         # width of tolerance around optimum

# Expected offspring count given compatibility score in [0,1]
function offspring_mean(score::Float64, model::Symbol)
    # All curves are tuned to peak near ~4 offspring so population maxima are comparable
    if model === :parabola
        peak = 4.0           # maximum expected offspring at optimum
        opt = 0.65           # compatibility giving the peak
        width = 0.70         # wider width so offspring stay viable over broader compatibility
        scale = 1 - ((score - opt) / width)^2
        base = peak * max(0.0, scale)
        return max(1.0, base)  # higher floor prevents collapse when compatibility is low
    elseif model === :asymptote
        max_children = 4.0   # horizontal asymptote
        k = 4.2              # steepness controls how fast it approaches max
        return max_children * (1 - exp(-k * score))
    elseif model === :mixed
        # Blend parabola and asymptote; peaks around 4 as well
        return 0.5 * offspring_mean(score, :parabola) + 0.5 * offspring_mean(score, :asymptote)
    else
        error("Unknown model: $model")
    end
end

# One simulation replicate for a chosen offspring curve
function run_replicate(model::Symbol; generations::Int=100, init_pop::Int=200, K::Int=400,
                       mutation_sd::Float64=0.05, rng=MersenneTwister(1))
    # initialize traits (0-1 range) with mild variation
    males = rand(rng, Uniform(0.3, 0.7), init_pop ÷ 2)
    females = rand(rng, Uniform(0.3, 0.7), init_pop - length(males))

    rows = Vector{NamedTuple}(undef, generations)

    for gen in 1:generations
        current_pop = length(males) + length(females)
        if min(length(males), length(females)) == 0
            # population collapsed; record zeros and continue
            rows[gen] = (model=String(model), rep=1, generation=gen, population=current_pop,
                         mean_trait=NaN, mean_male=NaN, mean_female=NaN,
                         mean_offspring_per_pair=0.0, mean_expected_offspring=0.0,
                         mean_score=NaN, carrying_factor=0.0)
            males = Float64[]; females = Float64[]
            continue
        end

        n_pairs = min(length(males), length(females))
        male_idx = randperm(rng, length(males))[1:n_pairs]
        female_idx = randperm(rng, length(females))[1:n_pairs]

        children_m = Float64[]
        children_f = Float64[]

        total_children = 0
        total_expected = 0.0
        total_score = 0.0
        total_scale = 0.0

        for i in 1:n_pairs
            mtrait = males[male_idx[i]]
            ftrait = females[female_idx[i]]
            score = compatibility(mtrait, ftrait)
            total_score += score

            base_mean = offspring_mean(score, model)
            # crowding feedback to keep population stable around K
            carrying = clamp(1 - current_pop / K, 0.2, 1.2)
            total_scale += carrying

            # MHC / relatedness penalty: ~10% lower success when traits are very similar
            related = abs(mtrait - ftrait) < REL_THRESHOLD
            mhc_factor = related ? (1 - MHC_PENALTY) : 1.0

            # Ovarian fluid boost scales up to 2x as compatibility rises past threshold
            ov_scale = clamp((score - OVARIAN_THRESHOLD) / (1 - OVARIAN_THRESHOLD), 0.0, 1.0)
            ovarian_factor = 1 + (OVARIAN_BOOST - 1) * ov_scale

            # Gestation-length stabilizing selection on female trait (peak near 0.5)
            gl_factor = exp(-((ftrait - GL_OPT) / GL_WIDTH)^2)

            λ = max(base_mean * carrying * mhc_factor * ovarian_factor * gl_factor, 0.0)
            total_expected += λ

            nchild = rand(rng, Poisson(λ))
            total_children += nchild

            for _ in 1:nchild
                trait = clamp((mtrait + ftrait) / 2 + randn(rng) * mutation_sd, 0.0, 1.0)
                if rand(rng) < 0.5
                    push!(children_m, trait)
                else
                    push!(children_f, trait)
                end
            end
        end

        # update population for next generation
        males = children_m
        females = children_f

        new_pop = length(males) + length(females)
        mean_trait = new_pop > 0 ? mean(vcat(males, females)) : NaN
        mean_m = isempty(males) ? NaN : mean(males)
        mean_f = isempty(females) ? NaN : mean(females)

        rows[gen] = (model=String(model), rep=1, generation=gen, population=new_pop,
                     mean_trait=mean_trait, mean_male=mean_m, mean_female=mean_f,
                     mean_offspring_per_pair=total_children / n_pairs,
                     mean_expected_offspring=total_expected / n_pairs,
                     mean_score=total_score / n_pairs,
                     carrying_factor=total_scale / n_pairs)
    end

    return DataFrame(rows)
end

function run_models(; generations=100, init_pop=200, K=400, mutation_sd=0.05, reps=5, seed=42)
    models = [:parabola, :asymptote, :mixed]
    all_runs = DataFrame()
    for (m_idx, model) in enumerate(models)
        for r in 1:reps
            rng = MersenneTwister(seed + 1000 * m_idx + r)
            df = run_replicate(model; generations=generations, init_pop=init_pop, K=K,
                               mutation_sd=mutation_sd, rng=rng)
            df.rep .= r
            append!(all_runs, df)
        end
    end
    return all_runs
end

function main()
    generations = parse(Int, get(ENV, "GENERATIONS", "100"))
    init_pop    = parse(Int, get(ENV, "INIT_POP", "200"))
    K           = parse(Int, get(ENV, "K", "400"))
    mutation_sd = parse(Float64, get(ENV, "MUTATION_SD", "0.05"))
    reps        = parse(Int, get(ENV, "REPS", "5"))

    results = run_models(generations=generations, init_pop=init_pop, K=K,
                         mutation_sd=mutation_sd, reps=reps)

    timestamp = Dates.format(now(), "yyyy-mm-dd_HHMMSS")
    outpath = joinpath(@__DIR__, "child_testing_results_$(timestamp).csv")
    CSV.write(outpath, results)
    println("Saved results -> " * outpath)
    println("Columns: ", names(results))
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
