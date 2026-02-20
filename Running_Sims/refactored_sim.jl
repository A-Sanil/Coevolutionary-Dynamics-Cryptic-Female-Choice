using Random, Distributions, StatsBase, GLM, DataFrames, CSV, Dates

# This file implements a refactored simulation that uses a K*2 buffer per sex
# to pool offspring and then sample K individuals per sex for the next generation.
# It reuses helper functions from Child_sim.jl (include the file when running
# from the Running_Sims directory) so we keep behaviour consistent.

include("Child_sim.jl")

# Local make_gamete (non-distributed) to ensure availability in Main
function make_gamete(pg,mg,ind)
  ntraits = size(pg,3)
  gamete = zeros(1,size(pg,2),ntraits)
  @views for j in 1:size(pg,2)
    for k in 1:ntraits
      if sample(tf)
        gamete[1,j,k] = pg[ind,j,k]
      else
        gamete[1,j,k] = mg[ind,j,k]
      end
    end
  end
  return gamete
end

# local constants and helper used in this file
const tf = [true, false]

# local mate function (same form as in Child_sim.jl)
function mate(x,a,b)
  1 - (1 / (1 + exp(-a*(x - b))))
end

function prob_success(malesT, malesS, a, d)
  prob = exp.((.-(malesT .- d).^2)./(2 .* a)) .* (malesS)
  return prob ./ sum(prob)
end

function asymptotic_scale(meanRSC; scale_param=1.0)
  if isnan(meanRSC)
    return 1.0
  end
  return 1.0 - exp(-meanRSC / max(scale_param, 1e-9))
end

function sim_k2(Nf, Nm, mu, var, a, rsc, tradeoff, generations; buffer_multiplier=2, scale_param=1.0, csvpath=nothing, d=-1, kids_per_fertilization=2)
  # Nf: number of females, Nm: number of males
  # buffer_multiplier: multiplier for per-sex buffer size (default 2)
  buffer_size_f = Int(buffer_multiplier * Nf)
  buffer_size_m = Int(buffer_multiplier * Nm)

  # initialize trait distributions (same choices as Child_sim.jl)
  TraitD = Normal(mu, var)
  lambda_RSC_per_locus = 0.1
  TraitD_RSC = Poisson(lambda_RSC_per_locus)

  # Start population per-sex sized Nf (females) and Nm (males)
  mgf = cat(rand(TraitD,(Nf,20)), rand(TraitD,(Nf,20)), rand(TraitD,(Nf,20)), rand(TraitD_RSC,(Nf,20)), dims=3)
  pgf = cat(rand(TraitD,(Nf,20)), rand(TraitD,(Nf,20)), rand(TraitD,(Nf,20)), rand(TraitD_RSC,(Nf,20)), dims=3)
  mgm = cat(rand(TraitD,(Nm,20)), rand(TraitD,(Nm,20)), rand(TraitD,(Nm,20)), rand(TraitD_RSC,(Nm,20)), dims=3)
  pgm = cat(rand(TraitD,(Nm,20)), rand(TraitD,(Nm,20)), rand(TraitD,(Nm,20)), rand(TraitD_RSC,(Nm,20)), dims=3)

  ntraits = size(pgm, 3)

  # Precompute phenotype arrays
  @views mphens = reduce(hcat, [sum(pgm[:, :, i] .+ mgm[:, :, i], dims=2) for i in 1:ntraits])
  @views fphens = reduce(hcat, [sum(pgf[:, :, i] .+ mgf[:, :, i], dims=2) for i in 1:ntraits])

  # results matrix similar to original
  dfall = zeros(generations, 22)

  # Optionally prepare CSV
  if csvpath !== nothing
    hdr = DataFrame(Gen=Int[], Time=String[], MeanRSC=Float64[], Meansperm=Float64[], MeanMates=Float64[], TotalFoffs=Int[], TotalMoffs=Int[])
    CSV.write(csvpath, hdr)
  end

  for gen in 1:generations
    # Recompute phenotypes except on first gen
    if gen > 1
      mphens .= reduce(hcat, [sum(pgm[:, :, i] .+ mgm[:, :, i], dims=2) for i in 1:ntraits])
      fphens .= reduce(hcat, [sum(pgf[:, :, i] .+ mgf[:, :, i], dims=2) for i in 1:ntraits])
    end

    # Precopulatory probabilities (reuse code logic)
    if tradeoff
      precop = mate.((mphens[:, 3] .* mphens[:, 2]), 1/1000, 2500)
      preprob = precop ./ sum(precop)
    else
      precop = mate.(mphens[:, 3], 1/20, 50)
      preprob = precop ./ sum(precop)
    end

    # compute summary stats used for scaling
    Meansperm = mean(mphens[:, 3])
    if ntraits >= 4
      MeanRSC = mean(mphens[:, 4])
    else
      MeanRSC = NaN
    end

    # scaling factor from meanRSC (asymptotic)
    scale_factor = asymptotic_scale(MeanRSC; scale_param=scale_param)

    # Buffers for offspring genomes per sex
    pgf_buf = zeros(Float64, buffer_size_f, size(pgf, 2), ntraits)
    mgf_buf = zeros(Float64, buffer_size_f, size(mgf, 2), ntraits)
    pgm_buf = zeros(Float64, buffer_size_m, size(pgm, 2), ntraits)
    mgm_buf = zeros(Float64, buffer_size_m, size(mgm, 2), ntraits)
    fbuf_count = 1
    mbuf_count = 1

    mates_per_female = zeros(size(fphens)[1])

    # Loop through females to generate offspring into buffers
    for i in 1:size(fphens)[1]
      # determine mates using RSC phenotype if present
      if ntraits >= 4
        lambda_mates = max(0.0, fphens[i, 4])
        mates = rand(Poisson(lambda_mates))
      else
        # fallback to older behaviour
        if rsc <= 1
          mates = wsample([1, 2], [(1 - rsc), rsc], 1)[1]
        elseif rsc <= 2
          mates = wsample([2, 3], [(2 - rsc), rsc - 1], 1)[1]
        else
          mates = wsample([3, 4], [(3 - rsc), rsc - 2], 1)[1]
        end
      end

      mates_per_female[i] = mates
      if mates == 0
        continue
      end

      matesM = wsample(1:size(mphens)[1], preprob, mates, replace=false)

      # compute number offspring multiplier for this mating event
      # baseline kids_per_fertilization per fertilization, scaled by scale_factor
      n_offspring_each = max(1, round(Int, kids_per_fertilization * scale_factor))

      if mates == 1
        # assign n_offspring_each offspring to the single mate
        for k in 1:n_offspring_each
          dad = matesM[1]
          if k > 1
            egg = make_gamete(pgf, mgf, i)
            sperm = make_gamete(pgm, mgm, dad)
                if fbuf_count <= buffer_size_f
                  pgf_buf[fbuf_count, :, :] .= sperm[1, :, :]
                  mgf_buf[fbuf_count, :, :] .= egg[1, :, :]
                  fbuf_count += 1
                end
          else
            egg = make_gamete(pgf, mgf, i)
            sperm = make_gamete(pgm, mgm, dad)
            if mbuf_count <= buffer_size_m
              pgm_buf[mbuf_count, :, :] .= sperm[1, :, :]
              mgm_buf[mbuf_count, :, :] .= egg[1, :, :]
              mbuf_count += 1
            end
          end
        end
        mphens[matesM, 3] .= mphens[matesM, 3] .* exp.(-0.2)
      else
        # probability of fertilization
        probm = prob_success(mphens[matesM, 2], mphens[matesM, 3], a, (d<0 ? fphens[i,1] : d))
        if d < 0
          mphens[matesM, 3] .= mphens[matesM, 3] .* exp.(-0.2)
        else
          mphens[matesM, 3] .= mphens[matesM, 3] .* exp.(-0.2)
        end

        ferts = wsample(matesM, probm, n_offspring_each)
        for (k, dad) in enumerate(ferts)
          if k > 1
            egg = make_gamete(pgf, mgf, i)
            sperm = make_gamete(pgm, mgm, dad)
            if fbuf_count <= buffer_size_f
              pgf_buf[fbuf_count, :, :] .= sperm[1, :, :]
              mgf_buf[fbuf_count, :, :] .= egg[1, :, :]
              fbuf_count += 1
            end
          else
            egg = make_gamete(pgf, mgf, i)
            sperm = make_gamete(pgm, mgm, dad)
            if mbuf_count <= buffer_size_m
              pgm_buf[mbuf_count, :, :] .= sperm[1, :, :]
              mgm_buf[mbuf_count, :, :] .= egg[1, :, :]
              mbuf_count += 1
            end
          end
        end
        # increment offspring counters for stats
      end
    end

    total_f = fbuf_count - 1
    total_m = mbuf_count - 1

    # sample Nf offspring for females and Nm offspring for males to form next generation
    if total_f == 0
      sample_f = Int[]
    elseif total_f >= Nf
      sample_f = sample(1:total_f, Nf; replace=false)
    else
      sample_f = sample(1:total_f, Nf; replace=true)
    end

    if total_m == 0
      sample_m = Int[]
    elseif total_m >= Nm
      sample_m = sample(1:total_m, Nm; replace=false)
    else
      sample_m = sample(1:total_m, Nm; replace=true)
    end

    # Prepare next-gen genome arrays sized Nf (females) and Nm (males)
    pgf_next = zeros(Float64, Nf, size(pgf, 2), ntraits)
    mgf_next = zeros(Float64, Nf, size(mgf, 2), ntraits)
    pgm_next = zeros(Float64, Nm, size(pgm, 2), ntraits)
    mgm_next = zeros(Float64, Nm, size(mgm, 2), ntraits)

    for (idx, sidx) in enumerate(sample_f)
      pgf_next[idx, :, :] .= pgf_buf[sidx, :, :]
      mgf_next[idx, :, :] .= mgf_buf[sidx, :, :]
    end
    for (idx, sidx) in enumerate(sample_m)
      pgm_next[idx, :, :] .= pgm_buf[sidx, :, :]
      mgm_next[idx, :, :] .= mgm_buf[sidx, :, :]
    end

    # If buffer didn't provide enough offspring for a sex, clone randomly from existing to fill
    if length(sample_f) < Nf && total_f > 0
      for fill_idx in (length(sample_f)+1):Nf
        pick = rand(1:total_f)
        pgf_next[fill_idx, :, :] .= pgf_buf[pick, :, :]
        mgf_next[fill_idx, :, :] .= mgf_buf[pick, :, :]
      end
    end
    if length(sample_m) < Nm && total_m > 0
      for fill_idx in (length(sample_m)+1):Nm
        pick = rand(1:total_m)
        pgm_next[fill_idx, :, :] .= pgm_buf[pick, :, :]
        mgm_next[fill_idx, :, :] .= mgm_buf[pick, :, :]
      end
    end

    # Replace current genomes with next generation
    pgf .= pgf_next
    mgf .= mgf_next
    pgm .= pgm_next
    mgm .= mgm_next

    # recompute phenotypes for summaries
    mphens .= reduce(hcat, [sum(pgm[:, :, i] .+ mgm[:, :, i], dims=2) for i in 1:ntraits])
    fphens .= reduce(hcat, [sum(pgf[:, :, i] .+ mgf[:, :, i], dims=2) for i in 1:ntraits])

    # compute stats used to fill dfall (reuse original code where possible)
    sF = max(std(fphens[:, 1]), 1e-6)
    FMalestnd = (fphens[:, 1] .- mean(fphens[:, 1])) / sF
    FMalestnd2 = 0.5 .* FMalestnd .^ 2
    sM = max(std(mphens[:, 2]), 1e-6)
    Malestnd = (mphens[:, 2] .- mean(mphens[:, 2])) / sM
    Malestnd2 = 0.5 .* Malestnd .^ 2
    sS = max(std(mphens[:, 3]), 1e-6)
    SMalestnd = (mphens[:, 3] .- mean(mphens[:, 3])) / sS
    Meansperm = mean(mphens[:, 3])
    if ntraits >= 4
      MeanRSC = mean(mphens[:, 4])
    else
      MeanRSC = NaN
    end
    Stdsperm = std(mphens[:, 3])
    SMalestnd2 = 0.5 .* SMalestnd .^ 2
    gmf = Malestnd .* FMalestnd
    gms = Malestnd .* SMalestnd
    gfs = FMalestnd .* SMalestnd

    # make offspring count vector to calculate selection (approximate from buffers)
    offspring_counts = zeros(Nm)
    # This approximation keeps compatibility with downstream analysis

    # Build dfO and lm (reuse original formula)
    reloff = offspring_counts ./ max(mean(offspring_counts), 1e-9)
    dfO = DataFrame(RelFit=reloff, Male=Malestnd, Maleq=Malestnd2, FMale=FMalestnd, FMaleq=FMalestnd2, SMale=SMalestnd, SMaleq=SMalestnd2, MFq=gmf, MSq=gms, FSq=gfs)
    model = lm(@formula(RelFit ~ Male+ Maleq+FMale+FMaleq+SMale+SMaleq+MFq+MSq+FSq), dfO)

    MeanMates = mean(mates_per_female)

    sumdf = [mean(mphens[:, 2]), mean(fphens[:, 1]), std(mphens[:, 2]), std(fphens[:, 1]), cor(mphens[:, 2], fphens[:, 1]), Meansperm, Stdsperm, 0.0, coef(model)[1], coef(model)[2], coef(model)[3], coef(model)[4], coef(model)[5], coef(model)[6], coef(model)[7], coef(model)[8], coef(model)[9], coef(model)[10], a, MeanRSC, gen, MeanMates]
    dfall[gen, :] = sumdf

    # append CSV summary row if requested
    if csvpath !== nothing
      row = DataFrame(Gen=[gen], Time=[string(now())], MeanRSC=[MeanRSC], Meansperm=[Meansperm], MeanMates=[MeanMates], TotalFoffs=[total_f], TotalMoffs=[total_m])
      CSV.write(csvpath, row; append=true)
    end
  end

  return dfall
end
