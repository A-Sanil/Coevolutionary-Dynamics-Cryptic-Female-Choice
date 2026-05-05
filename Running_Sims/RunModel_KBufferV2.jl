#Code to run the model for Kustra and Alonzo "The coevolutionary dynamics of cryptic female choice"
#Modified to include evolving RSC trait (4th trait)
#please send any questions to mkustra@ucsc.edu

#load up package for distributing work on many cpu's
using Distributed
using Dates
#add the number of procesess i.e. cores being used
# addprocs can be heavy on a laptop. It will only run if the environment
# variable NO_AUTO_ADDPROCS is not set. To skip adding workers for quick
# local runs set NO_AUTO_ADDPROCS=1 in your shell before invoking julia.
if !haskey(ENV, "NO_AUTO_ADDPROCS")
  cpu_threads = Sys.CPU_THREADS
  default_target = max(1, cpu_threads - 1)
  requested_target = try
    parse(Int, get(ENV, "N_WORKERS", string(default_target)))
  catch
    default_target
  end
  target_workers = clamp(requested_target, 1, max(1, cpu_threads - 1))
  workers_to_add = max(0, target_workers - nworkers())
  if workers_to_add > 0
    addprocs(workers_to_add)
  end
end

#load up packages across all cores
#the @everywhere tag executes the code across all cores
@everywhere using Random, Distributions, StatsBase, GLM, DataFrames,CSV,SharedArrays

#mutation distribution of alleles for 20 Loci runs
@everywhere const MTD=Normal(0,(4*0.25^2/40)^0.5)

#constant for true or false array to sample from
@everywhere const tf=[true,false]

#constant array for whether a mutation occurs
@everywhere const mutats=[0.005,1-0.005]

#probability of mating sucess eq.1 in main text
@everywhere function mate(x,a,b)
  1-(1/(1+exp(-a*(x-b))))
end

#generate population function with 4 traits
# TD = distribution for traits 1-3, TD_RSC = distribution for RSC (trait 4)
@everywhere function start_geno(TD, TD_RSC, n, l)
  #####initialization of simulations
  #make maternal genome females - 4 traits
  # First 3 traits use TD, RSC (trait 4) uses TD_RSC (can be negative)
  mgf=cat(rand(TD,(n,l)),rand(TD,(n,l)),rand(TD,(n,l)),rand(TD_RSC,(n,l)),dims=3)

  #NO CLAMPING - allow negative genotypes to evolve naturally

  #make paternal genome females - 4 traits
  pgf=cat(rand(TD,(n,l)),rand(TD,(n,l)),rand(TD,(n,l)),rand(TD_RSC,(n,l)),dims=3)

  #NO CLAMPING - allow negative genotypes to evolve naturally

  #make maternal genome males - 4 traits
  mgm=cat(rand(TD,(n,l)),rand(TD,(n,l)),rand(TD,(n,l)),rand(TD_RSC,(n,l)),dims=3)

  #NO CLAMPING - allow negative genotypes to evolve naturally

  #make paternal genome males - 4 traits
  pgm=cat(rand(TD,(n,l)),rand(TD,(n,l)),rand(TD,(n,l)),rand(TD_RSC,(n,l)),dims=3)

  #NO CLAMPING - allow negative genotypes to evolve naturally

  return(mgf,pgf,mgm,pgm)
end

#generate populations with different starting trait averages
@everywhere function start_genoD(TD,TD2,n,l)
  #####initialization of simulations
  #make maternal genome females - 4 traits
  mgf=cat(rand(TD2,(n,l)),rand(TD,(n,l)),rand(TD,(n,l)),rand(TD,(n,l)),dims=3)

  #convert negative genotypic values to 0
  mgf[mgf.< 0] .= 0

  #make paternal genome females - 4 traits
  pgf=cat(rand(TD2,(n,l)),rand(TD,(n,l)),rand(TD,(n,l)),rand(TD,(n,l)),dims=3)

  #convert negative genotypic values to 0
  pgf[pgf.< 0] .= 0

  #make maternal genome males - 4 traits
  mgm=cat(rand(TD2,(n,l)),rand(TD,(n,l)),rand(TD,(n,l)),rand(TD,(n,l)),dims=3)

  #convert negative genotypic values to 0
  mgm[mgm.< 0] .= 0

  #make paternal genome males - 4 traits
  pgm=cat(rand(TD2,(n,l)),rand(TD,(n,l)),rand(TD,(n,l)),rand(TD,(n,l)),dims=3)

  #convert negative genotypic values to 0
  pgm[pgm.< 0] .= 0

  return(mgf,pgf,mgm,pgm)
end

#Probability of fertilization sucess function (eq.9 in text)
# '.' makes function vectorized
@everywhere function prob_success(malesT,malesS,a,d)
  prob = exp.((.-(malesT .- d).^2)./(2 .* a)) .* (malesS)
  # guard against all-zero or non-finite weights
  # clamp negative values from numeric drift or negative trait values before sampling
  prob[prob .< 0] .= 0.0
  # replace non-finite entries with 0
  prob[.!isfinite.(prob)] .= 0.0
  s = sum(prob)
  if s == 0.0 || !isfinite(s)
    # fallback to uniform probabilities across males provided
    prob .= 1.0
    s = sum(prob)
  end
  return prob ./ s
end

#Probability of fertilization sucess for fair raffle
@everywhere function prob_successFR(malesS)
  prob = copy(malesS)
  prob[prob .< 0] .= 0.0
  s = sum(prob)
  if s == 0.0 || !isfinite(s)
    prob .= 1.0
    s = sum(prob)
  end
  return prob ./ s
end

# Helper: determine per-female offspring based on mate count and chosen function
# Arguments:
#  - mates: number of males female mated with (0..max_mates)
#  - max_offspring: maximum offspring value to enforce upper bound
#  - offspring_function: symbol controlling formula
#      :poisson -> current RSC-based Poisson rate
#      :logistic -> logistic saturation curve
#      :expdecay -> exponential decay with mates (diminishing returns)
#      :gaussian -> peaked optimum at mu mates
#  - r, c, mu, sigma: parameters controlling curve shape
#  - offspring_scale: scaling used by :poisson mode
#  - female_rsc: female RSC phenotype for :poisson mode
@everywhere function offspring_from_mates(mates::Int;
                                           max_offspring::Int=5,
                                           offspring_function::Symbol=:poisson,
                                           r::Float64=1.0,
                                           c::Float64=2.0,
                                           mu::Float64=2.0,
                                           sigma::Float64=1.0,
                                           offspring_scale::Float64=1.0,
                                           female_rsc::Float64=0.0)
  if offspring_function == :poisson
    # Existing dynamic behavior using RSC and offspring scale
    lambda_offspring = max(0.0, offspring_scale * max(0.0, female_rsc) * (1.0 + 0.25*mates))
    n = rand(Poisson(lambda_offspring))
  elseif offspring_function == :logistic
    # max_offspring/(1+e^{-r*(mates-c)})
    n = round(Int, max_offspring / (1 + exp(-r * (mates - c))))
  elseif offspring_function == :expdecay
    # max_offspring * e^{-r*mates}
    n = round(Int, max_offspring * exp(-r * mates))
  elseif offspring_function == :gaussian
    # max_offspring * exp(-(mates - mu)^2 / (2*sigma^2))
    n = round(Int, max_offspring * exp(-((mates - mu)^2) / (2 * sigma^2)))
  else
    error("Unsupported offspring_function: $offspring_function")
  end
  return clamp(n, 0, max_offspring)
end

#mutation function takes in a single allele
@everywhere function mutate(gene)
  if rand() < 0.005 # fast probability check, 0.5% chance
    gene += rand(MTD) # draw scalar from mutation distribution
    return gene < 0.0 ? 0.0 : gene # clamp to 0 if negative
  end
  return gene
end

#Enhanced mutation function for RSC trait
# RSC genotypes can be negative, so mutations are additive on the raw genotypic scale
@everywhere function mutate_rsc(gene)
  if rand() < 0.005 # fast probability check, 0.5% chance
    # Use mutation variance appropriate for RSC (higher variance trait)
    # Scale mutation variance relative to RSC trait variance (which is 3x standard)
    rsc_mut_sigma = (4*0.25^2/40)^0.5 * 1.5  # Higher mutation variance for RSC
    gene += rand(Normal(0, rsc_mut_sigma)) # add mutational effect (can be negative)
  end
  return gene
end

# Optimized in-place gamete generation
# writes directly to the pre-allocated target array to avoid memory allocations
@everywhere function make_gamete!(target, target_row, pg, mg, ind)
  nloci = size(pg, 2)
  ntraits = size(pg, 3)
  @inbounds for j in 1:nloci
    for k in 1:ntraits
      if rand(Bool) # fast optimized 50/50 boolean random
        target[target_row, j, k] = pg[ind, j, k]
      else
        target[target_row, j, k] = mg[ind, j, k]
      end
    end
  end
end

#simulation function with evolving RSC trait
# V2 adds optional dynamic offspring mode and true offspring totals.
@everywhere function sim(N,mu,var,a,rsc,tradeoff,generations,d=-1,K=N,maintain_sex_ratio=true,show_gui=false,dynamic_offspring_mode=false,offspring_mode=:poisson,offspring_scale=1.0,max_mates=5,max_offspring=5,offspring_r=1.0,offspring_c=2.0,offspring_mu=2.0,offspring_sigma=1.0,progress_interval=10,checkpoint_interval=0,repid=0)
  # sim() runs a full cohort of generation dynamics in K-buffer model V2.
  # Parameters:
  #  - N: initial male/female count each, so initial pop = 2*N
  #  - K: carrying capacity (default matches N), P = K total after cull
  #  - dynamic_offspring_mode: true to use offspring_from_mates formulas
  #  - offspring_mode: :poisson/:logistic/:expdecay/:gaussian
  #  - max_mates: mate cap per female (plate range 0..max_mates)
  #  - max_offspring: cap on target offspring per female
  #  - offspring_r, c, mu, sigma: shape parameters for logistic/gaussian/expdecay
  #  - rsc: used to set Poisson mate means and indirectly offspring in :poisson mode
  # return: dfall matrix (generations × 29 columns plus Rep appended later)
  # need to create deepcopies of all genomes to prevent overwriting.
  # Distribution for traits 1-3 (standard traits)
  TraitD=Normal(mu,var)
  
  # RSC initialized using Poisson distribution
  # Lambda per locus chosen so that sum across 20 loci has reasonable mean (~1-3)
  # Using lambda = 0.1 per locus gives mean ≈ 2 across 20 loci
  lambda_RSC_per_locus = 0.1
  TraitD_RSC = Poisson(lambda_RSC_per_locus)
  
  #Initialize population with 4 traits (female, male, sperm, RSC)
  # RSC genotypes initialized from Poisson distribution (sum of Poissons is Poisson)
  # After mutations (which are continuous), RSC can evolve continuously
  mgf_init, pgf_init, mgm_init, pgm_init = start_geno(TraitD, TraitD_RSC, N, 20)
  
  #Preallocating results
  # 29 columns = prior 26 + offspring_desired + offspring_realized + offspring_dropped
  dfall=zeros(generations,29)

  # Carrying capacity K should be at least N
  K = max(1, max(N, K))

  # offspring staging capacity (all cat-like arrays sized to at least 2*K and max offspring potential)
  bufsize = max(Int(2*K), max_offspring * N, max_mates * N)

  # fixed-capacity population arrays; active population can vary each generation
  ntraits = size(pgm_init,3)
  mgf = zeros(bufsize,20,ntraits)
  pgf = zeros(bufsize,20,ntraits)
  mgm = zeros(bufsize,20,ntraits)
  pgm = zeros(bufsize,20,ntraits)
  mphens = zeros(bufsize,ntraits)
  fphens = zeros(bufsize,ntraits)
  pgm_adults = zeros(bufsize,20,ntraits)
  mgm_adults = zeros(bufsize,20,ntraits)
  pgf_adults = zeros(bufsize,20,ntraits)
  mgf_adults = zeros(bufsize,20,ntraits)
  mgf[1:N,:,:] = mgf_init  
  pgf[1:N,:,:] = pgf_init
  mgm[1:N,:,:] = mgm_init
  pgm[1:N,:,:] = pgm_init

  # active counts (these can scale with offspring production)
  Nf_curr = N
  Nm_curr = N
  
  #maternal genome females next gen
  mgf2=zeros(bufsize,20,ntraits)
  #paternal genome females next gen
  pgf2=zeros(bufsize,20,ntraits)
  #maternal genome males next gen
  mgm2=zeros(bufsize,20,ntraits)
  #paternal genome males next gen
  pgm2=zeros(bufsize,20,ntraits)

  #Now for loop to simulate until specified generation.
  @inbounds @views for gen in 1:generations
    # Progress display
    if show_gui
      # Console progress bar updates every generation
      prog = Int(round(gen/generations * 50))
      bar = "█"^prog * "░"^(50-prog)
      percent = Int(round(gen/generations * 100))
      if myid() == 1
        print("\r[$bar] $percent% | Gen: $gen/$generations | Pop: $(Nm_curr + Nf_curr) (♂$(Nm_curr) ♀$(Nf_curr))  ")
        flush(stdout)
      end
    else
        # Silent mode: periodic updates controlled by progress_interval
        if gen % progress_interval == 0 || gen == 1
          if myid() == 1
            println("  Generation $gen/$generations (Pop: $(Nm_curr + Nf_curr))")
          end
        end
        # Checkpoint: write partial CSV for this replicate if requested
        if checkpoint_interval > 0 && (gen % checkpoint_interval == 0 || gen == generations)
          try
            cols = [:MeanMale, :MeanFemale, :SDMale, :SDFemale, :cor,
                    :MeanCount, :SDCount, :is, :int,
                    :BMale, :GMale, :BFemale, :GFemale, :BSperm, :GSperm,
                    :GMF, :GMS, :GFS, :a, :MeanRSC, :Generation, :MeanMates,
                    :population, :males, :females, :offspring,
                    :offspring_desired, :offspring_realized, :offspring_dropped]
            nrows = gen
            partial = DataFrame(dfall[1:nrows,1:29], cols)
            outroot = get(ENV, "BRC_OUTPUT_ROOT", "CSV data")
            runstamp = Dates.format(now(), "yyyy-mm-ddTHH-MM-SS")
            fname = joinpath(outroot, string("kbuffer_v2_rep", repid, "_progress_", runstamp, ".csv"))
            isdir(outroot) || mkpath(outroot)
            CSV.write(fname, partial)
          catch e
            @warn "checkpoint write failed" error=(e, catch_backtrace())
          end
        end
    end
    
    #note on indexing for the genotypes
    #pgm[row,column,other]
    #pgm[individual,loci,trait]

    #note on indexing for phenotypes
    #first column is female trait
    #second column is male trait
    #third column is sperm number
    #fourth column is RSC
    
  
    # recalculate phenotypes on active rows only using allocation-free loops
    @inbounds for i in 1:Nm_curr
      for t in 1:ntraits
        s = 0.0
        for l in 1:20
          s += pgm[i, l, t] + mgm[i, l, t]
        end
        mphens[i, t] = s
      end
    end
    @inbounds for i in 1:Nf_curr
      for t in 1:ntraits
        s = 0.0
        for l in 1:20
          s += pgf[i, l, t] + mgf[i, l, t]
        end
        fphens[i, t] = s
      end
    end

    ####Mating
    offspring = zeros(Nm_curr)
    
    #if tradeoff weight probability of precop sucess by both male phenotype and sperm number (eq.1 in text)
    if tradeoff
      precop= mate.((mphens[1:Nm_curr,3].*mphens[1:Nm_curr,2]),1/1000,2500)
      preprob=precop./sum(precop)
    #if not a tradeoff weight probability of precop success by sperm number (eq. 2 in text)
    else
      precop = mate.(mphens[1:Nm_curr,3],1/20,50)
      preprob=precop./sum(precop)
    end
    
    #need to standardize traits for selection analysis before sperm depletion
    #standardized female phenotypes
  sF = max(std(fphens[1:Nf_curr,1]), 1e-6)
  FMalestnd=(fphens[1:Nf_curr,1] .- mean(fphens[1:Nf_curr,1]))/sF

    #standardized female traits squared for gamma selection coeffients
    FMalestnd2=0.5 .* FMalestnd .^ 2

    #standardized male traits for selection analysis
  sM = max(std(mphens[1:Nm_curr,2]), 1e-6)
  Malestnd=(mphens[1:Nm_curr,2] .- mean(mphens[1:Nm_curr,2]))/sM

    #standardized male traits squared for gamma selection coeffients
    Malestnd2=0.5 .* Malestnd .^ 2

    #standardized sperm number for selection analysis
  sS = max(std(mphens[1:Nm_curr,3]), 1e-6)
  SMalestnd=(mphens[1:Nm_curr,3] .- mean(mphens[1:Nm_curr,3]))/sS

    #calculate mean sperm number to save for simulation output
    Meansperm=mean(mphens[1:Nm_curr,3])

    #calculate mean RSC to save for simulation output (if present)
    # RSC phenotypes are sum of alleles (additive model) with floor, so report mean directly
    if ntraits >= 4
      # Report mean RSC phenotype (additive model, positive values only)
      MeanRSC = mean(mphens[1:Nm_curr,4])
    else
      MeanRSC = NaN
    end

    #calculate standard deviation of sperm number to save for model output
  Stdsperm=std(mphens[1:Nm_curr,3])

    #standardized sperm number for gamma selection coeffients
    SMalestnd2=0.5 .* SMalestnd .^ 2

    # align lengths for selection model when active male/female counts differ
    nsel = min(length(Malestnd), length(FMalestnd), length(SMalestnd), length(offspring))
    if nsel == 0
      continue
    end

    Malestnd_sel = Malestnd[1:nsel]
    Malestnd2_sel = Malestnd2[1:nsel]
    FMalestnd_sel = FMalestnd[1:nsel]
    FMalestnd2_sel = FMalestnd2[1:nsel]
    SMalestnd_sel = SMalestnd[1:nsel]
    SMalestnd2_sel = SMalestnd2[1:nsel]
    offspring_sel = offspring[1:nsel]

    #gamma coeffient for male x female
    gmf=Malestnd_sel.*FMalestnd_sel

    #gamma coeffient for male x sperm number
    gms=Malestnd_sel.*SMalestnd_sel

    #gamma coeffient for female x sperm number
    gfs=FMalestnd_sel.*SMalestnd_sel

    #keep count of female offspring indexing
    fcount=1

    #keep count of male offspring indexing
    mcount=1
    
    #initialize array to track number of mates per female for this generation
    mates_per_female = zeros(Nf_curr)
    offspring_desired_total = 0

    #next part of code is to loop through all females to mate and reproduce
    for i in 1:Nf_curr
      #mates = number of males a female mates with
      #Use evolving RSC trait to determine number of mates
      #RSC phenotype directly serves as lambda (mean) for Poisson distribution
      if ntraits >= 4
        # Get female RSC phenotype (can be negative - natural evolution)
        female_rsc_phenotype = fphens[i,4]
        
        # Use RSC as lambda for Poisson - only clamp at point of use (Poisson requires non-negative)
        # This allows RSC to evolve freely, but ensures Poisson sampling works
        lambda_mates = max(0.0, female_rsc_phenotype)
        
        # Sample number of mates from Poisson distribution
        mates = rand(Poisson(lambda_mates))
        # enforce maximum mates (plate range 0..max_mates)
        mates = clamp(mates, 0, max_mates)
      else
        # Fallback to old static rsc-based sampling (shouldn't reach here)
        if rsc<=1
          mates=wsample([1,2],[(1-rsc),rsc],1)[1]
        elseif rsc<=2
          mates=wsample([2,3],[(2-rsc),rsc-1],1)[1]
        else
          mates=wsample([3,4],[(3-rsc),rsc-2],1)[1]
        end
      end
      
      # Record number of mates for this female
      mates_per_female[i] = mates
      
      # Skip reproduction if mates == 0 (RSC too low, no mating occurs)
      if mates == 0
        continue
      end
      
      #sample from male population to get males female mates with
      #weighted by precopulatory sucess calculated above
      sample_n = min(mates, Nm_curr)
      matesM=wsample(1:Nm_curr,preprob,sample_n,replace=false)

      # determine target offspring count for this female
      if dynamic_offspring_mode
        target_offspring = offspring_from_mates(mates;
                                                max_offspring=max_offspring,
                                                offspring_function=offspring_mode,
                                                r=offspring_r,
                                                c=offspring_c,
                                                mu=offspring_mu,
                                                sigma=offspring_sigma,
                                                offspring_scale=offspring_scale,
                                                female_rsc=fphens[i,4])
      else
        target_offspring = 2
      end
      offspring_desired_total += target_offspring

      # enforce available staging capacity before writing
      slots_left = max(0, bufsize - ((mcount - 1) + (fcount - 1)))
      realized_offspring = min(target_offspring, slots_left)

      if realized_offspring == 0
        continue
      end

      # LEVEL 1 (single-mate path): no sperm competition, one sire gets all realized offspring
      if mates==1
        dad = matesM[1]
        offspring[dad] += realized_offspring
        # Post-mating sperm depletion cost is still applied to the single sire
        mphens[matesM,3]=mphens[matesM,3].*exp.(-0.2)
        for _ in 1:realized_offspring
          if rand(Bool)
            make_gamete!(pgm2, mcount, pgm, mgm, dad) # sperm to paternal
            make_gamete!(mgm2, mcount, pgf, mgf, i)   # egg to maternal
            mcount+=1
          else
            make_gamete!(pgf2, fcount, pgm, mgm, dad) # sperm to paternal
            make_gamete!(mgf2, fcount, pgf, mgf, i)   # egg to maternal
            fcount+=1
          end
        end
      else
        # LEVEL 2 (multi-mate path): post-copulatory sperm competition among candidate sires
        # d < 0  -> non-CFC: female phenotype sets optimum in prob_success
        # d >= 0 -> CFC: fixed optimum d sets fertilization bias
        if d<0
          # Calculate fertilization weights for mates in non-CFC mode
          probm=prob_success(mphens[matesM,2],mphens[matesM,3],a,fphens[i,1])
          # Optional fair raffle baseline (sperm count only):
          # probm=prob_successFR(mphens[matesM,3])
        else
          # Calculate fertilization weights for mates in cryptic female choice mode
          probm=prob_success(mphens[matesM,2],mphens[matesM,3],a,d)
        end

        # Apply sperm depletion to all males that mated with this female.
        # This mirrors original RunModel behavior and ensures ejaculate cost is paid
        # in both non-CFC and CFC competition modes.
        mphens[matesM,3]=mphens[matesM,3].*exp.(-0.2)

        # LEVEL 3: draw sires for each realized offspring from sperm-competition weights
        ferts=wsample(matesM,probm,realized_offspring,replace=true)
        for dad in ferts
          offspring[dad]=offspring[dad]+1
          if rand(Bool)
            make_gamete!(pgm2, mcount, pgm, mgm, dad) # sperm to paternal
            make_gamete!(mgm2, mcount, pgf, mgf, i)   # egg to maternal
            mcount+=1
          else
            make_gamete!(pgf2, fcount, pgm, mgm, dad) # sperm to paternal
            make_gamete!(mgf2, fcount, pgf, mgf, i)   # egg to maternal
            fcount+=1
          end
        end
      end
    end
    #calculate opportunity for selection
    is=(sum((offspring_sel .- mean(offspring_sel)) .^ 2 ) ./ length(offspring_sel)) .* (1 ./ mean(offspring_sel) .^ 2)
    reloff=offspring_sel./mean(offspring_sel)
    #make data frame to calculate selection coeffients
    dfO=DataFrame(RelFit=reloff,Male=Malestnd_sel,Maleq=Malestnd2_sel,FMale=FMalestnd_sel,FMaleq=FMalestnd2_sel,SMale=SMalestnd_sel,SMaleq=SMalestnd2_sel,MFq=gmf,MSq=gms,FSq=gfs)
    #calculate selection coeffients
    model=lm(@formula(RelFit ~ Male+ Maleq+FMale+FMaleq+SMale+SMaleq+MFq+MSq+FSq),dfO)
    
    #calculate mean mates per female
    MeanMates = mean(mates_per_female)

    # current population counts
    females_now = Nf_curr
    males_now = Nm_curr
    population_now = females_now + males_now

    #put all model results together
    #mean male,mean female, std male, std female,cor,sperm count, sperm count std,is,int,beta,gamma,A,MeanRSC,a,gen,MeanMates,population,males,females,offspring
    # offspring will be calculated after generation transition
    ncor = min(Nm_curr, Nf_curr)
    cor_mf = ncor > 1 ? cor(mphens[1:ncor,2], fphens[1:ncor,1]) : NaN
    sumdf=[mean(mphens[1:Nm_curr,2]),mean(fphens[1:Nf_curr,1]),std(mphens[1:Nm_curr,2]),std(fphens[1:Nf_curr,1]),cor_mf,Meansperm,Stdsperm,is,coef(model)[1],coef(model)[2],coef(model)[3],coef(model)[4],coef(model)[5],coef(model)[6],coef(model)[7],coef(model)[8],coef(model)[9],coef(model)[10],a,MeanRSC,gen,MeanMates,population_now,males_now,females_now,0.0,0.0,0.0,0.0]
    # produced offspring counts this generation
    produced_females = fcount - 1
    produced_males = mcount - 1
    total_offspring = produced_females + produced_males
    
    # true offspring accounting
    desired_total = offspring_desired_total
    realized_total = total_offspring
    dropped_total = max(0, desired_total - realized_total)

    # store offspring count in sumdf
    sumdf[26] = total_offspring
    sumdf[27] = desired_total
    sumdf[28] = realized_total
    sumdf[29] = dropped_total
    dfall[gen,:]=sumdf
    #next generation

    # keep existing adults (including non-mating individuals) in the next generation
    # while staying within fixed 2*K capacity
    male_survivors = min(Nm_curr, max(0, bufsize - produced_males))
    female_survivors = min(Nf_curr, max(0, bufsize - produced_females))
    Nm_next = produced_males + male_survivors
    Nf_next = produced_females + female_survivors

    # snapshot current adults before writing next generation
    pgm_adults[1:Nm_curr, :, :] .= pgm[1:Nm_curr, :, :]
    mgm_adults[1:Nm_curr, :, :] .= mgm[1:Nm_curr, :, :]
    pgf_adults[1:Nf_curr, :, :] .= pgf[1:Nf_curr, :, :]
    mgf_adults[1:Nf_curr, :, :] .= mgf[1:Nf_curr, :, :]

    # Apply proper mutations to newly produced offspring using fast vectorized broadcasting.
    # Traits 1-3 use standard mutation (clamped at 0 to prevent negative phenotypic values).
    # Trait 4 (RSC) uses mutate_rsc (unclamped, allowing negative evolution).
    ntraits = size(pgm,3)
    if produced_males > 0
      @views pgm[1:produced_males, :, 1:3] .= mutate.(pgm2[1:produced_males, :, 1:3])
      @views mgm[1:produced_males, :, 1:3] .= mutate.(mgm2[1:produced_males, :, 1:3])
      if ntraits >= 4
        @views pgm[1:produced_males, :, 4:end] .= mutate_rsc.(pgm2[1:produced_males, :, 4:end])
        @views mgm[1:produced_males, :, 4:end] .= mutate_rsc.(mgm2[1:produced_males, :, 4:end])
      end
    end

    if produced_females > 0
      @views pgf[1:produced_females, :, 1:3] .= mutate.(pgf2[1:produced_females, :, 1:3])
      @views mgf[1:produced_females, :, 1:3] .= mutate.(mgf2[1:produced_females, :, 1:3])
      if ntraits >= 4
        @views pgf[1:produced_females, :, 4:end] .= mutate_rsc.(pgf2[1:produced_females, :, 4:end])
        @views mgf[1:produced_females, :, 4:end] .= mutate_rsc.(mgf2[1:produced_females, :, 4:end])
      end
    end

    # Append surviving adults unchanged after offspring block
    if male_survivors > 0
      pgm[(produced_males+1):Nm_next, :, :] .= pgm_adults[1:male_survivors, :, :]
      mgm[(produced_males+1):Nm_next, :, :] .= mgm_adults[1:male_survivors, :, :]
    end
    if female_survivors > 0
      pgf[(produced_females+1):Nf_next, :, :] .= pgf_adults[1:female_survivors, :, :]
      mgf[(produced_females+1):Nf_next, :, :] .= mgf_adults[1:female_survivors, :, :]
    end

    # Apply carrying capacity constraint: if population exceeds K, sample down to K
    total_pop = Nm_next + Nf_next
    if total_pop > K
      K_int = Int(floor(K))
      
      if maintain_sex_ratio
        # Maintain 50/50 sex ratio: sample K/2 males and K/2 females separately
        males_to_keep = Int(floor(K_int / 2))
        females_to_keep = K_int - males_to_keep
        
        # Sample males
        if males_to_keep < Nm_next
          keep_male_idx = sample(1:Nm_next, males_to_keep, replace=false)
          pgm[1:males_to_keep, :, :] .= pgm[keep_male_idx, :, :]
          mgm[1:males_to_keep, :, :] .= mgm[keep_male_idx, :, :]
          Nm_next = males_to_keep
        end
        
        # Sample females
        if females_to_keep < Nf_next
          keep_female_idx = sample(1:Nf_next, females_to_keep, replace=false)
          pgf[1:females_to_keep, :, :] .= pgf[keep_female_idx, :, :]
          mgf[1:females_to_keep, :, :] .= mgf[keep_female_idx, :, :]
          Nf_next = females_to_keep
        end
      else
        # Random sampling from combined pool (sex ratio can drift)
        # Create index mappings: males are indices 1:Nm_next, females are Nm_next+1:total_pop
        keep_indices = sample(1:total_pop, K_int, replace=false)
        
        # Separate into male and female indices
        male_indices = filter(i -> i <= Nm_next, keep_indices)
        female_indices = filter(i -> i > Nm_next, keep_indices) .- Nm_next
        
        males_to_keep = length(male_indices)
        females_to_keep = length(female_indices)
        
        # Copy selected individuals using .= operator (GC optimized)
        if males_to_keep > 0 && males_to_keep < Nm_next
          pgm[1:males_to_keep, :, :] .= pgm[male_indices, :, :]
          mgm[1:males_to_keep, :, :] .= mgm[male_indices, :, :]
        end
        
        if females_to_keep > 0 && females_to_keep < Nf_next
          pgf[1:females_to_keep, :, :] .= pgf[female_indices, :, :]
          mgf[1:females_to_keep, :, :] .= mgf[female_indices, :, :]
        end
        
        Nm_next = males_to_keep
        Nf_next = females_to_keep
      end
    end

    # update active counts for next generation
    Nm_curr = Nm_next
    Nf_curr = Nf_next
  end
  
  # Clear progress bar line if GUI was shown
  if show_gui
    println()  # New line after progress bar
  end
  
  return(dfall)
end

#function or run simulation so I can put it in a for loop below
@everywhere function runsim(reps,N,mu,var,a,rsc,tradeoff,gens,d=-1,K=N,maintain_sex_ratio=true,show_gui=false,dynamic_offspring_mode=false,offspring_mode=:poisson,offspring_scale=1.0,max_mates=5,max_offspring=5,offspring_r=1.0,offspring_c=2.0,offspring_mu=2.0,offspring_sigma=1.0,progress_interval=10,checkpoint_interval=0)
  resultsP=SharedArray{Float64}(reps*gens,30)
  @sync @distributed for i in 1:reps
    resultsP[(1+(i-1)*gens):(gens*i),1:29]=sim(N,mu,var,a,rsc,tradeoff,gens,d,K,maintain_sex_ratio,show_gui,dynamic_offspring_mode,offspring_mode,offspring_scale,max_mates,max_offspring,offspring_r,offspring_c,offspring_mu,offspring_sigma,progress_interval,checkpoint_interval,i)
    resultsP[(1+(i-1)*gens):(gens*i),30]=fill(i,gens)
  end
  return(resultsP)
end

# Single-threaded runner (copied/adapted from noeverywhere.jl)
function runsim_serial(reps,N,mu,var,a,rsc,tradeoff,gens,d=-1,K=N,maintain_sex_ratio=true,show_gui=false,dynamic_offspring_mode=false,offspring_mode=:poisson,offspring_scale=1.0,max_mates=5,max_offspring=5,offspring_r=1.0,offspring_c=2.0,offspring_mu=2.0,offspring_sigma=1.0,progress_interval=10,checkpoint_interval=0)
  resultsP=zeros(Float64, reps*gens, 30)
  for i in 1:reps
    println("Running replicate $i of $reps...")
    resultsP[(1+(i-1)*gens):(gens*i),1:29]=sim(N,mu,var,a,rsc,tradeoff,gens,d,K,maintain_sex_ratio,show_gui,dynamic_offspring_mode,offspring_mode,offspring_scale,max_mates,max_offspring,offspring_r,offspring_c,offspring_mu,offspring_sigma,progress_interval,checkpoint_interval,i)
    resultsP[(1+(i-1)*gens):(gens*i),30]=fill(i,gens)
  end
  return(resultsP)
end

#@everywhere block

@everywhere gens=30000
@everywhere mu=1.25
@everywhere var=(4*5^2/40)^0.5

# Main loop - rsc parameter is now ignored internally, but kept for backwards compatibility
# Guarded: only runs if RUN_FULL_SIMULATION environment variable is set
if get(ENV, "RUN_FULL_SIMULATION", "0") == "1"
  for j in [true,false]
    for k in [1,12.5,50]
      for l in [0.25,0.5,0.75,1]
        results=runsim(10,500,mu,var,k,l,j,gens)
        outdir = "Results_HighVar_20_1000_RSC"
        isdir(outdir) || mkdir(outdir)
        outfile = joinpath(outdir, string("HV_20_1000_",j,"_",k,"_",l,".csv"))
        data = DataFrame(results, [:MeanMale,:MeanFemale,:SDMale,:SDFemale,:cor,:MeanCount,:SDCount,:is,:int,:BMale,:GMale,:BFemale,:GFemale,:BSperm,:GSperm,:GMF,:GMS,:GFS,:a,:MeanRSC,:Generation,:MeanMates,:population,:males,:females,:offspring,:offspring_desired,:offspring_realized,:offspring_dropped,:Rep])
        CSV.write(outfile, data)
      end
    end
  end
end

# Test runs
# Guarded single-process test run (won't run unless you set environment variable RUN_SINGLE_TEST=1)
if get(ENV, "RUN_SINGLE_TEST", "0") == "1"
  println("Running single-process test simulation (guarded)...")
  # Test parameters (smaller values for testing)
  gens_test = 100
  mu_test = mu
  var_test = var

  # Single test run with one configuration (adjust values as needed)
  tradeoff_test = true
  a_test = 1
  rsc_test = 0.25
  results = runsim_serial(5, 100, mu_test, var_test, a_test, rsc_test, tradeoff_test, gens_test)
  data = DataFrame(results, [:MeanMale,:MeanFemale,:SDMale,:SDFemale,:cor,:MeanCount,:SDCount,:is,:int,:BMale,:GMale,:BFemale,:GFemale,:BSperm,:GSperm,:GMF,:GMS,:GFS,:a,:MeanRSC,:Generation,:MeanMates,:population,:males,:females,:offspring,:offspring_desired,:offspring_realized,:offspring_dropped,:Rep])
  timestamp = Dates.format(now(), "yyyy-mm-dd_HHMMSS")
  outfile = "test_simulation_results_$(timestamp).csv"
  CSV.write(outfile, data)
  println("Test simulation completed and saved to: ", outfile)
end

# Single run entrypoint for K-buffer model (guarded)
if get(ENV, "RUN_ONE_KBUFFER", "0") == "1"
  println("Running one K-buffer simulation...")
  N_one = 1000
  gens_one = 100
  reps_one = 1
  mu_one = mu
  var_one = var
  a_one = 1
  rsc_one = 0.25
  tradeoff_one = true

  results_one = runsim_serial(reps_one, N_one, mu_one, var_one, a_one, rsc_one, tradeoff_one, gens_one, -1, N_one^2, true, false, true, 1.0)
  data_one = DataFrame(results_one, [:MeanMale,:MeanFemale,:SDMale,:SDFemale,:cor,:MeanCount,:SDCount,:is,:int,:BMale,:GMale,:BFemale,:GFemale,:BSperm,:GSperm,:GMF,:GMS,:GFS,:a,:MeanRSC,:Generation,:MeanMates,:population,:males,:females,:offspring,:offspring_desired,:offspring_realized,:offspring_dropped,:Rep])

  runstamp = Dates.format(now(), "yyyy-mm-ddTHH-MM-SS")
  outdir = joinpath("CSV data", "run_" * runstamp)
  isdir(outdir) || mkpath(outdir)
  outfile = joinpath(outdir, "kbuffer_one_run_" * runstamp * ".csv")
  CSV.write(outfile, data_one)
  println("Saved one-run K-buffer CSV to: ", outfile)
end
