#Code to run the model for Kustra and Alonzo "The coevolutionary dynamics of cryptic female choice"
#Modified to include evolving RSC trait (4th trait)
#please send any questions to mkustra@ucsc.edu

#load up package for distributing work on many cpu's
using Distributed
#add the number of procesess i.e. cores being used
# addprocs can be heavy on a laptop. It will only run if the environment
# variable NO_AUTO_ADDPROCS is not set. To skip adding workers for quick
# local runs set NO_AUTO_ADDPROCS=1 in your shell before invoking julia.
if !haskey(ENV, "NO_AUTO_ADDPROCS")
  addprocs(23)
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

# Safely normalize a vector of weights; fall back to uniform if values are invalid
@everywhere function safe_normalize(w)
  if isempty(w)
    return w
  end
  if any(!isfinite, w)
    return fill(1/length(w), length(w))
  end
  total = sum(w)
  if total <= 0
    return fill(1/length(w), length(w))
  end
  return w ./ total
end

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
  prob=exp.((.-(malesT .- d).^2)./(2 .*a)) .* (malesS)
  return(prob./sum(prob))
end

#Probability of fertilization sucess for fair raffle
@everywhere function prob_successFR(malesS)
  return(malesS)./sum(malesS)
end

#mutation function takes in a single allele
@everywhere function mutate(gene)
  @views if wsample(tf,mutats,1)[1]#seeing if mutation happens
    mut=(rand(MTD,1))[1] #if mutation happens draw from mutation distribution
    gene=gene.+mut #add mutational effect
    # Note: RSC trait (trait 4) is allowed to have negative values, so no clamping here
    return(gene)
  else #if no mutation occurs just return same gene value
    return(gene)
  end
end

#Enhanced mutation function for RSC trait
# RSC genotypes can be negative, so mutations are additive on the raw genotypic scale
@everywhere function mutate_rsc(gene)
  @views if wsample(tf,mutats,1)[1]#seeing if mutation happens
    # Use mutation variance appropriate for RSC (higher variance trait)
    # Scale mutation variance relative to RSC trait variance (which is 3x standard)
    rsc_mut_sigma = (4*0.25^2/40)^0.5 * 1.5  # Higher mutation variance for RSC
    mut=(rand(Normal(0, rsc_mut_sigma),1))[1] #if mutation happens draw from mutation distribution
    gene=gene.+mut #add mutational effect (can be negative)
    return(gene)
  else #if no mutation occurs just return same gene value
    return(gene)
  end
end

#Function to make a gamete, i.e.
#sample a single allele per loci randomly from paternal and maternal copies
#pg = paternal genome, mg = maternal genome
#ind = index/ID of dad or mom
@everywhere function make_gamete(pg,mg,ind)
  #initialization of new gamete
  ntraits = size(pg,3)
  gamete=zeros(1,size(pg)[2],ntraits)
  @views for j in 1:size(pg)[2] #iterate through loci
    for k in 1:ntraits #iterate through traits
      if sample(tf) #sample true or false
        gamete[1,j,k]=pg[ind,j,k] #if true take from paternal
      else
        gamete[1,j,k]=mg[ind,j,k] #if false take from maternal
      end
    end
  end
  return gamete #return gamate
end

#simulation function with evolving RSC trait
@everywhere function sim(N,mu,var,a,rsc,tradeoff,generations,d=-1)
  #need to create a deepcopies of all genomes to prevent overwriting.
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
  mgf, pgf, mgm, pgm=start_geno(TraitD, TraitD_RSC, N, 20)
  
  #Preallocating results (22 columns: original 21 + MeanMates)
  dfall=zeros(generations,22)
  
  #Calculating phenotypes for all individuals by adding up paternal and maternal genomes
  #NO CLAMPING - allow natural evolution of all traits
  #first make male phenotype array
  ntraits = size(pgm,3)
  @views mphens=reduce(hcat,[sum(pgm[:,:,i],dims=2).+sum(mgm[:,:,i],dims=2) for i in 1:ntraits])
  
  #now make female phenotype array
  @views fphens=reduce(hcat,[sum(pgf[:,:,i],dims=2).+sum(mgf[:,:,i],dims=2) for i in 1:ntraits])
  
  ####Mating
  #allocate empty vector to keep track of offspring
  offspring=zeros(N)

  #maternal genome females next gen
  mgf2=zeros(N,20,ntraits)
  #paternal genome females next gen
  pgf2=zeros(N,20,ntraits)
  #maternal genome males next gen
  mgm2=zeros(N,20,ntraits)
  #paternal genome males next gen
  pgm2=zeros(N,20,ntraits)

  #Now for loop to simulate until specified generation.
  @inbounds @views for gen in 1:generations
    #note on indexing for the genotypes
    #pgm[row,column,other]
    #pgm[individual,loci,trait]

    #note on indexing for phenotypes
    #first column is female trait
    #second column is male trait
    #third column is sperm number
    #fourth column is RSC
    
    #Don't need to redo phenotypes if it is the first generation
    if gen == 1

    #otherwise recalculate phenoypes and reset offsrping to zero
    #.= reasigns variable without allocating more memory
    else
      #NO CLAMPING - recalculate phenotypes naturally
      mphens.=reduce(hcat,[sum(pgm[:,:,i]+mgm[:,:,i],dims=2) for i in 1:ntraits])
      
      #female phenotypes
      fphens.=reduce(hcat,[sum(pgf[:,:,i]+mgf[:,:,i],dims=2) for i in 1:ntraits])
      
      ####Mating
      offspring.=zeros(N)
    end
    
    #if tradeoff weight probability of precop sucess by both male phenotype and sperm number (eq.1 in text)
    if tradeoff
      precop= mate.((mphens[:,3].*mphens[:,2]),1/1000,2500)
      preprob=safe_normalize(precop)
    #if not a tradeoff weight probability of precop success by sperm number (eq. 2 in text)
    else
      precop = mate.(mphens[:,3],1/20,50)
      preprob=safe_normalize(precop)
    end
    
    #need to standardize traits for selection analysis before sperm depletion
    #standardized female phenotypes
  sF = max(std(fphens[:,1]), 1e-6)
  FMalestnd=(fphens[:,1] .- mean(fphens[:,1]))/sF

    #standardized female traits squared for gamma selection coeffients
    FMalestnd2=0.5 .* FMalestnd .^ 2

    #standardized male traits for selection analysis
  sM = max(std(mphens[:,2]), 1e-6)
  Malestnd=(mphens[:,2] .- mean(mphens[:,2]))/sM

    #standardized male traits squared for gamma selection coeffients
    Malestnd2=0.5 .* Malestnd .^ 2

    #standardized sperm number for selection analysis
  sS = max(std(mphens[:,3]), 1e-6)
  SMalestnd=(mphens[:,3] .- mean(mphens[:,3]))/sS

    #calculate mean sperm number to save for simulation output
    Meansperm=mean(mphens[:,3])

    #calculate mean RSC to save for simulation output (if present)
    # RSC phenotypes are sum of alleles (additive model) with floor, so report mean directly
    if ntraits >= 4
      # Report mean RSC phenotype (additive model, positive values only)
      MeanRSC = mean(mphens[:,4])
    else
      MeanRSC = NaN
    end

    #calculate standard deviation of sperm number to save for model output
  Stdsperm=std(mphens[:,3])

    #standardized sperm number for gamma selection coeffients
    SMalestnd2=0.5 .* SMalestnd .^ 2

    #gamma coeffient for male x female
    gmf=Malestnd.*FMalestnd

    #gamma coeffient for male x sperm number
    gms=Malestnd.*SMalestnd

    #gamma coeffient for female x sperm number
    gfs=FMalestnd.*SMalestnd

    #keep count of female offspring indexing
    fcount=1

    #keep count of male offspring indexing
    mcount=1
    
    #initialize array to track number of mates per female for this generation
    mates_per_female = zeros(size(fphens)[1])

    #next part of code is to loop through all females to mate and reproduce
    for i in 1:size(fphens)[1]
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
        # No clamping - allow RSC to evolve freely, including values < 1 or very large
        mates = rand(Poisson(lambda_mates))
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
      # Never try to sample more mates than are available
      mates = clamp(mates, 0, size(mphens, 1))
      
      # Record number of mates for this female
      mates_per_female[i] = mates
      
      # Skip reproduction if mates == 0 (RSC too low, no mating occurs)
      if mates == 0
        continue
      end
      
      #sample from male population to get males female mates with
      #weighted by precopulatory sucess calculated above
      matesM=wsample(1:size(mphens)[1],preprob,mates,replace=false)

      #if there is only one male no need to model risk of sperm competition
      if mates==1
        #add 2 offspring to index of male
        offspring[matesM[1]]=offspring[matesM[1]]+2
        #model sperm depletion of male after mating
        mphens[matesM,3]=mphens[matesM,3].*exp.(-0.2)
        #loop through number of offspring
        #could be extended beyond 2 which is why it is
        #in a for loop format
        for k in 1:2
          #get dad's index for make_gamete function
         dad=matesM[1]
         if k >1 #make females
           #make egg
           egg=make_gamete(pgf,mgf,i)
           #make sperm
           sperm=make_gamete(pgm,mgm,dad)
           #save paternal genome for female
           pgf2[fcount,:,:]=sperm
           #save maternal genome for female
           mgf2[fcount,:,:]=egg
           #updated female counter
           fcount+=1
         else #make males
           #make egg
           egg=make_gamete(pgf,mgf,i)
           #make sperm
           sperm=make_gamete(pgm,mgm,dad)
           #save paternal genome for male
           pgm2[mcount,:,:]=sperm
           #save maternal genome for male
           mgm2[mcount,:,:]=egg
           #increase male index by 1
           mcount+=1
         end
       end
      else
        #check if running selection analysis not cfc
        if d<0
          #calcualte probs of fertilization sucess for males
          probm=prob_success(mphens[matesM,2],mphens[matesM,3],a,fphens[i,1])
          #fair raffle uncomment line below 
          #probm=prob_successFR(mphens[matesM,3])
          #deplete ejaculation for males
         
        else
          #calculate probs of fertilization sucess for males for cryptic female choice
          probm=prob_success(mphens[matesM,2],mphens[matesM,3],a,d)
          #probm=prob_successFR(mphens[matesM,3])
          #deplete ejaculation for males
          mphens[matesM,3]=mphens[matesM,3].*exp.(-0.2)
        end
        #calculate who got fertilization sucess based on probs
        probm = safe_normalize(probm)
        ferts=wsample(matesM,probm,2)
        #specific female is i
        #specific male is in the vector ferts
          for k in 1:length(ferts)
          dad=ferts[k]
          if k >1 #make females
            #make egg
            egg=make_gamete(pgf,mgf,i)
            #make sperm
            sperm=make_gamete(pgm,mgm,dad)
            #add paternal genome  for female
            pgf2[fcount,:,:]=sperm
            #add maternal genome for female
            mgf2[fcount,:,:]=egg
            #increase index counter for female offspring
            fcount+=1
          else #make males
            #make egg
            egg=make_gamete(pgf,mgf,i)
            #make sperm
            sperm=make_gamete(pgm,mgm,dad)
            #make paternal genome for male
            pgm2[mcount,:,:]=sperm
            #make maternal genome for male
            mgm2[mcount,:,:]=egg
            #increase index counter for male offspring
            mcount+=1
          end
        end
        #next add offspring numbers to offpsring counter
        for c in ferts
          offspring[c]=offspring[c]+1
        end
      end
    end
    #calculate opportunity for selection
    is=(sum((offspring .- mean(offspring)) .^ 2 ) ./ length(offspring)) .* (1 ./ mean(offspring) .^ 2)
    reloff=offspring./mean(offspring)
    #make data frame to calculate selection coeffients
    dfO=DataFrame(RelFit=reloff,Male=Malestnd,Maleq=Malestnd2,FMale=FMalestnd,FMaleq=FMalestnd2,SMale=SMalestnd,SMaleq=SMalestnd2,MFq=gmf,MSq=gms,FSq=gfs)
    #calculate selection coeffients
    model=lm(@formula(RelFit ~ Male+ Maleq+FMale+FMaleq+SMale+SMaleq+MFq+MSq+FSq),dfO)
    
    #calculate mean mates per female
    MeanMates = mean(mates_per_female)

    #put all model results together
    #mean male,mean female, std male, std female,cor,sperm count, sperm count std,is,int,beta,gamma,A,MeanRSC,a,gen,MeanMates
    sumdf=[mean(mphens[:,2]),mean(fphens[:,1]),std(mphens[:,2]),std(fphens[:,1]),cor(mphens[:,2],fphens[:,1]),Meansperm,Stdsperm,is,coef(model)[1],coef(model)[2],coef(model)[3],coef(model)[4],coef(model)[5],coef(model)[6],coef(model)[7],coef(model)[8],coef(model)[9],coef(model)[10],a,MeanRSC,gen,MeanMates]
    dfall[gen,:]=sumdf
    #next generation

    # Apply mutations with specialized function for RSC trait
    ntraits = size(pgm,3)
    for i in 1:size(pgm,1)
      for j in 1:size(pgm,2)
        for k in 1:ntraits
          if k == 4  # RSC trait
            pgm[i,j,k] = mutate_rsc(pgm2[i,j,k])
            mgm[i,j,k] = mutate_rsc(mgm2[i,j,k])
            pgf[i,j,k] = mutate_rsc(pgf2[i,j,k])
            mgf[i,j,k] = mutate_rsc(mgf2[i,j,k])
          else  # Other traits (1-3): NO CLAMPING - allow natural evolution
            pgm[i,j,k] = mutate(pgm2[i,j,k])
            mgm[i,j,k] = mutate(mgm2[i,j,k])
            pgf[i,j,k] = mutate(pgf2[i,j,k])
            mgf[i,j,k] = mutate(mgf2[i,j,k])
          end
        end
      end
    end
  end
  return(dfall)
end

#function or run simulation so I can put it in a for loop below
@everywhere function runsim(reps,N,mu,var,a,rsc,tradeoff,gens,d=-1)
  resultsP=SharedArray{Float64}(reps*gens,23)
  @sync @distributed for i in 1:reps
    @async resultsP[(1+(i-1)*gens):(gens*i),1:22]=sim(N,mu,var,a,rsc,tradeoff,gens,d)
    @async resultsP[(1+(i-1)*gens):(gens*i),23]=fill(i,gens)
  end
  return(resultsP)
end

# Single-threaded runner (copied/adapted from noeverywhere.jl)
function runsim_serial(reps,N,mu,var,a,rsc,tradeoff,gens,d=-1)
  resultsP=zeros(Float64, reps*gens, 23)
  for i in 1:reps
    println("Running replicate $i of $reps...")
    resultsP[(1+(i-1)*gens):(gens*i),1:22]=sim(N,mu,var,a,rsc,tradeoff,gens,d)
    resultsP[(1+(i-1)*gens):(gens*i),23]=fill(i,gens)
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
        data = DataFrame(results, [:MeanMale,:MeanFemale,:SDMale,:SDFemale,:cor,:MeanCount,:SDCount,:is,:int,:BMale,:GMale,:BFemale,:GFemale,:BSperm,:GSperm,:GMF,:GMS,:GFS,:a,:MeanRSC,:Generation,:Rep])
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
  data = DataFrame(results, [:MeanMale,:MeanFemale,:SDMale,:SDFemale,:cor,:MeanCount,:SDCount,:is,:int,:BMale,:GMale,:BFemale,:GFemale,:BSperm,:GSperm,:GMF,:GMS,:GFS,:a,:MeanRSC,:Generation,:Rep])
  timestamp = Dates.format(now(), "yyyy-mm-dd_HHMMSS")
  outfile = "test_simulation_results_$(timestamp).csv"
  CSV.write(outfile, data)
  println("Test simulation completed and saved to: ", outfile)
end
