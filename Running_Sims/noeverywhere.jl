#Code to run the model for Kustra and Alonzo "The coevolutionary dynamics of cryptic female choice"
#Modified to include evolving RSC trait (4th trait)
#Single-threaded version (no parallel processing)
#please send any questions to mkustra@ucsc.edu

#load up packages
using Random, Distributions, StatsBase, GLM, DataFrames, CSV

#mutation distribution of alleles for 20 Loci runs
const MTD=Normal(0,(4*0.25^2/40)^0.5)

#constant for true or false array to sample from
const tf=[true,false]

#constant array for whether a mutation occurs
const mutats=[0.005,1-0.005]

#probability of mating sucess eq.1 in main text
function mate(x,a,b)
  1-(1/(1+exp(-a*(x-b))))
end

#generate population function with 4 traits
function start_geno(TD,n,l)
  #####initialization of simulations
  #make maternal genome females - 4 traits
  mgf=cat(rand(TD,(n,l)),rand(TD,(n,l)),rand(TD,(n,l)),rand(TD,(n,l)),dims=3)

  #convert negative genotypic values to 0
  mgf[mgf.< 0] .= 0

  #make paternal genome females - 4 traits
  pgf=cat(rand(TD,(n,l)),rand(TD,(n,l)),rand(TD,(n,l)),rand(TD,(n,l)),dims=3)

  #convert negative genotypic values to 0
  pgf[pgf.< 0] .= 0

  #make maternal genome males - 4 traits
  mgm=cat(rand(TD,(n,l)),rand(TD,(n,l)),rand(TD,(n,l)),rand(TD,(n,l)),dims=3)

  #convert negative genotypic values to 0
  mgm[mgm.< 0] .= 0

  #make paternal genome males - 4 traits
  pgm=cat(rand(TD,(n,l)),rand(TD,(n,l)),rand(TD,(n,l)),rand(TD,(n,l)),dims=3)

  #convert negative genotypic values to 0
  pgm[pgm.< 0] .= 0

  return(mgf,pgf,mgm,pgm)
end

#generate populations with different starting trait averages
function start_genoD(TD,TD2,n,l)
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
function prob_success(malesT,malesS,a,d)
  prob=exp.((.-(malesT .- d).^2)./(2 .*a)) .* (malesS)
  return(prob./sum(prob))
end

#Probability of fertilization sucess for fair raffle
function prob_successFR(malesS)
  return(malesS)./sum(malesS)
end

#mutation function takes in a single allele
function mutate(gene)
  @views if wsample(tf,mutats,1)[1]#seeing if mutation happens
    mut=(rand(MTD,1))[1] #if mutation happens draw from mutation distribution
    gene=gene.+mut #add mutational effect
    if gene<0 #if gene is now less than 0 make it 0
      return(0)
    else #otherwise return gene
      return(gene)
    end
  else #if no mutation occurs just return same gene value
    return(gene)
  end
end

#Function to make a gamete, i.e.
#sample a single allele per loci randomly from paternal and maternal copies
#pg = paternal genome, mg = maternal genome
#ind = index/ID of dad or mom
function make_gamete(pg,mg,ind)
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
function sim(N,mu,var,a,rsc,tradeoff,generations,d=-1)
  #need to create a deepcopies of all genomes to prevent overwriting.
  TraitD=Normal(mu,var)
  
  #Initialize population with 4 traits (female, male, sperm, RSC)
  mgf, pgf, mgm, pgm=start_geno(TraitD,N,20)
  
  #Initialize RSC trait (trait column 4) with mean=1, small variance
  ntraits = size(pgm,3)
  if ntraits >= 4
    rsc_init_sigma = 0.1  # Small starting variance around mean of 1
    pgf[:,:,4] .= rand(Normal(1, rsc_init_sigma), size(pgf[:,:,4]))
    mgf[:,:,4] .= rand(Normal(1, rsc_init_sigma), size(mgf[:,:,4]))
    pgm[:,:,4] .= rand(Normal(1, rsc_init_sigma), size(pgm[:,:,4]))
    mgm[:,:,4] .= rand(Normal(1, rsc_init_sigma), size(mgm[:,:,4]))
    
    #Enforce minimum RSC of 0.5
    pgf[pgf[:,:,4].<0.5, 4] .= 0.5
    mgf[mgf[:,:,4].<0.5, 4] .= 0.5
    pgm[pgm[:,:,4].<0.5, 4] .= 0.5
    mgm[mgm[:,:,4].<0.5, 4] .= 0.5
  end
  
  #Preallocating results
  dfall=zeros(generations,21)
  
  #Calculating phenotypes for all individuals by adding up paternal and maternal genomes
  #first make male phenotype array
  ntraits = size(pgm,3)
  @views mphens=reduce(hcat,[sum(pgm[:,:,i],dims=2).+sum(mgm[:,:,i],dims=2) for i in 1:ntraits])
  #only clamp the first three trait-columns (female trait, male trait, sperm number) to minimum 1
  for c in 1:min(3,size(mphens,2))
    mphens[mphens[:,c].<1,c] .= 1
  end
  #Enforce RSC minimum separately
  if ntraits >= 4
    mphens[mphens[:,4].<0.5, 4] .= 0.5
  end
  
  #now make female phenotype array
  @views fphens=reduce(hcat,[sum(pgf[:,:,i],dims=2).+sum(mgf[:,:,i],dims=2) for i in 1:ntraits])
  for c in 1:min(3,size(fphens,2))
    fphens[fphens[:,c].<1,c] .= 1
  end
  if ntraits >= 4
    fphens[fphens[:,4].<0.5, 4] .= 0.5
  end
  
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
      mphens.=reduce(hcat,[sum(pgm[:,:,i]+mgm[:,:,i],dims=2) for i in 1:ntraits])
      for c in 1:min(3,size(mphens,2))
        mphens[mphens[:,c].<1,c] .= 1
      end
      if ntraits >= 4
        mphens[mphens[:,4].<0.5, 4] .= 0.5
      end
      
      #female phenotypes
      fphens.=reduce(hcat,[sum(pgf[:,:,i]+mgf[:,:,i],dims=2) for i in 1:ntraits])
      for c in 1:min(3,size(fphens,2))
        fphens[fphens[:,c].<1,c] .= 1
      end
      if ntraits >= 4
        fphens[fphens[:,4].<0.5, 4] .= 0.5
      end
      
      ####Mating
      offspring.=zeros(N)
    end
    
    #if tradeoff weight probability of precop sucess by both male phenotype and sperm number (eq.1 in text)
    if tradeoff
      precop= mate.((mphens[:,3].*mphens[:,2]),1/1000,2500)
      preprob=precop./sum(precop)
    #if not a tradeoff weight probability of precop success by sperm number (eq. 2 in text)
    else
      precop = mate.(mphens[:,3],1/20,50)
      preprob=precop./sum(precop)
    end
    
    #need to standardize traits for selection analysis before sperm depletion
    #standardized female phenotypes
    FMalestnd=(fphens[:,1] .- mean(fphens[:,1]))/std(fphens[:,1])

    #standardized female traits squared for gamma selection coeffients
    FMalestnd2=0.5 .* FMalestnd .^ 2

    #standardized male traits for selection analysis
    Malestnd=(mphens[:,2] .- mean(mphens[:,2]))/std(mphens[:,2])

    #standardized male traits squared for gamma selection coeffients
    Malestnd2=0.5 .* Malestnd .^ 2

    #standardized sperm number for selection analysis
    SMalestnd=(mphens[:,3] .- mean(mphens[:,3]))/std(mphens[:,3])

    #calculate mean sperm number to save for simulation output
    Meansperm=mean(mphens[:,3])

    #calculate mean RSC to save for simulation output (if present)
    MeanRSC = ntraits >= 4 ? mean(mphens[:,4]) : NaN

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

    #next part of code is to loop through all females to mate and reproduce
    for i in 1:size(fphens)[1]
      #mates = number of males a female mates with
      #Use evolving RSC trait to determine number of mates via Poisson
      if ntraits >= 4
        # Average female's RSC and mean male RSC for lambda
        female_rsc = fphens[i,4]
        mean_male_rsc = mean(mphens[:,4])
        lambda_r = (female_rsc + mean_male_rsc) / 2
        
        # Draw from Poisson with averaged lambda
        mates_draw = rand(Poisson(lambda_r))
        
        # Clamp to reasonable range [1, 4]
        mates = clamp(mates_draw, 1, 4)
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
          mphens[matesM,3]=mphens[matesM,3].*exp.(-0.2)
        else
          #calculate probs of fertilization sucess for males for cryptic female choice
          probm=prob_success(mphens[matesM,2],mphens[matesM,3],a,d)
          #probm=prob_successFR(mphens[matesM,3])
          #deplete ejaculation for males
          mphens[matesM,3]=mphens[matesM,3].*exp.(-0.2)
        end
        #calculate who got fertilization sucess based on probs
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

    #put all model results together
    #mean male,mean female, std male, std female,cor,sperm count, sperm count std,is,int,beta,gamma,A,MeanRSC,a,gen
    sumdf=[mean(mphens[:,2]),mean(fphens[:,1]),std(mphens[:,2]),std(fphens[:,1]),cor(mphens[:,2],fphens[:,1]),Meansperm,Stdsperm,is,coef(model)[1],coef(model)[2],coef(model)[3],coef(model)[4],coef(model)[5],coef(model)[6],coef(model)[7],coef(model)[8],coef(model)[9],coef(model)[10],a,MeanRSC,gen]
    dfall[gen,:]=sumdf
    #next generation

    pgm .=mutate.(pgm2)
    mgm .=mutate.(mgm2)
    pgf.=mutate.(pgf2)
    mgf.=mutate.(mgf2)
  end
  return(dfall)
end

#function to run simulation - single threaded version
function runsim(reps,N,mu,var,a,rsc,tradeoff,gens,d=-1)
  resultsP=zeros(Float64, reps*gens, 22)
  for i in 1:reps
    println("Running replicate $i of $reps...")
    resultsP[(1+(i-1)*gens):(gens*i),1:21]=sim(N,mu,var,a,rsc,tradeoff,gens,d)
    resultsP[(1+(i-1)*gens):(gens*i),22]=fill(i,gens)
  end
  return(resultsP)
end

# Test parameters (smaller values for testing)
gens=100  # Reduced from 30000 to 100 for testing
mu=1.25
var=(4*5^2/40)^0.5

# Single test run with one configuration
println("Running test simulation...")
j = true  # Using tradeoff = true
k = 1     # Using a = 1
l = 0.25  # Using rsc = 0.25
results = runsim(10, 1000, mu, var, k, l, j, gens)  # Running 5 replicates
data = DataFrame(results, [:MeanMale,:MeanFemale,:SDMale,:SDFemale,:cor,:MeanCount,:SDCount,:is,:int,:BMale,:GMale,:BFemale,:GFemale,:BSperm,:GSperm,:GMF,:GMS,:GFS,:a,:MeanRSC,:Generation,:Rep])
CSV.write("test_simulation_results.csv", data)
println("Test simulation completed and saved to test_simulation_results.csv")

# Test runs (uncomment to use)
#results=runsim(50,500,mu,var,1,0.25,true,100)
#sim(500,mu,var,1,0.25,true,100,-1)