Description of code and data files affiliated with **Kustra and Alonzo "The coevolutionary dynamics of cryptic female choice"** . All code files (.jl or .R) are uploaded on Zenodo in **"Code.zip"** (or in github) and all data files are uploaded on dryad in **"Data.zip."** Dryad link: https://datadryad.org/stash/dataset/doi:10.7291/D1310W

To properly run code first unzip both **"Code.zip"** and **"Data.zip."** Then make sure the Code and Data files are in the following directories described below. This is most easily accomplished by moving all code files in the **[Analysis]{.underline}** directory of the unzipped **"Code.zip"** to the **[Analysis]{.underline}** directory of the unzipped **"Data.zip."** Then move **"App.R"** in the **[SI_Web_App]{.underline}** directory of the unzipped **"Code.zip"** to the **[SI_Web_App]{.underline}** directory of the unzipped **"Data.zip."** 

Below is a description of all the code and data files broken down by the directories in which code and data files should be organized in. 

**[Running_Sims:]{.underline}** Folder that contains all the files needed to run simulations. Only occurs in **"Code.zip"** and does not require any data.

-   **"RunModel.jl**": Julia script that was used for our simulations. Currently set up to run simulations at 20 loci and population size of 10000. Code is heavily commented and indicates what needs to be changed to run simulations with different loci and population sizes and starting variations.

-   **"Job_ibm.mpi"**: Example SLURM job script that was used to run simulations on super computer.

**[Analysis:]{.underline}** Folder that contains all R code to generate figures/summarize data. Occurs in both **"Code.zip"** and **"Data.zip."** Either move the **"Data"** folder from the **[Analysis:]{.underline}** folder in the unzipped **"Data.zip."** to the **[Analysis:]{.underline}** folder in  unzipped **"Code.zip"** or vice versa. 

-   **"Selection_Function_Graphing_SIFigure_3.R"**: R script to generate SI Figure 3 in the paper.

-   **"Figures.R**": R script to generate figures used in the paper. Uses files in the **[Data]{.underline}** folder

-   **"AdaptiveDynamics.jl**": Julia script used to solve the adaptive dynamics analytical model and generate SI Figure 1.

-   **[Data:]{.underline}** Folder with all the data from the model that were used in analysis/graphing.

    -   "**All_runs_shortened.rds**": Data file of simulation output thinned every 50 generations. There are 50 replicates/populations per parameter combination. This is used by "**Figures.R**". This data file was saved as a rds file to allow faster loading into R and for file compression. Here are the descriptions of the column variables:

        -   ***MeanMale***: Average male sperm trait value.

        -   ***MeanFemale***: Average cryptic female choice trait value.

        -   ***SDMale***: Standard deviation in the male sperm trait value.

        -   ***SDFemale***: Standard deviation in the cryptic female choice trait value.

        -   ***Cor***: Pearson correlation coefficient between male sperm trait and cryptic female choice trait within a population.

        -   ***MeanCount***: Average male sperm number.

        -   ***SDCount***: Standard deviation in the sperm number.

        -   ***Generation*:** Generation of the simulation.

        -   ***Rep*:** Population replicate of the simulation.

        -   ***Var*:** Starting phenotypic variation of that run. "HV" is large variation (SD = 10), "MV" is medium variation (SD = 5), and "LV" is low variation (SD = 2.5).

        -   ***Loci*:** Number of loci that determined each trait. Either "2" or "20".

        -   ***PopSize*:** Population size of the simulation run. Either is "1000" or "10,000".

        -   ***GMale*:** Gamma selection estimate of sperm trait.

        -   ***GFemale*:** Gamma selection estimate of cryptic choice trait.

        -   ***BSperm*:** Beta directional selection estimate of sperm number.

        -   ***GSperm*:** Gamma selection estimate of sperm number.

        -   ***GMF*:** Correlated gamma selection estimate between male and female trait.

        -   ***GMS*:** Correlated gamma selection estimate between sperm number and sperm trait.

        -   ***GFS*:** Correlated gamma selection estimate between sperm number and female cryptic choice trait.

        -   ***Tradeoff*:** Indicates whether a tradeoff between sperm number and male sperm trait existed. "true" means there was a tradeoff and "false" means there was no tradeoff.

        -   ***A*:** Strength of selection, *a*, parameter in the paper. "1.0" is strong selection, "12.5" is moderate selection, and "50.0" is weak selection.

        -   ***Level*:** Risk of sperm competition for that simulation run. Takes on one of the four values: 0.25, 0.5, 0.75, and 1.

    -   "**Compiled_Shortened_DataS_MiddleVar_L\_20.rds**": Data file of simulation output thinned every 50 generations for stabilizing selection with medium variation, population size of 10,000 and 20 loci per trait. There are 50 replicates/populations per parameter combination. This is used by "**Figures.R**". This data file was saved as a rds file to allow faster loading into R and for file compression. The column names and descriptions are the same as "**All_runs_shortened.rds.**"

    -   "**Last2000_MiddleVar_L\_20.rds**": Data file of simulation output for the last 2000 generations with no thinning with medium variation, population size of 10,000 and 20 loci per trait. There are 50 replicates/populations per parameter combination. This is used by "**Figures.R**". This data file was saved as a rds file to allow faster loading into R and for file compression. The column names and descriptions are the same as "**All_runs_shortened.rds.**"
   
    -   "**Compiled_Data_Shortened_FairRaffle.rds**": Data file of simulation output thinned every 50 generations for fair raffle with medium variation, population size of 10,000 and 20 loci per trait. There are 50 replicates/populations per parameter combination. This is used by "**Figures.R**". This data file was saved as a rds file to allow faster loading into R and for file compression. The column names and descriptions are the same as "**All_runs_shortened.rds**"  except it lacks selection estimates.
    
        -   "**Compiled_Shortened_Data_D.rds**": Data file of simulation output thinned every 50 generations for directional selection (females start much higher than males) with medium variation, population size of 10,000 and 20 loci per trait. There are 50 replicates/populations per parameter combination. This is used by "**Figures.R**". This data file was saved as a rds file to allow faster loading into R and for file compression. The column names and descriptions are the same as "**All_runs_shortened.rds**" except it has additional directional selection estimates. Specifically, ***BMale*:** is the beta directial selection estimate of sperm trait and ***BFemale*:** is the beta directional selection estimate of cryptic choice trait.
        
            -   "**Compiled_First1000_D.rds**": Data file of simulation output during the first 1000 generations for directional selection (females start much higher than males) with medium variation, population size of 10,000 and 20 loci per trait. There are 50 replicates/populations per parameter combination. This is used by "**Figures.R**". This data file was saved as a rds file to allow faster loading into R and for file compression. The column names and descriptions are the same as "**All_runs_shortened.rds**" except it has additional directional selection estimates. Specifically, ***BMale*:** is the beta directial selection estimate of sperm trait and ***BFemale*:** is the beta directional selection estimate of cryptic choice trait.
            
    -   "**SumRunsAll.rds**": Data file of summary statistics across all 50 replicates of a given parameter combination thinned every 50 generations. This is used by "**Figures.R**." This data file was saved as a rds to allow faster loading into R and take up less space. Here are the descriptions of the column variables:

        -   ***Generation*:** Generation of the simulation.

        -   ***Level*:** Risk of sperm competition for that simulation run. Takes on one of the four values: 0.25, 0.5, 0.75, and 1.

        -   ***A*:** Strength of selection, *a*, parameter in the paper. "1.0" is strong selection, "12.5" is moderate selection, and "50.0" is weak selection.

        -   ***Tradeoff*:** Indicates whether a tradeoff between sperm number and male sperm trait existed. "true" means there was a tradeoff and "false" means there was no tradeoff.

        -   ***Loci*:** Number of loci that determined each trait. Either "2" or "20".

        -   ***Var*:** Starting phenotypic variation of that run. "HV" is large variation (SD = 10), "MV" is medium variation (SD = 5), and "LV" is low variation (SD = 2.5).

        -   ***PopSize*:** Population size of the simulation run. Either is "1000" or "10,000".

        -   ***MeanMale_mean***: Average average male sperm trait value of a given population/simulation run at that generation.

        -   ***MeanFemale_mean***: Average average female choice trait value of a given population/simulation run at that generation.

        -   ***SDMale_mean***: Average standard deviation in the male sperm trait value of a given population/simulation run at that generation.

        -   ***SDFemale_mean***: Average standard deviation in the cryptic female choice trait value of a given population/simulation run at that generation.

        -   ***Cor_mean***: Average Pearson correlation coefficient between male sperm trait and cryptic female choice trait within a population at that generation.

        -   ***MeanCount_mean***: Average average sperm number of a given population/simulation run at that generation.

        -   ***SDCount_mean***: Average standard deviation in the sperm number of a given population/simulation run at that generation.

        -   ***Rep_mean*:** Average rep count (25.5), this column was just a sanity check that average was calculated properly.

        -   ***GMale_mean*:** Mean gamma selection estimate of sperm trait.

        -   ***GFemale_mean*:** Mean gamma selection estimate of cryptic choice trait.

        -   ***BSperm_mean*:** Mean beta directional selection estimate of sperm number.

        -   ***GSperm_mean*:** Mean gamma selection estimate of sperm number.

        -   ***GMF_mean*:** Mean correlated gamma selection estimate between male and female trait.

        -   ***GMS_mean*:** Mean correlated gamma selection estimate between sperm number and sperm trait.

        -   ***GFS_mean*:** Mean correlated gamma selection estimate between sperm number and female cryptic choice trait.

        -   ***CVCount_mean***: Average coefficient of variation of sperm number of a given population/simulation run at that generation.

        -   ***CVMale_mean***: Average coefficient of variation of male sperm trait value of a given population/simulation run at that generation.

        -   ***CVFemale_mean***: Average coefficient of variation of female choice trait value of a given population/simulation run at that generation.

        -   ***MeanMale_sd***: standard deviation of average male sperm trait value of a given population/simulation run at that generation.

        -   ***MeanFemale_sd***: standard deviation of average female choice trait value of a given population/simulation run at that generation.

        -   ***SDMale_sd***: standard deviation of standard deviation in the male sperm trait value of a given population/simulation run at that generation.

        -   ***SDFemale_sd:*** Standard deviation of standard deviation in the cryptic female choice trait value of a given population/simulation run at that generation.

        -   ***Cor_sd***: Standard deviation of Pearson correlation coefficient between male sperm trait and cryptic female choice trait within a population at that generation.

        -   ***MeanCount_sd***: Standard deviation of average sperm number of a given population/simulation run at that generation.

        -   ***SDCount_sd***: Standard deviation of standard deviation in the sperm number of a given population/simulation run at that generation.

        -   ***Rep_sd*:** Standard deviation of rep count (14.6), this column was just a sanity check that average was calculated properly.

        -   ***GMale_sd*:** Standard deviation of gamma selection estimate of sperm trait.

        -   ***GFemale_sd*:** Standard deviation of gamma selection estimate of cryptic choice trait.

        -   ***BSperm_sd*:** Standard deviation of beta directional selection estimate of sperm number.

        -   ***GSperm_sd:*** Standard deviation of gamma selection estimate of sperm number.

        -   ***GMF_sd*:** Standard deviation of correlated gamma selection estimate between male and female trait.

        -   ***GMS_sd*:** Standard deviation of correlated gamma selection estimate between sperm number and sperm trait.

        -   ***GFS_sd*:** Standard deviation of correlated gamma selection estimate between sperm number and female cryptic choice trait.

        -   ***CVCount_sd***: Standard deviation of coefficient of variation of sperm number of a given population/simulation run at that generation.

        -   ***CVMale_sd***: Standard deviation of coefficient of variation of male sperm trait value of a given population/simulation run at that generation.

        -   ***CVFemale_sd***: Standard deviation of coefficient of variation of female choice trait value of a given population/simulation run at that generation.

    -   "**All_sum_Last2000.rds**": Data file of summary statistics for model runs across all 50 replicates for the final 2,000 generations. This is used by "**Figures.R**." This data file was saved as a rds to allow faster loading into R and take up less space. Has same columns as **SumRunsAll.rds** with the addition of "median" summary statistics for each column (e.g., MeanCount_median,etc.).

    -   "**Lag_2000_all.rds**": Data file of evolutionary lag analysis performed on the last 2000 generations of all parameter combinations. This is used by "**Figures.R**." This data file was saved as a rds to allow faster loading into R and take up less space. Here are the descriptions of the column variables:

        -   ***Level*:** Risk of sperm competition for that simulation run. Takes on one of the four values: 0.25, 0.5, 0.75, and 1.

        -   ***A*:** Strength of selection, *a*, parameter in the paper. "1.0" is strong selection, "12.5" is moderate selection, and "50.0" is weak selection.

        -   ***Tradeoff*:** Indicates whether a tradeoff between sperm number and male sperm trait existed. "true" means there was a tradeoff and "false" means there was no tradeoff.

        -   ***PopSize*:** Population size of the simulation run. Either is "1000" or "10,000".

        -   ***Var*:** Starting phenotypic variation of that run. "HV" is large variation (SD = 10), "MV" is medium variation (SD = 5), and "LV" is low variation (SD = 2.5).

        -   ***Loci*:** Number of loci that determined each trait. Either "2" or "20".

        -   ***Rep:*** Population replicate of the simulation

        -   ***MinLag***: Generational lag that minimized error between male and female traits in the last 2000 generations.

        -   ***MinDiff***: Error with the generational lag that minimized the error.

        -   ***Lag0***: Error assuming no generational lag

**[SI_Web_App:]{.underline}** Folder that contains all the files and code needed to generate the supplemental web app.Occurs in both **"Code.zip"** and **"Data.zip."** Either move the **"data"** folder from the **[SI_Web_App:]{.underline}** folder in the unzipped **"Data.zip."** to the **[SI_Web_App:]{.underline}** folder in the unzipped **"Code.zip"** or vice versa. 

-   "**app.R**": The R code file that creates the supplemental web app. Code is not needed to run the web app but can be run locally if desired.

-   "**data**": Folder that contains R data file used in the supplemental web app: "run_data_shiny_shorter.rds".
