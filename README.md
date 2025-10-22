This repository contains the code and small helper files for the simulation study behind Kustra and Alonzo, "The coevolutionary dynamics of cryptic female choice".

Author / contributor notes
--------------------------
- Working under Dr. Mark Kustra (Kustra Lab). This branch contains development edits made to the original simulation framework to add an evolving RSC trait and make running small experiments easier on a local machine.

What I changed (high level)
---------------------------
- Added a fourth trait (RSC) to individual genotypes. RSC is stored in trait column 4 of the genotype arrays.
- RSC is initialized at the start of each simulation from a Normal(rsc, small_sigma) distribution so the population starts near the provided scalar rsc but has variation and can evolve.
- Per-female mate counts are now stochastic and depend on RSC: each female's number of mates is drawn from a Poisson with lambda = average(female_RSC, mean_male_RSC) and then clamped between 1 and 4. This models risk of sperm competition as an evolving trait rather than a fixed parameter.
- The code was refactored so small local runs do not spawn many workers by default. A single-process runner is provided for quick tests.

Quick run (single-process)
--------------------------
I added a small runner that avoids spawning distributed workers for quick debugging and testing.

- Runner: `Running_Sims/run_3_sims_single.jl`
- Output: `Running_Sims/run_3_sims_single_results.csv`

From the repository root (PowerShell):

```powershell
julia "Running_Sims/run_3_sims_single.jl"
```

Or from inside the `Running_Sims` folder:

```powershell
julia "run_3_sims_single.jl"
```

If you prefer to run under the VS Code Run/Debug panel, a launch configuration was added: `.vscode/launch.json` → "Run 3 sims (single)".

Dependencies
------------
Install these packages once in the Julia REPL:

```julia
import Pkg
Pkg.add(["Distributions","StatsBase","GLM","DataFrames","CSV"])
```

Notes about outputs and timing
------------------------------
- The single-process runner writes a CSV with 22 columns: the 21 columns produced by `sim()` plus a `Rep` column added by the runner.
- Long runs (large `N`, many `gens`, or many replicates) will be slow. Use the single-process runner for quick tests and reduce `N` / `gens` for debugging.

RSC details (technical)
-----------------------
- RSC is stored in genotype trait column 4 and is free to evolve (not clamped to 1).
- At start of each `sim(...)`, RSC allele values are initialized with `rand(Normal(rsc, rsc_init_sigma), ...)` where `rsc_init_sigma = max(1e-3, var/10)` so the starting variance scales with your trait variance parameter.
- Mating: for each female i, if RSC exists the code computes `lambda = max(0.001, (female_RSC + mean_male_RSC)/2)` then `mates = clamp(rand(Poisson(lambda)), 1, 4)`.

Figures and visuals
-------------------
Add your generated figures (PNG/PDF) to `Running_Sims/plots/` and reference them here. Example markdown to include an image:

```markdown
![Mean RSC across generations](Running_Sims/plots/mean_rsc_over_time.png)
```

Suggested figure captions
- Figure 1: Mean RSC (population average) across generations for three replicate runs.
- Figure 2: Distribution of mates-per-female in generation 1000 (histogram) showing Poisson-driven variance.

Troubleshooting
---------------
- If you see package errors, run the `Pkg.add(...)` commands above.
- If the script pauses in the VS Code debugger at `include(...)`, press Continue (▶) to let it run. The debugger sometimes pauses on first exceptions or breakpoints.
- If runs are too slow, reduce `N` and `gens` in `run_3_sims_single.jl` for testing.

Contact
-------
For questions about the model or the biological assumptions contact Dr. Mark Kustra (mkustra@ucsc.edu). For development questions about these local edits, reply in the repo issues or contact the developer.

