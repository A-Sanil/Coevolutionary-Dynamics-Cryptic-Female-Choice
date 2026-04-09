# RunModel_KBufferV2 Modifications

## 1. K and child buffer sizing
- `sim(..., K=N,...` now default.
- `K` is forced `max(N, K)` to avoid too-small capacities.
- `bufsize` set to `max(2*K, max_offspring*N, max_mates*N)`.
- This satisfies:
  - "Reduce size of K array to N"
  - "Storing children array max = 2K"
  - "max children array based on max_mates"

## 2. max mates and mate clamping
- New `max_mates` parameter (default 5).
- Mates from `Poisson(lambda_mates)` clamped by `clamp(mates, 0, max_mates)`.
- Ensures plates range from 0..4/5 depending on config.

## 3. offspring generation modes (Y offspring vs X mates)
- New `offspring_mode` parameter:
  - `:poisson` (existing RSC Poisson mode, with `offspring_scale`)
  - `:logistic` (`max/(1+exp(-r*(mates-c)))`)
  - `:expdecay` (`max*exp(-r*mates)`)
  - `:gaussian` (`max*exp(-((mates-mu)^2)/(2*sigma^2))`)
- Helper function: `offspring_from_mates(...)`.
- New parameters: `offspring_r`, `offspring_c`, `offspring_mu`, `offspring_sigma`, `max_offspring`.
- `dynamic_offspring_mode` must be `true` to enable this per-female calculation.

## 4. API updates
### `sim` signature
- Added parameters:
  - `offspring_mode`, `max_mates`, `max_offspring`, `offspring_r`, `offspring_c`, `offspring_mu`, `offspring_sigma`.

### `runsim`/`runsim_serial` signatures and internal calls
- Updated to forward new parameters to `sim`.

## 5. runner script updates (`run_kbuffer_simpleV2.jl`)
- Set defaults for new variables.
- Added param logging on start.
- Updated `runsim_serial` call to pass new parameters.

## 6. R graph script
- Added `offspring_mate_curves.R` to plot logistic, expdecay, gaussian.
- Usage example:
  - `Rscript Running_Sims/offspring_mate_curves.R max_offspring=5 r=1.4 c=3 mu=3 sigma=1 max_mates=10 outfile=plot.png`

## 7. runtime result (sample)
- For `N=40`, `generations=10`, `replicates=1`, measured 28.53s wall-clock.
- Peak RSS ~ 626 MB, user CPU 120.57s, sys 34.42s.

## 8. How to fix if issues
- If passing keywords, use positional args based on function signature in this version.
- If `offspring_mode` is invalid, it throws `Unsupported offspring_function`.
- To turn off new dynamic mode: `dynamic_offspring_mode=false`.
