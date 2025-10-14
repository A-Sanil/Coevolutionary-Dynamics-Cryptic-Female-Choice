Quick run instructions

Run the single-process 3-sim runner.

From PowerShell (repo root):

julia "Running_Sims/run_3_sims_single.jl"

Or in VS Code:
- Open Run and Debug (Ctrl+Shift+D)
- Choose "Run 3 sims (single)"
- Press the green play button (or F5)

Notes:
- The runner sets NO_AUTO_ADDPROCS=1 so the script will not spawn distributed workers.
- Results are saved to Running_Sims/run_3_sims_single_results.csv
