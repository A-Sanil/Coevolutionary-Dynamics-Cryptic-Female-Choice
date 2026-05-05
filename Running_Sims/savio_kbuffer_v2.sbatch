#!/bin/bash
#SBATCH --job-name=kbufferv2
#SBATCH --account=co_youraccount
#SBATCH --partition=savio3
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=08:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail

module load julia

# Use explicit username path as requested.
SCRATCH_ROOT="/global/scratch/users/asanil"

# Move to the simulation folder in scratch storage.
cd "$SCRATCH_ROOT/Coevolutionary-Dynamics-Cryptic-Female-Choice/Running_Sims"

# Keep outputs in scratch storage, not in home.
export BRC_OUTPUT_ROOT="$SCRATCH_ROOT/kbuffer_csv"
mkdir -p "$BRC_OUTPUT_ROOT" logs

# Run worker count should generally be cpus-per-task minus 1.
CPUS="${SLURM_CPUS_PER_TASK:-16}"
export N_WORKERS=$((CPUS>1 ? CPUS-1 : 1))

julia --project=@. run_kbuffer_simpleV2.jl
