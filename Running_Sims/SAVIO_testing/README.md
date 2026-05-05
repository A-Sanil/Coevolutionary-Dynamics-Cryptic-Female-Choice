SAVIO_testing

This folder contains a helper script to upload selected simulation files to Savio scratch and a short usage guide.

Files
- upload_to_savio.sh : rsync-based uploader. Edit env vars as needed.

How to use
1. Make the script executable:
   chmod +x upload_to_savio.sh

2. Optional: set environment variables to match your cluster and key:
   export SAVIO_USER=asanil
   export SAVIO_HOST=hpc.brc.berkeley.edu
   export KEYFILE=/path/to/your/private_key   # optional
   export SRC=.                                 # local path to Running_Sims

3. Run the uploader from `Running_Sims` (or set SRC):
   ./SAVIO_testing/upload_to_savio.sh

Notes
- By default the script uploads `run_kbuffer_simpleV2.jl`, `savio_kbuffer_v2.sbatch`, and `RunModel_KBufferV2.jl` if present. If none of those are found it falls back to syncing the entire `SRC` directory.
- Remote target default: `/global/scratch/users/<SAVIO_USER>/Coevolutionary-Dynamics-Cryptic-Female-Choice/Running_Sims/SAVIO_testing`.
- After upload, `sbatch -A YOUR_REAL_ACCOUNT savio_kbuffer_v2.sbatch` from the remote directory will submit the job (set account as needed).

Security
- If you supply `KEYFILE`, it will be used for SSH. Ensure the key has correct permissions (`chmod 600 keyfile`).
