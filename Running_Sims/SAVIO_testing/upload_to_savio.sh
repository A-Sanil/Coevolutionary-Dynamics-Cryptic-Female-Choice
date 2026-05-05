#!/usr/bin/env bash
set -euo pipefail

# upload_to_savio.sh
# Usage: ./upload_to_savio.sh
# Environment variables (optional):
#  SAVIO_USER - defaults to 'asanil'
#  SAVIO_HOST - defaults to 'hpc.brc.berkeley.edu'
#  KEYFILE    - optional path to your SSH private key
#  SRC        - local source directory to sync (default: current directory)
#  REMOTE_DIR - remote directory to place files (default: /global/scratch/users/$SAVIO_USER/Coevolutionary-Dynamics-Cryptic-Female-Choice/Running_Sims/SAVIO_testing)

SAVIO_USER=${SAVIO_USER:-asanil}
SAVIO_HOST=${SAVIO_HOST:-hpc.brc.berkeley.edu}
KEYFILE=${KEYFILE:-}
SRC=${SRC:-.}
REMOTE_DIR=${REMOTE_DIR:-/global/scratch/users/$SAVIO_USER/Coevolutionary-Dynamics-Cryptic-Female-Choice/Running_Sims/SAVIO_testing}

# Files to sync by default (relative to SRC)
DEFAULT_FILES=(run_kbuffer_simpleV2.jl savio_kbuffer_v2.sbatch RunModel_KBufferV2.jl)

# Create local staging dir list
FILES_TO_SYNC=()
for f in "${DEFAULT_FILES[@]}"; do
  if [ -e "$SRC/$f" ]; then
    FILES_TO_SYNC+=("$SRC/$f")
  fi
done

# If nothing found, fallback to syncing entire SRC folder
if [ ${#FILES_TO_SYNC[@]} -eq 0 ]; then
  echo "No primary files found in $SRC; will sync entire folder contents."
  FILES_TO_SYNC=("$SRC/")
else
  echo "Found ${#FILES_TO_SYNC[@]} files to sync: ${FILES_TO_SYNC[*]}"
fi

# Build rsync command
RSYNC_OPTS=(-avz --progress --delete)
if [ -n "$KEYFILE" ]; then
  SSH_CMD="ssh -i $KEYFILE -o StrictHostKeyChecking=accept-new"
else
  SSH_CMD="ssh -o StrictHostKeyChecking=accept-new"
fi

REMOTE_TARGET="$SAVIO_USER@$SAVIO_HOST:$REMOTE_DIR"

echo "Uploading to $REMOTE_TARGET"
ssh ${SAVIO_USER}@${SAVIO_HOST} "mkdir -p '$REMOTE_DIR'"

# Run rsync for each path
for path in "${FILES_TO_SYNC[@]}"; do
  if [ -d "$path" ]; then
    rsync "${RSYNC_OPTS[@]}" -e "$SSH_CMD" "$path" "$REMOTE_TARGET/"
  else
    rsync "${RSYNC_OPTS[@]}" -e "$SSH_CMD" "$path" "$REMOTE_TARGET/"
  fi
done

echo "Upload complete. Remote files are in: $REMOTE_DIR"

echo "To submit on Savio, ssh to the cluster and run (example):"
echo "  ssh ${SAVIO_USER}@${SAVIO_HOST}"
echo "  cd $REMOTE_DIR"
echo "  sbatch -A YOUR_REAL_ACCOUNT savio_kbuffer_v2.sbatch"
