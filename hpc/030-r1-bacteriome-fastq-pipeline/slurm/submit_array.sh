#!/bin/bash
#SBATCH --no-requeue
#SBATCH -t 23:55:00
#SBATCH -p cpu_p
#SBATCH -q cpu_normal

# Array-job wrapper: each task processes one SRA run from data/geng2024_sra_accessions.tsv.
# Resolve $RUN_ID by SLURM_ARRAY_TASK_ID, then dispatch to bin/030r1.slurm <subcmd>.
#
# Usage:
#   sbatch --array=1-341%20 -c 8  --mem 32G  slurm/submit_array.sh download
#   sbatch --array=1-341%10 -c 16 --mem 32G  slurm/submit_array.sh metaphlan4
#
# Override which TSV / which column to read with TSV / RUN_COL env vars.

set -euo pipefail

if [ -z "${VPIPEBIN:-}" ]; then
  [ -d "${HOME}/vpipe/bin" ] && export VPIPEBIN="${HOME}/vpipe/bin"
fi
source "${VPIPEBIN}/00-config.sh"
source "${VPIPEBIN}/01-functions.sh"

# When sbatch'ed, the script is copied to /var/spool/slurmd/... so BASH_SOURCE
# points there. Prefer SLURM_SUBMIT_DIR (the dir from which sbatch was run).
if [ -n "${SLURM_SUBMIT_DIR:-}" ] && [ -d "${SLURM_SUBMIT_DIR}/bin" ]; then
  BUNDLE_DIR="$SLURM_SUBMIT_DIR"
else
  BUNDLE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
fi
TSV="${TSV:-${BUNDLE_DIR}/data/geng2024_sra_accessions.tsv}"
SUBCMD=${1:?"first arg must be the bin/030r1.slurm subcommand (download|metaphlan4)"}

# 1-based array index → row of TSV (skipping header)
TASK_ID=${SLURM_ARRAY_TASK_ID:-1}
RUN_ID=$(awk -v i="$TASK_ID" 'NR==i+1 {print $1; exit}' "$TSV")
require_var "RUN_ID" "RUN_ID for array task $TASK_ID"

log_info "Array task $TASK_ID — dispatching $SUBCMD on RUN_ID=$RUN_ID"
bash "${BUNDLE_DIR}/bin/030r1.slurm" "$SUBCMD" "$RUN_ID"
