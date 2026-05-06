#!/bin/bash
#SBATCH --no-requeue
#SBATCH -t 23:55:00
#SBATCH -p cpu_p
#SBATCH -q cpu_normal

# Array-job wrapper for 070-r1. Two modes:
#   1. Per-run download:
#        sbatch --array=1-8 -c 4 --mem 8G slurm/submit_array.sh download
#   2. Per-condition step (flye/virsorter2/checkv/coverage):
#        sbatch --array=1-4 -c 16 --mem 64G slurm/submit_array.sh virsorter2

set -euo pipefail

if [ -z "${VPIPEBIN:-}" ]; then
  [ -d "${HOME}/vpipe/bin" ] && export VPIPEBIN="${HOME}/vpipe/bin"
fi
source "${VPIPEBIN}/00-config.sh"
source "${VPIPEBIN}/01-functions.sh"

if [ -n "${SLURM_SUBMIT_DIR:-}" ] && [ -d "${SLURM_SUBMIT_DIR}/bin" ]; then
  BUNDLE_DIR="$SLURM_SUBMIT_DIR"
else
  BUNDLE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
fi
TSV="${TSV:-${BUNDLE_DIR}/data/yahara2021_run_accessions.tsv}"
SUBCMD=${1:?"first arg must be subcommand: download|flye|virsorter2|checkv|coverage"}

TASK_ID=${SLURM_ARRAY_TASK_ID:-1}

case "$SUBCMD" in
  download)
    # iterate over runs (rows of TSV)
    KEY=$(awk -v i="$TASK_ID" 'NR==i+1 {print $1; exit}' "$TSV")
    ;;
  flye|virsorter2|checkv|coverage)
    # iterate over unique conditions (col 3 of TSV)
    KEY=$(awk -F'\t' 'NR>1 {print $3}' "$TSV" | sort -u | sed -n "${TASK_ID}p")
    ;;
  *)
    log_error "Unknown subcommand: $SUBCMD"; exit 1
    ;;
esac

require_var "KEY" "Array key for task $TASK_ID"
log_info "Array task $TASK_ID — dispatching $SUBCMD on KEY=$KEY"
bash "${BUNDLE_DIR}/bin/070r1.slurm" "$SUBCMD" "$KEY"
