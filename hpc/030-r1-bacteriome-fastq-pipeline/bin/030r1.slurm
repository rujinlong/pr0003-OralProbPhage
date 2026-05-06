#!/bin/bash
#SBATCH --no-requeue
#SBATCH -t 23:55:00
#SBATCH -p cpu_p
#SBATCH -q cpu_normal

# ============================================================================
# 030r1.slurm — Bacteriome fastq → MetaPhlAn4 species table (HPC, vpipe-style)
# ----------------------------------------------------------------------------
# Re-process the Geng 2024 cohorts' SRA fastq into per-sample species
# relative-abundance profiles, then merge into per-bioproject tables that
# match Geng's published Figure 01.Rdata format. Validation is the gate.
#
# Subcommands:
#   download <RUN_ID>           prefetch + fasterq-dump + pigz
#   metaphlan4 <RUN_ID>         MetaPhlAn4 species profile (uses DB_BIOBAKERY)
#   merge <BIOPROJECT>          merge per-run profiles -> wide + long TSV
#   validate <BIOPROJECT>       compare to Geng's Rdata (median-delta gate)
#
# Usage as a single job:
#   sbatch -c 8 --mem 32G bin/030r1.slurm download DRR285642
#   sbatch -c 16 --mem 32G bin/030r1.slurm metaphlan4 DRR285642
#
# Usage as a SLURM array over an accession TSV (see slurm/submit_array.sh):
#   sbatch --array=1-42%20 -c 16 --mem 32G slurm/submit_array.sh metaphlan4
# ============================================================================

set -euo pipefail

# Bootstrap VPIPEBIN — vpipe must be installed at ~/vpipe (or VPIPEBIN exported)
if [ -z "${VPIPEBIN:-}" ]; then
  if [ -d "${HOME}/vpipe/bin" ]; then
    export VPIPEBIN="${HOME}/vpipe/bin"
  else
    echo "FATAL: VPIPEBIN not set and ~/vpipe/bin not found" >&2; exit 1
  fi
fi
source "${VPIPEBIN}/00-config.sh"
source "${VPIPEBIN}/01-functions.sh"
log_start "$@"

SUBCMD=${1:-help}

# Bundle paths (relative to this script's directory)
BUNDLE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ACCESSIONS="${BUNDLE_DIR}/data/geng2024_sra_accessions.tsv"
OUTDIR="${OUTDIR:-${BUNDLE_DIR}/results}"
mkdir -p "$OUTDIR"

# --- compute-node env setup ----------------------------------------------------
# This HPC has no system-wide metaphlan/sra-tools. We activate the humann4 conda
# env (provides: metaphlan, merge_metaphlan_tables.py, Rscript+packages, python3
# w/ yaml) and wrap prefetch/fasterq-dump via apptainer + a project-local NCBI
# settings file (sra-tools 2.11 ships unconfigured).
if ! command -v metaphlan >/dev/null 2>&1; then
  source "${HOME}/miniconda3/etc/profile.d/conda.sh"
  conda activate humann4
fi
SRA_IMG="${SRA_IMG:-/ictstr01/project/genomics/jru/singularity/depot.galaxyproject.org-singularity-sra-tools-2.11.0--pl5321ha49a11a_3.img}"
NCBI_CFG="${NCBI_CFG:-${BUNDLE_DIR}/../../oralprob/.ncbi/user-settings.mkfg}"
SRA_RUN="apptainer exec -B /ictstr01 -B /localscratch -B /tmp --env VDB_CONFIG=${NCBI_CFG} ${SRA_IMG}"

case "$SUBCMD" in
  download)
    RUN_ID=${2:?"missing RUN_ID"}
    require_file "$ACCESSIONS" "Accession TSV"
    OUT_FQ="${OUTDIR}/fastq"
    mkdir -p "$OUT_FQ"
    cd "$OUT_FQ"

    if [ -f "${RUN_ID}_1.fastq.gz" ] || [ -f "${RUN_ID}.fastq.gz" ]; then
      log_info "Already downloaded: $RUN_ID — skipping"
      exit 0
    fi

    log_info "Prefetching $RUN_ID"
    $SRA_RUN prefetch --max-size 200g --output-directory . "$RUN_ID"
    log_info "fasterq-dump $RUN_ID"
    $SRA_RUN fasterq-dump --threads "$SLURM_CPUS_PER_TASK" --split-files --skip-technical \
                 --temp ./tmp "./${RUN_ID}/${RUN_ID}.sra"
    rm -rf "$RUN_ID" ./tmp
    for f in "${RUN_ID}"*.fastq; do
      [ -f "$f" ] && pigz -p "$SLURM_CPUS_PER_TASK" "$f"
    done
    log_success "Downloaded $RUN_ID -> $OUT_FQ"
    ;;

  metaphlan4)
    RUN_ID=${2:?"missing RUN_ID"}
    require_db "DB_BIOBAKERY"
    require_file "$ACCESSIONS" "Accession TSV"

    FQ_DIR="${OUTDIR}/fastq"
    if   [ -f "${FQ_DIR}/${RUN_ID}_1.fastq.gz" ] && [ -f "${FQ_DIR}/${RUN_ID}_2.fastq.gz" ]; then
      IN="${FQ_DIR}/${RUN_ID}_1.fastq.gz,${FQ_DIR}/${RUN_ID}_2.fastq.gz"
    elif [ -f "${FQ_DIR}/${RUN_ID}.fastq.gz" ]; then
      IN="${FQ_DIR}/${RUN_ID}.fastq.gz"
    else
      log_error "No fastq found for $RUN_ID under $FQ_DIR — run 'download' first"; exit 1
    fi

    BP=$(awk -v r="$RUN_ID" 'BEGIN{FS="\t"} NR>1 && $1==r {print $3; exit}' "$ACCESSIONS")
    : "${BP:=unknown}"
    OUT_BP="${OUTDIR}/metaphlan/${BP}"
    mkdir -p "$OUT_BP"
    cd "$OUT_BP"

    log_info "MetaPhlAn4: $RUN_ID (bioproject=$BP)"
    metaphlan "$IN" \
        --input_type fastq \
        --bowtie2db "${DB_BIOBAKERY}/humann4/metaphlan" \
        --index mpa_vOct22_CHOCOPhlAnSGB_202403 \
        --bowtie2out "${RUN_ID}.bowtie2.bz2" \
        --nproc "$SLURM_CPUS_PER_TASK" \
        --analysis_type rel_ab \
        --tax_lev s \
        -o "${RUN_ID}.metaphlan.tsv"
    log_success "Profile: $OUT_BP/${RUN_ID}.metaphlan.tsv"
    ;;

  merge)
    BP=${2:?"missing BIOPROJECT"}
    PROFILE_DIR="${OUTDIR}/metaphlan/${BP}"
    OUT_DIR="${OUTDIR}/species_abundance"
    mkdir -p "$OUT_DIR"
    require_dir "$PROFILE_DIR" "Per-run profile dir"

    profiles=( "${PROFILE_DIR}"/*.metaphlan.tsv )
    [ ${#profiles[@]} -eq 0 ] && { log_error "No *.metaphlan.tsv in $PROFILE_DIR"; exit 1; }

    cd "$OUT_DIR"
    log_info "Merging ${#profiles[@]} profiles for $BP"
    merge_metaphlan_tables.py "${profiles[@]}" > "${BP}.tsv"
    python3 "${BUNDLE_DIR}/bin/long_format.py" "${BP}.tsv" "${BP}_long.tsv" "$BP"
    log_success "Wrote ${OUT_DIR}/${BP}.tsv (+ _long.tsv)"
    ;;

  validate)
    BP=${2:?"missing BIOPROJECT"}
    LONG_TSV="${OUTDIR}/species_abundance/${BP}_long.tsv"
    require_file "$LONG_TSV" "Long-format profile TSV (run merge first)"
    Rscript "${BUNDLE_DIR}/bin/validate_against_geng.R" "$LONG_TSV" "$BP"
    ;;

  help|*)
    cat <<EOF
Usage:
  sbatch [-c N --mem MGiB] bin/030r1.slurm <subcmd> [args]

Subcommands:
  download <RUN_ID>        prefetch + fasterq-dump + pigz a single SRA run
  metaphlan4 <RUN_ID>      MetaPhlAn4 species profile (uses DB_BIOBAKERY/humann4/metaphlan)
  merge <BIOPROJECT>       merge per-run *.metaphlan.tsv into wide + long tables
  validate <BIOPROJECT>    diff our long table vs Geng's Figure 01.Rdata

Env:
  OUTDIR    output root (default: \$BUNDLE_DIR/results)
  VPIPEBIN  required, defaults to ~/vpipe/bin

Tools needed natively in PATH:
  prefetch, fasterq-dump, pigz, metaphlan, merge_metaphlan_tables.py
  (these are part of biobakery and sra-tools that vpipe assumes already)
EOF
    ;;
esac

log_end "$@"
