#!/bin/bash
#SBATCH --no-requeue
#SBATCH -t 23:55:00
#SBATCH -p cpu_p
#SBATCH -q cpu_normal

# ============================================================================
# 070r1.slurm — Virome fastq → vOTU catalog (HPC, vpipe-style)
# ----------------------------------------------------------------------------
# Re-process the Yahara 2021 oral phageome (PromethION + HiSeq) into a
# self-built vOTU catalog with CheckV quality flags and per-sample coverage,
# matching MIUViG / vClust conventions used elsewhere in the project.
#
# This script delegates to the existing vpipe SLURM scripts where possible:
#   - virsorter2 / checkv  -> ${VPIPEBIN}/anno_contig.slurm
#   - vOTU clustering      -> ${VPIPEBIN}/anno_contig.slurm vclust_votu
# It only adds the steps vpipe doesn't already have a wrapper for:
#   - SRA download
#   - Flye long-read assembly orchestration
#   - per-sample coverage against the catalog
#   - validation against Yahara's MOESM7 published table
#
# Subcommands:
#   download <RUN_ID>           prefetch + fasterq-dump + pigz (long or short)
#   flye <CONDITION>            Flye long-read assembly (uses COMMAND_PREFIX_FLYE)
#   virsorter2 <CONDITION>      call vpipe's virsorter2 + post-process
#   checkv <CONDITION>          call vpipe's checkv + filter to >=Medium quality
#   vclust                      vOTU dereplication across all conditions
#   coverage <CONDITION>        minimap2 short reads × catalog -> rpkm
#   validate                    compare catalog vs Yahara MOESM7 (ANI gate)
#
# Usage:
#   sbatch -c 4  --mem 8G       bin/070r1.slurm download DRR251800
#   sbatch -c 32 --mem 384G     bin/070r1.slurm flye t1
#   sbatch -c 16 --mem 64G      bin/070r1.slurm virsorter2 t1
#   sbatch -c 8  --mem 32G      bin/070r1.slurm checkv t1
#   sbatch -c 16 --mem 32G      bin/070r1.slurm vclust
#   sbatch -c 8  --mem 32G      bin/070r1.slurm coverage t1
# ============================================================================

set -euo pipefail

if [ -z "${VPIPEBIN:-}" ]; then
  [ -d "${HOME}/vpipe/bin" ] && export VPIPEBIN="${HOME}/vpipe/bin" \
    || { echo "FATAL: VPIPEBIN not set" >&2; exit 1; }
fi
source "${VPIPEBIN}/00-config.sh"
source "${VPIPEBIN}/01-functions.sh"
log_start "$@"

SUBCMD=${1:-help}
BUNDLE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ACCESSIONS="${BUNDLE_DIR}/data/yahara2021_run_accessions.tsv"
OUTDIR="${OUTDIR:-${BUNDLE_DIR}/results}"
mkdir -p "$OUTDIR"

# --- compute-node env setup ----------------------------------------------------
# Activate humann4 conda env for python3+yaml/openpyxl/Rscript+packages used by
# the validate step; and wrap prefetch/fasterq-dump in the sra-tools apptainer
# image + project-local NCBI settings file.
if ! command -v metaphlan >/dev/null 2>&1; then
  source "${HOME}/miniconda3/etc/profile.d/conda.sh"
  conda activate humann4
fi
SRA_IMG="${SRA_IMG:-/ictstr01/project/genomics/jru/singularity/depot.galaxyproject.org-singularity-sra-tools-2.11.0--pl5321ha49a11a_3.img}"
NCBI_CFG="${NCBI_CFG:-${BUNDLE_DIR}/../../oralprob/.ncbi/user-settings.mkfg}"
SRA_RUN="apptainer exec -B /ictstr01 -B /localscratch -B /tmp --env VDB_CONFIG=${NCBI_CFG} ${SRA_IMG}"

# Helper: list runs matching a column value in the accession TSV
runs_in_condition() {
  local cond=$1 platform=$2
  awk -v c="$cond" -v p="$platform" -F'\t' '
    NR>1 && $3==c && tolower($2) ~ tolower(p) { print $1 }
  ' "$ACCESSIONS"
}

case "$SUBCMD" in
  download)
    RUN_ID=${2:?"missing RUN_ID"}
    OUT_FQ="${OUTDIR}/fastq"
    mkdir -p "$OUT_FQ"; cd "$OUT_FQ"

    if compgen -G "${RUN_ID}*.fastq.gz" > /dev/null; then
      log_info "Already downloaded: $RUN_ID — skipping"; exit 0
    fi

    PLATFORM=$(awk -v r="$RUN_ID" -F'\t' 'NR>1 && $1==r {print tolower($2); exit}' "$ACCESSIONS")
    SPLIT_OPT="--split-files"
    if [[ "$PLATFORM" =~ promethion|gridion|minion|nanopore ]]; then
      SPLIT_OPT=""
      log_info "Long-read run detected ($PLATFORM)"
    fi

    log_info "Prefetching $RUN_ID"
    $SRA_RUN prefetch --max-size 200g --output-directory . "$RUN_ID"
    log_info "fasterq-dump $RUN_ID"
    $SRA_RUN fasterq-dump --threads "$SLURM_CPUS_PER_TASK" $SPLIT_OPT --skip-technical \
                 --temp ./tmp "./${RUN_ID}/${RUN_ID}.sra"
    rm -rf "$RUN_ID" ./tmp
    for f in "${RUN_ID}"*.fastq; do
      [ -f "$f" ] && pigz -p "$SLURM_CPUS_PER_TASK" "$f"
    done
    log_success "Downloaded $RUN_ID -> $OUT_FQ"
    ;;

  flye)
    CONDITION=${2:?"missing CONDITION"}
    OUT_ASM="${OUTDIR}/assembly/${CONDITION}"
    mkdir -p "$OUT_ASM"; cd "$OUT_ASM"

    LONG_RUNS=( $(runs_in_condition "$CONDITION" "promethion|gridion|minion|nanopore") )
    [ ${#LONG_RUNS[@]} -eq 0 ] && { log_error "No long-read runs for condition=$CONDITION"; exit 1; }
    READS=()
    for r in "${LONG_RUNS[@]}"; do
      f="${OUTDIR}/fastq/${r}.fastq.gz"
      require_file "$f" "Long-read fastq for $r"
      READS+=( "$f" )
    done

    FLYE_MODE=${FLYE_MODE:---nano-raw}
    log_info "Flye assembly: condition=$CONDITION, runs=${#READS[@]}, mode=$FLYE_MODE"
    $COMMAND_PREFIX_FLYE flye $FLYE_MODE "${READS[@]}" \
        --out-dir flye_out --meta --threads "$SLURM_CPUS_PER_TASK"
    cp flye_out/assembly.fasta .
    cp flye_out/assembly_info.txt .
    cp flye_out/flye.log .
    log_success "Assembly: $OUT_ASM/assembly.fasta"
    ;;

  virsorter2)
    CONDITION=${2:?"missing CONDITION"}
    ASM="${OUTDIR}/assembly/${CONDITION}/assembly.fasta"
    require_file "$ASM" "Flye assembly (run 'flye' first)"
    OUT_VC="${OUTDIR}/viral_calls/${CONDITION}"
    mkdir -p "$OUT_VC"; cd "$OUT_VC"

    log_info "Delegating to vpipe anno_contig.slurm virsorter2 (mode=phage)"
    bash "${VPIPEBIN}/anno_contig.slurm" virsorter2 phage "$ASM" 5000 0.5 "$CONDITION"
    log_success "VirSorter2 done: $OUT_VC/out_${CONDITION}_virsorter2/"
    ;;

  checkv)
    CONDITION=${2:?"missing CONDITION"}
    OUT_VC="${OUTDIR}/viral_calls/${CONDITION}"
    VS2_OUT="${OUT_VC}/out_${CONDITION}_virsorter2/final-viral-combined.fa"
    require_file "$VS2_OUT" "VirSorter2 output (run 'virsorter2' first)"
    cd "$OUT_VC"

    log_info "Delegating to vpipe anno_contig.slurm checkv"
    bash "${VPIPEBIN}/anno_contig.slurm" checkv "$VS2_OUT" "$CONDITION"

    # Filter to Medium-quality or better, write high_quality.fasta for vclust input
    QSUM="out_${CONDITION}_checkv/quality_summary.tsv"
    require_file "$QSUM" "CheckV quality_summary.tsv"
    awk -F'\t' 'NR>1 && ($8=="Complete" || $8=="High-quality" || $8=="Medium-quality"){print $1}' \
        "$QSUM" > "${CONDITION}_keep.txt"
    if [ -s "${CONDITION}_keep.txt" ]; then
      $COMMAND_PREFIX_VPFTAXA seqkit grep -f "${CONDITION}_keep.txt" "$VS2_OUT" \
          > "${CONDITION}_high_quality.fasta"
    else
      : > "${CONDITION}_high_quality.fasta"
    fi
    log_success "Filtered phages: $OUT_VC/${CONDITION}_high_quality.fasta ($(grep -c '^>' ${CONDITION}_high_quality.fasta || echo 0) seqs)"
    ;;

  vclust)
    OUT_VOTU="${OUTDIR}/vOTU"
    mkdir -p "$OUT_VOTU"; cd "$OUT_VOTU"

    # Concatenate all per-condition high-quality FASTAs
    shopt -s nullglob
    INPUTS=( "${OUTDIR}"/viral_calls/*/[A-Za-z0-9_]*_high_quality.fasta )
    shopt -u nullglob
    [ ${#INPUTS[@]} -eq 0 ] && { log_error "No *_high_quality.fasta — run 'checkv' for each condition first"; exit 1; }
    cat "${INPUTS[@]}" > all_phages.fasta
    log_info "Concatenated $(grep -c '^>' all_phages.fasta) phage contigs from ${#INPUTS[@]} conditions"

    log_info "Delegating to vpipe anno_contig.slurm vclust_votu"
    bash "${VPIPEBIN}/anno_contig.slurm" vclust_votu all_phages.fasta vOTU
    # Outputs: vOTU_votu.fasta, vOTU_votus.tsv, vOTU_clusters.tsv, ...
    log_success "vOTU catalog: $OUT_VOTU/vOTU_votu.fasta ($(grep -c '^>' vOTU_votu.fasta) vOTUs)"
    ;;

  coverage)
    CONDITION=${2:?"missing CONDITION"}
    CATALOG="${OUTDIR}/vOTU/vOTU_votu.fasta"
    require_file "$CATALOG" "vOTU catalog (run 'vclust' first)"
    OUT_COV="${OUTDIR}/coverage"
    mkdir -p "$OUT_COV"; cd "$OUT_COV"

    SHORT_RUNS=( $(runs_in_condition "$CONDITION" "hiseq|miseq|nextseq|novaseq|illumina") )
    [ ${#SHORT_RUNS[@]} -eq 0 ] && { log_error "No short-read runs for condition=$CONDITION"; exit 1; }

    # Build merged R1/R2 lists for coverm
    R1S=(); R2S=()
    for r in "${SHORT_RUNS[@]}"; do
      a="${OUTDIR}/fastq/${r}_1.fastq.gz"
      b="${OUTDIR}/fastq/${r}_2.fastq.gz"
      [ -f "$a" ] || { log_error "missing $a"; exit 1; }
      R1S+=( "$a" ); R2S+=( "$b" )
    done

    log_info "CoverM contig-level coverage: condition=$CONDITION, runs=${#R1S[@]}"
    $COMMAND_PREFIX_COVERM coverm contig \
        --reference "$CATALOG" \
        --coupled "${R1S[@]}" "${R2S[@]}" \
        --threads "$SLURM_CPUS_PER_TASK" \
        --methods rpkm mean trimmed_mean covered_fraction count length \
        --min-read-percent-identity 95 \
        --min-read-aligned-percent 75 \
        --output-format dense \
        > "${CONDITION}_vs_catalog.tsv"
    log_success "Coverage: $OUT_COV/${CONDITION}_vs_catalog.tsv"
    ;;

  validate)
    CATALOG="${OUTDIR}/vOTU/vOTU_votu.fasta"
    require_file "$CATALOG" "vOTU catalog"
    YAHARA_XLSX="${BUNDLE_DIR}/../../analyses/data/00-raw/d020-data-source-registry/yahara2021/MOESM7.xlsx"
    require_file "$YAHARA_XLSX" "Yahara MOESM7"
    log_info "ANI compare to Yahara MOESM7 (this requires you to first export their phage IDs to FASTA — see expected_output/README.md)"
    python3 "${BUNDLE_DIR}/bin/validate_against_yahara.py" \
        --catalog "$CATALOG" --yahara_xlsx "$YAHARA_XLSX" \
        --out "${OUTDIR}/validation/yahara_vs_ours.tsv"
    log_success "Validation: ${OUTDIR}/validation/yahara_vs_ours.tsv"
    ;;

  help|*)
    cat <<EOF
Usage:
  sbatch [-c N --mem MGiB] bin/070r1.slurm <subcmd> [args]

Subcommands:
  download <RUN_ID>          SRA prefetch + fasterq-dump (auto: long vs short)
  flye <CONDITION>           Flye --meta long-read assembly per condition
  virsorter2 <CONDITION>     phage identification (delegates to vpipe)
  checkv <CONDITION>         quality + filter (delegates to vpipe)
  vclust                     vOTU dereplication across all conditions (vpipe vclust_votu)
  coverage <CONDITION>       short-read coverage vs catalog (CoverM)
  validate                   ANI vs Yahara MOESM7

Env (override defaults):
  OUTDIR     output root (default: \$BUNDLE_DIR/results)
  FLYE_MODE  --nano-raw (default) | --nano-hq | --nano-corr
  VPIPEBIN   vpipe bin dir, defaults to ~/vpipe/bin
EOF
    ;;
esac

log_end "$@"
