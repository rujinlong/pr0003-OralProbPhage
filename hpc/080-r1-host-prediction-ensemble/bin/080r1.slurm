#!/bin/bash
#SBATCH --no-requeue
#SBATCH -t 23:55:00
#SBATCH -p cpu_p
#SBATCH -q cpu_normal

# ============================================================================
# 080r1.slurm — Phage host prediction ensemble (HPC, vpipe-style)
# ----------------------------------------------------------------------------
# Combines multiple host predictors over a vOTU FASTA:
#   1. iPHoP             — meta-predictor (BLAST + CRISPR + tRNA + composition)
#   2. Sylph vs IMG/VR   — fast ANI sketch -> match -> host transitivity
# Both delegate to vpipe's existing anno_contig.slurm subcommands. The
# ensemble step combines per-method outputs via weighted voting.
#
# Subcommands:
#   iphop <contigs> <prefix> [min_score=90]
#   sylph_imgvr <contigs> <prefix>
#   sylph_to_host <prefix>             # join sylph hits to IMG/VR host taxonomy
#   ensemble <prefix>                   # combine all method outputs
#
# Optional methods (require extra images NOT in stock vpipe; skip in MVP):
#   spacepharer <contigs> <prefix>      # CRISPR spacer host calls
#   rafah <contigs> <prefix>            # random-forest family-level
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
OUTDIR="${OUTDIR:-${BUNDLE_DIR}/results}"
mkdir -p "$OUTDIR"

case "$SUBCMD" in
  iphop)
    CONTIGS=${2:?"missing CONTIGS FASTA"}
    PREFIX=${3:?"missing PREFIX"}
    MINSCORE=${4:-90}
    require_file "$CONTIGS" "Phage FASTA"
    require_db "DB_IPHOP"
    OUT="${OUTDIR}/iphop"
    mkdir -p "$OUT"; cd "$OUT"

    log_info "Delegating to vpipe anno_contig.slurm iphop"
    bash "${VPIPEBIN}/anno_contig.slurm" iphop "$CONTIGS" "$PREFIX" "$MINSCORE"
    # Normalize: vOTU_id, host_genus, host_score, method
    SRC="out_${PREFIX}_iphop/Host_prediction_to_genus_m90.csv"
    require_file "$SRC" "iPHoP genus output"
    python3 - "$SRC" "${OUTDIR}/iphop/iphop.tsv" <<'PY'
import csv, sys
src, dst = sys.argv[1], sys.argv[2]
with open(src) as fh, open(dst, "w", newline="") as out:
    r = csv.DictReader(fh)
    w = csv.writer(out, delimiter="\t")
    w.writerow(["vOTU_id","host_genus","host_score","method"])
    for row in r:
        w.writerow([row["Virus"], row["Host genus"],
                    row["Confidence score"], "iphop"])
PY
    log_success "iPHoP normalized: ${OUTDIR}/iphop/iphop.tsv"
    ;;

  sylph_imgvr)
    CONTIGS=${2:?"missing CONTIGS FASTA"}
    PREFIX=${3:?"missing PREFIX"}
    require_file "$CONTIGS" "Phage FASTA"
    require_db "DB_SYLPH_IMGVR"
    OUT="${OUTDIR}/sylph"
    mkdir -p "$OUT"; cd "$OUT"

    log_info "Delegating to vpipe anno_contig.slurm sylph imgvr fa"
    bash "${VPIPEBIN}/anno_contig.slurm" sylph "$PREFIX" imgvr fa "$CONTIGS"
    # Output: ac_${PREFIX}_sylph_imgvr.tsv (raw sylph hits)
    log_success "Sylph IMG/VR raw: $OUT/ac_${PREFIX}_sylph_imgvr.tsv"
    ;;

  sylph_to_host)
    PREFIX=${2:?"missing PREFIX"}
    OUT="${OUTDIR}/sylph"
    SRC="${OUT}/ac_${PREFIX}_sylph_imgvr.tsv"
    require_file "$SRC" "Sylph IMG/VR output (run 'sylph_imgvr' first)"

    # IMG/VR metadata path — adjust if your DB layout differs
    IMGVR_META="${IMGVR_META:-${DBPATH}/virome_db/imgvr/IMGVR_all_Sequence_information-high_confidence.tsv}"
    require_file "$IMGVR_META" "IMG/VR metadata TSV (set IMGVR_META to override)"

    log_info "Joining sylph hits to IMG/VR host taxonomy"
    python3 "${BUNDLE_DIR}/bin/sylph_to_host.py" \
        --sylph "$SRC" --imgvr_meta "$IMGVR_META" \
        --out "${OUT}/sylph_imgvr.tsv"
    log_success "Sylph→host: $OUT/sylph_imgvr.tsv"
    ;;

  ensemble)
    PREFIX=${2:?"missing PREFIX"}
    OUT="${OUTDIR}/ensemble"
    mkdir -p "$OUT"

    IPHOP="${OUTDIR}/iphop/iphop.tsv"
    SYLPH="${OUTDIR}/sylph/sylph_imgvr.tsv"
    [ -f "$IPHOP" ] || { log_warn "Missing $IPHOP; ensemble will skip iphop method"; IPHOP=""; }
    [ -f "$SYLPH" ] || { log_warn "Missing $SYLPH; ensemble will skip sylph method"; SYLPH=""; }

    log_info "Combining methods into ensemble"
    python3 "${BUNDLE_DIR}/scripts/ensemble.py" \
        --iphop "$IPHOP" \
        --blast "$SYLPH" \
        --calibration "${BUNDLE_DIR}/scripts/method_weights.tsv" \
        --out_prefix "${OUT}/${PREFIX}_ensemble"
    log_success "Ensemble: ${OUT}/${PREFIX}_ensemble_predictions.tsv"
    ;;

  spacepharer|rafah)
    log_warn "$SUBCMD requires an apptainer image not in stock vpipe."
    log_warn "Either add a COMMAND_PREFIX_$(echo $SUBCMD | tr a-z A-Z) image and a wrapper to anno_contig.slurm,"
    log_warn "or skip this method (ensemble.py handles missing methods gracefully)."
    exit 2
    ;;

  help|*)
    cat <<EOF
Usage:
  sbatch [-c N --mem MGiB] bin/080r1.slurm <subcmd> [args]

Subcommands (all run on a vOTU FASTA):
  iphop <CONTIGS> <PREFIX> [MIN_SCORE=90]
  sylph_imgvr <CONTIGS> <PREFIX>
  sylph_to_host <PREFIX>
  ensemble <PREFIX>

Typical pipeline (after cd hpc/080-r1-host-prediction-ensemble):
  CONTIGS=../070-r1-virome-fastq-pipeline/results/vOTU/vOTU_votu.fasta
  sbatch -c 16 --mem 64G bin/080r1.slurm iphop "\$CONTIGS" oralprob 90
  sbatch -c 8  --mem 16G bin/080r1.slurm sylph_imgvr "\$CONTIGS" oralprob
  sbatch -c 2  --mem 8G  bin/080r1.slurm sylph_to_host oralprob
  sbatch -c 2  --mem 8G  bin/080r1.slurm ensemble oralprob

Env:
  OUTDIR        output root (default: \$BUNDLE_DIR/results)
  IMGVR_META    path to IMGVR_all_Sequence_information-high_confidence.tsv
EOF
    ;;
esac

log_end "$@"
