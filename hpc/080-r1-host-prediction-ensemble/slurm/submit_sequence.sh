#!/bin/bash
# Convenience driver: submit the 4-step host-prediction sequence with
# inter-job dependencies. Each step is a separate sbatch with its own
# resources, chained via afterok.
#
# Usage:
#   bash slurm/submit_sequence.sh <CONTIGS> <PREFIX> [MIN_IPHOP_SCORE=90]

set -euo pipefail

CONTIGS=${1:?"missing CONTIGS"}
PREFIX=${2:?"missing PREFIX"}
MIN_SCORE=${3:-90}
BUNDLE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
mkdir -p "$BUNDLE_DIR/logs"

J1=$(sbatch --parsable -c 16 --mem 64G -t 12:00:00 \
            -o "$BUNDLE_DIR/logs/iphop-%j.log" \
            "$BUNDLE_DIR/bin/080r1.slurm" iphop "$CONTIGS" "$PREFIX" "$MIN_SCORE")
echo "submitted iphop:           $J1"

J2=$(sbatch --parsable -c 8 --mem 32G -t 4:00:00 \
            -o "$BUNDLE_DIR/logs/sylph_imgvr-%j.log" \
            "$BUNDLE_DIR/bin/080r1.slurm" sylph_imgvr "$CONTIGS" "$PREFIX")
echo "submitted sylph_imgvr:     $J2"

J3=$(sbatch --parsable --dependency=afterok:$J2 -c 2 --mem 8G -t 1:00:00 \
            -o "$BUNDLE_DIR/logs/sylph_to_host-%j.log" \
            "$BUNDLE_DIR/bin/080r1.slurm" sylph_to_host "$PREFIX")
echo "submitted sylph_to_host:   $J3 (after $J2)"

J4=$(sbatch --parsable --dependency=afterok:$J1:$J3 -c 2 --mem 8G -t 1:00:00 \
            -o "$BUNDLE_DIR/logs/ensemble-%j.log" \
            "$BUNDLE_DIR/bin/080r1.slurm" ensemble "$PREFIX")
echo "submitted ensemble:        $J4 (after $J1, $J3)"

echo
echo "Watch progress:"
echo "  squeue -u \$USER -j ${J1},${J2},${J3},${J4}"
echo "Final output (after all complete):"
echo "  $BUNDLE_DIR/results/ensemble/${PREFIX}_ensemble_predictions.tsv"
