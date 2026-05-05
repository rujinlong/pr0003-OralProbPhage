#!/usr/bin/env bash
# Build a FASTA of the MVP source phage contigs from Yahara 2021's
# published BioProject. We need the actual sequences (MVP only stores IDs +
# host predictions). Two options:
#   1. If 070-r1 has run, the FASTA at analyses/data/070-m2a-virome-processing/r1/vOTU/vOTU_catalog.fasta
#      already covers most of these.
#   2. Otherwise, download Yahara's deposited assemblies from DDBJ/NCBI.

set -euo pipefail

OUT="$(dirname "$0")/yahara_mvp_contigs.fasta"
PROJECT="${PROJECT:-PRJDB10485}"

# Try the 070-r1 catalog first
R1_CATALOG="../../../analyses/data/070-m2a-virome-processing/r1/vOTU/vOTU_catalog.fasta"
if [ -f "$R1_CATALOG" ]; then
    echo "[extract_mvp_contigs] using 070-r1 catalog as MVP-equivalent source"
    cp "$R1_CATALOG" "$OUT"
    echo "wrote $OUT ($(grep -c '^>' "$OUT") contigs)"
    exit 0
fi

# Fall back: fetch Yahara assembly accessions
if ! command -v efetch >/dev/null 2>&1; then
    echo "Need NCBI efetch (entrez-direct) when 070-r1 catalog not available." >&2
    exit 1
fi

esearch -db assembly -query "$PROJECT" | \
    elink -target nuccore | \
    efetch -format fasta > "$OUT"

echo "wrote $OUT ($(grep -c '^>' "$OUT") contigs)"
