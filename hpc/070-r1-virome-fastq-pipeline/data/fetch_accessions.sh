#!/usr/bin/env bash
# Fetch the actual run accessions for Yahara 2021 (BioProject PRJDB10485)
# from SRA. Run this once on a node with internet access; output replaces
# data/yahara2021_run_accessions.tsv.
#
# Usage:
#   bash data/fetch_accessions.sh > data/yahara2021_run_accessions.tsv

set -euo pipefail

PROJECT="${PROJECT:-PRJDB10485}"

# Need NCBI's pysradb or efetch. Both work; pysradb is often less brittle.
if command -v pysradb >/dev/null 2>&1; then
    pysradb metadata "$PROJECT" --detailed | \
        awk -F'\t' 'NR==1{
                       for(i=1;i<=NF;i++) col[$i]=i
                       print "run_id\tplatform\tcondition\tnotes"
                       next
                    }
                    { printf "%s\t%s\t%s\t%s\n", $col["run_accession"], $col["instrument_model"], $col["library_name"], "" }'
elif command -v efetch >/dev/null 2>&1; then
    esearch -db sra -query "$PROJECT" | \
        efetch -format runinfo | \
        awk -F',' 'NR==1{print "run_id\tplatform\tcondition\tnotes"; next}
                   {print $1"\t"$19"\t"$25"\t"}'
else
    echo "Need pysradb (pip install pysradb) or NCBI efetch (entrez-direct)" >&2
    exit 1
fi
