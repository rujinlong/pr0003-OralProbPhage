# Inputs

## `yahara2021_run_accessions.tsv`

The 8 SRA runs of the Yahara 2021 long-read oral phageome (BioProject
PRJDB10485 — 4 PromethION + 4 HiSeq companion). Columns:

| col | description |
|-----|-------------|
| `run_id` | DRA / SRA accession |
| `platform` | "PromethION" / "HiSeq 2500" / etc. |
| `condition` | timepoint / donor label (used to group long+short reads) |
| `notes` | free-form; e.g. basecaller version |

The shipped TSV is a placeholder — run `data/fetch_accessions.sh` to populate it
from SRA. (Yahara 2021 uses BioProject PRJDB10485; if SRA returns different runs
than the placeholder DRR251800–DRR251807, trust SRA.)

## Reference databases (NOT shipped)

```
/scratch/$USER/dbs/checkv-db-v1.5     # ~6 GB
/scratch/$USER/dbs/virsorter2-db      # ~10 GB
```

Install:

```bash
checkv download_database /scratch/$USER/dbs/checkv-db-v1.5
virsorter setup -d /scratch/$USER/dbs/virsorter2-db -j 8
```
