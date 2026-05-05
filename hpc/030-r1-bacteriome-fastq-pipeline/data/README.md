# Inputs

## `geng2024_sra_accessions.tsv`

341 SRA runs from the Geng et al. 2024 8-cohort meta-analysis, extracted from
Table S2 of the paper's supplementary. Columns:

| col | description |
|-----|-------------|
| `run_id` | SRA / DRA / ERA run accession (DRR* / SRR* / ERR*) |
| `biosample_id` | NCBI BioSample accession |
| `bioproject_id` | NCBI BioProject (one of: PRJDB11203, PRJDB6966, PRJNA230363, PRJNA396840, PRJNA552294, PRJNA678453, PRJNA717815, PRJNA932553) |
| `group` | `Case` (periodontitis) or `Control` (healthy) |
| `country` | sample origin |
| `sex` | reported sex |
| `host_age` | reported age (numeric or NA) |

The 6 "modeling" cohorts have published per-sample abundance tables in
Geng's `Figure 01.Rdata`; the two validation cohorts (`PRJDB6966`, `PRJNA552294`)
do not, and are an MVP **gap** that 030-r1 closes.

## Reference databases (NOT shipped here)

You need a MetaPhlAn4 database (~25 GB). Download once on the cluster:

```bash
metaphlan --install --bowtie2db /scratch/$USER/dbs/metaphlan4
```

Then point `params.metaphlan_db` in `nextflow/nextflow.config` at that path.
