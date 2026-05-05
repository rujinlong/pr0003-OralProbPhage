# Expected output layout

After `bin/030r1.slurm` runs all four steps (`download → metaphlan4 → merge → validate`), the `results/` tree (or whatever `OUTDIR` is set to) looks like:

```
results/
├── fastq/                                   # raw fastq.gz, one per SRA run
│   ├── DRR285642_1.fastq.gz
│   └── ...
├── metaphlan/<bioproject>/<run_id>.metaphlan.tsv
├── metaphlan/<bioproject>/<run_id>.bowtie2.bz2
├── species_abundance/
│   ├── <BP>.tsv                             # wide, merge_metaphlan_tables.py output
│   ├── <BP>_long.tsv                        # long, DBI-ingestible (run_id, species, rel_abund)
│   └── <BP>_validation.tsv                  # per-(run, species) delta vs Geng
└── (no other artefacts — all intermediate work happens inline in $TMPDIR)
```

## Importing back into the contract layer

After `rsync -av <hpc>:~/oralprob/030-r1/results/species_abundance/ ./analyses/data/030-preprocessing-pipeline/species_abundance/`:

1. Update `analyses/030-preprocessing-pipeline.qmd` to read each `<BP>_long.tsv` and emit `geng_vs_ours.tsv` per cohort.
2. Optional 040-r2: switch its source from `Figure 01.Rdata` to our long-format TSVs.

## Validation gate

The 030-r1 ripple is "done" when `bin/validate_against_geng.R` reports, for **every** modeling cohort:

- median |delta| (per-run, per-species relative-abundance gap vs Geng's matrix) < 0.05
- < 10% of Geng's species are missing from our profile

If both hold, the MVP signal is robust to pipeline reproduction. The R script `quit(status=1)` if either fails, so SLURM exit codes propagate.

## Bonus: validating the two cohorts WITHOUT a published table

`PRJDB6966` and `PRJNA552294` have NO entry in Geng's `Figure 01.Rdata` (they were used as "validation cohorts" in the meta-analysis, not training). They DO have SRA fastq, so 030-r1 can still process them — there's just no Geng baseline to validate against. The `bin/validate_against_geng.R` script will skip the gate for these two and report "no Geng baseline available; emitting profile only".
