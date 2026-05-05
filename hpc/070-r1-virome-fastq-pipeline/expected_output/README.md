# Expected output layout

After all `bin/070r1.slurm` steps finish, `results/` (or `OUTDIR`) looks like:

```
results/
├── fastq/                                       # raw reads
│   ├── DRR251800.fastq.gz                       # PromethION (single-end)
│   ├── DRR251804_1.fastq.gz                     # HiSeq pair R1
│   └── DRR251804_2.fastq.gz                     # HiSeq pair R2
├── assembly/<condition>/
│   ├── assembly.fasta                           # Flye --meta output
│   ├── assembly_info.txt
│   ├── flye.log
│   └── flye_out/                                # full Flye dir (optional retain)
├── viral_calls/<condition>/
│   ├── out_<condition>_virsorter2/              # vpipe VirSorter2 raw output
│   ├── out_<condition>_checkv/                  # vpipe CheckV raw output
│   ├── <condition>_keep.txt                     # CheckV ≥Medium-quality contig IDs
│   └── <condition>_high_quality.fasta           # filtered phage contigs (input to vclust)
├── vOTU/
│   ├── all_phages.fasta                         # concatenation of per-condition high-quality FASTAs
│   ├── vOTU_votu.fasta                          # dereplicated catalog (95% ANI / 85% cov)
│   ├── vOTU_votus.tsv                           # per-vOTU representatives + cluster size
│   └── vOTU_clusters.tsv                        # full cluster membership
├── coverage/
│   └── <condition>_vs_catalog.tsv               # CoverM contig-level: rpkm / mean / trimmed_mean / breadth
└── validation/
    └── yahara_vs_ours.tsv                       # ID-level (default) or mash-ANI (with --yahara_fasta)
```

## Importing back into the contract layer

After `rsync` of `results/vOTU/` and `results/coverage/` to `analyses/data/070-m2a-virome-processing/r1/`:

1. Update `analyses/070-m2a-virome-processing.qmd` to read the self-built catalog + coverage tables when `params$source == "r1"`. Keep the MVP code path (read MOESM7) so the comparison persists.
2. The CoverM `rpkm` column gives real abundance values — switch `taxon_profile.abundance_kind` from `presence` to `rpkm`.
3. The 080-r1 host-prediction bundle takes `vOTU_votu.fasta` as its primary input.

## Validation gate

The ripple is "done" when one of the following holds (per `bin/validate_against_yahara.py` output):

- (default, ID-level) ≥ 80% of Yahara's MOESM7 phageIDs are recovered (per-sample sheet overlap with our catalog) at coarse-grained match
- (with `--yahara_fasta`) ≥ 50% of Yahara's contigs match a vOTU at ≥ 95% ANI

## Notes

- VirSorter2 default cuts categories 3 + 6 (low-confidence prophage / phage); CheckV "Medium-quality or better" is the secondary filter. Together they are typical for MIUViG-aligned vOTU catalogs.
- For Yahara FASTA-level validation, fetch their assemblies from DDBJ BioProject PRJDB10485 first, concatenate to one FASTA, then rerun `bin/070r1.slurm validate` with `YAHARA_FASTA=<path> bin/070r1.slurm validate` (the script reads that env var via the python helper).
