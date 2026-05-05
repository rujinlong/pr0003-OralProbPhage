# 070-r1 — Virome fastq → vOTU catalog (HPC, vpipe-style)

**Goal**: Re-process the Yahara et al. 2021 long-read PromethION saliva data into a self-built vOTU catalog with CheckV quality flags, MIUViG-conformant clustering (95% ANI / 85% coverage), and per-sample short-read coverage. Validate against the published `MOESM7.xlsx` contigs the MVP loads.

This is the first ripple of `070-m2a-virome-processing`.

## vpipe integration

All heavy steps delegate to `${VPIPEBIN}/anno_contig.slurm`:

| Step | vpipe subcommand | apptainer image | DB |
|------|------------------|-----------------|-----|
| Long-read assembly | (this bundle uses `$COMMAND_PREFIX_FLYE`) | `flye.sif` | — |
| Viral identification | `anno_contig.slurm virsorter2 phage ...` | `quay-virsorter2.img` | `DB_VIRSORTER2` |
| Quality flagging | `anno_contig.slurm checkv ...` | `quay-checkv.img` | `DB_CHECKV` |
| vOTU clustering | `anno_contig.slurm vclust_votu ...` | `jinlongru-cluster*.img` | — |
| Read coverage | (this bundle uses `$COMMAND_PREFIX_COVERM`) | `quay-coverm.img` | — |

Native PATH tools assumed: `prefetch` / `fasterq-dump` / `pigz`.

No conda env files, no Nextflow.

## Inputs (in this folder)

- `data/yahara2021_run_accessions.tsv` — placeholder (8 runs from BioProject PRJDB10485). Refresh from SRA via `bash data/fetch_accessions.sh > data/yahara2021_run_accessions.tsv`. Columns: `run_id`, `platform`, `condition` (=donor/timepoint label that groups long+short reads), `notes`.
- `bin/070r1.slurm` — vpipe-style driver, subcommands: `download | flye | virsorter2 | checkv | vclust | coverage | validate`.
- `bin/validate_against_yahara.py` — ID-level (default) or mash-ANI (with `--yahara_fasta`) cross-check.
- `slurm/submit_array.sh` — array submitter. For per-run tasks (download) it iterates over rows; for per-condition tasks (flye/virsorter2/checkv/coverage) it iterates over distinct `condition` values.

## Quick start on HPC

```bash
rsync -av hpc/070-r1-virome-fastq-pipeline/ <hpc>:~/oralprob/070-r1/
ssh <hpc>
cd ~/oralprob/070-r1
chmod +x bin/070r1.slurm slurm/submit_array.sh data/fetch_accessions.sh

# 0. (one-time) Refresh the accession TSV from SRA
bash data/fetch_accessions.sh > data/yahara2021_run_accessions.tsv

# 1. Download all 8 runs (4 long + 4 short)
sbatch --array=1-8 -c 4 --mem 8G slurm/submit_array.sh download

# 2. Flye long-read assembly per condition (4 conditions, very heavy)
sbatch --array=1-4 -c 32 --mem 384G -t 36:00:00 slurm/submit_array.sh flye

# 3. Viral identification + CheckV per condition
sbatch --array=1-4 -c 16 --mem 64G slurm/submit_array.sh virsorter2
sbatch --array=1-4 -c 8  --mem 32G slurm/submit_array.sh checkv

# 4. vOTU dereplication across all conditions (single job)
sbatch -c 16 --mem 32G bin/070r1.slurm vclust

# 5. Per-condition short-read coverage vs the catalog
sbatch --array=1-4 -c 8 --mem 32G slurm/submit_array.sh coverage

# 6. Validation
sbatch -c 2 --mem 8G bin/070r1.slurm validate
```

## Resource estimates

- Flye on PromethION reads (~30 GB / condition): ~1 day on 32 cores, 384 GB RAM
- VirSorter2: ~6 h / condition on 16 cores
- CheckV: ~2 h / condition
- vclust + CoverM: ~1 h each
- Total wall: ~3 days with `%4` concurrency

## Validation gate

`bin/validate_against_yahara.py` writes `results/validation/yahara_vs_ours.tsv`. The ripple is "done" when one of the following holds:

- (default, ID-only) ≥ 80% of Yahara's high-quality phage IDs (CheckV `Medium-quality` or better when re-checked) match a vOTU representative in our catalog by sample-level overlap, OR
- (with `--yahara_fasta` available) ≥ 50% of Yahara's contigs match a vOTU at ≥ 95% ANI / 85% mash containment.

## Importing back into the contract layer

After `rsync results/vOTU/` and `results/coverage/` to `analyses/data/070-m2a-virome-processing/r1/`:

1. Update `analyses/070-m2a-virome-processing.qmd` to read the self-built catalog (`vOTU_votu.fasta` + `vOTU_votus.tsv` + per-condition `coverage/*.tsv`) when `params$source == "r1"`.
2. The CoverM `rpkm` column gives real abundance values — use `abundance_kind = "rpkm"` instead of MVP's `presence`.
3. The 080-r1 host-prediction bundle takes the resulting `vOTU_votu.fasta` as its primary input.
