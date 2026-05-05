# 030-r1 — Bacteriome fastq → species table (HPC, vpipe-style)

**Goal**: Re-process the Geng et al. 2024 cohorts' raw shotgun fastq from SRA into per-sample species relative-abundance tables (MetaPhlAn4 + biobakery DB), then validate against the published `Figure 01.Rdata` matrix that the MVP loads in `040-m1a-data-processing.qmd`.

This is the first ripple of `030-preprocessing-pipeline`. The MVP loads Geng's published table directly; this pipeline produces our own table from raw reads, which (a) lets us verify the MVP signal isn't an artifact of upstream choices and (b) gives us a reproducible pipeline for cohorts that don't ship pre-computed tables (the two Geng validation cohorts `PRJDB6966` and `PRJNA552294` have no published per-sample table — they are first-class targets for this ripple).

## vpipe integration

This bundle assumes `vpipe` is installed at `~/vpipe/` on the HPC. It uses:

- `${VPIPEBIN}/00-config.sh` — exports `DB_*` (incl. `DB_BIOBAKERY`) and `COMMAND_PREFIX_*` apptainer wrappers.
- `${VPIPEBIN}/01-functions.sh` — `log_*`, `require_*`, `get_files_for_task` helpers.
- MetaPhlAn4 + `merge_metaphlan_tables.py` — assumed available natively in the user env (this is the same convention `${VPIPEBIN}/anno_community.slurm metaphlan4` follows).
- `prefetch` / `fasterq-dump` / `pigz` — natively in PATH (sra-tools).

No conda envs, no Nextflow. Per-sample work is via SLURM array jobs that delegate to `bin/030r1.slurm`.

## Inputs (in this folder)

- `data/geng2024_sra_accessions.tsv` — 341 SRA runs across 8 bioprojects (PRJDB11203, PRJDB6966, PRJNA230363, PRJNA396840, PRJNA552294, PRJNA678453, PRJNA717815, PRJNA932553). Columns: `run_id`, `biosample_id`, `bioproject_id`, `group`, `country`, `sex`, `host_age`.
- `bin/030r1.slurm` — vpipe-style driver with subcommands `download | metaphlan4 | merge | validate`.
- `bin/long_format.py` — wide → long conversion helper.
- `bin/validate_against_geng.R` — diff per-(run, species) against Geng's `Figure 01.Rdata`.
- `slurm/submit_array.sh` — array-job orchestrator (one task per SRA run).

## Quick start on HPC

```bash
# On the HPC (assumes ~/vpipe is set up, conda envs installed, dbs provisioned)
rsync -av hpc/030-r1-bacteriome-fastq-pipeline/ <hpc>:~/oralprob/030-r1/

ssh <hpc>
cd ~/oralprob/030-r1
chmod +x bin/030r1.slurm slurm/submit_array.sh

# 1. Download all 341 SRA runs as an array (8 GB / run, IO-bound)
sbatch --array=1-341%20 -c 8 --mem 16G slurm/submit_array.sh download

# 2. Profile each run with MetaPhlAn4 (after download array completes)
sbatch --array=1-341%20 -c 16 --mem 32G slurm/submit_array.sh metaphlan4

# 3. Merge per-bioproject (one sbatch per BP — ~9 of them; quick)
for bp in PRJDB11203 PRJDB6966 PRJNA230363 PRJNA396840 PRJNA552294 \
          PRJNA678453 PRJNA717815 PRJNA932553; do
  sbatch -c 2 --mem 8G bin/030r1.slurm merge $bp
done

# 4. Validate one cohort against Geng's published matrix
sbatch -c 2 --mem 8G bin/030r1.slurm validate PRJDB11203
```

## Resource estimates

- Download: ~30 min / run, ~10 GB SRA storage / run, network-bound
- MetaPhlAn4: ~2 CPU-h / run on 16 cores
- Total wall: ~1 day with `%20` array concurrency

## Validation gate

A bioproject passes 030-r1 when `bin/validate_against_geng.R` reports:

- median |delta| (per-run-per-species relative-abundance gap vs Geng) < 0.05
- < 10% of Geng's species are missing from our profile

If both hold, the MVP signal is robust to pipeline reproduction. The script writes `results/species_abundance/<BP>_validation.tsv` per (run, species) for diagnostics.

## Importing back into the contract layer

After the MetaPhlAn tables come back to your local repo as `analyses/data/030-preprocessing-pipeline/species_abundance/<BP>_long.tsv`:

1. Update `analyses/030-preprocessing-pipeline.qmd` to read the long TSVs and emit `geng_vs_ours.tsv` per cohort.
2. Optional 040-r2: switch its source from `Figure 01.Rdata` to our long-format TSVs. This makes the bacteriome substrate fully reproducible from raw fastq.

## What's deliberately NOT here

- ❌ `env/*.yml` conda envs — vpipe's apptainer images are the source of truth.
- ❌ Nextflow DSL2 modules — vpipe-style SLURM is the primary path. If you want Nextflow at scale (>1000 runs), wrap `bin/030r1.slurm` in a single Nextflow process — Nextflow becomes a thin orchestrator that just `sbatch`es this script.
