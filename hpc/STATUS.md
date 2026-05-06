# OralProbPhage HPC ripple status

## 2026-05-05T22:35 — prerequisites: PARTIAL ✓ / PATH-tools blocker

Verified on `hpc-build01.scidom.de` (login) and `cpusrv32.scidom.de` (compute, via diagnostic jobs 36057360/36057362/36057366):

### OK
- vpipe at `~/vpipe/bin/` — config + 01-functions + anno_contig.slurm all present.
- DB paths exported correctly when `source ~/vpipe/bin/00-config.sh` runs at top level of a `#!/bin/bash` script: `DB_BIOBAKERY`, `DB_VIRSORTER2`, `DB_CHECKV`, `DB_IPHOP`, `DB_SYLPH_IMGVR` all present.
- Apptainer 1.4.2 + all needed singularity images (flye.sif, virsorter2, checkv, iphop, sylph, sra-tools) under `/ictstr01/project/genomics/jru/singularity/`.
- SLURM partition `cpu_p` healthy; max walltime 3-00:00:00.
- `/ictstr01/scratch/users/jinlong.ru/` is writable on compute, but per user direction we use `~/project/pr0003-OralProbPhage/oralprob/` (lustre /ictstr01) instead.

### Original blockers (resolved by patches below)
- `prefetch`/`fasterq-dump`/`metaphlan`/`merge_metaphlan_tables.py` not in compute-node PATH.
- sra-tools 2.11 image ships unconfigured (`vdb-config --interactive` prompt).

## 2026-05-06T06:35 — patches applied (per user direction)

Per user: "use env humann4 (you might need to run `ici` first). put to `~/project/pr0003-OralProbPhage/oralprob`. Refer to `~/vpipe/bin/anno_community.slurm` for humann4 usage."

### Bundle script edits (in this repo only — vpipe untouched)
- `hpc/030-r1-bacteriome-fastq-pipeline/bin/030r1.slurm`
- `hpc/070-r1-virome-fastq-pipeline/bin/070r1.slurm`
  - Added at top of script after `mkdir -p "$OUTDIR"`:
    ```bash
    if ! command -v metaphlan >/dev/null 2>&1; then
      source "${HOME}/miniconda3/etc/profile.d/conda.sh"
      conda activate humann4
    fi
    SRA_IMG="${SRA_IMG:-/ictstr01/.../singularity/depot.galaxyproject.org-singularity-sra-tools-2.11.0--pl5321ha49a11a_3.img}"
    NCBI_CFG="${NCBI_CFG:-${BUNDLE_DIR}/../../oralprob/.ncbi/user-settings.mkfg}"
    SRA_RUN="apptainer exec -B /ictstr01 -B /localscratch -B /tmp --env VDB_CONFIG=${NCBI_CFG} ${SRA_IMG}"
    ```
  - Replaced `prefetch ...` with `$SRA_RUN prefetch ...` and same for `fasterq-dump`.
- `hpc/030-r1-bacteriome-fastq-pipeline/slurm/submit_array.sh`
- `hpc/070-r1-virome-fastq-pipeline/slurm/submit_array.sh`
  - Switched `BUNDLE_DIR` resolution to prefer `$SLURM_SUBMIT_DIR` (sbatch copies the script to `/var/spool/slurmd/job<id>/slurm_script`, breaking `BASH_SOURCE[0]`).

### One-time setup (also in `oralprob/`)
- `oralprob/.ncbi/user-settings.mkfg` — project-local sra-toolkit config (sets `/config/default = "true"` + a stub GUID; pointed at via `VDB_CONFIG` env var into the apptainer container).

### NOT in this repo (consider committing)
- The above patches live on the HPC at `/ictstr01/project_copy/.../pr0003-OralProbPhage/`. Once you `rsync` back, please commit them so the laptop-side bundle and HPC stay in sync. Or replace this with a cleaner setup later (e.g. ship a wrapper image that bundles sra-tools + humann4).

### Side effect to clean up
- `vdb-config --restore-defaults` (run during diagnosis before I picked the project-local config approach) created `~/.ncbi/user-settings.mkfg`. Harmless but it's outside this project. You may delete it: `rm -rf ~/.ncbi`.

## 2026-05-06T06:40 — 030-r1 smoke download PASS (52 s)

DRR285661 downloaded successfully (193M paired fastq.gz), conda+apptainer+sra-tools wrapper verified end-to-end.

## 2026-05-06T06:43 — 030-r1 full chain submitted

```
download=36058449   array=1-341%20   c=8  mem=16G  t=8h    (running 20 parallel)
metaphlan=36058450  array=1-341%20   c=16 mem=32G  t=12h   (afterok download)
merge=36058474      single           c=2  mem=8G   t=2h    (afterok metaphlan)
validate=36058475   single           c=2  mem=8G   t=1h    (afterok merge — 6 modeling cohorts)
```

Tracker: `oralprob/jobs.tracker`. Logs: `oralprob/logs/030r1-*-<JID>_<task>.{out,err}`. Estimated wall: ~1 day per the brief.

## 2026-05-06T06:35 — 070-r1 BLOCKED on data discrepancy

`hpc/070-r1-virome-fastq-pipeline/data/fetch_accessions.sh` (run with `pysradb` from the `py3` conda env) returned **188 runs from `PRJDB10485`, all `Illumina MiSeq`**. The bundle was designed against an expected 8 runs (4 PromethION long-read + 4 paired Illumina, grouped by `t1`/`t2`/`t3`/`t4`). The placeholder accessions DRR251800–DRR251803 (marked PromethION) **do not appear in `PRJDB10485` at all** when queried via pysradb.

Concretely:
- `pysradb metadata PRJDB10485 --detailed | awk -F'\t' 'NR>1 {print $19}' | sort -u` → only `Illumina MiSeq`.
- `pysradb metadata DRR251800 ... DRR251807` → all return empty result rows (no metadata).

Implications:
- Without long-reads, **Flye is moot** for this bundle. The whole long-read assembly path doesn't apply.
- Library names are all `-` (no t1/t2/t3/t4 grouping), so the per-condition logic in 070-r1 doesn't have natural keys.

Original placeholder TSV is preserved at `data/yahara2021_run_accessions.tsv.placeholder.bak`. The freshly-fetched (188-row) TSV is **not** written to disk yet — I held back so you can decide.

### Question for you
Three ways forward:

1. **Different BioProject for the long-reads** — if the Yahara 2021 PromethION data is under a sister project (not PRJDB10485), tell me the project ID and I'll re-fetch and re-shape.
2. **Skip 070-r1's Flye path entirely** — use the 188 Illumina runs only, switch to a short-read assembly (SPAdes/MEGAHIT) on a chosen subset (e.g. ≤4 samples for the MVP). This requires a small bundle change (replace Flye with SPAdes call).
3. **Use Option B fallback (MVP contigs)** — `bash data/extract_mvp_contigs.sh` fetches Yahara assemblies from DDBJ; we skip 070-r1 entirely and just feed `data/yahara_mvp_contigs.fasta` directly to 080-r1. Faster and matches the brief's "Option B".

`080-r1` is held until you decide.

## 2026-05-06T06:35 — 080-r1 deferred

Deferred until 070-r1 emits a vOTU FASTA (or you choose Option B above).
