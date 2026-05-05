# HPC Agent Brief — OralProbPhage ripple execution

You are Claude Code running on an HPC login node. Your job is to execute three
"ripple" bioinformatics pipelines for the **OralProbPhage** project (oral
probiotic + phage co-therapy data layer). Each pipeline is a self-contained
bundle under `hpc/`; submit its SLURM jobs in the right order, monitor them,
run the validation gate, and report results.

**Do NOT modify pipeline code unless validation fails AND root cause is in the
script itself.** The bundles were designed and reviewed on the laptop side;
your job is execution, not redesign.

---

## TL;DR

| Bundle | Wall time (est) | Key resources | Depends on |
|---|---|---|---|
| `030-r1-bacteriome-fastq-pipeline` | ~1 day | 341 SRA runs × MetaPhlAn4 | nothing |
| `070-r1-virome-fastq-pipeline` | ~3 days | 4 Flye assemblies (32 cpu, 384 GB RAM, 36 h each) | nothing |
| `080-r1-host-prediction-ensemble` | ~12 h | iPHoP (DB read 120 GB) + Sylph IMG/VR | 070-r1 catalog (or MVP fallback) |

030-r1 and 070-r1 are independent — submit both in parallel. 080-r1 starts
after 070-r1 emits `vOTU/vOTU_votu.fasta` (or use the MVP-fallback FASTA if
you decide to start 080 early).

---

## Step 0 — Prerequisites verification (BEFORE any sbatch)

Run these checks and STOP if anything fails. Don't try to fix vpipe itself —
ping the user via the report file.

```bash
# 1. vpipe is installed and config sourceable
test -f "${HOME}/vpipe/bin/00-config.sh" && echo "vpipe OK" || echo "FAIL: no vpipe"
test -f "${HOME}/vpipe/bin/01-functions.sh" && echo "functions OK" || echo "FAIL"
test -f "${HOME}/vpipe/bin/anno_contig.slurm" && echo "anno_contig OK" || echo "FAIL"

# 2. Source the config and verify key DB vars exist
export VPIPEBIN="${HOME}/vpipe/bin"
source "${VPIPEBIN}/00-config.sh" 2>/dev/null
for v in DB_BIOBAKERY DB_VIRSORTER2 DB_CHECKV DB_IPHOP DB_SYLPH_IMGVR; do
  if [ -e "${!v:-}" ]; then echo "$v OK ($(printf '%s' "${!v}"))"
  else echo "FAIL: $v not set or path missing"; fi
done

# 3. Native tools in PATH (vpipe expects these)
for cmd in prefetch fasterq-dump pigz metaphlan merge_metaphlan_tables.py sbatch squeue sacct; do
  command -v "$cmd" >/dev/null && echo "$cmd OK" || echo "FAIL: $cmd not in PATH"
done

# 4. SLURM partition sanity
sinfo -h -p cpu_p --format='%P %a %D %T' 2>&1 | head -3
```

If any of these fail, write a one-line summary into `hpc/STATUS.md` like:

```
[BLOCKED] 2026-05-05T12:00 — DB_IPHOP missing at $DBPATH/viroprofiler/iphop/Jun_2025_pub_rw. Run `iphop download` before retrying 080-r1.
```

Then stop. Don't try to install databases — that's a multi-hour + curated
step the user owns.

---

## Step 1 — Pick a scratch directory

`$HOME` typically has quotas; outputs (especially fastq + Flye work dirs) go to
scratch:

```bash
export OUTDIR_BASE="/scratch/$USER/oralprob"
mkdir -p "$OUTDIR_BASE"/{030-r1,070-r1,080-r1}
```

Each bundle's SLURM script honours `$OUTDIR` to override `bundle/results/`.
Pass it via `--export=ALL,OUTDIR=...` on sbatch, or `export` it before
calling `bin/<bundle>.slurm`. Examples below.

---

## Step 2 — Execute the three ripples

### 2A. 030-r1 — bacteriome (parallel with 2B)

Read `hpc/030-r1-bacteriome-fastq-pipeline/README.md` for full context. Then:

```bash
cd "${REPO_ROOT}/hpc/030-r1-bacteriome-fastq-pipeline"
chmod +x bin/030r1.slurm slurm/submit_array.sh

OUTDIR="$OUTDIR_BASE/030-r1"
N_RUNS=$(($(wc -l < data/geng2024_sra_accessions.tsv) - 1))   # =341

# Chain the entire pipeline with afterok dependencies. Capture JobIDs.
JOB_DL=$(OUTDIR=$OUTDIR sbatch --parsable \
    --array=1-${N_RUNS}%20 -c 8 --mem 16G \
    slurm/submit_array.sh download)
echo "030-r1 download array: $JOB_DL"

JOB_MP=$(OUTDIR=$OUTDIR sbatch --parsable \
    --dependency=afterok:$JOB_DL \
    --array=1-${N_RUNS}%20 -c 16 --mem 32G \
    slurm/submit_array.sh metaphlan4)
echo "030-r1 metaphlan array: $JOB_MP"

# Merge per bioproject (one job that loops over BPs)
JOB_MERGE=$(OUTDIR=$OUTDIR sbatch --parsable \
    --dependency=afterok:$JOB_MP -c 2 --mem 8G -t 2:00:00 \
    --wrap='set -e; cd '"$PWD"'; for bp in PRJDB11203 PRJDB6966 PRJNA230363 PRJNA396840 PRJNA552294 PRJNA678453 PRJNA717815 PRJNA932553; do bash bin/030r1.slurm merge $bp; done')
echo "030-r1 merge: $JOB_MERGE"

# Validate the 6 modeling cohorts (the 2 validation cohorts have no Geng baseline)
JOB_VAL=$(OUTDIR=$OUTDIR sbatch --parsable \
    --dependency=afterok:$JOB_MERGE -c 2 --mem 8G -t 1:00:00 \
    --wrap='set -e; cd '"$PWD"'; for bp in PRJDB11203 PRJNA230363 PRJNA396840 PRJNA678453 PRJNA717815 PRJNA932553; do bash bin/030r1.slurm validate $bp || true; done')
echo "030-r1 validate: $JOB_VAL"

# Save the chain in a tracker
mkdir -p "$OUTDIR"
echo "030-r1 download=$JOB_DL metaphlan=$JOB_MP merge=$JOB_MERGE validate=$JOB_VAL" \
  >> "$OUTDIR_BASE/jobs.tracker"
```

### 2B. 070-r1 — virome (parallel with 2A)

```bash
cd "${REPO_ROOT}/hpc/070-r1-virome-fastq-pipeline"
chmod +x bin/070r1.slurm slurm/submit_array.sh data/fetch_accessions.sh

OUTDIR="$OUTDIR_BASE/070-r1"

# 0. Refresh the accession TSV from SRA (placeholder accessions ship with the bundle)
if [ -f data/yahara2021_run_accessions.tsv ] && \
   awk 'NR>1 && $4 ~ /placeholder/' data/yahara2021_run_accessions.tsv | grep -q .; then
  echo "Refreshing accession TSV from SRA..."
  bash data/fetch_accessions.sh > data/yahara2021_run_accessions.tsv.new && \
    mv data/yahara2021_run_accessions.tsv.new data/yahara2021_run_accessions.tsv
fi
N_RUNS=$(($(wc -l < data/yahara2021_run_accessions.tsv) - 1))      # ~8
N_COND=$(awk -F'\t' 'NR>1 {print $3}' data/yahara2021_run_accessions.tsv | sort -u | wc -l)  # ~4

# Chain: download -> flye -> virsorter2 -> checkv -> vclust -> coverage -> validate
JOB_DL=$(OUTDIR=$OUTDIR sbatch --parsable \
    --array=1-${N_RUNS} -c 4 --mem 8G -t 8:00:00 \
    slurm/submit_array.sh download)

JOB_FLYE=$(OUTDIR=$OUTDIR sbatch --parsable \
    --dependency=afterok:$JOB_DL \
    --array=1-${N_COND} -c 32 --mem 384G -t 36:00:00 -p cpu_p -q cpu_normal \
    slurm/submit_array.sh flye)

JOB_VS2=$(OUTDIR=$OUTDIR sbatch --parsable \
    --dependency=afterok:$JOB_FLYE \
    --array=1-${N_COND} -c 16 --mem 64G -t 12:00:00 \
    slurm/submit_array.sh virsorter2)

JOB_CV=$(OUTDIR=$OUTDIR sbatch --parsable \
    --dependency=afterok:$JOB_VS2 \
    --array=1-${N_COND} -c 8 --mem 32G -t 6:00:00 \
    slurm/submit_array.sh checkv)

JOB_VCLUST=$(OUTDIR=$OUTDIR sbatch --parsable \
    --dependency=afterok:$JOB_CV \
    -c 16 --mem 32G -t 6:00:00 \
    bin/070r1.slurm vclust)

JOB_COV=$(OUTDIR=$OUTDIR sbatch --parsable \
    --dependency=afterok:$JOB_VCLUST \
    --array=1-${N_COND} -c 8 --mem 32G -t 6:00:00 \
    slurm/submit_array.sh coverage)

JOB_VAL70=$(OUTDIR=$OUTDIR sbatch --parsable \
    --dependency=afterok:$JOB_COV -c 2 --mem 8G -t 1:00:00 \
    bin/070r1.slurm validate)

echo "070-r1 download=$JOB_DL flye=$JOB_FLYE virsorter2=$JOB_VS2 checkv=$JOB_CV vclust=$JOB_VCLUST coverage=$JOB_COV validate=$JOB_VAL70" \
  >> "$OUTDIR_BASE/jobs.tracker"
```

### 2C. 080-r1 — host prediction (depends on 070-r1)

The bundle ships a `slurm/submit_sequence.sh` that builds the 4-step chain
itself. You only need to wait for 070-r1's `vOTU_votu.fasta` to exist, then
launch.

```bash
cd "${REPO_ROOT}/hpc/080-r1-host-prediction-ensemble"
chmod +x bin/080r1.slurm bin/sylph_to_host.py scripts/ensemble.py slurm/submit_sequence.sh

OUTDIR="$OUTDIR_BASE/080-r1"
CONTIGS_R1="$OUTDIR_BASE/070-r1/vOTU/vOTU_votu.fasta"

# Option A (recommended): launch 080 with --dependency on 070's vclust job, so
# everything is queued upfront. Read $JOB_VCLUST from jobs.tracker.
JOB_VCLUST=$(awk '/^070-r1/ {for(i=1;i<=NF;i++) if($i ~ /^vclust=/) {split($i,a,"="); print a[2]}}' \
    "$OUTDIR_BASE/jobs.tracker")
echo "Will start 080-r1 after 070-r1 vclust ($JOB_VCLUST)"

# We can't use submit_sequence.sh directly (it doesn't take an --dependency
# arg), so submit the chain manually with the upstream wait baked into iphop+sylph.
J_IPHOP=$(OUTDIR=$OUTDIR sbatch --parsable \
    --dependency=afterok:$JOB_VCLUST \
    -c 16 --mem 64G -t 12:00:00 \
    bin/080r1.slurm iphop "$CONTIGS_R1" oralprob 90)

J_SYLPH=$(OUTDIR=$OUTDIR sbatch --parsable \
    --dependency=afterok:$JOB_VCLUST \
    -c 8 --mem 32G -t 4:00:00 \
    bin/080r1.slurm sylph_imgvr "$CONTIGS_R1" oralprob)

J_S2H=$(OUTDIR=$OUTDIR sbatch --parsable \
    --dependency=afterok:$J_SYLPH \
    -c 2 --mem 8G -t 1:00:00 \
    bin/080r1.slurm sylph_to_host oralprob)

J_ENS=$(OUTDIR=$OUTDIR sbatch --parsable \
    --dependency=afterok:$J_IPHOP:$J_S2H \
    -c 2 --mem 8G -t 1:00:00 \
    bin/080r1.slurm ensemble oralprob)

echo "080-r1 iphop=$J_IPHOP sylph=$J_SYLPH sylph_to_host=$J_S2H ensemble=$J_ENS" \
  >> "$OUTDIR_BASE/jobs.tracker"
```

**Option B**: If 070-r1 catalog isn't available and you want to start 080-r1
early on the MVP-fallback FASTA, run `bash data/extract_mvp_contigs.sh` first
(fetches Yahara assemblies from DDBJ if 070-r1 absent), then point `CONTIGS`
at `data/yahara_mvp_contigs.fasta` and use `slurm/submit_sequence.sh` directly:

```bash
bash slurm/submit_sequence.sh data/yahara_mvp_contigs.fasta oralprob_mvp 90
```

---

## Step 3 — Monitoring without burning context

**Don't sleep-poll.** SLURM takes hours to days. Use this status helper, and
only re-check at long intervals (every ~30 min, or after a known long step).

```bash
# Quick status of every JobID in the tracker
status_all() {
  awk '{for(i=2;i<=NF;i++) {split($i,a,"="); if(a[2]) print a[2]}}' \
      "$OUTDIR_BASE/jobs.tracker" | sort -u | while read jid; do
    sacct -j "$jid" -X --format=JobID%14,JobName%20,State%12,Elapsed,ExitCode -n -P \
        | head -1
  done
}
status_all
```

For long waits (e.g. Flye 36 h), prefer:

- **Background monitor** — submit one extra job with `--dependency=afterany:$LAST_JOB --wrap='echo done'` and just wait on its sacct state. SLURM does the work; you don't.
- **Status batch** — write a `status.sh` script that prints a one-line summary, and run it on demand (`bash status.sh`) when you check back.

Avoid `while ...; do squeue ...; sleep 60; done` loops in the foreground —
they burn context and don't give the user useful feedback.

---

## Step 4 — Validation (the actual gate)

A SLURM job exiting with code 0 means "the command ran". The ripple is only
"done" when validation passes. Each bundle's last job IS the validation
script:

| Ripple | Validation pass criterion | Read |
|---|---|---|
| 030-r1 | `validate_against_geng.R` reports median |delta| < 0.05 AND < 10% missing species, **for ≥ 4 of 6 modeling cohorts** | `030-r1-bacteriome-fastq-pipeline/expected_output/README.md` |
| 070-r1 | `validate_against_yahara.py` recovers ≥ 80% of Yahara phageIDs (ID-level), OR ≥ 50% at 95% ANI if `YAHARA_FASTA` available | `070-r1-virome-fastq-pipeline/expected_output/README.md` |
| 080-r1 | ≥ 50% of input vOTUs get a host call; iPHoP/sylph Jaccard ≥ 0.2; top hosts include Streptococcus / Veillonella / Prevotella / Fusobacterium / Porphyromonas / Treponema | `080-r1-host-prediction-ensemble/expected_output/README.md` |

Validation outputs land at:

- `030-r1: $OUTDIR/030-r1/species_abundance/<BP>_validation.tsv`
- `070-r1: $OUTDIR/070-r1/validation/yahara_vs_ours.tsv`
- `080-r1: $OUTDIR/080-r1/ensemble/oralprob_ensemble_predictions.tsv` + `_method_pairwise_overlap.tsv`

If validation fails: don't blindly retry. Investigate.

---

## Step 5 — Failure handling

When a SLURM job exits non-zero, OR validation fails, follow this triage:

1. **Find the failing job's log**:
   ```bash
   sacct -j $JOBID -X --format=JobID,State,ExitCode,DerivedExitCode,Reason -n -P
   # Logs default to slurm-<jobid>.out / .err in the directory you sbatch'ed from
   ls -lt slurm-${JOBID}* logs/*${JOBID}*  2>/dev/null | head
   ```

2. **Classify the error** by reading the last 50 lines of stderr:
   - **OOM** (`oom-kill`, `Out of memory`): bump `--mem` and resubmit only the failed array tasks
   - **Timeout** (`DUE TO TIME LIMIT`): bump `-t` and resubmit
   - **Database missing**: stop, write blocker to `STATUS.md`, ask the user
   - **Network / SRA flake** (download step only): re-submit the failed array tasks. SRA throttles routinely
   - **Tool error** (e.g. CheckV crash on a malformed contig): isolate the contig, decide whether to skip or fix

3. **Re-submit only the failed tasks** (for arrays):
   ```bash
   FAILED=$(sacct -j $JOBID -X --state=F -n -P --format=JobID | grep -oE '_[0-9]+$' | tr -d _ | tr '\n' ',' | sed 's/,$//')
   echo "Re-submitting failed array tasks: $FAILED"
   sbatch --array=$FAILED -c 16 --mem 64G slurm/submit_array.sh metaphlan4   # adjust resource as needed
   ```

4. **Never modify** `${VPIPEBIN}/...` files. If a bug is in vpipe, write it to `STATUS.md` and stop.

5. **Modifying the bundle's own scripts** is allowed only if the failure is clearly inside `bin/<bundle>.slurm` (e.g. a typo in a path resolution, an argument forwarding bug). Edit, commit a note in `STATUS.md` describing the fix, and resume.

---

## Step 6 — Reporting

Write progress to `hpc/STATUS.md` at:

- Step 0 finish (prereq verification result)
- Each ripple's submission (job IDs)
- Each validation pass/fail
- Any blocker requiring user action

Format:

```markdown
# OralProbPhage HPC ripple status

## 2026-05-05T14:30 — prerequisites verified
All vpipe paths + DBs present. Scratch at /scratch/jru/oralprob/.

## 2026-05-05T14:35 — 030-r1 submitted
download=12345678, metaphlan=12345679, merge=12345680, validate=12345681

## 2026-05-05T14:36 — 070-r1 submitted
download=..., flye=..., virsorter2=..., checkv=..., vclust=..., coverage=..., validate=...

## 2026-05-06T08:12 — 030-r1 download array: 339/341 ok, 2 SRA timeouts
Re-submitted DRR285661 + DRR251000 — array 12345699.

## 2026-05-08T03:00 — 070-r1 PASS
Yahara recovery 87.4% (gate: 80%). vOTU catalog: 4,213 vOTUs from 4 samples.

## 2026-05-08T15:00 — 080-r1 PASS
1,847 / 4,213 vOTUs got host calls (44%). Top hosts: Streptococcus (511), Veillonella (203), ...
```

When all three ripples PASS:

1. Append a final summary block to `STATUS.md`
2. Tar the deliverables for download:

```bash
tar czf "$OUTDIR_BASE/oralprob_ripple_results.tar.gz" \
    "$OUTDIR_BASE/030-r1/species_abundance" \
    "$OUTDIR_BASE/070-r1/vOTU" "$OUTDIR_BASE/070-r1/coverage" "$OUTDIR_BASE/070-r1/validation" \
    "$OUTDIR_BASE/080-r1/iphop" "$OUTDIR_BASE/080-r1/sylph" "$OUTDIR_BASE/080-r1/ensemble"
ls -lh "$OUTDIR_BASE/oralprob_ripple_results.tar.gz"
echo "READY FOR DOWNLOAD: $OUTDIR_BASE/oralprob_ripple_results.tar.gz" >> hpc/STATUS.md
```

The user will `rsync` this back to the laptop and re-import into the
contract layer. You don't need to do that step.

---

## Quick reference

```bash
# Bundle paths
hpc/030-r1-bacteriome-fastq-pipeline/    # MetaPhlAn4 species profiling
hpc/070-r1-virome-fastq-pipeline/        # Flye + VirSorter2 + CheckV + vclust
hpc/080-r1-host-prediction-ensemble/     # iPHoP + Sylph IMG/VR + ensemble

# vpipe integration
${VPIPEBIN}/anno_contig.slurm     virsorter2|checkv|iphop|sylph|vclust_votu  ...
${VPIPEBIN}/00-config.sh          DB_*, COMMAND_PREFIX_*
${VPIPEBIN}/01-functions.sh       log_*, require_*

# Resource cheat sheet (per-task)
#   download    : -c 4-8   --mem 8-16G    -t 8h
#   metaphlan4  : -c 16    --mem 32G      -t 4h
#   flye        : -c 32    --mem 384G     -t 36h
#   virsorter2  : -c 16    --mem 64G      -t 12h
#   checkv      : -c 8     --mem 32G      -t 6h
#   vclust      : -c 16    --mem 32G      -t 6h
#   coverage    : -c 8     --mem 32G      -t 6h
#   iphop       : -c 16    --mem 64G      -t 12h
#   sylph_imgvr : -c 8     --mem 32G      -t 4h
#   ensemble    : -c 2     --mem 8G       -t 1h

# Useful one-liners
sacct -j $JOBID -X --format=JobID,JobName,State,Elapsed,ExitCode -n -P
squeue -u $USER -t RUNNING,PENDING --format='%i %j %T %M %l %r'
scancel <JOBID>           # if you really need to cancel
scancel -t PD -u $USER    # cancel only pending jobs
```

---

## What you should NOT do

- Don't `pip install`, `conda install`, or otherwise touch the vpipe environment.
- Don't run jobs on the login node — always `sbatch`.
- Don't write to `$HOME` for big outputs (use `$OUTDIR_BASE` on scratch).
- Don't push to GitHub — there's no remote auth set up here.
- Don't `git commit` results back into the project — the user will do that on the laptop.
- Don't run anything labeled "Optional" in the READMEs unless explicitly told to (spacepharer, RaFAH, mash-ANI validation with full Yahara FASTA).
- Don't poll-loop in the foreground for hours. Submit, sleep on dependencies, check back at long intervals.
