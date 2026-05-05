# 080-r1 — Phage host prediction ensemble (HPC, vpipe-style)

**Goal**: Predict the bacterial host(s) of each oral phage with two independent methods, then combine via weighted voting. The MVP uses only Yahara 2021's published CAT taxonomy assignments — a single-method baseline. This ripple replaces it with **iPHoP + Sylph→IMG/VR**, both already provisioned in vpipe.

## vpipe integration

| Step | vpipe subcommand | apptainer image | DB |
|------|------------------|-----------------|-----|
| iPHoP host prediction | `anno_contig.slurm iphop ...` | `quay.io-biocontainers-iphop-1.4.1...img` | `DB_IPHOP` |
| Sylph IMG/VR sketch query | `anno_contig.slurm sylph imgvr fa ...` | `jinlongru-genomesearch...img` | `DB_SYLPH_IMGVR` |
| IMG/VR host taxonomy join | `bin/sylph_to_host.py` (native python3) | — | IMG/VR metadata TSV |
| Ensemble voting | `scripts/ensemble.py` | — | — |

Optional methods (NOT in stock vpipe — would need a new image + a new wrapper subcommand in `${VPIPEBIN}/anno_contig.slurm`):

- **spacepharer** — CRISPR-spacer-only host calls; high precision, low recall
- **RaFAH** — random forest, family-level

Skip them in MVP-r1; the ensemble script handles missing methods gracefully.

## Inputs

The pipeline runs on a vOTU FASTA. Two sources, in priority order:

1. **From 070-r1** (recommended): `../070-r1-virome-fastq-pipeline/results/vOTU/vOTU_votu.fasta`. Run 070-r1 first, then point 080-r1 at its output.
2. **MVP fallback**: `data/yahara_mvp_contigs.fasta`. Build via `data/extract_mvp_contigs.sh` — falls back to the 070-r1 catalog if present, otherwise fetches Yahara assemblies from DDBJ.

The IMG/VR metadata TSV (host-taxonomy column) is needed at:

```
${IMGVR_META:-${DBPATH}/virome_db/imgvr/IMGVR_all_Sequence_information-high_confidence.tsv}
```

If your IMG/VR install puts it elsewhere, set `IMGVR_META` in the environment.

## Quick start on HPC

```bash
rsync -av hpc/080-r1-host-prediction-ensemble/ <hpc>:~/oralprob/080-r1/
ssh <hpc>
cd ~/oralprob/080-r1
chmod +x bin/080r1.slurm slurm/submit_sequence.sh data/extract_mvp_contigs.sh

# Stage the input vOTU FASTA
CONTIGS=~/oralprob/070-r1/results/vOTU/vOTU_votu.fasta   # or the MVP fallback

# Submit the chained 4-step pipeline (jobs auto-depend via afterok)
bash slurm/submit_sequence.sh "$CONTIGS" oralprob 90
```

That submits four jobs (iphop, sylph_imgvr, sylph_to_host, ensemble) with the right `--dependency=afterok:` chains. Logs go to `logs/`. Final output:

```
results/ensemble/oralprob_ensemble_predictions.tsv
```

## Manual stepwise (if you want to gate each step)

```bash
sbatch -c 16 --mem 64G  bin/080r1.slurm iphop "$CONTIGS" oralprob 90
sbatch -c 8  --mem 32G  bin/080r1.slurm sylph_imgvr "$CONTIGS" oralprob
# wait for sylph_imgvr to finish before:
sbatch -c 2  --mem 8G   bin/080r1.slurm sylph_to_host oralprob
# wait for both iphop AND sylph_to_host:
sbatch -c 2  --mem 8G   bin/080r1.slurm ensemble oralprob
```

## Resource estimates

- iPHoP: ~6 h on 16 cores, 64 GB RAM (DB read ~120 GB, IO-heavy)
- Sylph IMG/VR query: ~30 min on 8 cores
- Sylph→host + ensemble: < 10 min each
- Total wall: ~1 day if iphop + sylph run in parallel

## Outputs (per `bin/080r1.slurm` step)

```
results/
├── iphop/
│   ├── out_<prefix>_iphop/                       # vpipe iphop raw output
│   └── iphop.tsv                                 # normalized: vOTU_id, host_genus, host_score, method
├── sylph/
│   ├── ac_<prefix>_sylph_imgvr.tsv               # vpipe sylph raw output
│   └── sylph_imgvr.tsv                           # normalized + host taxonomy joined
└── ensemble/
    ├── <prefix>_ensemble_predictions.tsv         # vOTU → top-1 host + score
    ├── <prefix>_ensemble_per_vOTU_evidence.tsv   # long-format evidence
    └── <prefix>_ensemble_method_pairwise_overlap.tsv
```

## Importing back into the contract layer

After `rsync results/ analyses/data/080-m2a-host-prediction/r1/`:

1. Update `analyses/080-m2a-host-prediction.qmd` to read each method's TSV and the ensemble TSV. Write into `phage_host_link` with:
   - `evidence_kind = "iphop" | "sylph_imgvr" | "ensemble"` (one row per method per vOTU; the schema already supports multiple rows per vOTU)
   - `confidence` = method's host_score (clipped to [0, 1])
   - `support_json` = the per-row metadata (raw method score, JSON of the `support_json` column from `_ensemble_predictions.tsv`)
2. Module 100 should prefer `evidence_kind = "ensemble"` rows when computing `compatibility_score`.

## Validation gate

The ripple is "done" when:

- ≥ 50% of input vOTUs receive at least one host call
- iPHoP and Sylph→host method-pairwise Jaccard ≥ 0.2 at vOTU-host level (i.e. they aren't completely uncorrelated)
- Top hosts in the ensemble include the canonical oral genera: Streptococcus, Veillonella, Prevotella, Fusobacterium, Porphyromonas, Treponema (sanity check against the MVP CAT-only output, which already shows these dominate)

## Adding spacepharer / RaFAH later

The ensemble script already accepts arbitrary number of method TSVs (`--iphop`, `--blast` (sylph here), or any extra `--<name>` you wire in). To add spacepharer:

1. Build / pull an apptainer image for spacepharer (e.g. via mmseqs2 base or a spacepharer-specific build)
2. Add a `COMMAND_PREFIX_SPACEPHARER` line to `${VPIPEBIN}/00-config.sh`
3. Add a `spacepharer)` subcommand to `${VPIPEBIN}/anno_contig.slurm`
4. Add a corresponding subcommand to this bundle's `bin/080r1.slurm`
5. Update the ensemble call to include the new TSV

Same recipe for RaFAH.
