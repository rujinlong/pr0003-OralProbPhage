# Expected output layout

After all four steps finish (via `slurm/submit_sequence.sh` or manual sbatch chain):

```
results/
├── iphop/
│   ├── out_<prefix>_iphop/                       # raw vpipe iphop output (Host_prediction_to_genus_m90.csv etc.)
│   └── iphop.tsv                                 # normalized: vOTU_id, host_genus, host_score, method
├── sylph/
│   ├── ac_<prefix>_sylph_imgvr.tsv               # raw vpipe sylph output (ANI hits to IMG/VR)
│   └── sylph_imgvr.tsv                           # normalized + IMG/VR host taxonomy joined
└── ensemble/
    ├── <prefix>_ensemble_predictions.tsv         # per-vOTU top-1 host + ensemble_score
    ├── <prefix>_ensemble_per_vOTU_evidence.tsv   # long-format evidence (one row per (vOTU, method))
    └── <prefix>_ensemble_method_pairwise_overlap.tsv
```

## Importing back into the contract layer

After `rsync results/ analyses/data/080-m2a-host-prediction/r1/`:

1. Update `analyses/080-m2a-host-prediction.qmd` to read each TSV and write `phage_host_link` rows with:
   - `evidence_kind = "iphop" | "sylph_imgvr" | "ensemble"`
   - `confidence` = method's `host_score` clipped to [0,1] (or `ensemble_score` for the ensemble row)
   - `support_json` = the `support_json` column from `_ensemble_predictions.tsv` for ensemble rows; raw row JSON for per-method rows
2. Module 100 should prefer `evidence_kind = "ensemble"` rows when computing `compatibility_score`.

## Validation gate

The 080-r1 ripple is "done" when:

- ≥ 50% of input vOTUs receive at least one host call (across all methods)
- iPHoP and Sylph→host method-pairwise Jaccard ≥ 0.2 at vOTU-host level (per `_method_pairwise_overlap.tsv`)
- Top hosts include the canonical oral genera: Streptococcus, Veillonella, Prevotella, Fusobacterium, Porphyromonas, Treponema (sanity check against MVP CAT-only output, which already shows these dominate)

## Adding spacepharer / RaFAH later

The ensemble script accepts arbitrary methods (`--iphop`, `--blast` (sylph here), `--rafah`, `--spacepharer`). To add one of the optional methods:

1. Build / pull an apptainer image
2. Add `COMMAND_PREFIX_<NAME>` to `${VPIPEBIN}/00-config.sh`
3. Add a wrapper subcommand in `${VPIPEBIN}/anno_contig.slurm`
4. Add a delegating subcommand in this bundle's `bin/080r1.slurm`
5. Pass the new TSV to `scripts/ensemble.py` in the `ensemble` step
