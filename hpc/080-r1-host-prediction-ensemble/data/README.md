# Inputs

## Phage FASTA

The pipeline expects one of these files:

- `data/vOTU_catalog.fasta` — typically the 070-r1 catalog. Drop in via:
   ```bash
   cp ../../analyses/data/070-m2a-virome-processing/r1/vOTU/vOTU_catalog.fasta data/
   ```
- `data/yahara_mvp_contigs.fasta` — MVP source. Build via:
   ```bash
   bash data/extract_mvp_contigs.sh
   ```
   This script re-extracts the actual VirSorter contigs from Yahara's
   MOESM7 sheet IDs against the 070-r1 assembly (or the published Yahara
   GenBank accessions). If 070-r1 hasn't run yet, the script fetches
   contigs from the BioProject's deposited assemblies on NCBI.

## Reference databases (NOT shipped — install once)

```
/scratch/$USER/dbs/iphop-db        # ~120 GB
/scratch/$USER/dbs/imgvr-v4        # ~30 GB
/scratch/$USER/dbs/spacepharer/    # CRISPR spacers from oral bacterial genomes (eHOMD)
/scratch/$USER/dbs/rafah-models/   # RaFAH random-forest models (~5 GB)
```

Install instructions:

```bash
# iPHoP database (≈120 GB)
iphop download --db_dir /scratch/$USER/dbs/iphop-db

# IMG/VR v4 (request access from JGI; place fasta + metadata under)
mkdir -p /scratch/$USER/dbs/imgvr-v4
# manual download from https://genome.jgi.doe.gov/portal/IMG_VR/IMG_VR.home.html

# spacepharer DB — build from eHOMD genomes
spacepharer createsetdb /scratch/$USER/dbs/eHOMD_genomes /scratch/$USER/dbs/spacepharer/eHOMD nucl --batchsize 4

# RaFAH — pretrained models
git clone https://github.com/felipehcoutinho/RaFAH.git /scratch/$USER/dbs/rafah-models
```
