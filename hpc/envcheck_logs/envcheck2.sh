#!/bin/bash
#SBATCH -p cpu_p
#SBATCH -q cpu_normal
#SBATCH -c 1
#SBATCH --mem 1G
#SBATCH -t 00:05:00
#SBATCH -J envcheck2
#SBATCH -o /ictstr01/project_copy/genomics/jru/project/pr0003-OralProbPhage/hpc/envcheck_logs/envcheck2-%j.out
#SBATCH -e /ictstr01/project_copy/genomics/jru/project/pr0003-OralProbPhage/hpc/envcheck_logs/envcheck2-%j.err

set +u
echo "=== node: $(hostname) ==="
echo "=== shell ==="
echo "BASH_VERSION=${BASH_VERSION:-<unset>}"
echo "SHELL=$SHELL"
echo "=== PATH ==="
echo "$PATH"
echo "=== sourcing vpipe (real bash script) ==="
export VPIPEBIN="${HOME}/vpipe/bin"
source "${VPIPEBIN}/00-config.sh" 2>&1 | tail -5
echo "DB_BIOBAKERY=${DB_BIOBAKERY:-<unset>}"
echo "DB_VIRSORTER2=${DB_VIRSORTER2:-<unset>}"
echo "COMMAND_PREFIX_FLYE=${COMMAND_PREFIX_FLYE:-<unset>}"
echo "=== conda activate humann4 ==="
source ${HOME}/miniconda3/etc/profile.d/conda.sh
conda activate humann4
which metaphlan prefetch fasterq-dump merge_metaphlan_tables.py 2>&1 | head -10
echo "=== sra-tools apptainer image ==="
SRA_IMG=/ictstr01/project/genomics/jru/singularity/depot.galaxyproject.org-singularity-sra-tools-2.11.0--pl5321ha49a11a_3.img
[ -e "$SRA_IMG" ] && echo "sra-tools image OK" || echo "MISSING sra-tools image"
apptainer exec "$SRA_IMG" prefetch --version 2>&1 | head -2
apptainer exec "$SRA_IMG" fasterq-dump --version 2>&1 | head -2
echo "=== check we can write to /ictstr01/scratch/users/$USER ==="
mkdir -p /ictstr01/scratch/users/$USER/oralprob_test && rmdir /ictstr01/scratch/users/$USER/oralprob_test && echo "ictstr01 scratch writable" || echo "ictstr01 scratch NOT writable"
echo "=== done ==="
