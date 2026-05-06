#!/bin/bash
#SBATCH -p cpu_p
#SBATCH -q cpu_normal
#SBATCH -c 1
#SBATCH --mem 1G
#SBATCH -t 00:05:00
#SBATCH -J envcheck3
#SBATCH -o /ictstr01/project_copy/genomics/jru/project/pr0003-OralProbPhage/hpc/envcheck_logs/envcheck3-%j.out
#SBATCH -e /ictstr01/project_copy/genomics/jru/project/pr0003-OralProbPhage/hpc/envcheck_logs/envcheck3-%j.err

set +u
echo "=== node $(hostname) ==="
echo "=== python3 yaml on compute ==="
which python3
python3 --version
python3 -c 'import yaml; print("yaml OK", yaml.__version__)' 2>&1
echo "=== run yaml_to_env directly on compute ==="
python3 ${HOME}/vpipe/bin/yaml_to_env.sh 2>&1 | head -20
echo "=== with PYTHONPATH inspection ==="
python3 -c 'import sys; print("\n".join(sys.path))'
echo "=== source vpipe with full stderr captured ==="
export VPIPEBIN="${HOME}/vpipe/bin"
( source "${VPIPEBIN}/00-config.sh" ) >/tmp/vp_stdout.$$ 2>/tmp/vp_stderr.$$
echo "--- stdout ---"
cat /tmp/vp_stdout.$$
echo "--- stderr ---"
cat /tmp/vp_stderr.$$
rm -f /tmp/vp_stdout.$$ /tmp/vp_stderr.$$
echo "=== now actually re-source and check vars ==="
source "${VPIPEBIN}/00-config.sh"
echo "DBPATH=${DBPATH:-<unset>}"
echo "DB_BIOBAKERY=${DB_BIOBAKERY:-<unset>}"
echo "COMMAND_PREFIX_FLYE=${COMMAND_PREFIX_FLYE:-<unset>}"
echo "=== _VPIPE_CONF_DIR ==="
echo "_VPIPE_CONF_DIR=${_VPIPE_CONF_DIR:-<unset>}"
ls "${_VPIPE_CONF_DIR}/database.yml" 2>&1
echo "=== done ==="
