#!/bin/bash
# Verification arm of the Perlmutter scan: runs after the Fortran, JAX-CPU
# and JAX-GPU batches (submit_all.sh sets the afterok dependency) and writes
# the statistical comparison under $SCRATCH/aao_rad_scan/verify/.
#SBATCH --job-name=aao_verify
#SBATCH --account=m3792
#SBATCH --output=slurm_%j.out
#SBATCH --constraint=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=04:00:00

set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=${AAO_SCAN_ROOT:-${SCRATCH:?SCRATCH is not set}/aao_rad_scan}
PY=${AAO_PYTHON:-${AAO_VENV:-$SCRATCH/aao_venv}/bin/python}

"$PY" "$HERE/verify_statistics.py" --root "$ROOT"
