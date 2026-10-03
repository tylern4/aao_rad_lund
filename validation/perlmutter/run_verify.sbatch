#!/bin/bash
# Verification arm of the Perlmutter scan: runs after the Fortran, JAX-CPU
# and JAX-GPU batches (submit_all.sh sets the afterok dependency) and writes
# the statistical comparison under $SCRATCH/aao_rad_scan/verify/.

# Account m3792's default path lands in the `debug` QoS, whose wall limit is 30
# minutes -- too short for a 192-run arm.  `shared` allows 2 days on one node,
# which is what these jobs need; AAO_QOS overrides it.
#SBATCH --qos=shared
#SBATCH --job-name=aao_verify
#SBATCH --account=m3792
#SBATCH --output=/global/u1/t/tylern/aao_%j.out
#SBATCH --constraint=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=2-00:00:00

set -euo pipefail
# Locate the checkout: $0 is not the script's own path when the job arrives
# through the nersc API, so take the first candidate that really is this repo.
REPO=""
for cand in "${AAO_REPO:-}" \
            "$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)/../.." \
            "${SLURM_SUBMIT_DIR:-}" \
            "$HOME/aao_rad_lund"; do
  if [ -n "$cand" ] && [ -f "$cand/pyproject.toml" ] && [ -d "$cand/parms" ]; then
    REPO=$(cd "$cand" && pwd)
    break
  fi
done
if [ -z "$REPO" ]; then
  echo "cannot find the aao_rad_lund checkout -- set AAO_REPO" >&2
  exit 1
fi
HERE="$REPO/validation/perlmutter"
ROOT=${AAO_SCAN_ROOT:-${SCRATCH:?SCRATCH is not set}/aao_rad_scan}
PY=${AAO_PYTHON:-${AAO_VENV:-$SCRATCH/aao_venv}/bin/python}

"$PY" "$HERE/verify_statistics.py" --root "$ROOT"
