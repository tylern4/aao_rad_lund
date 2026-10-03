#!/bin/bash
# Fortran arm of the Perlmutter scan: build (if needed), generate the grid,
# then run every manifest job on one PM-CPU node, one run per core.
#
# Submit through submit_all.sh (it sets --account and the log paths), or on
# its own with:
#   sbatch --account=<project> -o $SCRATCH/aao_rad_scan/logs/fortran_%j.out \
#          run_fortran.sbatch
#SBATCH --job-name=aao_fortran
#SBATCH --constraint=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=08:00:00

set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$HERE/../.." && pwd)
ROOT=${AAO_SCAN_ROOT:-${SCRATCH:?SCRATCH is not set}/aao_rad_scan}
PY=${AAO_PYTHON:-python3}
mkdir -p "$ROOT/logs"

module load gcc
if [ ! -x "$REPO/bin/aao_rad_lund" ]; then
  echo "== building aao_rad (gfortran)"
  make -C "$REPO"
fi

"$PY" "$HERE/make_grid.py" --root "$ROOT"

# One run per core; each run gets its own CWD (the MAID tables are opened
# relative to CWD) and a symlinked spp_tbl.  AAO_TRIALS=0 keeps the 0.5 GB
# trial dump off while the event n-tuple is still written.
"$PY" "$HERE/run_batch.py" --code fortran --root "$ROOT" --repo "$REPO" \
  --jobs "$(nproc)" --pin-cpu

# The Fortran seeds myran from unixtime (1 s resolution): guard against two
# runs of one configuration starting in the same second.
"$PY" "$HERE/run_batch.py" --code fortran --root "$ROOT" --repo "$REPO" \
  --check-collisions
