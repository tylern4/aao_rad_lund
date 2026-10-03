#!/bin/bash
# Fortran arm of the Perlmutter scan: build (if needed), generate the grid,
# then run every manifest job on one PM-CPU node, one run per core.
#
# Submit through submit_all.sh (it sets --account and the log paths), or on
# its own with:
#   sbatch --account=m3792 -o $SCRATCH/aao_rad_scan/logs/fortran_%j.out \
#          run_fortran.sbatch
#SBATCH --job-name=aao_fortran
#SBATCH --constraint=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=12:00:00

set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$HERE/../.." && pwd)
ROOT=${AAO_SCAN_ROOT:-${SCRATCH:?SCRATCH is not set}/aao_rad_scan}
PY=${AAO_PYTHON:-${AAO_VENV:-$SCRATCH/aao_venv}/bin/python}
mkdir -p "$ROOT/logs"

# The system python3 here is 3.6, which cannot even parse these scripts, so the
# venv python submit_all.sh built is used for the driver and the grid.
if [ ! -x "$PY" ]; then
  echo "venv python not found at $PY -- run submit_all.sh first" >&2
  exit 1
fi

# gcc-native is already on PATH on the compute nodes; `module load gcc` fails
# there, so it is only a fallback.
command -v gfortran >/dev/null || module load gcc || true
if [ ! -x "$REPO/bin/aao_rad_lund" ]; then
  echo "== building aao_rad (gfortran)"
  command -v gfortran >/dev/null || {
    echo "no gfortran on PATH -- load a gcc/craype compiler module" >&2
    exit 1
  }
  make -C "$REPO"
fi

"$PY" "$HERE/make_grid.py" --root "$ROOT"

# Each run gets its own CWD with a symlinked spp_tbl (the MAID tables are opened
# relative to CWD) and AAO_TRIALS=0, which keeps the 0.5 GB trial dump off while
# the event n-tuple is still written -- both handled by run_batch.py.
#
# The original is serial Fortran (no OpenMP), so one process per core is the
# only way to fill the node: 128 single-core workers.
"$PY" "$HERE/run_batch.py" --code fortran --root "$ROOT" --repo "$REPO" \
  --jobs "$(nproc)" --cpus-per-job 1

# The Fortran seeds myran from unixtime (1 s resolution): guard against two
# runs of one configuration starting in the same second.
"$PY" "$HERE/run_batch.py" --code fortran --root "$ROOT" --repo "$REPO" \
  --check-collisions
