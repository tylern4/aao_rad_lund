#!/bin/bash
# Fortran arm of the Perlmutter scan: build (if needed), generate the grid,
# then run every manifest job on one PM-CPU node, one run per core.
#
# Submit through submit_all.sh (it sets --account and the log paths), or on
# its own with:
#   sbatch --account=m3792 -o $SCRATCH/aao_rad_scan/logs/fortran_%j.out \
#          run_fortran.sbatch

# Account m3792's default path lands in the `debug` QoS, whose wall limit is 30
# minutes -- too short for a 192-run arm.  `shared` allows 2 days on one node,
# which is what these jobs need; AAO_QOS overrides it.
#SBATCH --qos=shared
#SBATCH --job-name=aao_fortran
#SBATCH --account=m3792
#SBATCH --output=/global/u1/t/tylern/aao_%j.out
#SBATCH --constraint=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=2-00:00:00

set -euo pipefail
# Locate the checkout.  $0 is not the script's own path when the job arrives
# through the nersc API -- that runs the script with a working directory of $HOME
# -- so deriving REPO from it silently yields a directory with no src/ in it and
# every run then fails with ModuleNotFoundError.  Take the first candidate that
# actually looks like this repository instead.
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
