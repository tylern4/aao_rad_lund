#!/bin/bash
# JAX-on-CPU arm of the Perlmutter scan: every manifest job from one process,
# with XLA's CPU thread pool spanning the whole node.  Runs are independent
# 20k-event batches, so the unit of parallelism is the run; one process with
# all 128 threads keeps one copy of the tables and one compilation cache hot
# instead of 128 of each.
#
# Needs the venv built by submit_all.sh (AAO_VENV, default $SCRATCH/aao_venv).
#SBATCH --job-name=aao_py_cpu
#SBATCH --account=m3792
#SBATCH --output=/global/u1/t/tylern/aao_%j.out
#SBATCH --constraint=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=12:00:00

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

if [ ! -x "$PY" ]; then
  echo "venv python not found at $PY -- run submit_all.sh first" >&2
  exit 1
fi

# A PM-CPU node is 128 cores / 256 threads, so the mask has one entry per hardware
# thread.  XLA sizes its pool from that mask, and 256 threads on 128 cores
# oversubscribes; AAO_CPU_THREADS caps it (default: one thread per core, i.e.
# half the mask).  calibrate.sbatch measures whether threads help at all here.
AFFINITY=$(nproc)
THREADS=${AAO_CPU_THREADS:-$(( AFFINITY / 2 ))}
echo "affinity $AFFINITY threads, using $THREADS"
export OMP_NUM_THREADS="$THREADS"
export JAX_NUM_THREADS="$THREADS"
unset XLA_FLAGS

"$PY" "$HERE/make_grid.py" --root "$ROOT"

# --cpus-per-job defaults to 0, which leaves the process the whole allocation so
# XLA sizes its pool from the node's affinity mask.  AAO_CPU_JOBS exists to try a
# processes x threads split without editing this file.
"$PY" "$HERE/run_batch.py" --code py_cpu --root "$ROOT" --repo "$REPO" \
  --python "$PY" --jobs "${AAO_CPU_JOBS:-1}"
