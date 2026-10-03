#!/bin/bash
# JAX-on-CPU arm of the Perlmutter scan: every manifest job from one process,
# with XLA's CPU thread pool spanning the whole node.  Runs are independent
# 20k-event batches, so the unit of parallelism is the run; one process with
# all 128 threads keeps one copy of the tables and one compilation cache hot
# instead of 128 of each.
#
# Needs the venv built by submit_all.sh (AAO_VENV, default $SCRATCH/aao_venv).
#SBATCH --job-name=aao_py_cpu
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

if [ ! -x "$PY" ]; then
  echo "venv python not found at $PY -- run submit_all.sh first" >&2
  exit 1
fi

"$PY" "$HERE/make_grid.py" --root "$ROOT"

# --cpus-per-job defaults to 0, which leaves the single process the whole
# allocation so XLA sizes its pool from the node's affinity mask.
# AAO_CPU_JOBS exists to try a processes x threads split without editing this
# file.
"$PY" "$HERE/run_batch.py" --code py_cpu --root "$ROOT" --repo "$REPO" \
  --python "$PY" --jobs "${AAO_CPU_JOBS:-1}"
