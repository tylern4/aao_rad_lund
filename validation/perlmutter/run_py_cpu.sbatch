#!/bin/bash
# JAX-on-CPU arm of the Perlmutter scan: every manifest job on one PM-CPU
# node, one single-threaded worker per core (taskset keeps each jax process
# from spawning a full-node thread pool).
#
# Needs the venv built by submit_all.sh (AAO_VENV, default $SCRATCH/aao_venv).
#SBATCH --job-name=aao_py_cpu
#SBATCH --constraint=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=04:00:00

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

"$PY" "$HERE/run_batch.py" --code py_cpu --root "$ROOT" --repo "$REPO" \
  --python "$PY" --jobs "$(nproc)" --pin-cpu
