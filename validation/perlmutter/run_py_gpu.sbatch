#!/bin/bash
# JAX-on-GPU arm of the Perlmutter scan: one process per A100 (4 on a PM-GPU
# node), each pinned to its device and given a quarter of the host cores for the
# host side of the work.
#
# Needs the venv built by submit_all.sh (jax[cuda12]).
#SBATCH --job-name=aao_py_gpu
#SBATCH --constraint=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --gpus=4
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

# JAX recompiles per configuration (the beam energy enters as a constant);
# a persistent cache lets the CPU and GPU arms share compile time.
export JAX_COMPILATION_CACHE_DIR="$ROOT/jax_cache"

# One worker per GPU, each given a single device by run_batch.py; without this
# every process would try to reserve most of every visible device.
export XLA_PYTHON_CLIENT_PREALLOCATE=false

GPUS=4
"$PY" "$HERE/make_grid.py" --root "$ROOT"

# One process per GPU, each restricted to a quarter of the host cores so the
# four processes do not oversubscribe the node's threads.
"$PY" "$HERE/run_batch.py" --code py_gpu --root "$ROOT" --repo "$REPO" \
  --python "$PY" --gpus "$GPUS" --jobs "$GPUS" --cpus-per-job "$(( $(nproc) / GPUS ))"
