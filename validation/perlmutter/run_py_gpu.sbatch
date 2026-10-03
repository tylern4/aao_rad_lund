#!/bin/bash
# JAX-on-GPU arm of the Perlmutter scan: one process per A100 (4 on a PM-GPU
# node), each pinned to its device and given a quarter of the host cores for the
# host side of the work.
#
# Needs the venv built by submit_all.sh (jax[cuda12]).
#
# A PM-GPU node has 64 cores / 128 threads and 4 A100s, and Slurm counts cores
# here: asking for all 4 GPUs means the whole node, so 64 cores is the matching
# request.  (A PM-CPU node is 128 cores / 256 threads, which is why the CPU arms
# ask for 128.)
#SBATCH --job-name=aao_py_gpu
#SBATCH --account=m3792
#SBATCH --output=/global/u1/t/tylern/aao_%j.out
#SBATCH --constraint=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --gpus=4
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
