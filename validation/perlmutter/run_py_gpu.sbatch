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

# Account m3792's default path lands in the `debug` QoS, whose wall limit is 30
# minutes -- too short for a 192-run arm.  `shared` allows 2 days on one node,
# which is what these jobs need; AAO_QOS overrides it.
#SBATCH --qos=shared
#SBATCH --job-name=aao_py_gpu
#SBATCH --account=m3792
#SBATCH --output=/pscratch/sd/t/tylern/aao_rad_scan/logs/aao_py_gpu_%j.out
# Slurm does not expand $SCRATCH here -- --output is read at submit time, before any
# of this script's own setup -- so the path is spelled out.  It cannot be $HOME or
# /tmp either: $HOME is 74% full, and /tmp is not shared with the compute nodes.
#SBATCH --constraint=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --gpus=4
#SBATCH --time=2-00:00:00

set -euo pipefail
# Locate the checkout.  Neither $0 nor $HOME works when the job arrives through
# the nersc API: Slurm is handed /var/spool/slurmd/job<id>/slurm_script as the
# script, and the job runs with $HOME=/global/homes/t/tylern but a working
# directory of /global/u1/t/tylern.  Deriving REPO from either gives a directory
# with no src/ in it, and every run then dies with ModuleNotFoundError.  Try the
# plausible anchors and accept the first that really is this repository.  The
# test is for src/aao_rad, not parms/: parms/ is gitignored and only exists after
# submit_all.sh unpacks parms.tar.gz, so requiring it would reject a perfectly
# good checkout on a fresh clone.
REPO=""
for cand in "${AAO_REPO:-}" \
            "$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)/../.." \
            "${SLURM_SUBMIT_DIR:-}" \
            "${SCRATCH:-}/aao_rad_lund" \
            "$PWD/aao_rad_lund"; do
  if [ -n "$cand" ] && [ -f "$cand/pyproject.toml" ] && [ -d "$cand/src/aao_rad" ]; then
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
