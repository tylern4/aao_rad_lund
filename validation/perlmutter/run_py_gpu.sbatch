#!/bin/bash
# JAX-on-GPU arm of the Perlmutter scan: one process per A100, each pinned to its
# device and given an equal share of the node's host threads for the host side of
# the work.
#
# Needs the venv built by submit_all.sh (jax[cuda12]).

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

# A PM-GPU node has 4 A100s and 128 hardware threads, and the gpu_shared_ss11
# QoS charges 32 cores per GPU.  The earlier --cpus-per-task=64 was refused for
# exactly that reason ("requires you to request 32.0 cores per GPU, job
# requested 64 cores for 4 GPUs"), and the whole node is precisely the 128 that
# satisfies it.  So this takes the node whole, all 4 GPUs.

#SBATCH --constraint=gpu
#SBATCH --exclusive
#SBATCH --gpus=4
#SBATCH --time=2-00:00:00
# No --nodes/--ntasks/--cpus-per-task on purpose.  A probe with no directives at
# all (debug 59300338) came back with 128 CPUs, 229 GB and
# CUDA_VISIBLE_DEVICES=0,1,2,3 -- the whole node.  Asking for less is what left
# the previous attempt with 112 GB, which 2 processes exhausted (59269609,
# OUT_OF_MEMORY, 95 of 192 runs killed).  --exclusive keeps a co-tenant off the
# memory as well as off the cores.

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

# Must match --gpus above; override both together with AAO_GPUS if the node shape
# ever changes.
GPUS=${AAO_GPUS:-4}
"$PY" "$HERE/make_grid.py" --root "$ROOT"

# One process per GPU, each restricted to an equal share of the node's affinity
# entries so the processes do not oversubscribe it.  Here, unlike in the CPU
# arms, a run does have use for host threads -- it is feeding a device rather
# than saturating a core -- so the node's threads are split between the GPUs
# rather than given one run each.  nproc is read at run time, so this holds
# whether or not the node is hyperthreaded.
"$PY" "$HERE/run_batch.py" --code py_gpu --root "$ROOT" --repo "$REPO" \
  --python "$PY" --gpus "$GPUS" --jobs "$GPUS" \
  --cpus-per-job "$(( $(nproc) / GPUS ))"
