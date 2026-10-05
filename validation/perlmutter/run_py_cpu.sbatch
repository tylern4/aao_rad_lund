#!/bin/bash
# JAX-on-CPU arm of the Perlmutter scan: every manifest job as its own
# process, one hardware thread apiece, all of them in flight.  Runs are
# independent 20k-event batches and a run gains almost nothing from extra
# threads (see the layout note below), so the node is filled with runs rather
# than with threads inside them.
#
# Needs the venv built by submit_all.sh (AAO_VENV, default $SCRATCH/aao_venv).

# `regular`, not `shared`: both allow 2 days, but `shared` caps a job at one
# node's worth of shared access while `regular` is the project's own allocation
# and has no per-job TRES cap.  AAO_QOS overrides it.
#SBATCH --qos=regular
#SBATCH --job-name=aao_py_cpu
#SBATCH --account=m3792
#SBATCH --output=/pscratch/sd/t/tylern/aao_rad_scan/logs/aao_py_cpu_%j.out
# Slurm does not expand $SCRATCH here -- --output is read at submit time, before any
# of this script's own setup -- so the path is spelled out.  It cannot be $HOME or
# /tmp either: $HOME is 74% full, and /tmp is not shared with the compute nodes.

#SBATCH --constraint=cpu
#SBATCH --exclusive
#SBATCH --time=2-00:00:00
# No --nodes/--ntasks/--cpus-per-task on purpose.  Slurm sizes a step from those,
# and --cpus-per-task=128 does *not* mean a whole node: a probe with no
# directives at all (debug 59300339) came back with 256 CPUs and 487 GB, so the
# 128 an explicit request bought was 64 cores' worth.  The old 128 processes x 2
# threads did not fit in that either: it is 256 slots on 64 cores, so every core
# ran two processes and four threads.  Asking for nothing gets the whole node,
# and --exclusive keeps a co-tenant off it.  How many hardware threads that is in
# practice is read from nproc at run time, so nothing here depends on whether the
# node is hyperthreaded.

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
# Every one of these jobs can be killed at its wall limit, and Python block-buffers
# stdout when it is not a terminal.  A killed job then loses its own progress
# reporting -- the '[42/46] cfg_091_s0: ok (102s)' lines that are the only record of
# which runs finished and which were still in flight.  Unbuffer it for every
# invocation below, including the `python -m aao_rad.cli` runs: run_batch.py
# captures each run's output to its own out.txt, but a run killed mid-flight loses
# whatever its own stdout had not yet flushed, which is exactly the progress lines
# that distinguish a slow run from a stalled one.
export PYTHONUNBUFFERED=1
mkdir -p "$ROOT/logs"

if [ ! -x "$PY" ]; then
  echo "venv python not found at $PY -- run submit_all.sh first" >&2
  exit 1
fi

# A PM-CPU node hands out 256 hardware threads, which is 128 cores.  calibrate.sbatch
# measured what a run does with more than one thread: 35.0 s at 1 thread, 28.8 at
# 2, 20.2 at 4, and 8.4 at 64.  The per-run cost in *cores* therefore rises as
# 35.0, 57.6, 80.8, ... core-seconds -- a second thread buys 1.2x of work for 2x the
# cores -- so runs are the unit of parallelism and each gets exactly one core.
#
# One worker per *core*, not per hardware thread.  Asking for $(nproc) put 256
# runs on 128 cores, which is both oversubscription the calibration says is pure
# loss and 384 GB of resident JAX processes on a 503 GB node.  taskset_prefix cuts
# the mask from the real core list, so nproc/2 workers land on 128 distinct cores.
# 192 runs over 128 workers is 1.5 waves either way, so this costs no wall time.
echo "affinity $(($(nproc) / 2)) cores, one core per run, --cpus-per-job 1"
# An inherited XLA_FLAGS could pin the child's intra-op pool to a thread count
# that contradicts the slice it is given below.
unset XLA_FLAGS

"$PY" "$HERE/make_grid.py" --root "$ROOT"

# AAO_CPU_JOBS and AAO_CPU_CPUS override the layout, which is also how a fresh
# calibration is applied without editing this file.
"$PY" "$HERE/run_batch.py" --code py_cpu --root "$ROOT" --repo "$REPO" \
  --python "$PY" --jobs "${AAO_CPU_JOBS:-$(($(nproc) / 2))}" \
  --cpus-per-job "${AAO_CPU_CPUS:-1}"
