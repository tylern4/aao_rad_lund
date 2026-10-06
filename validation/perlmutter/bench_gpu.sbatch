#!/bin/bash
# Speedup benchmark, JAX-on-GPU arm: the same four configurations run twice on
# one PM-GPU node, once with an empty XLA cache and once with the cache the
# first pass filled.
#
# Why two passes.  The port's wall clock splits into a loop (a few tens of
# seconds) and everything before it, and the question is how much of that
# preamble is compilation rather than fixed interpreter cost.  The split is
# real: generate.py starts its own clock at t_start (generate.py:863) and
# reports the two apart as `startup` and `wall time` in bench_report.py.
# run_batch.py's warm-up runs one throwaway process before the fan-out
# (run_batch.py:446), but _make_sampler and _make_stepper bake the kinematics
# and the batch size into the compiled program (generate.py:742-755), so a
# warm-up of runs[0] warms exactly one of the 96 configurations and the other
# 95 compile themselves inside the fan-out.  Pass 1 starts from an empty XLA
# cache and pass 2 from the one pass 1 filled: startup falling from pass 1 to
# pass 2 is the compile, startup not moving is fixed interpreter/import/table
# cost that no cache will help.
#
# The ~457 s of startup this script was originally written to explain was not
# a measurement of either.  It was birth-to-mtime, and on Lustre `stat -c %W`
# reports the *server's* clock while mtime is stamped by the node -- a ~490 s
# skew charged to every run, which is why the first submission read 490.0 s
# for eight runs that had each taken about twenty seconds.  bench_report.py
# now takes the wall clock from the driver's own timing.json
# (run_batch.execute), so this is the first submission whose two passes
# compare a real compile against a real compile.
#
# Why the full 20,000-event quota.  This arm is the one whose loop time is
# compared against the Fortran's, and at production length the loop is 20-70 s
# -- long enough that a 0.1 s print resolution is a fraction of a percent.  The
# Fortran arm (bench_fortran.sbatch) cannot use the same quota: 20,000 events of
# cfg_064 is 5e8 trials and two hours of serial work, so it runs a reduced
# quota and the speedup is reported as a ratio of *rates*, each arm measured at
# a workload it can finish inside the 30 minutes the debug QoS allows.  No
# quota exists at which both arms run long enough to time directly -- that is
# the speedup.

#SBATCH --qos=debug
#SBATCH --job-name=aao_bench_gpu
#SBATCH --account=m3792
#SBATCH --output=/pscratch/sd/t/tylern/aao_rad_scan/logs/aao_bench_gpu_%j.out
# Slurm does not expand $SCRATCH in --output, which is read at submit time, so
# the path is spelled out.  Not $HOME (74% full) and not /tmp (not shared with
# the compute nodes).
#
# `debug`, not `gpu_debug`.  A `gpu_debug` QoS exists and sacctmgr reports its
# MaxWall as 00:30:00, but every request against it is refused outright --
# `sbatch --test-only` comes back with "Job request does not match any supported
# policy" for all of gres=gpu:1 and gres=gpu:4, with and without
# --constraint=gpu, with and without --exclusive, and for the bare QoS on its
# own.  `debug` accepts the same request unchanged, so the GPU arm runs there.
# `gpu_interactive` and `gpu_regular` refuse it too; `regular` accepts it but
# its predicted start was six days out.
#SBATCH --constraint=gpu
#SBATCH --gpus=4
#SBATCH --exclusive
# The debug QoS refuses anything longer (QOSMaxWallDurationPerJobLimit).  The
# first submission's two passes took 145 s of job time in all; the rest is
# queueing slack.
#SBATCH --time=00:30:00

set -euo pipefail

# Locate the checkout.  Neither $0 nor $HOME works when the job arrives through
# the nersc API: Slurm is handed /var/spool/slurmd/job<id>/slurm_script as the
# script, and the job runs with $HOME=/global/homes/t/tylern but a working
# directory of /global/u1/t/tylern.  Deriving REPO from either gives a directory
# with no src/ in it, and every run then dies with ModuleNotFoundError.
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
ROOT=${AAO_BENCH_ROOT:-${SCRATCH:?SCRATCH is not set}/aao_bench}
BENCH="$ROOT/full"
PY=${AAO_PYTHON:-${AAO_VENV:-$SCRATCH/aao_venv}/bin/python}
# Python block-buffers stdout when it is not a terminal, and every one of these
# jobs can be stopped at its wall limit -- unbuffer so the log keeps them.
export PYTHONUNBUFFERED=1
mkdir -p "$ROOT/logs"

if [ ! -x "$PY" ]; then
  echo "venv python not found at $PY -- run submit_all.sh first" >&2
  exit 1
fi

# Four configurations, one wave across the node's four A100s: beam energies 2,
# 6, 10 and 12 GeV, and both values of delta, which is the axis that decides how
# expensive a configuration is (median trials per event 58,566 at delta 0.005
# against 9,118 at delta 0.05 -- the tighter explicit-photon cut costs ~6x).
CFG=${AAO_BENCH_CFGS:-cfg_000,cfg_034,cfg_064,cfg_082}
GPUS=${AAO_GPUS:-4}
CACHE=$ROOT/xla_cache

# Empty on purpose.  A cache left by an earlier attempt would turn the cold
# measurement into a warm one and the whole job would report a startup that
# does not exist anywhere in production.
rm -rf "$CACHE"
mkdir -p "$CACHE"
export AAO_JAX_CACHE_DIR="$CACHE"
export JAX_COMPILATION_CACHE_DIR="$CACHE"
# One worker per GPU, each handed a single device by run_batch.py; without this
# every process would try to reserve most of every visible device.
export XLA_PYTHON_CLIENT_PREALLOCATE=false

"$PY" "$HERE/make_grid.py" --root "$BENCH"

archive() {
  local tag=$1 dest="$BENCH/timing/$1"
  rm -rf "$dest"
  for d in "$BENCH"/py_gpu/*/; do
    [ -f "$d/out.txt" ] || continue
    mkdir -p "$dest/$(basename "$d")"
    cp -f "$d/out.txt" "$dest/$(basename "$d")/"
    # timing.json has to travel with out.txt.  It holds the driver's absolute
    # start and end epochs for *this* pass, and pass 2 overwrites both files in
    # the live directory -- so a report run against the archive without it would
    # fall back to the filesystem clock, which on Lustre is the server's clock
    # and ~490 s behind the node's.  cp is safe with it because the epochs inside
    # are absolute; nothing here depends on the copy's own timestamps.
    if [ -f "$d/timing.json" ]; then
      cp -f "$d/timing.json" "$dest/$(basename "$d")/"
    fi
  done
}

run_pass() {
  # run_batch.py's --only-cfg implies a re-run (select_pending, run_batch.py:481),
  # and clear_stale_outputs() then removes the previous out.txt (run_batch.py:543)
  # -- so pass 2 re-runs rather than reporting "0 to run" and exiting 0 having
  # measured nothing, and pass 1's numbers have to be archived before it.
  "$PY" "$HERE/run_batch.py" --code py_gpu --root "$BENCH" --repo "$REPO" \
    --python "$PY" --gpus "$GPUS" --jobs "$GPUS" \
    --cpus-per-job "$(( $(nproc) / GPUS ))" \
    --only-cfg "$CFG" --warmup-events 0
}

echo "nproc: $(nproc)   gpus: $GPUS   cache: $CACHE"

echo
echo "== pass 1 -- cold cache (started $(date -u +%H:%M:%S))"
t0=$SECONDS
run_pass
archive pass_cold
echo "pass 1 took $(( SECONDS - t0 ))s; cache entries: $(ls "$CACHE" | wc -l)"
"$PY" "$HERE/bench_report.py" --label "pass 1: cold XLA cache" \
  --concurrency "$GPUS" "$BENCH/timing/pass_cold"

echo
echo "== pass 2 -- same cache, same configurations (started $(date -u +%H:%M:%S))"
t0=$SECONDS
run_pass
echo "pass 2 took $(( SECONDS - t0 ))s; cache entries: $(ls "$CACHE" | wc -l)"
"$PY" "$HERE/bench_report.py" --label "pass 2: warm XLA cache" \
  --concurrency "$GPUS" "$BENCH/py_gpu"
archive pass_warm

echo
echo "== elapsed for the job: $SECONDS s"
