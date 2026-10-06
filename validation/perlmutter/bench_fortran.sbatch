#!/bin/bash
# Speedup benchmark, Fortran arm: four configurations on one PM-CPU node at a
# per-configuration event quota sized to burn ~600 s of serial work each.
#
# Why a reduced quota.  A debug job runs 30 minutes, no longer
# (QOSMaxWallDurationPerJobLimit), and 20,000 events of cfg_034 is 1.6e8 trials.
# Measured in this scan's own Fortran arm that is 23,792 s -- seven hours -- so
# a production-length reference cannot be timed in a test at all.  The quota is
# therefore set per configuration rather than uniformly, because the four span
# a factor of 16 in trials per second: cfg_034 manages 6,599 trials/s where
# cfg_082 manages 108,907, and one shared budget would either finish cfg_082
# before the clock has resolution or strand cfg_034 at the wall limit.
#
# Each quota is (600 s) x (that configuration's measured production rate) /
# (its trials per event), in bench_events.csv, which puts every run at ~600 s
# *under production contention* -- the 192-process, 128-core layout the real arm
# uses.  Four runs on 128 cores have no contention at all, so the expectation is
# that they finish sooner; that direction is the safe one, since a benchmark
# that ran longer than the budget would be killed by the QoS rather than
# reported.
#
# What it buys.  bench_report.py prints each run's trials per second, which is
# what the speedup is a ratio of.  The port arm (bench_gpu.sbatch) cannot run
# the same quota and time it: it is roughly two orders of magnitude faster, so
# at 600 s of Fortran work its loop would be six seconds -- below the resolution
# worth reporting.  Both arms are therefore measured at a workload each can
# finish comfortably, and the comparison is rate against rate, each arm's own
# events per trial included.

#SBATCH --qos=debug
#SBATCH --job-name=aao_bench_fort
#SBATCH --account=m3792
#SBATCH --output=/pscratch/sd/t/tylern/aao_rad_scan/logs/aao_bench_fort_%j.out
# Slurm does not expand $SCRATCH in --output, which is read at submit time, so
# the path is spelled out.  Not $HOME (74% full) and not /tmp (not shared with
# the compute nodes).
#SBATCH --constraint=cpu
# A whole node, so the four runs land on four separate cores and the rate they
# report is one uncontended core's -- which is what the node-hours arithmetic
# multiplies out.  Under `shared` --exclusive is refused outright.
#SBATCH --exclusive
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
BENCH="$ROOT/reduced"
PY=${AAO_PYTHON:-${AAO_VENV:-$SCRATCH/aao_venv}/bin/python}
# Python block-buffers stdout when it is not a terminal, and every one of these
# jobs can be stopped at its wall limit -- unbuffer so the log keeps them.
export PYTHONUNBUFFERED=1
mkdir -p "$ROOT/logs"

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

# The same four configurations the GPU arm runs: beam energies 2, 6, 10 and
# 12 GeV, both values of delta -- the axis that sets how expensive a
# configuration is (median trials per event 58,566 at delta 0.005 against
# 9,118 at delta 0.05).
CFG=${AAO_BENCH_CFGS:-cfg_000,cfg_034,cfg_064,cfg_082}

# Written into <root>/grid/ before make_grid runs so it is picked up as the
# default override (make_grid.py:252-256) and survives any later regeneration
# of this root.  Passing it only as --n-events-override would set the quota for
# this invocation and be lost the next time anything ran make_grid --root here,
# which is precisely the failure that flag's help text warns about.
mkdir -p "$BENCH/grid"
cp "$HERE/bench_events.csv" "$BENCH/grid/n_events_override.csv"
"$PY" "$HERE/make_grid.py" --root "$BENCH"

# One worker per hardware thread, as production is launched; only four of them
# are in flight here, so each has a core to itself.  taskset_prefix does the
# pinning (run_batch.py:194).
echo "fortran: $(nproc) workers available, four runs to do"
"$PY" "$HERE/run_batch.py" --code fortran --root "$BENCH" --repo "$REPO" \
  --python "$PY" --jobs "$(nproc)" --cpus-per-job 1 \
  --only-cfg "$CFG" --warmup-events 0

# 128 physical cores per PM-CPU node, so one run per core is the packing the
# node-hours are divided by; production in fact launched $(nproc) processes,
# two per core, which does not double their throughput.
"$PY" "$HERE/bench_report.py" --label "fortran, reduced quota" \
  --concurrency 128 "$BENCH/fortran"

echo
echo "== elapsed for the job: $SECONDS s"
