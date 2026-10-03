#!/bin/bash
# How does a single run scale with threads on this node, and which node layout
# finishes the 192-run grid first?
#
# Every run is one fixed 65,536-trial vectorised pass (generate.stream sizes its
# chunk as max(min(chunk_events, n_events), batch_size), so a run's cost is set
# by how many trials it takes to reach n_events accepted ones, and is the same
# whatever it keeps).  That makes this a clean measurement: time real runs at
# several CPU counts, then project the wall time of the candidate processes x
# threads layouts.  Thread counts come from taskset, because XLA sizes its pool
# from sched_getaffinity.
#
# Runs are SHORT here -- AAO_CALIB_EVENTS, 300 by default, not the grid's 20,000 --
# and the projection scales by the ratio.  This is not a rounding of the
# measurement: the 1-thread point is the whole reason the script exists, and at
# the grid's 20,000 events a single-threaded run takes hours.  A 300-event run at
# 1 thread takes tens of seconds, so all eight thread counts fit in the debug
# queue's 30 minutes; the first attempt used 20,000 and burned its whole limit on
# the warmup plus one label.  Compilation is a fixed few seconds against tens of
# seconds of work, so dropping the first repeat still leaves a clean steady-state
# number at every point.
#
# Each thread count is measured with REPEATS sequential runs *inside one
# process*, and only the repeats after the first are used.  The first pays for
# jax's compilation; a process that compiles once and then runs 192 batches pays
# it once, so timing fresh processes would charge every layout for a cost only
# the first process really pays -- and would penalise exactly the many-thread
# layouts whose runs are short.
#
# Results are appended to $ROOT/calib_results.txt as each label finishes, and the
# projection reads that file rather than in-memory state, so a job killed by its
# wall clock still leaves usable numbers instead of nothing at all.
#
# Run it before the scan (it takes a few minutes) and set run_py_cpu.sbatch's
# --jobs/--cpus-per-job to whatever wins here:
#
#   nersc submit $SCRATCH/aao_rad_lund/validation/perlmutter/calibrate.sbatch
#
# For short experiments this is `debug` rather than `shared`, deliberately
# inverted from the scan arms.  The shared CPU pool is often saturated enough
# that a 128-core request waits hours, while debug is a 30-minute queue that
# starts in seconds and still gets the whole node -- the first calibration run
# this way took 4s to start and got 128 cores.  The measurement has to fit in
# those 30 minutes, which the short runs above are sized to do.
#SBATCH --qos=debug
#SBATCH --job-name=aao_calib
#SBATCH --account=m3792
# All real work and all logs live under $SCRATCH.  Slurm does not expand $SCRATCH
# in --output (it is read at submit time, before any of this script's setup), so
# the path is spelled out; $HOME is no alternative since it is 74% full, and /tmp
# is not shared with the compute nodes.  The API's own default of
# /global/u1/<user> is a symlink into $HOME and lands there for the same reason.
#
# Do not move or relocate this file while the job runs.  $SCRATCH and $HOME are
# different filesystems, so "moving" a log is a copy plus an unlink: Slurm holds
# the original inode open and keeps writing to it, which is now unlinked, and the
# copy stops growing.  The per-label logs under $ROOT are written by the children
# themselves and are unaffected -- read those if this one goes quiet.
#SBATCH --output=/pscratch/sd/t/tylern/aao_rad_scan/logs/aao_calib_%j.out
#SBATCH --constraint=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=00:28:00

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
export PYTHONPATH="$REPO/src${PYTHONPATH:+:$PYTHONPATH}"

# jax is installed from jax[cuda12], so on a CPU node its CUDA plugin raises
# CUDA_ERROR_NO_DEVICE while initialising.  That is noise -- it falls back to CPU
# -- but it costs a few seconds per process and fills the per-label logs with a
# traceback, so pin the platform and keep the measurement clean.
export JAX_PLATFORMS=cpu

if [ ! -x "$PY" ]; then
  echo "venv python not found at $PY -- run submit_all.sh first" >&2
  exit 1
fi
mkdir -p "$ROOT" "$ROOT/logs"

# The production event count comes from the grid itself, so the projection stays
# honest if make_grid.py's N_EVENTS is ever changed.
PROD_EVENTS=$("$PY" -c 'import sys; sys.path.insert(0, sys.argv[1]); import make_grid; print(make_grid.N_EVENTS)' "$HERE")
CALIB_EVENTS=${AAO_CALIB_EVENTS:-300}
REPEATS=${AAO_CALIB_REPEATS:-2}
THREAD_COUNTS=${AAO_CALIB_THREADS:-"1 2 4 8 16 32 64 128"}
LABEL_TIMEOUT=${AAO_CALIB_LABEL_TIMEOUT:-600}
RESULTS="$ROOT/calib_results.txt"
# Start fresh: a stale file from an earlier attempt with different settings would
# be read back as if it belonged to this one.
rm -f "$RESULTS"

# A real configuration from the grid, so the measurement includes the physics.
CARD=$("$PY" - <<'EOF'
import csv, os
root = os.environ.get("AAO_SCAN_ROOT") or os.path.join(os.environ["SCRATCH"], "aao_rad_scan")
rows = list(csv.DictReader(open(os.path.join(root, "grid", "manifest.csv"))))
print(os.path.join(root, "grid", rows[0]["cfg_id"] + ".txt"))
EOF
)
NCPU=$("$PY" -c 'import os; print(len(os.sched_getaffinity(0)))')
FIRST=$("$PY" -c 'import os; print(min(os.sched_getaffinity(0)))')
echo "card   $CARD"
echo "node   $NCPU cpus in the affinity mask, lowest is $FIRST"
echo "timing $CALIB_EVENTS events x $REPEATS per thread count, grid uses $PROD_EVENTS"
echo "scale  ${THREAD_COUNTS}"
echo

# $1 = cpu count, $2 = label.  Prints the per-run seconds, space separated, or
# the word TIMEOUT if the label exceeded LABEL_TIMEOUT.
one_process() {
  local cpus
  # Take a slice of the real affinity mask, striding rather than taking the
  # lowest n: a PM-CPU node's 128 cores are 256 threads, and the lowest 128 of
  # those are one hardware thread per core with the sibling threads left idle.
  cpus=$("$PY" - "$1" <<'EOF'
import os, sys
count = int(sys.argv[1])
mask = sorted(os.sched_getaffinity(0))
if count >= len(mask):
    picked = mask
else:
    picked = mask[:: max(1, len(mask) // count)][:count]
print(",".join(str(c) for c in picked))
EOF
)
  if runs=$(timeout "$LABEL_TIMEOUT" taskset -c "$cpus" "$PY" - \
      "$CARD" "$REPEATS" "$ROOT" "$2" "$REPO" "$CALIB_EVENTS" \
      2> "$ROOT/calib_$2.log" <<'EOF'
import sys, time
from pathlib import Path

from aao_rad.cli import main

card, repeats, root, label, repo, events = sys.argv[1:7]
times = []
for i in range(int(repeats)):
    out = str(Path(root) / f"calib_{label}.npz")
    t0 = time.perf_counter()
    main([
        "--legacy-input", card,
        "--parms", str(Path(repo) / "parms"),
        "--n-events", events,
        "--seed", str(4242 + i),
        "--ek-sampling", "fortran",
        "--format", "npz",
        "-o", out,
    ])
    times.append(time.perf_counter() - t0)
print(" ".join(f"{t:.2f}" for t in times))
EOF
  ); then
    printf '%s' "$runs"
  else
    printf 'TIMEOUT'
  fi
}

# Warm the page cache on the MAID tables and pay the compilation once: every
# process reads them, so timing a cold first run would measure the filesystem
# rather than the layout.
echo "warmup: $(one_process "$NCPU" warmup)"
echo

for k in $THREAD_COUNTS; do
  [ "$k" -le "$NCPU" ] || { echo "threads $k: skipped, node has $NCPU"; continue; }
  runs=$(one_process "$k" "t$k")
  echo "threads $k: $runs"
  # Record as we go.  A label that timed out is written out as such and dropped
  # by the projection rather than guessed at: its true cost is unknown, and only
  # knowing "at least this slow" is not enough to interpolate from.
  printf 't%s\t%s\n' "$k" "${runs// /,}" >> "$RESULTS"
done

echo
"$PY" - "$NCPU" "$CALIB_EVENTS" "$PROD_EVENTS" "$RESULTS" <<'EOF'
import math, sys

ncpu = int(sys.argv[1])
calib_events, prod_events = int(sys.argv[2]), int(sys.argv[3])
results_path = sys.argv[4]

measured = {}
for line in open(results_path):
    line = line.strip()
    if not line:
        continue
    label, runs = line.split("\t")
    if runs == "TIMEOUT":
        print(f"  (dropped {label}: exceeded the per-label timeout)")
        continue
    measured[int(label[1:])] = [float(x) for x in runs.split(",")]

if not measured:
    sys.exit("no usable measurements")

# steady state: the repeats after the compile.  min, because a run interrupted
# by another tenant should not be mistaken for the layout's cost.
steady = {k: min(v[1:]) if len(v) > 1 else v[0] for k, v in measured.items()}
scale = prod_events / calib_events


def interp(k):
    """Steady-state seconds for k threads, log-interpolated between points."""
    ks = sorted(steady)
    if k <= ks[0]:
        return steady[ks[0]] * ks[0] / max(k, 1)
    if k >= ks[-1]:
        return steady[ks[-1]]
    lo = max(x for x in ks if x <= k)
    hi = min(x for x in ks if x >= k)
    if lo == hi:
        return steady[lo]
    f = (math.log(k) - math.log(lo)) / (math.log(hi) - math.log(lo))
    return steady[lo] * (steady[hi] / steady[lo]) ** f


runs = 192
one = steady.get(1)
if one is None:
    print(f"no 1-thread measurement: no single-thread baseline to compare "
          f"against, and no extrapolated cost below {min(steady)} threads")
    print(f"scaling {calib_events}-event runs to the grid's {prod_events}")
else:
    print(f"steady-state run at {calib_events} events, 1 thread: {one:.1f}s")
    print(f"scaled to the grid's {prod_events} events: {one * scale:.1f}s")
print()
print(f"{'layout':>14} {'threads':>8} {'waves':>6} {'projected':>10} {'vs 1x1':>8}")
best = None
for procs in (1, 2, 4, 8, 16, 32, 64, 128):
    if procs > ncpu:
        break
    threads = max(1, ncpu // procs)
    total = math.ceil(runs / procs) * interp(threads) * scale
    gain = f"{one * runs * scale / total:>7.1f}x" if one else f"{'-':>8}"
    print(f"{procs:>6} x {threads:<4}   {threads:>8} {math.ceil(runs / procs):>6} "
          f"{total / 60:>9.1f}m {gain}")
    if best is None or total < best[0]:
        best = (total, procs, threads)
print()
print(f"best layout: {best[1]} processes x {best[2]} threads "
      f"({best[0] / 60:.1f} min for {runs} runs)")
print(f"set run_py_cpu.sbatch: --jobs {best[1]} --cpus-per-job {best[2]}")
EOF
