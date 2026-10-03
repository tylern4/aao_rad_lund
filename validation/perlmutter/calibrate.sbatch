#!/bin/bash
# How does a single run scale with threads on this node, and which node layout
# finishes the 192-run grid first?
#
# Every run is one fixed 65,536-trial vectorised pass (generate.stream sizes its
# chunk as max(min(chunk_events, n_events), batch_size), so a run costs the same
# whatever it keeps), which makes this a clean measurement: time real runs at
# several CPU counts, then project the wall time of the candidate processes x
# threads layouts.  Thread counts come from taskset, because XLA sizes its pool
# from sched_getaffinity.
#
# Each thread count is measured with REPEATS sequential runs *inside one
# process*, and only the repeats after the first are used.  The first pays for
# jax's compilation; a process that compiles once and then runs 192 batches pays
# it once, so timing fresh processes would charge every layout for a cost only
# the first process really pays -- and would penalise exactly the many-thread
# layouts whose runs are short.
#
# Run it before the scan (it takes a few minutes) and set run_py_cpu.sbatch's
# --jobs/--cpus-per-job to whatever wins here:
#
#   nersc submit $SCRATCH/aao_rad_lund/validation/perlmutter/calibrate.sbatch

# `debug` rather than `shared`, deliberately inverted from the scan arms.  The
# shared CPU pool is often saturated enough that a 128-core request waits hours,
# while debug is a 30-minute queue that starts in seconds and still gets the whole
# node -- the first calibration run this way took 4s to start and got 128 cores.
# The measurement only has to fit in 30 minutes, which it does at REPEATS=2.
#SBATCH --qos=debug
#SBATCH --job-name=aao_calib
#SBATCH --account=m3792
#SBATCH --output=/global/u1/t/tylern/aao_%j.out
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
# plausible anchors and accept the first that really is this repository.
REPO=""
for cand in "${AAO_REPO:-}" \
            "$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)/../.." \
            "${SLURM_SUBMIT_DIR:-}" \
            "${SCRATCH:-}/aao_rad_lund" \
            "$PWD/aao_rad_lund"; do
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
export PYTHONPATH="$REPO/src${PYTHONPATH:+:$PYTHONPATH}"

if [ ! -x "$PY" ]; then
  echo "venv python not found at $PY -- run submit_all.sh first" >&2
  exit 1
fi
mkdir -p "$ROOT"

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
REPEATS=${AAO_CALIB_REPEATS:-2}
echo "card   $CARD"
echo "node   $NCPU cpus in the affinity mask, lowest is $FIRST, $REPEATS runs per thread count"
echo

# $1 = cpu count, $2 = label.  Prints the per-run seconds, space separated.
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
  taskset -c "$cpus" "$PY" - "$CARD" "$REPEATS" "$ROOT" "$2" "$REPO" \
    2> "$ROOT/calib_$2.log" <<'EOF'
import sys, time
from pathlib import Path

from aao_rad.cli import main

card, repeats, root, label, repo = sys.argv[1:6]
times = []
for i in range(int(repeats)):
    out = str(Path(root) / f"calib_{label}.npz")
    t0 = time.perf_counter()
    main([
        "--legacy-input", card,
        "--parms", str(Path(repo) / "parms"),
        "--n-events", "20000",
        "--seed", str(4242 + i),
        "--ek-sampling", "fortran",
        "--format", "npz",
        "-o", out,
    ])
    times.append(time.perf_counter() - t0)
print(" ".join(f"{t:.2f}" for t in times))
EOF
}

# Warm the page cache on the MAID tables: every process reads them, so timing a
# cold first run would measure the filesystem rather than the layout.
echo "warmup: $(one_process "$NCPU" warmup)"
echo

TIMES=""
for k in 1 2 4 8 16 32 64 128; do
  [ "$k" -le "$NCPU" ] || break
  runs=$(one_process "$k" "t$k")
  echo "threads $k: $runs"
  # keep every repeat; the projection drops the first (compilation)
  TIMES="$TIMES $k:${runs// /,}"
done

echo
"$PY" - "$NCPU" "$TIMES" <<'EOF'
import math, sys

ncpu = int(sys.argv[1])
measured = {}
for field in sys.argv[2].split():
    k, runs = field.split(":")
    measured[int(k)] = [float(x) for x in runs.split(",")]

# steady state: the repeats after the compile.  min, because a run interrupted
# by another tenant should not be mistaken for the layout's cost.
steady = {k: min(v[1:]) if len(v) > 1 else v[0] for k, v in measured.items()}
base = steady.get(1) or min(steady.values())


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
print(f"steady-state run at 1 thread: {base:.1f}s")
print(f"{'layout':>14} {'threads':>8} {'waves':>6} {'projected':>10} {'vs 1x1':>8}")
best = None
for procs in (1, 2, 4, 8, 16, 32, 64, 128):
    if procs > ncpu:
        break
    threads = max(1, ncpu // procs)
    total = math.ceil(runs / procs) * interp(threads)
    print(f"{procs:>6} x {threads:<4}   {threads:>8} {math.ceil(runs / procs):>6} "
          f"{total / 60:>9.1f}m {base * runs / total:>7.1f}x")
    if best is None or total < best[0]:
        best = (total, procs, threads)
print()
print(f"best layout: {best[1]} processes x {best[2]} threads "
      f"({best[0] / 60:.1f} min for {runs} runs)")
print(f"set run_py_cpu.sbatch: --jobs {best[1]} --cpus-per-job {best[2]}")
EOF