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
# For short experiments this is `debug` rather than the `regular` the scan arms
# use, deliberately.  The regular CPU pool is often saturated enough that a
# whole-node request waits a long time, while debug is a 30-minute queue that
# starts in seconds and still gets the whole node -- the first calibration run
# this way took 4s to start and got 128 cores, and the corrected run took a few
# minutes to get all 256 threads to itself.  The measurement has to fit in those
# 30 minutes, which the short runs above are sized to do.
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
#SBATCH --exclusive
#SBATCH --time=00:28:00
# No --nodes/--ntasks/--cpus-per-task on purpose, for the reason the scan arms
# give: --cpus-per-task=128 is half a PM-CPU node, not a whole one.  Asking for
# nothing gets the node whole and --exclusive keeps a co-tenant off it, which
# matters more here than anywhere -- a competing run would land in the middle of
# a timing point.  The earlier calibration was measured on 128 of the node's 256
# hardware threads, and its projection assumed that was the whole node.

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
# ...and how many of those entries are distinct cores.  The two differ wherever
# the node is hyperthreaded, and the difference is the whole ball game: a run
# given N threads by taskset_prefix gets them on N *different* cores, so the
# node's capacity is its core count, not its thread count.
CORES=$("$PY" -c 'import os, pathlib; mask=os.sched_getaffinity(0); keys=set()
for c in mask:
    try:
        raw=pathlib.Path(f"/sys/devices/system/cpu/cpu{c}/topology/thread_siblings_list").read_text()
    except OSError:
        keys.add(c); continue
    sib=set()
    for part in raw.strip().split(","):
        part=part.strip()
        if not part: continue
        if "-" in part:
            lo,hi=part.split("-"); sib.update(range(int(lo),int(hi)+1))
        else: sib.add(int(part))
    keys.add(min(sib) if sib else c)
print(len(keys))')
FIRST=$("$PY" -c 'import os; print(min(os.sched_getaffinity(0)))')
echo "card   $CARD"
echo "node   $NCPU threads in the affinity mask ($CORES cores), lowest is $FIRST"
echo "timing $CALIB_EVENTS events x $REPEATS per thread count, grid uses $PROD_EVENTS"
echo "scale  ${THREAD_COUNTS}"
echo

# $1 = cpu count, $2 = label.  Prints the per-run seconds, space separated, or
# the word TIMEOUT if the label exceeded LABEL_TIMEOUT.
one_process() {
  local cpus
  # "N threads" has to mean what --cpus-per-job N means to the arm being
  # calibrated, and that is N *different cores*: taskset_prefix hands CPUs out
  # round-robin over cores precisely so that no worker is given the other
  # hyperthread of a core another worker already has.  Striding the mask does
  # not do that.  On a 256-entry mask a stride of 128 picks entries 0 and 128,
  # which are the two threads of core 0, so the previous table's "2 threads"
  # measured one core's two hyperthreads and its "4 threads" measured two cores
  # -- a different quantity from the one the arm would use, and the reason its
  # projection came out wrong.  Group by core and take one thread per core.
  cpus=$("$PY" - "$1" "$HERE" <<'EOF'
import os, sys

sys.path.insert(0, sys.argv[2])
from run_batch import core_groups

groups = core_groups(sorted(os.sched_getaffinity(0)))
picked = [g[0] for g in groups[: int(sys.argv[1])]]
if not picked:
    sys.exit("no cores in the affinity mask")
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
echo "warmup: $(one_process "$CORES" warmup)"
echo

for k in $THREAD_COUNTS; do
  # Cores, not threads: one_process takes one hardware thread from each of N
  # distinct cores, which is what taskset_prefix does and therefore what the arm
  # under test will actually do.  A count above the core count would start
  # handing out a second thread of a core already in the slice -- measuring
  # hyperthreading, which the arms never ask for.
  [ "$k" -le "$CORES" ] || { echo "threads $k: skipped, node has $CORES cores"; continue; }
  runs=$(one_process "$k" "t$k")
  echo "threads $k: $runs"
  # Record as we go.  A label that timed out is written out as such and dropped
  # by the projection rather than guessed at: its true cost is unknown, and only
  # knowing "at least this slow" is not enough to interpolate from.
  printf 't%s\t%s\n' "$k" "${runs// /,}" >> "$RESULTS"
done

echo
"$PY" - "$NCPU" "$CORES" "$CALIB_EVENTS" "$PROD_EVENTS" "$RESULTS" <<'EOF'
import math, sys

ncpu, cores = int(sys.argv[1]), int(sys.argv[2])
calib_events, prod_events = int(sys.argv[3]), int(sys.argv[4])
results_path = sys.argv[5]

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
else:
    print(f"steady-state run at {calib_events} events, 1 thread: {one:.1f}s")
    print(f"scaled to the grid's {prod_events} events: {one * scale:.1f}s")
print(f"node: {ncpu} threads, {cores} cores")
print()


def wall(procs, threads):
    """Seconds to finish `runs` runs from `procs` processes of `threads` each.

    Two effects, and a layout needs both.  Waves: there are only `runs` of them,
    so a layout with fewer processes than runs pays for the last partial wave.
    Contention: `procs * threads` is more slots than the node has cores, and
    every run in the overfull layout then gets a fraction of the core it asked
    for.  Ignoring contention is what made the previous projection recommend 128
    processes x 2 threads -- 256 slots on 128 cores, twice oversubscribed -- as
    though the extra slots were free.
    """
    waves = math.ceil(runs / procs)
    return waves * interp(threads) * scale * max(1.0, procs * threads / cores)


# Every measured thread count is a candidate width, not just the one that fills
# the node at a given process count: which width wins is not obvious in advance.
# One row per process count, showing the width that suits it best.
widths = sorted(k for k in steady if k <= cores)
print(f"{'procs':>6} {'threads':>8} {'waves':>6} {'slots':>6} {'projected':>10}")
rows = []
for procs in (1, 2, 4, 8, 16, 32, 64, 128, runs):
    if procs > runs and procs > cores:
        continue
    total, threads = min((wall(procs, t), t) for t in widths)
    waves = math.ceil(runs / procs)
    print(f"{procs:>6} x {threads:<4}   {waves:>6} {procs * threads:>6} "
          f"{total / 60:>9.1f}m")
    rows.append((total, procs, threads))
rows.sort()
print()
best_total, best_procs, best_threads = rows[0]
print(f"best layout: {best_procs} processes x {best_threads} threads "
      f"({best_total / 60:.1f} min for {runs} runs on {cores} cores)")
if best_procs >= runs:
    print(f"that is every run in flight at once.  Workers beyond the run count are")
    print(f"idle, so run_py_cpu.sbatch's default of AAO_CPU_JOBS=$(nproc) already")
    print(f"gives it -- set AAO_CPU_CPUS={best_threads} and leave the count alone.")
else:
    print(f"set run_py_cpu.sbatch: AAO_CPU_JOBS={best_procs} "
          f"AAO_CPU_CPUS={best_threads}")
EOF
