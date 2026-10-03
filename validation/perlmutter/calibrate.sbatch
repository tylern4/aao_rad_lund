#!/bin/bash
# How does a single run scale with threads on this node, and which node layout
# finishes the 192-run grid first?
#
# Every run is one fixed 65,536-trial vectorised pass (generate.stream sizes its
# chunk as max(min(chunk_events, n_events), batch_size), so a run costs the same
# whatever it keeps), which makes this a clean measurement: time one real run at
# several CPU counts, then project the wall time for the candidate processes x
# threads layouts.  Thread counts come from taskset, because XLA sizes its pool
# from sched_getaffinity.
#
# Run it before the scan (it takes a few minutes) and set run_py_cpu.sbatch's
# --jobs/--cpus-per-job to whatever wins here.
#
#   sbatch --account=m3792 -o $SCRATCH/aao_rad_scan/logs/calib_%j.out \
#          calibrate.sbatch
#SBATCH --job-name=aao_calib
#SBATCH --constraint=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=00:30:00

set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$HERE/../.." && pwd)
ROOT=${AAO_SCAN_ROOT:-${SCRATCH:?SCRATCH is not set}/aao_rad_scan}
PY=${AAO_PYTHON:-${AAO_VENV:-$SCRATCH/aao_venv}/bin/python}
export PYTHONPATH="$REPO/src${PYTHONPATH:+:$PYTHONPATH}"

if [ ! -x "$PY" ]; then
  echo "venv python not found at $PY -- run submit_all.sh first" >&2
  exit 1
fi

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
echo "card  $CARD"
echo "node  $NCPU cpus, lowest is $FIRST"
echo

one_run() {  # $1 = cpu count, $2 = label; prints seconds
  local cpus t0 t1
  cpus=$("$PY" -c "print(','.join(str($FIRST + i) for i in range($1)))")
  t0=$(date +%s)
  taskset -c "$cpus" "$PY" -m aao_rad.cli --legacy-input "$CARD" \
    --parms "$REPO/parms" --n-events 20000 --seed 4242 \
    --ek-sampling fortran --format npz -o "$ROOT/calib_$2.npz" \
    > "$ROOT/calib_$2.log" 2>&1
  t1=$(date +%s)
  echo $((t1 - t0))
}

# Warm up first: the jax compile is a fixed cost that must not land in a
# thread-scaling number.
echo "warmup (compiles): $(one_run "$NCPU" warmup)s"
echo

TIMES=""
for k in 1 2 4 8 16 32 64 128; do
  [ "$k" -le "$NCPU" ] || break
  s=$(one_run "$k" "t$k")
  echo "threads $k: $s s"
  TIMES="$TIMES $k $s"
done

echo
"$PY" - "$NCPU" "$TIMES" <<'EOF'
import math, sys

ncpu = int(sys.argv[1])
nums = [int(x) for x in sys.argv[2].split()]
measured = dict(zip(nums[::2], nums[1::2]))
base = measured.get(1) or min(measured.values())


def interp(k):
    """Seconds for k threads, log-interpolated between measured points."""
    ks = sorted(measured)
    if k <= ks[0]:
        return measured[ks[0]] * ks[0] / max(k, 1)
    if k >= ks[-1]:
        return measured[ks[-1]]
    lo = max(x for x in ks if x <= k)
    hi = min(x for x in ks if x >= k)
    if lo == hi:
        return measured[lo]
    f = (math.log(k) - math.log(lo)) / (math.log(hi) - math.log(lo))
    return measured[lo] * (measured[hi] / measured[lo]) ** f


runs = 192
print(f"one run, 1 thread: {base}s")
print(f"{'layout':>16} {'threads':>8} {'waves':>6} {'projected':>10} {'speedup':>8}")
best = None
for procs in (1, 2, 4, 8, 16, 32, 64, 128):
    if procs > ncpu:
        break
    threads = max(1, ncpu // procs)
    per = interp(threads)
    waves = math.ceil(runs / procs)
    total = waves * per
    speed = base * waves / total if total else 0
    print(f"{procs:>6} x {threads:<3} cpu {threads:>8} {waves:>6} {total / 60:>9.1f}m {speed:>7.2f}x")
    if best is None or total < best[0]:
        best = (total, procs, threads)
print()
print(f"best layout: {best[1]} processes x {best[2]} threads "
      f"({best[0] / 60:.1f} min for {runs} runs)")
print(f"set run_py_cpu.sbatch: --jobs {best[1]} --cpus-per-job {best[2]}")
EOF