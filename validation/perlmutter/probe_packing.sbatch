#!/bin/bash
# Packing probe: does one Fortran process per physical core keep the rate the
# speedup and the node-hours are computed from?
#
# bench_fortran.sbatch put four runs on a 128-core node.  Four runs never
# contend with each other or with the node's memory system, so what it measured
# is a single idle core -- and bench_report.py then divides that wall by 128 to
# say what a run costs.  Nothing has tested the division.
#
# Production is the evidence that it might not hold: the arm ran across three
# cancelled jobs (14:27 + 5:04 + 10:44 = 30 node-hours) with 192 processes on
# 128 cores, two per core, and never once reported a clean completion.  That
# history cannot be used to correct the rate, though.  Each job restarted the
# runs the previous one left unfinished, so its wall clocks mix whatever
# contention existed with thrown-away work and an hour-long queue gap -- there
# is no measurement in there to subtract.
#
# So this runs the whole 192-run grid on one node with the same quota
# make_probe_quota.py sized for every configuration, one process per physical
# core, and times it with the driver's own clock (run_batch.execute).  The
# comparison is against the eight-run bench on the same node type:
#
#   ratio = median projected 20,000-event wall here / the same there
#
# which is 1.0 if packing a node costs nothing, and is exactly the factor the
# node-hours arithmetic has been missing.  It is deliberately 128 in flight
# rather than production's 192: the node-hours this reports are for one run per
# core, so that is the layout that has to be measured.  Oversubscribing to 192
# would measure a packing nobody should use and put an unlabelled number
# against a labelled one.
#
# The quota matters as much as the layout.  A uniform event quota would run
# delta 0.005 (median 58,566 trials per event) six times longer than delta 0.05
# (9,118), so the result would be decided by whichever configurations happen to
# be expensive.  Every configuration is sized to the same 20 s of uncontended
# work instead -- about 1,600 s of node time for the whole grid if packing costs
# nothing, and still inside the 30-minute debug cap if it costs 40x.

#SBATCH --qos=debug
#SBATCH --job-name=aao_pack128
#SBATCH --account=m3792
#SBATCH --output=/pscratch/sd/t/tylern/aao_rad_scan/logs/aao_pack128_%j.out
# Slurm does not expand $SCRATCH in --output, which is read at submit time, so
# the path is spelled out.  Not $HOME (74% full) and not /tmp (not shared with
# the compute nodes).
#SBATCH --constraint=cpu
# A whole node: 128 runs, one per physical core, is the point of the exercise.
# Under `shared` --exclusive is refused outright.
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
PROBE="$ROOT/pack128"
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

# The production scan, for each configuration's trials per event.  Sized from
# out.txt rather than a stored table so the probe follows the scan it is meant
# to reproduce, including the 23 configurations whose counter wraps.
SCAN=${AAO_SCAN_ROOT:-$SCRATCH/aao_rad_scan}

echo "== sizing 192 runs at 20s of uncontended work each"
"$PY" "$HERE/make_probe_quota.py" --scan "$SCAN" --out "$PROBE/grid/n_events_override.csv"

"$PY" "$HERE/make_grid.py" --root "$PROBE"

# One process per physical core, which is what the node-hours divide by.  No
# --only-cfg: the root is fresh, so select_pending takes all 192 runs.
echo "pack: $(nproc) hardware threads, 128 in flight"
rc=0
"$PY" "$HERE/run_batch.py" --code fortran --root "$PROBE" --repo "$REPO" \
  --python "$PY" --jobs 128 --cpus-per-job 1 --warmup-events 0 || rc=$?

# Report even when the driver failed or the wall limit cut the job short: runs
# that finished carry their own timing.json, and a partial answer beats an
# empty log.
mkdir -p "$ROOT/timing"
PROBE_JSON="$ROOT/timing/pack128.json"
REF_JSON="$ROOT/timing/bench8.json"
"$PY" "$HERE/bench_report.py" --label "fortran, 192 runs at one per core" \
  --concurrency 128 --json "$PROBE_JSON" "$PROBE/fortran" || true

# The eight-run bench, re-read from its own timing.json, as the uncontended
# baseline for the ratio below.
"$PY" "$HERE/bench_report.py" --label "reference, 8 runs on the same node type" \
  --concurrency 128 --json "$REF_JSON" "$ROOT/reduced/fortran" || true

"$PY" - "$PROBE_JSON" "$REF_JSON" <<'PY'
import json
import sys
from pathlib import Path

probe_p, ref_p = (Path(p) for p in sys.argv[1:3])
if not probe_p.is_file() or not ref_p.is_file():
    print("\n== packing penalty: no summary to compare")
    raise SystemExit(0)
probe, ref = json.loads(probe_p.read_text()), json.loads(ref_p.read_text())
print("\n== packing penalty (median wall projected to 20,000 events)")
for tag, s in (("8 runs, uncontended", ref), ("192 runs, one per core", probe)):
    proj = s.get("median_proj")
    print(f"  {tag:<24} {proj:.0f} s" if proj else f"  {tag:<24} n/a")
a, b = ref.get("median_proj"), probe.get("median_proj")
if a and b:
    print(f"  ratio {b / a:.2f} x   (1.0 means packing a node costs nothing)")
for tag, s in (("one per core (128)", ref), ("this layout (128)", probe)):
    nh = s.get("node_hours")
    print(f"  node-h per run, {tag:<20} {nh:.6f}" if nh else f"  node-h per run, {tag:<20} n/a")
if (nh := probe.get("node_hours")):
    print(f"  extrapolated 192-run grid at this rate: {nh * 192:.2f} node-h")
PY

echo
echo "== elapsed for the job: $SECONDS s"
exit $rc
