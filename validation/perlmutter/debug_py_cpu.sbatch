#!/bin/bash
# 30-minute smoke test for the JAX-CPU arm.
#
# Job 59322722 put 256 JAX processes on one 128-core node and every one of them
# ended up blocked in a futex at 0% CPU with the node 99.4% idle, having burned
# about two CPU-minutes each in six hours.  A whole node to find that out is a
# bad way to find it out, so this runs the same launch path at the same worker
# count for half an hour and reports what it managed to do.
#
# It uses its own scan root so it cannot contaminate the production arm's
# results: these runs are killed at the wall limit part-way through, and a
# half-written out.npz left next to the real ones is exactly the kind of thing
# that quietly poisons a later verification.

# `debug`, which is a short queue.  The production arms are `regular`.
#SBATCH --qos=debug
#SBATCH --job-name=aao_py_cpu_smoke
#SBATCH --account=m3792
#SBATCH --output=/pscratch/sd/t/tylern/aao_rad_scan/logs/aao_py_cpu_smoke_%j.out
# Slurm does not expand $SCRATCH in --output: it is read at submit time, before
# any of this script's setup runs.  It cannot be $HOME (/tmp is not shared with
# the compute nodes and $HOME is 74% full).

#SBATCH --constraint=cpu
#SBATCH --time=00:30:00
# No --nodes/--ntasks/--cpus-per-task on purpose: a probe with no directives at
# all (debug 59300339) came back with 256 CPUs and 487 GB, so asking for a size
# here would buy half the machine.  The worker count below is read from nproc at
# run time, so this does not have to know whether the node is hyperthreaded.

set -uo pipefail

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
ROOT=${AAO_SMOKE_ROOT:-${SCRATCH:?SCRATCH is not set}/aao_rad_scan_smoke}
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

# The fan-out gets killed at the wall limit, so the report has to survive that.
summarise() {
  local d="$ROOT/py_cpu"
  echo
  echo "== smoke summary at $(date +%T)"
  echo "  launched (out.txt exists): $(ls "$d"/*/out.txt 2>/dev/null | wc -l)"
  echo "  finished (sigma (MC)):    $(grep -l 'sigma (MC)' "$d"/*/out.txt 2>/dev/null | wc -l)"
  echo "  moving (progress line):   $(grep -l 'events after' "$d"/*/out.txt 2>/dev/null | wc -l)"
  echo "  past the pilot:           $(grep -l 'weight maximum' "$d"/*/out.txt 2>/dev/null | wc -l)"
  echo "  hit a slow-op alarm:      $(grep -l 'slow_operation_alarm' "$d"/*/out.txt 2>/dev/null | wc -l)"
  echo "  errored:                  $(grep -l -E 'Traceback|error:' "$d"/*/out.txt 2>/dev/null | wc -l)"
  echo "  XLA cache: $(ls "$ROOT/jax_cache" 2>/dev/null | wc -l) entries, $(du -sh "$ROOT/jax_cache" 2>/dev/null | cut -f1)"
  echo "  furthest run: $(grep -h 'events after' "$d"/*/out.txt 2>/dev/null | tail -1)"
  echo "  per-run CPU seconds, top 5:"
  for p in $(pgrep -f 'aao_rad.cli' 2>/dev/null | head -5); do
    awk '{printf "    pid %s  cpu %.1fs  %s\n", $1, ($14+$15)/100, $3}' "/proc/$p/stat" 2>/dev/null
  done
}
trap summarise EXIT

command -v gfortran >/dev/null || module load gcc || true
unset XLA_FLAGS

"$PY" "$HERE/make_grid.py" --root "$ROOT"

# One worker per core rather than per hardware thread, matching run_py_cpu.sbatch:
# 256 runs on 128 cores is both the oversubscription the calibration calls a loss
# and 384 GB of resident JAX on a 503 GB node.  AAO_SMOKE_JOBS overrides it.
WORKERS=${AAO_SMOKE_JOBS:-$(($(nproc) / 2))}
echo "smoke: $WORKERS workers (one core each), $(nproc) threads visible"
echo "smoke: wall limit $(scontrol show job "$SLURM_JOB_ID" 2>/dev/null | grep -o 'TimeLimit=[^ ]*')"

"$PY" "$HERE/run_batch.py" --code py_cpu --root "$ROOT" --repo "$REPO" \
  --python "$PY" --jobs "$WORKERS" --cpus-per-job 1 \
  --limit "${AAO_SMOKE_LIMIT:-128}"