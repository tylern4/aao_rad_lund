#!/bin/bash
# How long does the sampling loop's module take to compile, and what drives it?
#
# The JAX-CPU arm does not fail because of contention.  A single process on an
# otherwise idle node, with 0.9 of a core to itself, also stalls: 25.7 s of CPU
# in total, five of its six threads parked in futex_wait_queue, and nothing in
# the XLA cache.  Every run gets *past* the pilot -- "weight maximum estimated
# as ..." prints for all of them -- and no run ever emits a sampling-loop
# progress line, so the stall is at the first call of ``_make_stepper``'s
# jitted step rather than in the physics.
#
# The step writes accepted events into a ``(chunk_events, n_columns)`` buffer at
# a position carried across iterations, i.e. a dynamic update slice over a
# buffer that is 500,000 events wide by default.  That is the shape XLA's CPU
# backend is most likely to choke on, and its cost should scale with the buffer
# width.  So: run the same configuration at several --chunk-events values, one
# process, one core, and time how far each gets.
#
# Answer to look for: a chunk width small enough that the loop emits progress
# lines within the cap.  Anything that reaches "4000/4000" has compiled the
# step and is merely slow, which is a different fix from one that never gets
# there at all.

#SBATCH --qos=debug
#SBATCH --job-name=aao_jit_probe
#SBATCH --account=m3792
#SBATCH --output=/pscratch/sd/t/tylern/aao_rad_scan/logs/aao_jit_probe_%j.out
# Slurm does not expand $SCRATCH in --output.  It cannot be $HOME either: /tmp is
# not shared with the compute nodes and $HOME is 74% full.

#SBATCH --constraint=cpu
#SBATCH --time=00:25:00
# One core, on purpose: this is measuring a single compile, and a second busy
# process next to it is the variable being eliminated.

set -uo pipefail

REPO=""
for cand in "${AAO_REPO:-}" \
            "$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)/../.." \
            "${SLURM_SUBMIT_DIR:-}" \
            "${SCRATCH:-}/aao_rad_lund"; do
  if [ -n "$cand" ] && [ -f "$cand/pyproject.toml" ] && [ -d "$cand/src/aao_rad" ]; then
    REPO=$(cd "$cand" && pwd)
    break
  fi
done
if [ -z "$REPO" ]; then
  echo "cannot find the aao_rad_lund checkout -- set AAO_REPO" >&2
  exit 1
fi

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
WORK=${AAO_PROBE_ROOT:-$SCRATCH/aao_rad_scan_smoke}/jit_probe
mkdir -p "$WORK"

export PYTHONPATH="$REPO/src${PYTHONPATH:+:$PYTHONPATH}"
# Same as the batch driver gives a --cpus-per-job 1 worker, and the same cache it
# writes, so a compile that succeeds here is the one the arm will reuse.
export OMP_NUM_THREADS=1 JAX_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export JAX_PLATFORMS=cpu
export AAO_JAX_CACHE_DIR=${AAO_PROBE_ROOT:-$SCRATCH/aao_rad_scan_smoke}/jax_cache
unset XLA_FLAGS

# Any card from the generated grid will do -- the width being varied here is a
# CLI knob, not something the card carries.
CARD=$(ls "$SCRATCH"/aao_rad_scan_smoke/grid/cfg_*.txt 2>/dev/null | head -1)
if [ ! -f "$CARD" ]; then
  echo "no run card to work from -- run debug_py_cpu.sbatch once first" >&2
  exit 1
fi
echo "probe card: $CARD"
echo "nproc: $(nproc)   cache: $AAO_JAX_CACHE_DIR"

# CAP seconds per case, so one pathological width cannot eat the whole job.
CAP=${AAO_PROBE_CAP:-180}
EVENTS=${AAO_PROBE_EVENTS:-4000}

for chunk in ${AAO_PROBE_CHUNKS:-500000 100000 20000 4000}; do
  d="$WORK/chunk_$chunk"
  rm -rf "$d"; mkdir -p "$d"
  cp "$CARD" "$d/run_card.txt"
  echo
  echo "== chunk-events $chunk (buffer $(($chunk * 32 * 4 / 1048576)) MiB), cap ${CAP}s"
  start=$(date +%s)
  ( cd "$d" && timeout "$CAP" "$PY" -m aao_rad.cli \
      --legacy-input run_card.txt \
      --parms "$REPO/parms" \
      --n-events "$EVENTS" --seed 8675309 \
      --ek-sampling fortran --format npz \
      --chunk-events "$chunk" --progress 500 \
      -o "$d/out.npz" ) > "$d/out.txt" 2>&1
  rc=$?
  took=$(( $(date +%s) - start ))
  pilot=$(grep -c 'weight maximum' "$d/out.txt")
  lines=$(grep -c 'events after' "$d/out.txt")
  last=$(grep -h 'events after' "$d/out.txt" | tail -1)
  printf "   exit %s in %ss  past-pilot %s  progress-lines %s\n" "$rc" "$took" "$pilot" "$lines"
  [ -n "$last" ] && echo "   last: $last"
  if [ "$rc" -eq 124 ]; then
    echo "   -> hit the ${CAP}s cap without reaching the loop"
  elif [ "$lines" -gt 0 ]; then
    echo "   -> compiled the step and sampled; this width works"
  fi
done

echo
echo "== summary: chunk width vs time to first progress line"
grep -l . "$WORK"/*/out.txt 2>/dev/null | while read f; do
  printf "  %-28s %s\n" "$(basename "$(dirname "$f")")" \
    "$(grep -m1 'events after' "$f" || echo 'never reached the loop')"
done