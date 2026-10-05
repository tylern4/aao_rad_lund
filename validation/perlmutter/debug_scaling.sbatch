#!/bin/bash
# How does the port scale with the *number of processes*, which nothing measured?
#
# What is known so far:
#   - one worker, uncontended: a 20,000-event run reaches 16000/20000 in 12
#     minutes and would finish in ~14 core-minutes (job 59342401).  It works.
#   - 128 workers on one node: every run clears the pilot, none emits a single
#     sampling-loop progress line, each banks ~14 CPU-seconds in 30 minutes of
#     wall, and the node sits 99% idle (job 59339555).
#   - the compiles are cheap.  Four chunk widths all finished 4,000 events in
#     53-59 s including compilation (job 59342635), so the 7m32s "very slow
#     compile" in the 128-worker run is contention, not a pathological module --
#     and --chunk-events, the suspect, is irrelevant.
#
# calibrate.sbatch measured the port scaling with *threads inside one process*
# (1 to 64).  Nothing ever measured it scaling with processes, and 128 XLA CPU
# runtimes on one node is the configuration that failed.  This walks the process
# count and reports what each level actually achieved, so the arm can be launched
# at a concurrency known to work rather than the largest one the node offers.
#
# The last case repeats the worst level with the persistent cache switched off.
# Every run in production shares one cache directory on $SCRATCH, and XLA takes
# file locks in it, so a shared filesystem lock is a plausible reason 128
# processes would serialise on compiles while the node looks idle.

#SBATCH --qos=debug
#SBATCH --job-name=aao_scaling
#SBATCH --account=m3792
#SBATCH --output=/pscratch/sd/t/tylern/aao_rad_scan/logs/aao_scaling_%j.out
# Slurm does not expand $SCRATCH in --output, and it cannot be $HOME (/tmp is not
# shared with the compute nodes and $HOME is 74% full).

#SBATCH --constraint=cpu
# The debug QoS refuses anything longer (QOSMaxWallDurationPerJobLimit), so this
# is 30 minutes, and the per-case cap below is set so the five cases fit inside
# it even if every one of them hangs to the cap.
#SBATCH --time=00:30:00

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
ROOT=${AAO_SCALING_ROOT:-$SCRATCH/aao_rad_scan_smoke}/scaling
CARD=$(ls "$SCRATCH"/aao_rad_scan_smoke/grid/cfg_*.txt 2>/dev/null | head -1)
if [ -z "$CARD" ]; then
  echo "no run card to work from -- run debug_py_cpu.sbatch once first" >&2
  exit 1
fi

export PYTHONPATH="$REPO/src${PYTHONPATH:+:$PYTHONPATH}"
export OMP_NUM_THREADS=1 JAX_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export JAX_PLATFORMS=cpu
unset XLA_FLAGS

echo "nproc: $(nproc)   threads visible"
echo "card: $CARD"
EVENTS=${AAO_SCALING_EVENTS:-500}

# run_case <label> <n> <cache-dir|none>
run_case() {
  local label=$1 n=$2 cachedir=$3
  local d="$ROOT/$label"
  rm -rf "$d"; mkdir -p "$d/cache" "$d/runs"

  if [ "$cachedir" = "none" ]; then
    unset AAO_JAX_CACHE_DIR
  else
    export AAO_JAX_CACHE_DIR="$cachedir"
  fi

  echo
  echo "== $label: $n concurrent processes, events $EVENTS, cache ${cachedir##*/}"
  local start=$SECONDS
  local pids=""
  for i in $(seq 1 "$n"); do
    mkdir -p "$d/runs/r$i"
    cp "$CARD" "$d/runs/r$i/run_card.txt"
    ( cd "$d/runs/r$i" && exec timeout "${CASE_CAP:-180}" "$PY" -m aao_rad.cli \
        --legacy-input run_card.txt --parms "$REPO/parms" \
        --n-events "$EVENTS" --seed $((8675309 + i)) \
        --ek-sampling fortran --format npz \
        --progress "$EVENTS" -o "$d/runs/r$i/out.npz" ) \
      > "$d/runs/r$i/out.txt" 2>&1 &
    pids="$pids $!"
  done

  local ok=0 bad=0
  for p in $pids; do
    if wait "$p"; then ok=$((ok + 1)); else bad=$((bad + 1)); fi
  done
  local took=$(( SECONDS - start ))
  local progressed
  progressed=$(grep -l 'events after' "$d"/runs/*/out.txt 2>/dev/null | wc -l)
  local core_s=0
  for f in "$d"/runs/*/out.txt; do
    [ -f "$f" ] || continue
    local t
    t=$(grep -o 'sigma (MC)' "$f" 2>/dev/null | wc -l)
    core_s=$((core_s + t))
  done
  printf "   %d/%d finished in %ss   reached the loop: %s   cache entries: %s\n" \
    "$ok" "$n" "$took" "$progressed" \
    "$(ls "$cachedir" 2>/dev/null | wc -l)"
}

# Shared cache directory, as production does.  Each level gets a fresh one so
# every case pays its own compiles, which is the situation a fresh arm is in.
FRESH="$ROOT/caches"
for n in ${AAO_SCALING_LEVELS:-1 16 64 128}; do
  run_case "shared_n$n" "$n" "$FRESH/n$n"
done
# The worst level again, with no persistent cache at all: separates "too many XLA
# runtimes" from "128 processes fighting over file locks in one cache".
run_case "nocache_n128" 128 none

echo
echo "== summary (finished / launched, wall seconds)"
for d in "$ROOT"/*/; do
  [ -d "$d" ] || continue
  printf "  %-16s %s\n" "$(basename "$d")" \
    "$(ls "$d"/runs/r*/out.npz 2>/dev/null | wc -l) npz written"
done