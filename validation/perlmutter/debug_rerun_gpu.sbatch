#!/bin/bash
# 30-minute smoke test for the reduced-quota rerun, JAX-on-GPU arm.
#
# The reduced quotas are sized at 0.9 * 2**31 trials -- about 1.9e9 -- which the
# serial Fortran needs 2.2 hours to walk and this arm needs about five minutes.
# So this is the arm that can prove a reduced-quota run *completes* inside the
# debug window: card taken, 15,080 events produced, cross section printed, and
# an n-tuple whose record count matches the quota the manifest now declares.
#
# That last part is the one worth having.  A quota mismatch between the card and
# the manifest would show up later as a run that looks short against its
# denominator, and that is a much more expensive thing to find in verification
# than here.
#
# `debug`, not `gpu_debug`.  A gpu_debug QoS exists and sacctmgr reports MaxWall
# 00:30:00, but every request against it is refused outright -- `sbatch
# --test-only` comes back with "Job request does not match any supported policy"
# for all of gres=gpu:1 and gres=gpu:4, with and without --constraint=gpu, with
# and without --exclusive, and for the bare QoS on its own.  `debug` accepts the
# same request unchanged and NERSC remaps it to gpu_debug with the m3792_g
# account on submission.  This header is the one bench_gpu.sbatch ran under.

#SBATCH --qos=debug
#SBATCH --job-name=aao_rerun_gpu
#SBATCH --account=m3792
#SBATCH --output=/pscratch/sd/t/tylern/aao_rad_scan/logs/aao_rerun_gpu_%j.out
# Slurm does not expand $SCRATCH in --output, which is read at submit time, so
# the path is spelled out.  Not $HOME (74% full) and not /tmp (not shared with
# the compute nodes).

#SBATCH --constraint=gpu
#SBATCH --gpus=4
#SBATCH --exclusive
#SBATCH --time=00:30:00

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
ROOT=${AAO_SCAN_ROOT:-${SCRATCH:?SCRATCH is not set}/aao_rad_scan}
PY=${AAO_PYTHON:-${AAO_VENV:-$SCRATCH/aao_venv}/bin/python}
export PYTHONUNBUFFERED=1
mkdir -p "$ROOT/logs"

if [ ! -x "$PY" ]; then
  echo "venv python not found at $PY -- run submit_all.sh first" >&2
  exit 1
fi

# JAX recompiles per configuration (the beam energy enters as a constant); a
# persistent cache lets the CPU and GPU arms share compile time.
export JAX_COMPILATION_CACHE_DIR="$ROOT/jax_cache"
# One worker per GPU, each given a single device by run_batch.py; without this
# every process would try to reserve most of every visible device.
export XLA_PYTHON_CLIENT_PREALLOCATE=false
unset XLA_FLAGS

summarise() {
  local d="$ROOT/py_gpu"
  echo
  echo "== rerun smoke (py_gpu) summary at $(date +%T)"
  for cfg in cfg_040 cfg_088; do
    echo "  $cfg quota: $(tail -2 "$ROOT/grid/$cfg.txt" 2>/dev/null | head -1) events"
  done
  echo "  launched (out.txt exists): $(ls "$d"/*/out.txt 2>/dev/null | wc -l)"
  echo "  past the pilot:            $(grep -l 'weight maximum' "$d"/*/out.txt 2>/dev/null | wc -l)"
  echo "  finished (sigma (MC)):     $(grep -l 'sigma (MC)' "$d"/*/out.txt 2>/dev/null | wc -l)"
  echo "  errored:                   $(grep -l -E 'Traceback|error:' "$d"/*/out.txt 2>/dev/null | wc -l)"
  echo "  XLA cache: $(ls "$ROOT/jax_cache" 2>/dev/null | wc -l) entries"
  echo "  furthest run: $(grep -h 'events after' "$d"/*/out.txt 2>/dev/null | tail -1)"
}
trap summarise EXIT

GPUS=${AAO_GPUS:-4}
"$PY" "$HERE/make_grid.py" --root "$ROOT"

SMOKE_CFGS="$ROOT/grid/smoke_configs.txt"
printf '%s\n' ${AAO_SMOKE_CFGS:-cfg_040 cfg_088} > "$SMOKE_CFGS"
echo "smoke: $(tr '\n' ' ' < "$SMOKE_CFGS")"
echo "smoke: $GPUS workers, one per GPU, $(( $(nproc) / GPUS )) host threads each"

"$PY" "$HERE/run_batch.py" --code py_gpu --root "$ROOT" --repo "$REPO" \
  --python "$PY" --gpus "$GPUS" --jobs "$GPUS" \
  --cpus-per-job "$(( $(nproc) / GPUS ))" \
  --only-cfg "@$SMOKE_CFGS" --limit "${AAO_SMOKE_LIMIT:-4}"
