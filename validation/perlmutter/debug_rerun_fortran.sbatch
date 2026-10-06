#!/bin/bash
# 30-minute smoke test for the reduced-quota rerun, Fortran arm.
#
# plan_ref_rerun.py cut 24 configurations' event quotas so integer*4 ntries stays
# inside 90% of 2**31 -- cfg_040 at 15,080 events against the original 20,000.
# That quota reaches the binary only through the run card, and build_job copies
# the card fresh for exactly this reason: a copy left at 20,000 wraps the counter
# all over again and prints a cross section that looks no different from the
# broken one.  The only proof it arrived is a run that starts on it.
#
# These configurations cannot finish here.  Every reduced quota is the same
# 1.9e9 trials (0.9 * 2**31), and at ~240,000 trials/s that is about 2.2 hours of
# serial work against a 30-minute limit.  So this job proves the card, the
# dispatch and that the binary runs on the new nmax -- not that a run completes.
# Completion at the reduced quota is the GPU arm's job: it does the same 1.9e9
# trials in a few minutes (debug_rerun_gpu.sbatch).
#
# It runs against the production root on purpose.  These are configurations
# whose Fortran reference is unusable anyway, run_batch clears their stale
# outputs before it starts, and the production rerun re-runs whatever this
# leaves part-written.
#
# The collision check is left to the production job: it md5-compares the two
# seeds of a configuration, and here neither seed has finished.

#SBATCH --qos=debug
#SBATCH --job-name=aao_rerun_fortran
#SBATCH --account=m3792
#SBATCH --output=/pscratch/sd/t/tylern/aao_rad_scan/logs/aao_rerun_fortran_%j.out
# Slurm does not expand $SCRATCH in --output: it is read at submit time, before
# any of this script's setup runs.  Not $HOME (74% full) and not /tmp (not
# shared with the compute nodes).

#SBATCH --constraint=cpu
#SBATCH --exclusive
#SBATCH --time=00:30:00
# No --nodes/--ntasks/--cpus-per-task on purpose.  A probe with no directives at
# all (debug 59300339) came back with 256 CPUs and 487 GB, so asking for a size
# here would buy half the machine.  The worker count is read from nproc at run
# time, exactly as the production arm reads it.

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
# The job may be killed at the wall limit and Python block-buffers stdout when
# it is not a terminal, so unbuffer every invocation: the progress lines in each
# out.txt are the only record of how far a run got.
export PYTHONUNBUFFERED=1
mkdir -p "$ROOT/logs"

if [ ! -x "$PY" ]; then
  echo "venv python not found at $PY -- run submit_all.sh first" >&2
  exit 1
fi

command -v gfortran >/dev/null || module load gcc || true
if [ ! -x "$REPO/bin/aao_rad_lund" ]; then
  echo "== building aao_rad (gfortran)"
  make -C "$REPO"
fi

# The kill at 30 minutes is expected, so the report has to survive it.
summarise() {
  local d="$ROOT/fortran"
  echo
  echo "== rerun smoke (fortran) summary at $(date +%T)"
  for cfg in cfg_040 cfg_088; do
    echo "  $cfg quota: $(tail -2 "$ROOT/grid/$cfg.txt" 2>/dev/null | head -1) events"
  done
  echo "  launched (out.txt exists):  $(ls "$d"/*/out.txt 2>/dev/null | wc -l)"
  echo "  checkpointed:               $(grep -l 'ntries, nevent' "$d"/*/out.txt 2>/dev/null | wc -l)"
  echo "  finished (cross section):   $(grep -l 'Integrated cross section' "$d"/*/out.txt 2>/dev/null | wc -l)"
  echo "  errored:                    $(grep -l -E 'Fortran runtime error|SIGFPE|Traceback' "$d"/*/out.txt 2>/dev/null | wc -l)"
  echo "  furthest checkpoint: $(grep -h 'ntries, nevent' "$d"/*/out.txt 2>/dev/null | tail -1)"
}
trap summarise EXIT

"$PY" "$HERE/make_grid.py" --root "$ROOT"

# A file rather than a comma list, because this is the form the production job
# uses and the one worth exercising: --only-cfg @path.
SMOKE_CFGS="$ROOT/grid/smoke_configs.txt"
printf '%s\n' ${AAO_SMOKE_CFGS:-cfg_040 cfg_088} > "$SMOKE_CFGS"
echo "smoke: $(tr '\n' ' ' < "$SMOKE_CFGS")"
echo "smoke: $(nproc) workers, one hardware thread each"

"$PY" "$HERE/run_batch.py" --code fortran --root "$ROOT" --repo "$REPO" \
  --jobs "$(nproc)" --cpus-per-job 1 \
  --only-cfg "@$SMOKE_CFGS" --limit "${AAO_SMOKE_LIMIT:-4}"
