#!/bin/bash
# Verification arm of the Perlmutter scan: runs after the Fortran, JAX-CPU
# and JAX-GPU batches (submit_all.sh sets the afterok dependency) and writes
# the statistical comparison under $SCRATCH/aao_rad_scan/verify/.

# `regular`, not `shared`, like the three scan arms; AAO_QOS overrides it.
#SBATCH --qos=regular
#SBATCH --job-name=aao_verify
#SBATCH --account=m3792
#SBATCH --output=/pscratch/sd/t/tylern/aao_rad_scan/logs/aao_verify_%j.out
# Slurm does not expand $SCRATCH here -- --output is read at submit time, before any
# of this script's own setup -- so the path is spelled out.  It cannot be $HOME or
# /tmp either: $HOME is 74% full, and /tmp is not shared with the compute nodes.
# Unlike the scan arms this one does ask for a size: it is a single process
# comparing 192 configurations' worth of n-tuples, so it wants a handful of
# cores, and taking a whole node for it would idle 124 of them.
#SBATCH --constraint=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=2-00:00:00

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
# Every one of these jobs can be killed at its wall limit, and Python block-buffers
# stdout when it is not a terminal.  A killed job then loses its own progress
# reporting -- the '[42/46] cfg_091_s0: ok (102s)' lines that are the only record of
# which runs finished and which were still in flight.  Unbuffer it for every
# invocation below, including the `python -m aao_rad.cli` runs: run_batch.py
# captures each run's output to its own out.txt, but a run killed mid-flight loses
# whatever its own stdout had not yet flushed, which is exactly the progress lines
# that distinguish a slow run from a stalled one.
export PYTHONUNBUFFERED=1

"$PY" "$HERE/verify_statistics.py" --root "$ROOT"
