#!/bin/bash
# Submit the whole Perlmutter scan: python venv (one-time), Fortran build,
# the grid, three batch jobs (Fortran, JAX-CPU, JAX-GPU), then a verification
# job that runs once all three finish.
#
# Usage (on a Perlmutter login node, from a clone of this repo):
#   export AAO_ACCOUNT=<your NERSC project account>
#   validation/perlmutter/submit_all.sh
#
# Environment:
#   AAO_ACCOUNT      NERSC project account (required)
#   AAO_SCAN_ROOT    scan root (default: $SCRATCH/aao_rad_scan)
#   AAO_VENV         python venv (default: $SCRATCH/aao_venv)
#
# Everything is written under $SCRATCH -- /tmp is not shared between
# Perlmutter's login and compute nodes and is routinely scrubbed.

set -euo pipefail
: "${AAO_ACCOUNT:?export AAO_ACCOUNT=<your NERSC project account>}"
: "${SCRATCH:?SCRATCH is not set -- run this on Perlmutter}"
HERE=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$HERE/../.." && pwd)
ROOT=${AAO_SCAN_ROOT:-$SCRATCH/aao_rad_scan}
VENV=${AAO_VENV:-$SCRATCH/aao_venv}
PYBIN="$VENV/bin/python"
mkdir -p "$ROOT/logs"

# 1. python with jax (CPU + GPU in one venv; the CUDA wheels are large, and
#    compute nodes have no internet, so this happens on the login node once).
if [ ! -x "$PYBIN" ]; then
  echo "== creating venv at $VENV (one-time; downloads ~2 GB of CUDA wheels)"
  module load python
  python3 -m venv "$VENV"
  "$VENV/bin/pip" install --upgrade pip
  "$VENV/bin/pip" install "jax[cuda12]" numpy scipy matplotlib
fi

# 2. the Fortran binary, with the repo's own Makefile (gfortran).
if [ ! -x "$REPO/bin/aao_rad_lund" ]; then
  echo "== building aao_rad with gfortran"
  module load gcc
  make -C "$REPO"
fi

# 3. the scan grid: 96 configurations x 2 seeds (pure stdlib).
module load python
python3 "$HERE/make_grid.py" --root "$ROOT"

# 4. the three arms run in parallel; verification runs after all of them.
sub() { sbatch --account="$AAO_ACCOUNT" "$@"; }
FORT=$(sub -o "$ROOT/logs/fortran_%j.out" "$HERE/run_fortran.sbatch" | awk '{print $4}')
CPU=$(sub -o "$ROOT/logs/py_cpu_%j.out" "$HERE/run_py_cpu.sbatch" | awk '{print $4}')
GPU=$(sub -o "$ROOT/logs/py_gpu_%j.out" "$HERE/run_py_gpu.sbatch" | awk '{print $4}')
VER=$(sub -o "$ROOT/logs/verify_%j.out" \
      -d afterok:"$FORT":"$CPU":"$GPU" "$HERE/run_verify.sbatch" | awk '{print $4}')

echo "submitted:"
echo "  fortran  $FORT   -> $ROOT/logs/fortran_$FORT.out"
echo "  py_cpu   $CPU    -> $ROOT/logs/py_cpu_$CPU.out"
echo "  py_gpu   $GPU    -> $ROOT/logs/py_gpu_$GPU.out"
echo "  verify   $VER    -> $ROOT/logs/verify_$VER.out (after the three above)"
echo "results:     $ROOT/verify/"
