#!/bin/bash
# Submit the whole Perlmutter scan: python venv (one-time), Fortran build,
# the grid, three batch jobs (Fortran, JAX-CPU, JAX-GPU), then a verification
# job that runs once all three finish.
#
# Usage (on a Perlmutter login node, from a clone of this repo):
#   validation/perlmutter/submit_all.sh
#
# Environment:
#   AAO_ACCOUNT      NERSC project account (default: m3792)
#   AAO_SCAN_ROOT    scan root (default: $SCRATCH/aao_rad_scan)
#   AAO_VENV         python venv (default: $SCRATCH/aao_venv)
#
# Everything is written under $SCRATCH -- /tmp is not shared between
# Perlmutter's login and compute nodes and is routinely scrubbed.

set -euo pipefail
: "${AAO_ACCOUNT:=m3792}"
: "${SCRATCH:?SCRATCH is not set -- run this on Perlmutter}"
HERE=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$HERE/../.." && pwd)
ROOT=${AAO_SCAN_ROOT:-$SCRATCH/aao_rad_scan}
VENV=${AAO_VENV:-$SCRATCH/aao_venv}
PYBIN="$VENV/bin/python"
mkdir -p "$ROOT/logs"

# 1. python with jax (CPU + GPU in one venv; the CUDA wheels are large, and
#    compute nodes have no internet, so this happens on the login node once).
#
#    Perlmutter's system python3 is 3.6 -- too old for the port (which needs
#    >= 3.10) and for jax -- and `module load python` does not reliably put an
#    interpreter on PATH on the login nodes, so the venv is built with uv
#    against a managed CPython.  Interpreter, uv cache and venv all go under
#    $SCRATCH: it is the only filesystem the compute nodes share with the login
#    node, and $HOME is small and slow.
UV=$(command -v uv || true)
if [ ! -x "$PYBIN" ]; then
  echo "== creating venv at $VENV (one-time; downloads ~2 GB of CUDA wheels)"
  export UV_PYTHON_INSTALL_DIR="$SCRATCH/uv/python"
  export UV_CACHE_DIR="$SCRATCH/uv/cache"
  mkdir -p "$UV_PYTHON_INSTALL_DIR" "$UV_CACHE_DIR" "$VENV"
  if [ -n "$UV" ]; then
    "$UV" venv --python 3.13 "$VENV"
    "$UV" pip install --python "$PYBIN" "jax[cuda12]" numpy scipy matplotlib
  else
    module load python || true
    python3 -m venv "$VENV" || python3 -m venv --without-pip "$VENV"
    "$VENV/bin/pip" install --upgrade pip
    "$VENV/bin/pip" install "jax[cuda12]" numpy scipy matplotlib
  fi
  "$PYBIN" - <<'EOF'
import sys
assert sys.version_info >= (3, 10), f"venv python too old: {sys.version}"
import jax, numpy, scipy  # noqa: F401
print(f"venv ok: python {sys.version.split()[0]}, jax {jax.__version__}")
EOF
fi

# 2. the Fortran binary, with the repo's own Makefile (gfortran).
if [ ! -x "$REPO/bin/aao_rad_lund" ]; then
  echo "== building aao_rad with gfortran"
  module load gcc
  make -C "$REPO"
fi

# 3. the scan grid: 96 configurations x 2 seeds (pure stdlib).
"$PYBIN" "$HERE/make_grid.py" --root "$ROOT"

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
