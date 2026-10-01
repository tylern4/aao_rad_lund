"""Validate the Python table reader against the Fortran reference reader.

``dump_table`` (Fortran, see ``dump_table.f90``) reproduces the read sequence of
``src/read_sf_file.f90`` for a single grid point and prints the 62 structure
functions at single precision.  This script drives it over a sample of grid
points and requires exact equality with what :mod:`aao_rad.table` parses.

Run from the repository root::

    python validation/compare_table_to_fortran.py

Exits non-zero if any amplitude disagrees.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
DUMP = HERE / "dump_table"
DEFAULT_TBL = HERE.parent / "parms" / "spp_tbl" / "maid07-PNpi.tbl"


def run_fortran(tbl: Path, iq2: int, iw: int) -> np.ndarray:
    """Return the 62 amplitudes the Fortran reader sees at 1-based (iq2, iw)."""
    out = subprocess.run(
        [str(DUMP), str(tbl), str(iq2), str(iw)],
        capture_output=True,
        text=True,
        check=True,
    ).stdout
    vals = np.empty(62)
    seen = 0
    for line in out.splitlines():
        toks = line.split()
        # Fortran writes '(a,i3,1x,es24.16)', i.e. "SF   1   1.3727999...E+01"
        if len(toks) == 3 and toks[0] == "SF" and toks[1].isdigit():
            vals[int(toks[1]) - 1] = float(toks[2])
            seen += 1
    if seen != 62:
        raise RuntimeError(
            f"fortran dump returned {seen}/62 amplitudes for ({iq2},{iw})"
        )
    return vals


def main(argv: list[str]) -> int:
    if not DUMP.is_file():
        print(f"error: {DUMP} not built.\n"
              f"  cd validation && gfortran -O2 -o dump_table dump_table.f90",
              file=sys.stderr)
        return 2

    from aao_rad.table import load_table

    # Usage: compare_table_to_fortran.py [--save] [TBL] [CHANNEL]
    save = "--save" in argv
    rest = [a for a in argv[1:] if a != "--save"]
    tbl = Path(rest[0]) if rest else DEFAULT_TBL
    channel = int(rest[1]) if len(rest) > 1 else 3
    table = load_table(channel, parms_dir=tbl.parent.parent)

    # Grid points to check: both corners, the centre, and a spread of interior
    # points (including the upper-right corner where interpolation is hardest).
    nq2, nw = table.q2.size, table.w.size
    points = [
        (1, 1), (1, nw), (nq2, 1), (nq2, nw),
        (nq2 // 2, nw // 2), (nq2 // 2, 1), (1, nw // 2),
        (nq2, nw // 2), (nq2 // 2, nw),
        (nq2 // 4, nw // 3), (3 * nq2 // 4, 2 * nw // 3),
    ]

    print(f"table   : {tbl}")
    print(f"grid    : {nq2} x {nw}  Q2={table.q2[0]:g}..{table.q2[-1]:g} "
          f"W={table.w[0]:g}..{table.w[-1]:g}")
    print()

    worst = 0.0
    worst_at = ""
    n_checked = 0
    for iq2, iw in points:
        if not (1 <= iq2 <= nq2 and 1 <= iw <= nw):
            continue
        ref = run_fortran(tbl, iq2, iw)
        got = table.amps[:, iq2 - 1, iw - 1]
        n_checked += 1

        # The Fortran stores SF as REAL (single precision); the table file only
        # carries ~5 significant digits, so float32 is the common denominator.
        ref32 = ref.astype(np.float32).astype(np.float64)
        diff = np.abs(got - ref32)
        scale = np.maximum(np.abs(ref32), 1e-30)
        rel = diff / scale
        k = int(np.argmax(rel))
        if rel[k] > worst:
            worst = float(rel[k])
            worst_at = f"Q2 index {iq2}, W index {iw}, SF{k + 1}: " \
                       f"fortran={ref32[k]:.8e} python={got[k]:.8e}"
        flag = "ok " if rel[k] <= 1e-6 else "BAD"
        print(f"  {flag} ({iq2:3d},{iw:3d})  max|rel| = {rel[k]:.3e}  "
              f"(SF{k + 1}: {ref32[k]: .6e})")

    print()
    print(f"checked {n_checked} grid points x 62 amplitudes")
    print(f"worst relative difference: {worst:.3e}")
    print(f"  at {worst_at}")

    if save:
        # Save the reference amplitudes *together with* the indices that were
        # actually populated, so tests/test_table.py compares exactly the
        # points the Fortran was asked about rather than guessing.
        ref = np.zeros((62, nq2, nw), dtype=np.float32)
        checked = np.zeros((0, 2), dtype=np.int32)
        for iq2, iw in points:
            if 1 <= iq2 <= nq2 and 1 <= iw <= nw:
                ref[:, iq2 - 1, iw - 1] = run_fortran(tbl, iq2, iw).astype(np.float32)
                checked = np.append(checked, [[iq2 - 1, iw - 1]], axis=0)
        out = HERE / "fortran_table.npz"
        np.savez(out, amps=ref, points=checked)
        # Drop the .npy an earlier revision of this script used to write.
        (HERE / "fortran_table.npy").unlink(missing_ok=True)
        print(f"saved reference to {out} "
              f"({len(checked)} of {nq2 * nw} grid points populated)")

    if worst <= 1e-6:
        print("PASS: python table reader agrees with the Fortran reference")
        return 0
    print("FAIL: disagreement exceeds float32 round-off")
    return 1


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
