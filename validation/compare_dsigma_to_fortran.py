#!/usr/bin/env python3
"""Compare the Python response functions against the original Fortran.

``dump_dsigma.f90`` calls the original ``dsigma()`` on a list of kinematic
points; this script generates those points, runs the Fortran (if the binary is
present), and evaluates the port's :func:`aao_rad.xsection.response_functions`
at the same points.  The two are then compared component by component.

This isolates the MAID response from every sampling and integration detail,
so a mismatch here is unambiguously a porting error in the response itself.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402

from aao_rad.amplitudes import multipole_amplitudes  # noqa: E402,F401
from aao_rad.generate import build_grid  # noqa: E402
from aao_rad.interpolate import InterpolationGrid  # noqa: E402
from aao_rad.xsection import Response, response_functions  # noqa: E402

FIELDS = ("sigma0", "sigma_t", "sigma_l", "sigma_tt", "sigma_lt", "sigma_ltp", "asym_p")

#: dsigma()'s positional output order is (sig0, sigu, sigt, sigl, sigi, sigip,
#: asym_p).  maid_lee.f90 fills those variables from sigma_0, sigma_t, sigma_tt,
#: sigma_l, sigma_lt, sigma_ltp -- note that the *names* ``sigt``/``sigl`` hold
#: sigma_tt/sigma_l, i.e. the name order and the value order differ.
FORT_COLUMNS = ("sigma_0", "sigma_t", "sigma_tt", "sigma_l", "sigma_lt", "sigma_ltp", "asym_p")


def sample_points(n: int, seed: int, channel: int) -> np.ndarray:
    """Random *physical* (theta, q2, W, cos(theta_cm), phi_cm, ehel) points.

    The point has to satisfy the on-shell conditions the generator guarantees:
    ``W > m_p + m_pi`` and ``nu_cm = (W^2 - m_p^2 - Q^2) / (2W) > 0`` (otherwise
    the hadronic CM kinematics are imaginary), plus ``Q^2 < 5`` and
    ``W < 2`` so the lookup stays inside the table.
    """
    rng = np.random.default_rng(seed)
    m_p, m_pi = 0.9382720813, (0.1349766 if channel == 1 else 0.13957039)

    theta = np.deg2rad(rng.uniform(5.0, 60.0, n))
    # Q^2 and W are drawn independently and then repaired, exactly the region the
    # generator actually visits: nu_cm > 0 is what makes the amplitude real.
    q2 = rng.uniform(0.05, 4.5, n)
    w_max = np.sqrt(q2 + m_p**2 + 1.4 * m_p)  # nu_cm >~ 0.7 GeV
    w_hi = np.minimum(w_max, 1.98)
    w = 1.10 + rng.uniform(0.0, 1.0, n) * (w_hi - 1.10)
    # Keep only points above the pion threshold; resample deterministically.
    good = (w > m_p + m_pi) & (w**2 - m_p**2 > q2) & (w < 2.0)
    for _ in range(64):
        if good.all():
            break
        theta = np.where(good, theta, np.deg2rad(rng.uniform(5.0, 60.0)))
        q2 = np.where(good, q2, rng.uniform(0.05, 4.5))
        w_hi = np.where(good, w_hi, np.minimum(np.sqrt(q2 + m_p**2 + 1.4 * m_p), 1.98))
        w = np.where(good, w, 1.10 + rng.uniform(0.0, 1.0) * (w_hi - 1.10))
        good = (w > m_p + m_pi) & (w**2 - m_p**2 > q2) & (w < 2.0)
    if not good.all():
        raise SystemExit("could not generate physical points")
    cscm = rng.uniform(-1.0, 1.0, n)
    phicm = rng.uniform(0.0, 360.0, n)
    ehel = rng.choice([-1, 0, 1], n)
    return np.column_stack([theta, q2, w, cscm, phicm, ehel])


def run_fortran(binary: Path, pts: np.ndarray, channel: int, cwd: Path) -> np.ndarray:
    lines = []
    for theta, q2, w, cscm, phicm, ehel in pts:
        lines.append(f"{theta:.10f} {q2:.10f} {w:.10f} {cscm:.10f} {phicm:.10f} "
                     f"{int(ehel)} 7 {channel} 0")
    proc = subprocess.run(
        [str(binary)],
        input="\n".join(lines) + "\n",
        capture_output=True,
        text=True,
        cwd=str(cwd),
    )
    if proc.returncode != 0:
        raise SystemExit(f"dump_dsigma failed ({proc.returncode}):\n{proc.stderr[-2000:]}")
    rows = [r.split() for r in proc.stdout.splitlines() if len(r.split()) == 36]
    if len(rows) != len(pts):
        raise SystemExit(
            f"dump_dsigma returned {len(rows)} rows for {len(pts)} points; "
            f"stderr:\n{proc.stderr[-2000:]}"
        )
    return np.array([[float(x) for x in r] for r in rows])


def port_response(
    grid: InterpolationGrid, pts: np.ndarray, channel: int, *, stages: bool = False
):
    """Response functions (and optionally the intermediate amplitudes) from the port."""
    theta, q2, w, cscm, phicm, ehel = (pts[:, i] for i in range(6))
    m_pi = 0.1349766 if channel == 1 else 0.13957039
    # dsigma() computes epsilon internally from Mp = 0.93827; mirror that.
    nu = 0.5 * (w**2 + q2 - 0.93827**2) / 0.93827
    eps = 1.0 / (1.0 + 2.0 * (1.0 + nu * nu / q2) * np.tan(0.5 * theta) ** 2)

    r: Response = response_functions(
        grid, q2, w, cscm, np.deg2rad(phicm), eps, ehel.astype(np.int32), m_pi,
        scheme="linear",
    )
    # dsigma()'s positional output order: sig0, sigu, sigt, sigl, sigi, sigip, asym_p
    out = np.column_stack([
        np.asarray(r.sigma0), np.asarray(r.sigma_t), np.asarray(r.sigma_tt),
        np.asarray(r.sigma_l), np.asarray(r.sigma_lt), np.asarray(r.sigma_ltp),
        np.asarray(r.asym_p),
    ])
    if not stages:
        return out

    from aao_rad.amplitudes import (
        cgln_amplitudes,
        helicity_amplitudes,
        legendre_polynomials,
    )
    from aao_rad.xsection import _center_of_mass

    w_c = np.clip(w, 1.1, 2.0)
    q2_c = np.minimum(q2, 5.0)
    p_pi_cm, qv_cm, nu_cm, _fkt = _center_of_mass(w_c, q2_c, m_pi)
    amps = grid(q2_c, w_c, scheme="linear")
    sp, sm, ep, em, mp, mm = multipole_amplitudes(amps, nu_cm, qv_cm)
    pol = legendre_polynomials(cscm)
    ff = cgln_amplitudes(pol, sp, sm, ep, em, mp, mm)
    hh = helicity_amplitudes(*ff, cscm)
    return out, {"ff": [np.asarray(a) for a in ff], "hh": [np.asarray(a) for a in hh]}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n", type=int, default=4000)
    ap.add_argument("--seed", type=int, default=20260930)
    ap.add_argument("--channel", type=int, default=3)
    ap.add_argument("--binary", type=Path, default=ROOT / "build" / "dump_dsigma")
    ap.add_argument("--scheme", default="linear", choices=["linear", "spline"])
    ap.add_argument(
        "--save", type=Path, default=ROOT / "validation" / "dsigma_comparison.npz",
        help="where to cache the Fortran reference",
    )
    args = ap.parse_args()

    pts = sample_points(args.n, args.seed, args.channel)

    if args.save.exists():
        ref = np.load(args.save)
        if ref["points"].shape == pts.shape and np.allclose(ref["points"], pts):
            ref_cols = ref["fortran"]
            print(f"loaded cached Fortran reference from {args.save}")
        else:
            ref_cols = None
    else:
        ref_cols = None

    if ref_cols is None:
        if not args.binary.exists():
            raise SystemExit(
                f"{args.binary} not found; build it with\n"
                "  cd build && gfortran -O2 -g ../validation/dump_dsigma.f90 "
                "libaao_libs.a -o dump_dsigma"
            )
        ref_cols = run_fortran(args.binary, pts, args.channel, ROOT / "parms")
        np.savez_compressed(args.save, points=pts, fortran=ref_cols)
        print(f"wrote Fortran reference to {args.save}")

    grid = build_grid(args.channel, scheme=args.scheme, parms_dir=str(ROOT / "parms"))
    got, stages = port_response(grid, pts, args.channel, stages=True)

    print()
    print(f"n = {pts.shape[0]}   scheme = {args.scheme}   channel = {args.channel}")
    print()

    # The Fortran writes Re/Im as separate real columns; rebuild complex pairs.
    def as_complex(start: int, count: int) -> list[np.ndarray]:
        return [
            ref_cols[:, start + 2 * k] + 1j * ref_cols[:, start + 2 * k + 1]
            for k in range(count)
        ]

    worst_all = 0.0

    def report(title, names, refs, gots) -> None:
        nonlocal worst_all
        print(title)
        print(f"{'quantity':>14} {'max |rel diff|':>15} {'median |rel diff|':>18} "
              f"{'max |abs diff|':>15}")
        w = 0.0
        for name, a, b in zip(names, refs, gots, strict=False):
            a = np.asarray(a)
            b = np.asarray(b, dtype=a.dtype)
            rel = np.abs(a - b) / np.maximum(np.abs(a), 1e-30)
            w = max(w, float(rel.max()))
            print(f"{name:>14} {rel.max():15.3e} {np.median(rel):18.3e} "
                  f"{np.abs(a - b).max():15.3e}")
        print(f"{'worst':>14} {w:15.3e}")
        print()
        worst_all = max(worst_all, w)

    report(
        "CGLN amplitudes ff1..ff6",
        ("ff1", "ff2", "ff3", "ff4", "ff5", "ff6"),
        as_complex(12, 6), stages["ff"],
    )
    report(
        "helicity amplitudes hh1..hh6",
        ("hh1", "hh2", "hh3", "hh4", "hh5", "hh6"),
        as_complex(24, 6), stages["hh"],
    )
    report(
        "response functions (dsigma output order)",
        FORT_COLUMNS,
        [ref_cols[:, 5 + i] for i in range(7)], [got[:, i] for i in range(7)],
    )

    print(f"worst relative difference overall: {worst_all:.3e}")
    return 0 if worst_all < 1e-4 else 1


if __name__ == "__main__":
    raise SystemExit(main())
