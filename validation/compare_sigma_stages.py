#!/usr/bin/env python3
"""Stage-by-stage comparison of ``sigma()`` between the Fortran and the port.

``validation/dump_sigma.f90`` sets up aao_rad's ``/ALPHA/`` and ``/radcal/``
COMMON blocks for one trial point, calls the *original* ``sigma()``, and dumps
every intermediate (qq, mf2, ffac, gfac, f, g, sig_r, sigf, epeps, the seven
response functions).  This script feeds it the radiative-branch rows of the
instrumented n-tuple and lines each Fortran intermediate up against the port's
``motsa_sigma`` recomputed in numpy, so a disagreement is localised to a single
factor instead of showing up as a ratio on the final cross section.

Only rows with ``ntp(17) = ek >= delta`` are usable.  The original jumps over
the ``sigr_raw``/``soft_ek`` assignments for the soft branch, so on soft rows
those two validation columns hold values left over from the previous *hard*
trial.  The soft branch is therefore *not* covered here; it is compared
through the full event weight in ``compare_sigma_points.py``.

Usage
-----
    validation/build_dump_sigma.sh
    cd /tmp/frun && python3 "$OLDPWD/validation/compare_sigma_stages.py" \\
        --ntuple /tmp/frun/aao_rad.ntuple --driver "$OLDPWD/build/dump_sigma"
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

import jax.numpy as jnp  # noqa: E402
import numpy as np  # noqa: E402

from aao_rad.config import GeneratorConfig  # noqa: E402
from aao_rad.constants import ALPHA, M_E, M_N, M_PIP, PI  # noqa: E402
from aao_rad.generate import build_grid  # noqa: E402
from aao_rad.motsa import _epsilon, _fg_amplitudes, motsa_sigma  # noqa: E402
from aao_rad.xsection import response_functions  # noqa: E402

WG = M_N + M_PIP + 5e-4  # aao_rad.f90: WG
wg2 = WG * WG

#: Output columns of dump_sigma, in order (32 total).
COLS = (
    "es ep q2 ek cstk phik csthcm phicm ehel sigma sigma_re "
    "qq qsq mf2 epw ffac gfac f g sig_r sigf_pre epeps sigf nu "
    "sig0 sigu sigt sigl sigi sigip asym_p eps_d"
).split()
#: dsigma's positional output order differs from its variable names
#: (maid_lee.f90:78-83): sigu carries sigma_t and sigt carries sigma_tt.
RESPONSE = {"sig0": "sigma0", "sigu": "sigma_t", "sigt": "sigma_tt",
            "sigl": "sigma_l", "sigi": "sigma_lt", "sigip": "sigma_ltp"}


def _is_numeric(line: str) -> bool:
    try:
        [float(t) for t in line.split()]
    except ValueError:
        return False
    return True


def run_driver(driver: Path, rows: np.ndarray, cwd: Path) -> np.ndarray:
    """Feed ``rows`` (es, ep_pre, q2, ek, cstk, phik_deg, csthcm, phicm, ehel)
    to the Fortran driver and return its output columns."""
    payload = "\n".join(" ".join(f"{v:.17e}" for v in r) for r in rows)
    out = subprocess.run(
        [str(driver)], input=payload, capture_output=True, text=True,
        cwd=cwd, check=True,
    )
    # maid_lee writes a two-line banner ("theory_opt,channel_opt,...") before
    # the data, and the data lines are 33 fixed-width columns; keep only lines
    # of exactly that width that also parse as numbers.
    data = np.array([[float(t) for t in line.split()]
                     for line in out.stdout.splitlines()
                     if len(line.split()) == len(COLS) and _is_numeric(line)])
    if data.shape[0] != rows.shape[0]:
        raise SystemExit(
            f"driver returned {data.shape[0]} rows for {rows.shape[0]} inputs"
        )
    return {c: data[:, i] for i, c in enumerate(COLS)}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ntuple", type=Path, default=Path("/tmp/frun/aao_rad.ntuple"))
    ap.add_argument("--driver", type=Path, default=ROOT / "build" / "dump_sigma")
    ap.add_argument("--cwd", type=Path, default=Path("/tmp/frun"),
                    help="must contain spp_tbl/ (maid_lee reads the table "
                         "relative to the working directory)")
    ap.add_argument("--delta", type=float, default=0.005)
    ap.add_argument("--max-rows", type=int, default=400)
    ap.add_argument("--channel", type=int, default=3)
    ap.add_argument("--scheme", default="linear")
    args = ap.parse_args()

    if not args.driver.exists():
        raise SystemExit(f"{args.driver} missing -- run validation/build_dump_sigma.sh")
    if not (args.cwd / "spp_tbl").exists() and "CLAS_PARMS" not in __import__("os").environ:
        raise SystemExit(f"{args.cwd} has no spp_tbl/ and CLAS_PARMS is unset")

    d = np.loadtxt(args.ntuple, ndmin=2)
    hard = d[:, 16] >= args.delta
    idx = np.where(hard)[0][: args.max_rows]
    print(f"{args.ntuple}: {len(d)} rows, {hard.sum()} radiative "
          f"(ek >= delta), using {len(idx)}")

    # 0-based columns: ntp(1..32) occupy 0..31, then the validation columns
    # weight, sigr_max, ntries, mcall, sigr_raw, ek_sigma, cstk, csthcm, phicm,
    # phik, ehel, es, ep, qsq, epw at 32..46.
    es = d[idx, 0]
    theta = np.deg2rad(d[idx, 2])
    q2 = d[idx, 29]
    ek = d[idx, 16]
    sigr_raw = d[idx, 36]
    cstk = d[idx, 38]
    csthcm = d[idx, 39]
    phicm = d[idx, 40]
    phik = d[idx, 41]
    ehel = d[idx, 42].astype(int)

    # sigma() saw the pre-exit-bremsstrahlung electron energy.  The recorded
    # theta is unaffected by that loss, and es is the post-entrance value, so
    # sin^2(theta_e/2) = Q^2/(4 E_s E_p) inverts it exactly.
    ep_pre = q2 / (4.0 * es * np.sin(theta / 2.0) ** 2)

    rows = np.column_stack([es, ep_pre, q2, ek, cstk, phik, csthcm, phicm,
                            ehel.astype(float)])
    fr = run_driver(args.driver, rows, args.cwd)

    print("\n-- sanity: driver vs the n-tuple's own sigma() record")
    r = fr["sigma"] / np.where(sigr_raw > 0, sigr_raw, 1.0)
    ok = sigr_raw > 0
    print(f"   n={ok.sum()}  ratio median {np.median(r[ok]):.6f}  "
          f"max|dev| {np.abs(r[ok] - 1).max():.2e}")
    r2 = fr["sigma_re"] / fr["sigma"]
    print(f"   replica inside the driver reproduces sigma(): max|dev| "
          f"{np.abs(r2 - 1).max():.2e}")

    # --- the port, at the same points ------------------------------------
    cfg = GeneratorConfig(channel=args.channel, w_max="clamp")
    grid = build_grid(cfg.channel, scheme=args.scheme,
                      parms_dir=str(ROOT / "parms"))

    cst0 = np.cos(theta)
    ps = np.sqrt(es**2 - M_E**2)
    pp = np.sqrt(ep_pre**2 - M_E**2)
    u0 = es - ep_pre + M_N
    pu = np.sqrt(ps**2 + pp**2 - 2 * ps * pp * cst0)
    uu = u0**2 - pu**2
    csths = (ps - pp * cst0) / pu
    csthp = (ps * cst0 - pp) / pu
    snths = np.sqrt(1 - csths**2)
    snthp = np.sqrt(1 - csthp**2)
    sntk = np.sin(np.arccos(cstk))
    sdotk = es * ek - ps * ek * cstk * csths - ps * ek * sntk * snths * np.cos(phik)
    pdotk = ep_pre * ek - pp * ek * cstk * csthp - pp * ek * sntk * snthp * np.cos(phik)

    sig, asym = motsa_sigma(
        grid,
        es=jnp.asarray(es), ep=jnp.asarray(ep_pre), th0=jnp.asarray(theta),
        cst0=jnp.asarray(cst0), ps=jnp.asarray(ps), pp=jnp.asarray(pp),
        sdotk=jnp.asarray(sdotk), pdotk=jnp.asarray(pdotk),
        u0=jnp.asarray(u0), pu=jnp.asarray(pu), uu=jnp.asarray(uu),
        cstk=jnp.asarray(cstk), ek=jnp.asarray(ek),
        csthcm=jnp.asarray(csthcm), phi=jnp.deg2rad(jnp.asarray(phicm)),
        e_hel=jnp.asarray(ehel),
        wg=M_N + M_PIP + 5e-4, m_pi=M_PIP, scheme=args.scheme,
    )
    sig = np.asarray(sig)

    print("\n-- motsa_sigma vs sigma()")
    live = sig > 1e-30
    rr = sig[live] / fr["sigma"][live]
    print(f"   n={live.sum()}  ratio median {np.median(rr):.6f}  "
          f"p1 {np.percentile(rr, 1):.6f}  p99 {np.percentile(rr, 99):.6f}  "
          f"max|dev| {np.abs(rr - 1).max():.3e}")

    # --- where does the difference come from? ----------------------------
    def show(name, port, ref, *, positive: bool = True):
        # Everything here is positive by construction except ``qq``, and the
        # original substitutes the 0.1e-30 sentinel for a failed ffac/gfac
        # test -- so require ref > 0 by default and let the sentinel drop out.
        ok = np.isfinite(port) & np.isfinite(ref)
        good = ok & ((ref > 0) & (port > 0) if positive else (np.abs(ref) > 0))
        if not good.any():
            print(f"   {name:<10s} no usable points")
            return
        q = port[good] / ref[good]
        print(f"   {name:<10s} n={good.sum():5d}  ratio median "
              f"{np.median(q):12.6f}  max|dev| {np.abs(q - 1).max():.3e}")

    print("\n-- geometry (Fortran internals, port recomputed in numpy)")
    mf2 = uu - 2 * ek * (u0 - pu * cstk)
    show("mf2", mf2, fr["mf2"])
    show("epw", np.sqrt(np.maximum(mf2, 0)), fr["epw"])
    # aao_rad.f90:1244 defines qq as the *negative* invariant mass squared;
    # qsq = -qq at line 1268.  Compare against qq itself, not qsq.
    show("qq", 2 * M_E**2 - 2 * es * ep_pre + 2 * ps * pp * cst0
         - 2 * ek * (es - ep_pre) + 2 * ek * pu * cstk, fr["qq"], positive=False)

    sdotk_s = np.where(sdotk == 0, 1.0, sdotk)
    pdotk_s = np.where(pdotk == 0, 1.0, pdotk)
    qq = 2 * M_E**2 - 2 * es * ep_pre + 2 * ps * pp * cst0\
        - 2 * ek * (es - ep_pre) + 2 * ek * pu * cstk
    qsq = -qq
    sp = es * ep_pre - ps * pp * cst0
    ffac = (-(M_E / pdotk_s) ** 2 * (2 * es * (ep_pre + ek) + qq / 2)
            - (M_E / sdotk_s) ** 2 * (2 * ep_pre * (es - ek) + qq / 2)
            - 2.0
            + 2.0 / sdotk_s / pdotk_s * (M_E**2 * (sp - ek**2)
                                         + sp * (2 * es * ep_pre - sp + ek * (es - ep_pre)))
            + (2 * (es * ep_pre + es * ek + ep_pre**2) + qq / 2 - sp - M_E**2) / pdotk_s
            - (2 * (es * ep_pre - ep_pre * ek + es**2) + qq / 2 - sp - M_E**2) / sdotk_s)
    gfac = (M_E**2 * (2 * M_E**2 + qq) * (1 / pdotk_s**2 + 1 / sdotk_s**2) + 4.0
            + 4 * sp * (sp - 2 * M_E**2) / pdotk_s / sdotk_s
            + (2 * sp + 2 * M_E**2 - qq) * (1 / pdotk_s - 1 / sdotk_s))
    show("ffac", ffac, fr["ffac"])
    show("gfac", gfac, fr["gfac"])

    sig_r = ((ALPHA**3 / (2 * PI * qq) ** 2) / M_N) * (ep_pre / es) * ek
    show("sig_r", sig_r, fr["sig_r"])

    # --- the amplitude block sigma() builds from dsigma's response ---------
    print("\n-- amplitudes (sigma() lines 1319-1360)")
    cfg = GeneratorConfig(channel=args.channel, w_max="clamp")
    grid = build_grid(cfg.channel, scheme=args.scheme, parms_dir=str(ROOT / "parms"))
    w_sq = jnp.asarray(np.where((mf2 > wg2) & (qq < 0.0), mf2, 0.0))
    qsq = -jnp.asarray(qq)
    epw = jnp.sqrt(jnp.maximum(w_sq, 0.0))
    cst0j = jnp.asarray(cst0)
    eps_r = _epsilon(jnp.asarray(es), jnp.asarray(ep_pre), cst0j,
                    (w_sq - M_N**2 + qsq) / (2.0 * M_N), qsq)
    resp = response_functions(grid, qsq, epw, jnp.asarray(csthcm),
                              jnp.deg2rad(jnp.asarray(phicm)), eps_r,
                              jnp.asarray(ehel), M_PIP, scheme=args.scheme)
    f_p, g_p, nu_p = _fg_amplitudes(resp, w_sq, qsq, M_PIP)
    show("nu", nu_p, fr["nu"])
    show("kfac", (w_sq - M_N**2) / 2.0 / M_N,
         (np.asarray(fr["mf2"]) - M_N**2) / 2.0 / M_N)
    show("epeps", _epsilon(jnp.asarray(es), jnp.asarray(ep_pre), cst0j,
                           jnp.asarray(es - ep_pre), qsq), fr["epeps"])
    show("eps_d", eps_r, fr["eps_d"])
    show("f", f_p, fr["f"])
    show("g", g_p, fr["g"])
    for name, attr in RESPONSE.items():
        show(attr, getattr(resp, attr), fr[name])
    show("asym_p", resp.asym_p, fr["asym_p"], positive=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
