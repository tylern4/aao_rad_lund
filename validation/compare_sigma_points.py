#!/usr/bin/env python3
"""Pointwise check of the integrand against the instrumented Fortran.

``aao_rad.f90`` was instrumented to write, next to every accepted event, both
the raw ``sigma()`` return value and the complete trial kinematics.  This
script re-evaluates the port at exactly those points.

Because the comparison is at fixed kinematics, no sampling, no importance
weight and no cross-section normalisation enters.  Four factors of the event
weight are checked separately, so a disagreement can be localised:

``rad``
    ``sigma()`` for the radiative branch (``ek > delta``), i.e. the Mo & Tsai
    integrand against the MAID07 amplitudes.  Already cross-checked at the
    amplitude level by ``compare_dsigma_to_fortran.py``.
``soft``
    ``signr * (1 + delta_r) * exp(delta_inf)`` for the soft branch
    (``ek < delta``), the ``0 < ek < delta`` integral Mo & Tsai folded into a
    single weight.
``jacob``
    ``exp(k_exp * ek) / k_exp / (2 E_s E_p) * Q^2``.
``mcfac * mpfac``
    The five-region importance-sampling normalisation (aao_rad.f90:693-729).
    Reconstructed from ``cos(theta_k)`` and checked against the discrete set of
    values the Fortran can produce.

Two recording subtleties this script has to undo:

* ``ntp(17) = ek`` is the drawn photon energy, unrescaled, and is what both
  ``sigma()`` and the soft branch saw.  The separate ``soft_ek`` validation
  column is written only by the hard branch (the soft branch jumps over
  aao_rad.f90:896-898), so it is stale on soft rows and must not be used.
* ``ntp(2) = ep`` is *post*-exit-bremsstrahlung; ``sigma()`` was called with the
  pre-exit value, which is recovered from the recorded angle and ``Q^2`` via
  ``sin^2(theta_e/2) = Q^2 / (4 E_s E_p)``.
* ``ntp(19) = phik`` is ``phi_k`` after the shift at aao_rad.f90:1091, i.e.
  ``phi_k + pi``.  Column 42 holds the unshifted value that ``sigma()`` saw.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

import jax.numpy as jnp  # noqa: E402
import numpy as np  # noqa: E402

from aao_rad.config import GeneratorConfig  # noqa: E402
from aao_rad.constants import M_E, M_N  # noqa: E402
from aao_rad.generate import build_grid, build_kinematics  # noqa: E402
from aao_rad.motsa import motsa_sigma, non_radiative_sigma  # noqa: E402

#: Column index (0-based) into ``aao_rad.ntuple``.
NTP_COLUMNS = (
    "es ep theta epw wreal ppx ppy ppz eprot ppix ppiy ppiz epi csthcm phicm "
    "mm2 ek cstk phik qvecx qvecz q0 cst0 ekx eky ekz vx vy vz qsq ehel asym_p"
).split()
NTP_EXTRA = ("weight sigr_max ntries mcall sigr_raw ek_sigma cstk csthcm phicm "
             "phik ehel es ep qsq epw").split()
assert len(NTP_COLUMNS) == 32 and len(NTP_EXTRA) == 15


def read_ntuple(path: Path, max_rows: int | None = None) -> dict[str, np.ndarray]:
    raw = np.loadtxt(path, ndmin=2)
    need = len(NTP_COLUMNS) + len(NTP_EXTRA)
    if raw.shape[1] < need:
        raise SystemExit(
            f"{path}: {raw.shape[1]} columns, need {need}.  Rebuild "
            "src/aao_rad.f90 with the validation instrumentation."
        )
    if max_rows:
        raw = raw[:max_rows]
    out = {n: raw[:, i].astype(np.float64) for i, n in enumerate(NTP_COLUMNS)}
    out.update(
        {n: raw[:, len(NTP_COLUMNS) + j].astype(np.float64)
         for j, n in enumerate(NTP_EXTRA)}
    )
    return out


def _stats(label: str, port: np.ndarray, ref: np.ndarray, mask: np.ndarray) -> None:
    ok = mask & (ref > 0) & np.isfinite(port)
    if not ok.any():
        print(f"  {label:<12s} no usable points")
        return
    r = port[ok] / ref[ok]
    print(f"  {label:<12s} n={ok.sum():6d}  ratio median {np.median(r):9.5f}  "
          f"mean {r.mean():9.5f}  p1 {np.percentile(r, 1):9.5f}  "
          f"p99 {np.percentile(r, 99):9.5f}  "
          f"max|dev| {np.abs(r - 1).max():.2e}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ntuple", type=Path, default=Path("/tmp/frun/aao_rad.ntuple"))
    ap.add_argument("--channel", type=int, default=3)
    ap.add_argument("--delta", type=float, default=0.005)
    ap.add_argument("--mm-cut", type=float, default=0.2)
    ap.add_argument("--k-exp", type=float, default=5.0)
    ap.add_argument("--csrng", type=float, default=0.04,
                    help="aao_rad.f90:205 hardcodes csrng = .04")
    ap.add_argument("--regions", type=float, nargs=4,
                    default=(0.20, 0.12, 0.20, 0.20),
                    help="reg1..reg4 from the run card (cumulative)")
    ap.add_argument("--max-rows", type=int, default=0)
    ap.add_argument("--scheme", default="linear")
    ap.add_argument("--dump-bad", type=int, default=6)
    args = ap.parse_args()

    d = read_ntuple(args.ntuple, args.max_rows or None)
    n = len(d["es"])
    print(f"{args.ntuple}: {n} accepted events\n")

    cfg = GeneratorConfig(channel=args.channel, min_photon_energy=args.delta,
                          missing_mass_cut=args.mm_cut, w_max="clamp")
    grid = build_grid(cfg.channel, scheme=args.scheme, parms_dir=str(ROOT / "parms"))
    kin = build_kinematics(cfg)

    es = d["es"]
    q2 = d["qsq"]
    th0 = np.deg2rad(d["theta"])
    cst0 = np.cos(th0)
    s2 = np.sin(th0 / 2.0) ** 2
    # sigma() saw the pre-exit-bremsstrahlung scattered-electron energy.
    ep = q2 / (4.0 * es * s2)
    ps = np.sqrt(np.maximum(es**2 - M_E**2, 0.0))
    pp = np.sqrt(np.maximum(ep**2 - M_E**2, 0.0))

    u0 = es - ep + M_N
    pu = np.sqrt(np.maximum(ps**2 + pp**2 - 2.0 * ps * pp * cst0, 0.0))
    uu = u0**2 - pu**2
    csths = (ps - pp * cst0) / np.where(pu > 0, pu, 1.0)
    csthp = (ps * cst0 - pp) / np.where(pu > 0, pu, 1.0)
    snths = np.sqrt(np.maximum(1.0 - csths**2, 0.0))
    snthp = np.sqrt(np.maximum(1.0 - csthp**2, 0.0))

    cstk = np.clip(d["cstk"], -1.0, 1.0)
    sntk = np.sin(np.arccos(cstk))
    ek = d["ek"]
    # column 42 is already in radians (the driver records it that way); only
    # phicm is in degrees.
    phik = d["phik"]

    sdotk = es * ek - ps * ek * cstk * csths - ps * ek * sntk * snths * np.cos(phik)
    pdotk = ep * ek - pp * ek * cstk * csthp - pp * ek * sntk * snthp * np.cos(phik)

    common = dict(
        es=jnp.asarray(es), ep=jnp.asarray(ep), th0=jnp.asarray(th0),
        cst0=jnp.asarray(cst0), ps=jnp.asarray(ps), pp=jnp.asarray(pp),
        csthcm=jnp.asarray(d["csthcm"]), phi=jnp.deg2rad(jnp.asarray(d["phicm"])),
        e_hel=jnp.asarray(d["ehel"].astype(np.int32)),
        scheme=args.scheme,
    )
    rad, _ = motsa_sigma(
        grid, sdotk=jnp.asarray(sdotk), pdotk=jnp.asarray(pdotk),
        u0=jnp.asarray(u0), pu=jnp.asarray(pu), uu=jnp.asarray(uu),
        cstk=jnp.asarray(cstk), ek=jnp.asarray(ek),
        wg=kin.wg, m_pi=kin.m_pi, **common,
    )
    w_sq = M_N**2 + 2.0 * M_N * (es - ep) - q2        # aao_rad.f90:616
    soft, _ = non_radiative_sigma(
        grid, q2=jnp.asarray(q2), w_sq=jnp.asarray(w_sq),
        delta=args.delta, m_pi=kin.m_pi, **common,
    )
    # non_radiative_sigma returns the *averaged* value (already divided by
    # delta * 4 * pi); sigr1 is what aao_rad.f90:875 computes.
    sigr1 = np.asarray(soft) * args.delta * 4.0 * np.pi

    rad = np.asarray(rad)
    ref_raw = d["sigr_raw"]
    ref_wt = d["weight"]

    # The original calls sigma() with the *drawn* ek even on the soft branch;
    # only its return value is discarded.  So the branch is decided by the
    # drawn energy, which the n-tuple kept as ntp(17).
    hard = d["ek"] >= args.delta

    print("radiative branch, sigma() vs aao_rad.f90:sigr_raw")
    _stats("rad", rad, ref_raw, hard)

# jacob and the five-region normalisation
    jac = np.exp(args.k_exp * ek) / args.k_exp / (2.0 * es * ep) * q2**2

    # aao_rad.f90:657-731.  cstk1 >= cstk2 after the swap at line 668, the band
    # is clamped twice (lines 670-671), csrngb = csrng/40 capped at csrnge/5, and
    # delphi is computed then *unconditionally* reset to pi/9 at line 700 -- the
    # geometric value above it is dead code.  Regions 1-4 all draw phik
    # uniformly in [-pi/18, pi/18]; region 5 draws phik over the full circle but
    # redraws while it is inside that window *and* inside a band (line 721), so
    # a wide-phik point inside a band is a region-5 draw, not region 3 or 4.
    cstk1, cstk2 = csths, csthp
    swap = cstk1 < cstk2
    cstk1 = np.where(swap, csthp, csths)
    cstk2 = np.where(swap, csths, csthp)
    csrnge = np.minimum.reduce(
        [np.full(n, args.csrng), 1.0 - cstk1, 0.5 * (cstk1 - cstk2)]
    )
    csrngb = np.minimum(args.csrng / 40.0, csrnge / 5.0)
    delphi = np.pi / 9.0
    reg1, reg2, reg3, reg4 = np.cumsum(np.asarray(args.regions))

    inband = np.abs(phik) < delphi / 2.0
    d1 = np.abs(cstk - cstk1)
    d2 = np.abs(cstk - cstk2)
    r1 = inband & (d1 < csrngb)
    r2 = inband & ~r1 & (d2 < csrngb)
    r3 = inband & ~r1 & ~r2 & (d1 < csrnge)
    r4 = inband & ~r1 & ~r2 & ~r3 & (d2 < csrnge)
    r5 = ~(r1 | r2 | r3 | r4)

    mcfac = np.empty(n)
    mpfac = np.empty(n)
    mcfac[r1] = csrngb[r1] / reg1
    mcfac[r2] = csrngb[r2] / (reg2 - reg1)
    mcfac[r3] = (csrnge - csrngb)[r3] / (reg3 - reg2)
    mcfac[r4] = (csrnge - csrngb)[r4] / (reg4 - reg3)
    mcfac[r5] = (1.0 - csrnge[r5] * delphi / np.pi) / (1.0 - reg4)
    mpfac[r5] = 1.0
    mpfac[~r5] = delphi / 2.0 / np.pi

    sigma_used = np.where(hard, ref_raw, sigr1 / args.delta / 4.0 / np.pi)
    pred = mcfac * mpfac * sigma_used * jac
    print("\nevent weight, mcfac*mpfac*sigma*jacob vs aao_rad.f90:sigr")
    _stats("weight", pred, ref_wt, np.ones(n, bool))

    # The same weight with the port's own cross section, evaluated on exactly
    # the same accepted events.  Both sigma estimators are a mean of the weight
    # over their *own* trials, so a residual gap in the totals after this ratio
    # is a difference in the trial density, not in the physics -- whereas this
    # ratio itself is purely the cross-section bias.
    port_sigma = np.where(hard, rad, np.asarray(soft))
    pred_port = mcfac * mpfac * port_sigma * jac
    ok = np.isfinite(pred_port) & (pred_port > 0) & (ref_wt > 0)
    def _sum(mask):
        if not (ok & mask).any():
            return float("nan")
        return pred_port[ok & mask].sum() / ref_wt[ok & mask].sum()
    print(f"\nweight with the port's cross section, on the same {ok.sum()} events")
    print(f"  sum(pred)/sum(f77) = {_sum(np.ones(n, bool)):.6f}"
          f"   [soft rows {_sum(~hard):.6f}, hard rows {_sum(hard):.6f}]")

    # The reconstruction of the region is discrete; report which values appear.
    print("\nreconstructed mcfac*mpfac values (rounded):")
    vals, counts = np.unique(np.round((mcfac * mpfac)[hard], 3), return_counts=True)
    for v, cnt in zip(vals, counts, strict=False):
        print(f"  {v:12.7g}  x{cnt}")
    exp_vals = {
        "region 1": csrngb[0] / reg1 * delphi / 2.0 / np.pi,
        "region 2": csrngb[0] / (reg2 - reg1) * delphi / 2.0 / np.pi,
        "region 3": (csrnge[0] - csrngb[0]) / (reg3 - reg2) * delphi / 2.0 / np.pi,
        "region 4": (csrnge[0] - csrngb[0]) / (reg4 - reg3) * delphi / 2.0 / np.pi,
        "region 5": (1.0 - csrnge[0] * delphi / np.pi) / (1.0 - reg4),
    }
    for k, v in exp_vals.items():
        print(f"  expect {k:9s} {v:12.7g}")
    if args.dump_bad:
        both = (ref_raw > 0) & (rad > 0)
        r = np.where(both, rad / np.where(ref_raw > 0, ref_raw, 1.0), 1.0)
        bad = hard & both & ~np.isclose(r, 1.0, rtol=1e-3)
        print(f"\nradiative points off by >0.1%: {bad.sum()} / {hard.sum()}")
        for i in np.where(bad)[0][: args.dump_bad]:
            print(f"  es={es[i]:8.4f} ep={ep[i]:8.4f} q2={q2[i]:7.4f} "
                  f"ek={ek[i]:9.6f} mf2={uu[i] - 2 * ek[i] * (u0[i] - pu[i] * cstk[i]):8.5f} "
                  f"cstk={cstk[i]:8.4f} fortran={ref_raw[i]:13.6e} "
                  f"port={rad[i]:13.6e} ratio={r[i]:8.4f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
