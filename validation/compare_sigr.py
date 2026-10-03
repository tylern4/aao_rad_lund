#!/usr/bin/env python3
"""Attribute a cross-section mismatch to a region of phase space.

The trial weight factorises as

    weight = sigr_raw * mcfac * mpfac * jacob

and the Fortran trial dump carries ``mcfac``, ``mpfac`` and the sampled
variables that ``jacob`` is a closed form of.  So the raw cross section is
*recoverable from every dumped row*::

    sigr_raw = weight / (mcfac * mpfac * jacob)
    jacob    = exp(k_exp * ek) / k_exp / (2 * es * ep) * q2**2

That matters because it removes the last need to instrument the Fortran to
compare the physics.  It also double-checks the geometry reconstruction row by
row: ``sigr_raw`` is only a sensible positive cross section if the factors were
really separable in the order the Fortran applied them.

Note on naming: column 10 of the dump is the Fortran variable ``sigr`` *after*
it has been multiplied by ``jacob`` and by ``mcfac * mpfac``, so it holds the
full weight, not the cross section (``aao_rad.f90:930-943``).  This script calls
it ``weight`` and recovers the cross section from it.

A single global ratio of means cannot say *where* a disagreement lives, and the
two candidate explanations need opposite fixes:

* a ratio **flat** across the sampled volume is a normalisation error -- a wrong
  constant, not a wrong function -- and no amount of shape work will fix it;
* a ratio that **varies** is a shape error: a wrong term, a wrong index, or a
  wrong interpolation somewhere in the W/Q^2/omega response functions, and it is
  only fixable once located.

So the script bins ``sum(sigr_python) / sum(sigr_fortran)`` along each sampled
variable separately and reports, per bin, the ratio and that bin's share of the
total signed discrepancy.  The share column is what says where to look.

Requires the instrumented Fortran run that produces ``aao_rad.trials``
(see the Validation section of ``README.md``).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

import jax  # noqa: E402

from aao_rad import build_grid, build_kinematics  # noqa: E402
from aao_rad.config import GeneratorConfig  # noqa: E402
from aao_rad.generate import draw_kinematics, draw_photon, integrand  # noqa: E402

#: Column names in the Fortran trial dump (aao_rad.f90:943, unit 14).  Column 10
#: is the Fortran's ``sigr``, which by then is the full trial weight.
F77_COLUMNS = (
    "es", "q2", "ep", "ek", "cstk", "phik", "intreg", "mcfac", "mpfac", "weight",
)

#: (label, column, lo, hi) -- the conditioning variables to bin along.
BIN_ALONG = (
    ("es (GeV)", "es", 3.9, 4.25),
    ("q2 (GeV^2)", "q2", 0.2, 1.9),
    ("ep (GeV)", "ep", 1.59, 2.91),
    ("ek (GeV)", "ek", 0.0, 0.60),
    ("cstk", "cstk", -1.0, 1.0),
    ("abs(phik)", "_aphik", 0.0, 3.15),
)


def recover_sigr(col: dict, k_exp: float) -> tuple[np.ndarray, np.ndarray]:
    """Undo the geometry factors on every dumped row.

    Returns the recovered cross section and the weight it came from.  The
    recovery is exact in algebra and in float64, so it also acts as a check on
    the factorisation: ``sigr_raw`` should be strictly positive and of order
    unity in these units, not spread over many decades.
    """
    jacob = (
        np.exp(k_exp * col["ek"]) / k_exp
        / (2.0 * col["es"] * col["ep"]) * col["q2"] ** 2
    )
    geom = col["mcfac"] * col["mpfac"] * jacob
    with np.errstate(divide="ignore", invalid="ignore"):
        sigr_raw = col["weight"] / geom
    return sigr_raw, geom


def _binned(x: np.ndarray, v: np.ndarray, edges: np.ndarray, how: str) -> np.ndarray:
    """Aggregate ``v`` over bins of ``x``.

    Percentiles cannot go through ``bincount``, so they are computed by grouping
    once and slicing, which is fine at the few-million-row scale this runs at.
    """
    if how == "sum":
        idx = np.clip(np.digitize(x, edges) - 1, 0, len(edges) - 2)
        return np.bincount(idx, weights=v, minlength=len(edges) - 1)
    q = 90.0 if how == "p90" else 50.0
    out = np.full(len(edges) - 1, np.nan)
    idx = np.clip(np.digitize(x, edges) - 1, 0, len(edges) - 2)
    order = np.argsort(idx, kind="stable")
    bounds = np.searchsorted(idx[order], np.arange(len(edges)))
    vals = v[order]
    for i in range(len(edges) - 1):
        sel = vals[bounds[i]:bounds[i + 1]]
        if sel.size:
            out[i] = np.percentile(sel, q)
    return out


def report(label: str, x_f: np.ndarray, w_f: np.ndarray,
           x_p: np.ndarray, w_p: np.ndarray, lo: float, hi: float,
           n_bins: int) -> None:
    """Binned comparison along one conditioning variable.

    Every statistic is a *mean* over the bin, never a raw sum.  The two samples
    hold very different numbers of rows -- the Fortran dump is a 1-in-191
    stride sample while the port's is contiguous -- so a sum ratio would just
    report the row-count ratio and look like a result while carrying no
    information about the physics at all.

    Three statistics per bin, because ``sigr_raw`` spans ten decades and they
    answer different questions:

    ``mean`` ratio
        What the cross-section integral sees, weighted by how much of the
        sampled volume each bin holds.  Heavily tail-influenced, so it is the
        number that must match but the noisiest one to read.
    ``p90`` ratio
        The 90th percentile within the bin.  Robust to the tail, so a mismatch
        here is a genuine *shape* disagreement rather than a handful of extreme
        events landing differently in two independent samples.
    ``median`` ratio
        The same idea further into the bulk.  If the mean moves while the median
        does not, the disagreement lives in the tail of the integrand and the
        place to look is the large-``omega`` behaviour; if all three move
        together, it is a normalisation.
    """
    edges = np.linspace(lo, hi, n_bins + 1)
    sum_f = _binned(x_f, w_f, edges, "sum")
    sum_p = _binned(x_p, w_p, edges, "sum")
    cnt_f = _binned(x_f, np.ones_like(w_f), edges, "sum")
    cnt_p = _binned(x_p, np.ones_like(w_p), edges, "sum")
    with np.errstate(divide="ignore", invalid="ignore"):
        mean_f = sum_f / np.where(cnt_f > 0, cnt_f, 1.0)
        mean_p = sum_p / np.where(cnt_p > 0, cnt_p, 1.0)
    p90_f = _binned(x_f, w_f, edges, "p90")
    p90_p = _binned(x_p, w_p, edges, "p90")
    med_f = _binned(x_f, w_f, edges, "median")
    med_p = _binned(x_p, w_p, edges, "median")

    ratio = np.where(mean_f > 0.0, mean_p / np.where(mean_f > 0.0, mean_f, 1.0), np.nan)
    p90r = np.where(p90_f > 0.0, p90_p / np.where(p90_f > 0.0, p90_f, 1.0), np.nan)
    medr = np.where(med_f > 0.0, med_p / np.where(med_f > 0.0, med_f, 1.0), np.nan)

    # Decompose the global mean difference into per-bin contributions, weighted
    # by the share of the *reference* population each bin holds.  That is what
    # answers "which part of the sampled space is the missing weight in", and
    # unlike a plain ratio it is not allowed to hide a large deficit in a small
    # bin behind a large surplus in a big one.
    share_f = cnt_f / max(cnt_f.sum(), 1.0)
    contrib = (mean_p - mean_f) * share_f
    tot = contrib.sum()

    print()
    print(f"sigr_raw along {label}   (row counts: fortran {cnt_f.sum():,.0f}, "
          f"python {cnt_p.sum():,.0f})")
    print(f"{'bin':>4} {'range':>20} {'n_fr':>8} {'mean_fr':>11} {'mean_py':>11} "
          f"{'mean r':>8} {'p90 r':>8} {'med r':>8} {'share':>7}")
    print("-" * 92)
    for i in range(len(edges) - 1):
        def fmt(a: float) -> str:
            return "       -" if not np.isfinite(a) else f"{a:8.4f}"
        print(f"{i:4d} {edges[i]:9.4g}..{edges[i + 1]:<9.4g} {cnt_f[i]:8.0f} "
              f"{mean_f[i]:11.4e} {mean_p[i]:11.4e} {fmt(ratio[i])} "
              f"{fmt(p90r[i])} {fmt(medr[i])} "
              f"{(contrib[i] / tot if tot != 0 else 0.0):+7.3f}")
    for name, arr in (("mean", ratio), ("p90", p90r), ("median", medr)):
        good = arr[np.isfinite(arr)]
        if good.size:
            print(f"  {name:>6} ratio min {good.min():.4f}  max {good.max():.4f}  "
                  f"spread {good.max() - good.min():.4f}  "
                  f"geometric mean {np.exp(np.mean(np.log(good))):.6f}")
    print(f"  global mean ratio {w_p.mean() / w_f.mean():.6f}   "
          f"rows are unequal, so this is a mean ratio and not sum_p/sum_f")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--trials", type=Path, default=Path("/tmp/frun/aao_rad.trials"))
    ap.add_argument("--n", type=int, default=1 << 23, help="port trials")
    ap.add_argument("--bins", type=int, default=12)
    ap.add_argument("--seed", type=int, default=20260902)
    ap.add_argument("--beam-energy", type=float, default=4.244)
    ap.add_argument("--q2-min", type=float, default=0.2)
    ap.add_argument("--q2-max", type=float, default=1.9)
    ap.add_argument("--ep-min", type=float, default=1.6)
    ap.add_argument("--ep-max", type=float, default=2.9)
    ap.add_argument("--delta", type=float, default=0.005)
    ap.add_argument("--k-exp", type=float, default=5.0)
    ap.add_argument("--mm-cut", type=float, default=0.2)
    args = ap.parse_args()

    if not args.trials.exists():
        print(f"{args.trials} not found; run the instrumented aao_rad first "
              f"(see the Validation section of README.md)", file=sys.stderr)
        return 1

    f77 = np.loadtxt(args.trials, ndmin=2)
    if f77.shape[1] != len(F77_COLUMNS):
        print(f"{args.trials}: expected {len(F77_COLUMNS)} columns, "
              f"got {f77.shape[1]}", file=sys.stderr)
        return 1
    col = {name: f77[:, i] for i, name in enumerate(F77_COLUMNS)}
    col["_aphik"] = np.abs(col["phik"])
    print(f"fortran trial dump: {f77.shape[0]:,} rows from {args.trials}")

    # ---------------------------------------------------- recover the physics
    col["sigr"], geom = recover_sigr(col, args.k_exp)
    finite = np.isfinite(col["sigr"]) & (col["sigr"] > 0.0)
    print(f"recovered sigr_raw from the dump: {finite.sum():,}/{col['sigr'].size:,} "
          f"positive and finite")
    q = np.percentile(col["sigr"][finite], [0.1, 1, 50, 99, 99.9])
    print(f"  recovered sigr_raw quantiles 0.1/1/50/99/99.9%: "
          f"{q[0]:.4g} {q[1]:.4g} {q[2]:.4g} {q[3]:.4g} {q[4]:.4g}")
    if finite.sum() != col["sigr"].size:
        print(f"  !! {col['sigr'].size - finite.sum():,} rows did not recover -- "
              f"the weight does not factorise as assumed")

    cfg = GeneratorConfig(
        beam_energy=args.beam_energy,
        q2_min=args.q2_min,
        q2_max=args.q2_max,
        ep_min=args.ep_min,
        ep_max=args.ep_max,
        min_photon_energy=args.delta,
        k_exp=args.k_exp,
        missing_mass_cut=args.mm_cut,
        n_events=1,
        batch_size=max(args.n, 1024),
        seed=args.seed,
        w_max="clamp",
        ek_sampling="fortran",
    )
    grid = build_grid(cfg.channel, scheme=cfg.interp_scheme,
                      parms_dir=str(ROOT / "parms"))
    kin = build_kinematics(cfg)

    @jax.jit
    def chunk(kv, ph):
        return integrand(grid, kin, kv, ph, cfg.interp_scheme)

    k1, k2 = jax.random.split(jax.random.PRNGKey(args.seed))
    kv = draw_kinematics(k1, kin, args.n)
    ph = draw_photon(k2, kin, kv, args.n)
    integ = chunk(kv, ph)
    # Condition on reaching the weight stage: the dump is written after the
    # missing-mass cut, so that is the population the decomposition is on.
    ok = np.asarray(integ.ok)
    print(f"port trials       : {args.n:,}  ({ok.sum():,} reach the weight stage)")
    print()

    port = {
        "es": np.asarray(kv["es"]), "q2": np.asarray(kv["q2"]),
        "ep": np.asarray(kv["ep"]), "ek": np.asarray(ph["ek"]),
        "cstk": np.asarray(ph["cstk"]),
        "_aphik": np.abs(np.asarray(ph["phik"])),
    }
    for name in ("es", "q2", "ep", "ek", "cstk", "_aphik"):
        port[name] = port[name][ok]
    port["sigr"] = np.asarray(integ.sigr)[ok]

    soft_f = col["intreg"] == 6
    # The port's soft branch is the low-omega one, selected on ek < delta
    # exactly as the Fortran selects it (aao_rad.f90, soft section).
    soft_p = port["ek"] < args.delta

    print("sigr_raw: overall and by radiative branch")
    print(f"{'':<22} {'fortran':>13} {'port':>13} {'ratio':>10} "
          f"{'n_fr':>10} {'n_py':>10}")
    print("-" * 82)
    for label, mf, mp in (
        ("all", np.ones_like(soft_f, bool), np.ones_like(soft_p, bool)),
        ("hard (ek >= delta)", ~soft_f, ~soft_p),
        ("soft (ek < delta)", soft_f, soft_p),
    ):
        sf, sp = col["sigr"][mf], port["sigr"][mp]
        print(f"{label:<22} {sf.mean():13.6e} {sp.mean():13.6e} "
              f"{sp.sum() / sf.sum():10.6f} {mf.sum():10,d} {mp.sum():10,d}")
    print()

    for label, name, lo, hi in BIN_ALONG:
        report(label, col[name], col["sigr"], port[name], port["sigr"], lo, hi,
               args.bins)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
