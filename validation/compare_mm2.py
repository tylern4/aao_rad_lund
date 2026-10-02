#!/usr/bin/env python3
"""Compare the missing-mass cut between the port and the Fortran.

The cut ``|mm2 - mm_exp| <= mm_cut`` acts on ``ek`` only through ``mm2``, and
``mm2`` appears in no accepted-event column, so a shifted or narrowed window is
invisible in the output n-tuple.  It shows up instead as a shifted ``ek``
distribution among the trials that survive, which is why ``ek`` was the one
sampled variable whose conditional distribution disagreed while ``es``, ``q2``,
``ep``, ``cstk`` and ``phik`` all sat at the noise floor.

Both sides are compared *after* their own cut, which is what makes this a clean
test rather than a circular one: if the two ``mm2`` reconstructions agree, then
both windows are the same window, and the survivors are drawn from the same
restricted prior, so the two ``mm2`` distributions must match.  A residual shift
therefore cannot be explained away as a selection effect -- it has to be a
different ``mm2`` or a different ``mm_exp``.

Three things are reported:

* the survivor ``mm2`` and ``w_real`` distributions, binned;
* the edges of each window as actually realised (the outermost populated bins),
  since a ``mm_exp`` offset would show up there first;
* the acceptance profile ``A(ek)``, formed by dividing the survivor density by
  the *proposal* density.  The proposal is ``exp(-k_exp*ek)`` truncated at
  ``ek_max``, so the port can supply it exactly, and the ratio of the two
  profiles is the ``ek`` discrepancy without the proposal in the way.

Requires the instrumented Fortran run (``aao_rad.trials`` with the 12-column
layout -- see ``validation/README.md``).
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
from aao_rad.generate import (  # noqa: E402
    draw_kinematics,
    draw_photon,
    integrand,
)
from aao_rad.kinematics import hadronic_final_state  # noqa: E402

F77_COLUMNS = (
    "es", "q2", "ep", "ek", "cstk", "phik", "intreg", "mcfac", "mpfac",
    "weight", "mm2", "wreal",
)


def window_edges(values: np.ndarray, edges: np.ndarray, frac: float = 0.002) -> None:
    """Report where the populated range of ``values`` actually sits.

    A difference in ``mm_exp`` moves the whole window, and shows up in the
    outermost populated bins long before it is visible in a KS statistic, so it
    is worth printing the realised edges directly.
    """
    counts, _ = np.histogram(values, bins=edges)
    nz = np.nonzero(counts)[0]
    if nz.size == 0:
        print("    (no populated bins)")
        return
    lo_i, hi_i = int(nz[0]), int(nz[-1])
    total = counts.sum()
    # Trim bins holding less than ``frac`` of the sample from each end, so a
    # single stray event does not define the edge.
    while lo_i < hi_i and counts[lo_i] < frac * total:
        lo_i += 1
    while hi_i > lo_i and counts[hi_i] < frac * total:
        hi_i -= 1
    print(f"    populated [{edges[lo_i]:.6f}, {edges[hi_i + 1]:.6f}]  "
          f"outer bins kept for display: {edges[nz[0]]:.6f} .. {edges[nz[-1] + 1]:.6f}")


def compare_dist(label: str, a: np.ndarray, b: np.ndarray, lo: float, hi: float,
                 n_bins: int, quantiles: bool = True) -> None:
    """Binned comparison of two survivor samples of the same variable."""
    edges = np.linspace(lo, hi, n_bins + 1)
    ha, _ = np.histogram(a, bins=edges)
    hb, _ = np.histogram(b, bins=edges)
    fa = ha / max(ha.sum(), 1)
    fb = hb / max(hb.sum(), 1)
    ratio = np.where(fb > 0, fa / np.where(fb > 0, fb, 1.0), np.nan)
    print()
    print(f"survivor {label}:  python {a.size:,}  fortran {b.size:,}")
    print(f"{'bin':>4} {'range':>22} {'n_fr':>8} {'n_py':>8} "
          f"{'frac_fr':>11} {'frac_py':>11} {'ratio':>8}")
    print("-" * 78)
    for i in range(n_bins):
        r = "      -" if not np.isfinite(ratio[i]) else f"{ratio[i]:8.4f}"
        print(f"{i:4d} {edges[i]:10.5f}..{edges[i + 1]:<10.5f} {hb[i]:8d} {ha[i]:8d} "
              f"{fb[i]:11.5e} {fa[i]:11.5e} {r}")
    good = ratio[np.isfinite(ratio)]
    if good.size:
        print(f"  ratio min {good.min():.4f}  max {good.max():.4f}  "
              f"spread {good.max() - good.min():.4f}")
    if quantiles and a.size and b.size:
        qs = [1, 5, 25, 50, 75, 95, 99]
        qa = np.percentile(a, qs)
        qb = np.percentile(b, qs)
        print(f"  {'quantile':>10} " + " ".join(f"{q:>9g}%" for q in qs))
        print(f"  {'fortran':>10} " + " ".join(f"{v:12.6f}" for v in qb))
        print(f"  {'python':>10} " + " ".join(f"{v:12.6f}" for v in qa))
        print(f"  {'py/fr':>10} " + " ".join(f"{v:11.6f} " for v in qa / qb))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--trials", type=Path,
                    default=Path("/tmp/frun_mm2/aao_rad.trials"))
    ap.add_argument("--n", type=int, default=1 << 23, help="port trials")
    ap.add_argument("--bins", type=int, default=20)
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
              f"(see validation/README.md)", file=sys.stderr)
        return 1

    f77 = np.loadtxt(args.trials, ndmin=2)
    if f77.shape[1] != len(F77_COLUMNS):
        print(f"{args.trials}: expected {len(F77_COLUMNS)} columns, got "
              f"{f77.shape[1]} -- rebuild aao_rad with the current instrumented "
              f"source", file=sys.stderr)
        return 1
    col = {name: f77[:, i] for i, name in enumerate(F77_COLUMNS)}
    print(f"fortran trial dump: {f77.shape[0]:,} rows from {args.trials}")

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
    ok = np.asarray(integ.ok)
    print(f"port trials       : {args.n:,}  ({ok.sum():,} reach the weight stage)")

    # The port's mm2/w_real on every trial, not just the survivors, so the
    # proposal can be recovered below.
    fs = hadronic_final_state(
        kin.e_beam, kv["es"], kv["ep"], kv["th0"], ph["ek"], ph["cstk"],
        ph["phik"], ph["csthcm"], ph["phicm_deg"], kin.m_pi, kin.is_pi0,
    )
    py_mm2 = np.asarray(fs.mm2)[ok]
    py_wreal = np.asarray(fs.w_real)[ok]
    py_ek = np.asarray(ph["ek"])[ok]
    print(f"port m_exp = {float(kin.m_exp):.6f}   mm_cut = {args.mm_cut}   "
          f"window [{float(kin.m_exp) - args.mm_cut:.6f}, "
          f"{float(kin.m_exp) + args.mm_cut:.6f}]")
    print(f"fortran mm_exp = {0.938 ** 2:.6f} (mp = .938, aao_rad.f90:277)   "
          f"window [{0.938 ** 2 - args.mm_cut:.6f}, "
          f"{0.938 ** 2 + args.mm_cut:.6f}]")

    # ------------------------------------------------------------- survivors
    compare_dist("mm2", py_mm2, col["mm2"], 0.6, 1.35, args.bins)
    compare_dist("w_real", py_wreal, col["wreal"], 1.0, 2.5, args.bins)

    print()
    print("realised window edges")
    print("  mm2   python:")
    window_edges(py_mm2, np.linspace(0.6, 1.35, args.bins + 1))
    print("  mm2   fortran:")
    window_edges(col["mm2"], np.linspace(0.6, 1.35, args.bins + 1))
    print("  wreal python:")
    window_edges(py_wreal, np.linspace(1.0, 2.5, args.bins + 1))
    print("  wreal fortran:")
    window_edges(col["wreal"], np.linspace(1.0, 2.5, args.bins + 1))

    # ------------------------------------------- acceptance profile in ek
    # The proposal is exp(-k_exp*ek) truncated at ek_max, so on the *drawn*
    # (pre-cut) population its density is known analytically.  Dividing the
    # survivor density by it recovers the acceptance profile, and the ratio of
    # the two profiles is the ek discrepancy with the proposal divided out.
    n_draw = args.n
    surv_f = np.histogram(col["ek"], bins=np.linspace(0, 0.6, 25))[0] / col["ek"].size
    surv_p = np.histogram(py_ek, bins=np.linspace(0, 0.6, 25))[0] / py_ek.size
    # Marginal proposal density, from the port's own draws before the cut.  Only
    # the *shape* is needed -- the constant normalises away in the ratio.
    prop = np.histogram(np.asarray(ph["ek"]), bins=np.linspace(0, 0.6, 25))[0] / n_draw
    with np.errstate(divide="ignore", invalid="ignore"):
        acc_f = surv_f / prop
        acc_p = surv_p / prop
        ratio = acc_p / acc_f
    print()
    print("acceptance profile A(ek)  (survivor density / proposal density)")
    edges = np.linspace(0, 0.6, 25)
    print(f"{'range':>18} {'prop_rel':>10} {'A_fr':>10} {'A_py':>10} {'A_py/A_fr':>10}")
    print("-" * 62)
    pmax = prop[prop > 0].max() if (prop > 0).any() else 1.0
    for i in range(len(edges) - 1):
        if prop[i] <= 0:
            continue
        r = "        -" if not np.isfinite(ratio[i]) else f"{ratio[i]:10.4f}"
        print(f"{edges[i]:8.3f}..{edges[i + 1]:<8.3f} {prop[i] / pmax:10.4f} "
              f"{acc_f[i]:10.5f} {acc_p[i]:10.5f} {r}")
    good = ratio[np.isfinite(ratio)]
    if good.size:
        print(f"  A_py/A_fr min {good.min():.4f}  max {good.max():.4f}  "
              f"spread {good.max() - good.min():.4f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
