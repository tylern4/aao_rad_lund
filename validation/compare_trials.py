#!/usr/bin/env python3
"""Compare the port's *trial* measure against the Fortran's.

The cross section is a product of two independent things::

    sigma = P(reach the weight stage) * mean(weight | reached) * phase_space

Comparing only the final ``sigma`` cannot tell a sampling-density mismatch from
a cross-section mismatch, and comparing only the accepted-event n-tuple cannot
tell them from each other at all, because the accepted events are importance
sampled *by the weight*.  ``aao_rad.f90`` therefore also dumps a 1-in-N
subsample of the trials that reach the weight stage (``unit 14,
``aao_rad.trials``), which is free of any importance sampling and directly
comparable with the port's own trials.

Two numbers carry the whole comparison:

``mean(weight | reached)``
    The trial-averaged integrand, evaluated on identically sampled points.  A
    disagreement here is a cross-section bug.
``P(reach the weight stage)``
    The fraction of *all* trials -- including those rejected before the weight
    is ever computed -- that survive every cut.  A disagreement here is a
    sampling bug.

The script reports both, their product, and the conditional distributions of
the sampled variables so a mismatch can be localised.

Requires ``aao_rad.trials`` next to the run, written by the instrumented
``aao_rad`` (see ``src/aao_rad.f90`` at the ``write(14, ...)`` statement).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

import jax  # noqa: E402
import jax.numpy as jnp  # noqa: E402
import numpy as np  # noqa: E402

from aao_rad.config import GeneratorConfig  # noqa: E402
from aao_rad.generate import (  # noqa: E402
    build_grid,
    build_kinematics,
    draw_kinematics,
    draw_photon,
    integrand,
)

#: ``(name, lo, hi)`` for every variable whose conditional distribution is
#: compared.  The ranges are the full kinematic window, not the observed one.
VARIABLES = (
    ("es", 3.80, 4.25),
    ("q2", 0.19, 1.91),
    ("ep", 1.59, 2.91),
    ("ek", 0.00, 0.60),
    ("cstk", -1.0, 1.0),
    ("phik", -3.15, 3.15),
)

#: The ten columns this script needs, in dump order.  The dump has since grown
#: two more (mm2, wreal, for ``compare_mm2.py``), so only the leading ten are
#: read and older ten-column dumps still work.
F77_COLUMNS = ("es", "q2", "ep", "ek", "cstk", "phik", "intreg", "mcfac", "mpfac", "sigr")
F77_MIN_COLUMNS = len(F77_COLUMNS)


def ks_distance(a: np.ndarray, b: np.ndarray, lo: float, hi: float, bins: int = 60):
    """Two-sample KS distance and the CDF-difference integral on a grid.

    The integral is the L1 distance between the two ECDFs.  For a correctly
    matched pair it is ~1/n; an order of magnitude larger means the shapes
    genuinely differ rather than just the noise.

    Also returns two z-scores.  ``z_diff`` says how significant the difference
    of the means is.  ``ref_z`` says how well the *reference* mean is resolved
    from zero, which is what decides whether a ratio of the two means means
    anything at all: a signed ``phi_k`` has expectation zero, so its ratio is
    the quotient of two noise values and can read 0.4 while the two samples
    agree perfectly.
    """
    edges = np.linspace(lo, hi, bins + 1)
    a_hist, _ = np.histogram(a, bins=edges)
    b_hist, _ = np.histogram(b, bins=edges)
    a_cdf = np.concatenate([[0.0], np.cumsum(a_hist) / a.size])
    b_cdf = np.concatenate([[0.0], np.cumsum(b_hist) / b.size])
    ks = float(np.max(np.abs(a_cdf - b_cdf)))
    l1 = float(np.abs(a_cdf - b_cdf).sum() / bins)
    a_mean = float(a.mean())
    b_mean = float(b.mean())
    a_se = float(a.std() / np.sqrt(a.size))
    b_se = float(b.std() / np.sqrt(b.size))
    z_diff = (b_mean - a_mean) / float(np.hypot(a_se, b_se)) if (a_se or b_se) else 0.0
    ref_z = a_mean / a_se if a_se > 0.0 else 0.0
    return ks, l1, a_mean, b_mean, z_diff, ref_z


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--trials", type=Path, default=Path("/tmp/frun/aao_rad.trials"))
    ap.add_argument("--fortran-ntries", type=float, default=None,
                    help="trial count from the Fortran run; if given, the "
                         "Fortran fraction reaching the weight stage is derived "
                         "from it instead of from --fortran-sigma")
    ap.add_argument("--fortran-sigma", type=float, default=None,
                    help="sig_sum printed by the Fortran run, in micro-barns")
    ap.add_argument("--tdump", type=int, default=97,
                    help="stride of the Fortran trial dump (tdump in "
                         "aao_rad.f90); needed to scale dumped rows back up "
                         "to a trial count when using --fortran-ntries")
    ap.add_argument("--n", type=int, default=1 << 24, help="port trials")
    ap.add_argument("--seed", type=int, default=20260902)
    ap.add_argument("--beam-energy", type=float, default=4.244)
    ap.add_argument("--q2-min", type=float, default=0.2)
    ap.add_argument("--q2-max", type=float, default=1.9)
    ap.add_argument("--ep-min", type=float, default=1.6)
    ap.add_argument("--ep-max", type=float, default=2.9)
    ap.add_argument("--delta", type=float, default=0.005)
    ap.add_argument("--k-exp", type=float, default=5.0)
    ap.add_argument("--mm-cut", type=float, default=0.2)
    ap.add_argument("--w-max", default="clamp")
    ap.add_argument("--ek-sampling", default="fortran", choices=["fortran", "truncated"])
    args = ap.parse_args()

    if not args.trials.exists():
        print(f"{args.trials} not found; run the instrumented aao_rad first", file=sys.stderr)
        return 1

    f77 = np.loadtxt(args.trials, ndmin=2)
    if f77.shape[1] < F77_MIN_COLUMNS:
        print(f"{args.trials}: expected at least {F77_MIN_COLUMNS} columns, "
              f"got {f77.shape[1]}", file=sys.stderr)
        return 1
    f77 = f77[:, :F77_MIN_COLUMNS]
    print(f"fortran trial dump: {f77.shape[0]:,} rows from {args.trials}")

    w_max = "clamp" if args.w_max == "clamp" else float(args.w_max)
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
        w_max=w_max,
        ek_sampling=args.ek_sampling,
    )
    grid = build_grid(cfg.channel, scheme=cfg.interp_scheme, parms_dir=str(ROOT / "parms"))
    kin = build_kinematics(cfg)

    @jax.jit
    def chunk(kv, ph):
        return integrand(grid, kin, kv, ph, cfg.interp_scheme)

    k1, k2 = jax.random.split(jax.random.PRNGKey(args.seed))
    kv = draw_kinematics(k1, kin, args.n)
    ph = draw_photon(k2, kin, kv, args.n)
    integ = chunk(kv, ph)
    print(f"port trials       : {args.n:,}")
    print(f"phase space       : {(4.0 * np.pi) ** 2 * 2.0 * np.pi * float(kin.uq2_range) * float(kin.ep_range):.6g}")
    print()

    # ---------------------------------------------------------------- weights
    f77_w = f77[:, F77_COLUMNS.index("sigr")]
    port_ok = np.asarray(integ.ok)
    port_w = np.asarray(jnp.where(port_ok, integ.weight, 0.0))
    port_mean_ok = float(port_w[port_ok].mean())

    print("mean(weight | reached the weight stage)")
    print(f"  fortran         : {f77_w.mean():.8e}")
    print(f"  port            : {port_mean_ok:.8e}")
    print(f"  ratio           : {port_mean_ok / f77_w.mean():.6f}")
    print()

    port_ok_frac = float(port_ok.mean())
    print("P(reach the weight stage)")
    if args.fortran_sigma is not None or args.fortran_ntries is not None:
        ps = (4.0 * np.pi) ** 2 * 2.0 * np.pi * float(kin.uq2_range) * float(kin.ep_range)
        if args.fortran_sigma is not None:
            f77_sigma_over_ps = args.fortran_sigma / ps
            # mean_dumped = sig_tot / n_dumped and sig_tot/ntries = sigma/PS, so
            # n_dumped/ntries -- the fraction of trials that reach the weight
            # stage -- is (sigma/PS)/mean_dumped.
            f77_ok_frac = f77_sigma_over_ps / f77_w.mean()
            f77_sigma = args.fortran_sigma
        else:
            # The dump is written when ``mod(ntries, tdump) .eq. 0``, i.e. it is
            # a 1-in-tdump stride over *all* trials, and only the ones that then
            # survive the missing-mass cut produce a row.  So the rows have to be
            # multiplied back up by tdump before they can be divided by the total
            # trial count -- dividing them directly understates the acceptance by
            # exactly that factor.
            f77_ok_frac = f77.shape[0] * args.tdump / args.fortran_ntries
            f77_sigma = f77_ok_frac * f77_w.mean() * ps
        print(f"  fortran         : {f77_ok_frac:.6f}"
              + (f"   [{f77.shape[0]:,} / {args.fortran_ntries:,.0f}]"
                 if args.fortran_ntries else ""))
        print(f"  port            : {port_ok_frac:.6f}")
        print(f"  ratio           : {port_ok_frac / f77_ok_frac:.6f}")
        print()
        print("sigma = P(reached) * mean(weight | reached) * PS")
        print(f"  fortran         : {f77_sigma:.8g} micro-barn")
        port_sigma = port_ok_frac * port_mean_ok * ps
        print(f"  port            : {port_sigma:.8g} micro-barn")
        print(f"  ratio           : {port_sigma / f77_sigma:.6f}")
        print(f"  from P(reached) : {(port_ok_frac / f77_ok_frac) - 1.0:+.4%}")
        print(f"  from mean(w)    : {(port_mean_ok / f77_w.mean()) - 1.0:+.4%}")
        print()

    # --------------------------------------------------------------- regions
    # aao_rad.f90:769 overwrites intreg with 6 for every soft trial, so the
    # region mix of the hard rows is the only comparable one.
    f77_hard = f77[:, 6] != 6
    csran = np.asarray(ph["csran"])
    reg_edges = np.asarray(kin.reg)
    port_reg = np.searchsorted(reg_edges, csran, side="right") + 1
    port_hard = np.asarray(ph["ek"]) >= args.delta

    print("region mix of the hard rows")
    print(f"  {'region':>7} {'fortran':>10} {'port':>10}")
    for r in range(1, 6):
        a = float((f77[f77_hard, 6] == r).mean())
        b = float((port_reg[port_hard & port_ok] == r).mean())
        print(f"  {r:>7} {a:>10.5f} {b:>10.5f}   ratio {b / a:.5f}")
    print()

    # ---------------------------------------------------- variable shapes
    print("conditional distributions (both rows already conditioned on reaching")
    print("the weight stage), KS distance and L1 ECDF gap")
    print(f"  {'variable':>8} {'KS':>9} {'L1':>9} {'fortran mean':>14} "
          f"{'port mean':>12} {'ratio':>8} {'dz':>7}")
    for i, (name, lo, hi) in enumerate(VARIABLES):
        f_vals = f77[:, i]
        p_vals = np.asarray(
            {"es": kv["es"], "q2": kv["q2"], "ep": kv["ep"], "ek": ph["ek"],
             "cstk": ph["cstk"], "phik": ph["phik"]}[name][port_ok]
        )
        ks, l1, a_mean, b_mean, z, ref_z = ks_distance(f_vals, p_vals, lo, hi)
        # Only quote a ratio when the reference mean is itself resolved from
        # zero; otherwise print "-", because the ratio of two means that are
        # both consistent with zero is not a comparison.
        ratio = f"{b_mean / a_mean:8.5f}" if abs(ref_z) > 5.0 and a_mean else "       -"
        print(f"  {name:>8} {ks:>9.5f} {l1:>9.5f} {a_mean:>14.6g} {b_mean:>12.6g} "
              f"{ratio} {z:>+7.2f}")
    print()

    print("noise floor for the KS distances above: 1.7/sqrt(N) with")
    print(f"  N_fortran = {f77.shape[0]:,}  ->  {1.7 / np.sqrt(f77.shape[0]):.5f}")
    print(f"  N_port    = {port_ok.sum():,}  ->  {1.7 / np.sqrt(max(port_ok.sum(), 1)):.5f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
