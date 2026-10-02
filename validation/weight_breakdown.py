#!/usr/bin/env python3
"""Decompose the port's mean event weight by importance-sampling region.

Prints the contribution of the soft (``ek < delta``) and hard branches of the
integrand, and of each ``cos(theta_k)`` region, to the mean weight.  Because
the cross section is ``mean_weight * phase_space``, a factor that only shows up
in one region is a normalisation bug in that region.

The integrand has a heavy tail -- ``sig_r`` in ``motsa_sigma`` carries a
``1/Q^2`` factor and ``Q^2 -> 0`` lies inside the requested window -- so the
mean weight converges far more slowly than a binomial count would suggest.
The script therefore runs several *independent* batches and reports the
between-batch spread, which is the honest error estimate.  The naive Poisson
error ``sigma / sqrt(n_trials)`` is printed as well, to show how optimistic it
is.
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


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--beam-energy", type=float, default=4.244)
    ap.add_argument("--q2-min", type=float, default=0.2)
    ap.add_argument("--q2-max", type=float, default=1.9)
    ap.add_argument("--ep-min", type=float, default=1.6)
    ap.add_argument("--ep-max", type=float, default=2.9)
    ap.add_argument("--delta", type=float, default=0.005)
    ap.add_argument("--k-exp", type=float, default=5.0)
    ap.add_argument("--mm-cut", type=float, default=0.2)
    ap.add_argument("--w-max", default="clamp")
    ap.add_argument("--ek-sampling", default="fortran")
    ap.add_argument("--n", type=int, default=1 << 22, help="trials per batch")
    ap.add_argument("--repeats", type=int, default=8, help="independent batches")
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument(
        "--fortran-sigma", type=float, default=None,
        help="cross section printed by the Fortran run, in micro-barns",
    )
    args = ap.parse_args()

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

    ps_all = (4.0 * np.pi) ** 2 * 2.0 * np.pi * float(kin.uq2_range) * float(kin.ep_range)
    print(f"trials per batch   : {args.n}   batches: {args.repeats}")
    print(f"phase space        : {ps_all:.6g}")
    print()

    sigmas, soft_share, hard_share, ok_frac = [], [], [], []
    for r in range(args.repeats):
        k1, k2 = jax.random.split(jax.random.PRNGKey(args.seed + 7919 * r))
        kv = draw_kinematics(k1, kin, args.n)
        ph = draw_photon(k2, kin, kv, args.n)
        integ = chunk(kv, ph)
        o = integ.ok
        w = np.asarray(jnp.where(o, integ.weight, 0.0))
        ek = np.asarray(ph["ek"])
        soft = ek < args.delta
        sigmas.append(w.mean() * ps_all)
        soft_share.append(w[soft].sum() / w.sum())
        hard_share.append(w[~soft].sum() / w.sum())
        ok_frac.append(float(np.asarray(o).mean()))
        print(f"  batch {r}: sigma = {sigmas[-1]:.6g} micro-barn")

    s = np.array(sigmas)
    total = w.size * args.repeats
    print()
    print(f"accepted fraction  : {np.mean(ok_frac) * 100:.4f}%")
    print(f"soft  (ek<delta)   : {np.mean(soft_share) * 100:.2f}% of the cross section")
    print(f"hard  (ek>=delta)  : {np.mean(hard_share) * 100:.2f}% of the cross section")
    print()
    print(f"sigma              : {s.mean():.6g} micro-barn")
    if args.repeats > 1:
        print(f"between-batch sd   : {s.std(ddof=1):.3g}")
        print(f"standard error     : {s.std(ddof=1) / np.sqrt(args.repeats):.3g}")
    print(f"(naive Poisson)    : {s.mean() * np.sqrt(1.0 / total):.3g}"
          f"   [{total:.3g} trials]")

    if args.fortran_sigma is not None:
        ratio = s.mean() / args.fortran_sigma
        print()
        print(f"fortran            : {args.fortran_sigma:.6g} micro-barn")
        print(f"ratio python/fortran: {ratio:.4f}")
        se = args.fortran_sigma * s.std(ddof=1) / s.mean() / np.sqrt(args.repeats) \
            if args.repeats > 1 else 0.0
        print(f"deviation          : {(ratio - 1.0) * s.mean() / max(se, 1e-30):+.1f} sigma"
              f"   (ratio 1.00 = agreement)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
