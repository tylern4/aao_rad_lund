#!/usr/bin/env python3
"""Compute the port's integrated cross section for a given run card.

The original prints its own answer in the log::

    Integrated cross section (MC, numerical) =  1.53775308E-02  mu-barns

which is the *mean-weight* estimator ``sum(weight) * phase_space / n_trials``
(``aao_rad.f90:1142``).  That estimator needs no importance-sampling ceiling
and is directly comparable with :attr:`GenerationStats.sigma_accepted`, which is
the same quantity written a different way.

This script runs the port under an explicit run card so the number can be
compared with a Fortran run made from the same card.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402

from aao_rad import EventGenerator, GeneratorConfig, build_grid  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--beam-energy", type=float, default=4.244)
    ap.add_argument("--q2-min", type=float, default=0.2)
    ap.add_argument("--q2-max", type=float, default=1.9)
    ap.add_argument("--ep-min", type=float, default=1.6)
    ap.add_argument("--ep-max", type=float, default=2.9)
    ap.add_argument("--mm-cut", type=float, default=0.2)
    ap.add_argument("--delta", type=float, default=0.005)
    ap.add_argument("--k-exp", type=float, default=5.0)
    ap.add_argument("--target-length", type=float, default=5.0)
    ap.add_argument("--n-events", type=int, default=20000)
    ap.add_argument("--batch-size", type=int, default=1 << 17)
    ap.add_argument("--seed", type=int, default=None)
    ap.add_argument("--w-max", default="clamp", help="'clamp' or a number or 'table'")
    ap.add_argument("--ek-sampling", default="fortran")
    ap.add_argument(
        "--fortran-sigma",
        type=float,
        default=None,
        help="the value printed by the Fortran run, in micro-barns",
    )
    args = ap.parse_args()

    w_max: float | str | None
    if args.w_max == "table":
        w_max = None
    elif args.w_max == "clamp":
        w_max = "clamp"
    else:
        w_max = float(args.w_max)

    cfg = GeneratorConfig(
        beam_energy=args.beam_energy,
        q2_min=args.q2_min,
        q2_max=args.q2_max,
        ep_min=args.ep_min,
        ep_max=args.ep_max,
        missing_mass_cut=args.mm_cut,
        min_photon_energy=args.delta,
        k_exp=args.k_exp,
        target_length_cm=args.target_length,
        n_events=args.n_events,
        batch_size=args.batch_size,
        seed=args.seed,
        w_max=w_max,
        ek_sampling=args.ek_sampling,
    )

    gen = EventGenerator(
        build_grid(cfg.channel, scheme=cfg.interp_scheme, parms_dir=str(ROOT / "parms"))
    )
    _, stats = gen.generate(cfg)

    print()
    print(stats)
    print(f"phase space        : {stats.phase_space:.8g}")
    print(f"mean weight        : {stats.mean_weight:.8g}")
    # Poisson error on the mean-weight estimator: sigma * sqrt(1/n_trials).
    err = stats.sigma_mc * np.sqrt(1.0 / max(stats.n_trials, 1))
    print(f"sigma (mean-weight): {stats.sigma_mc:.6g} +- {err:.3g} micro-barn")

    if args.fortran_sigma is not None:
        ref = args.fortran_sigma
        ratio = stats.sigma_mc / ref
        dev = (ratio - 1.0) * np.sqrt(stats.n_trials)
        print()
        print(f"fortran            : {ref:.6g} micro-barn")
        print(f"ratio python/fortran: {ratio:.5f}")
        print(f"deviation (sigma)    : {dev:+.2f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
