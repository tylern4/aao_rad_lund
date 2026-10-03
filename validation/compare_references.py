#!/usr/bin/env python3
"""Compare two Fortran n-tuples run from the same card, column by column.

Whenever the instrumented ``aao_rad.f90`` is changed -- a bug fix, a new dump
column, a switch like ``AAO_TRIALS`` -- the question is whether the change did
what it claimed and nothing else.  Two runs of the *same* binary from the same
run card answer that, because they are independent samples of the same physics
with different ``myran`` seeds (``aao_rad.f90:373`` seeds from ``unixtime``, so
row-for-row equality is impossible and distributional agreement is the only
available criterion).

So each column is compared by two-sample Kolmogorov-Smirnov, against the
reference noise floor ``1.36 * sqrt(1/N_a + 1/N_b)``.  A column that a change
was *not* supposed to touch must sit at the floor.  A column it was supposed to
fix must move clear of it.

This is a much sharper instrument than the port comparison for validating a
reference edit: both sides are the same code modulo the edit, so the difference
is attributable to the edit alone, with no port involved.

Usage
-----
    python validation/compare_references.py \\
        --before /tmp/frun/aao_rad.ntuple \\
        --after  /tmp/frun_fixed/aao_rad.ntuple \\
        --expect-changed asym_p
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "validation"))

from compare_distributions import NTP_COLUMNS, NTP_EXTRA, load_fortran_ntuple  # noqa: E402

#: Run-level bookkeeping rather than observables, so a KS distance on them is
#: meaningless: ``ntries`` is a monotone counter, ``sigr_max`` is the scanned
#: importance-sampling ceiling (a single number per run, sampled), and ``mcall``
#: is the accept/reject outcome.  Pass ``--include-run-state`` to look anyway.
RUN_STATE = ("ntries", "sigr_max", "mcall")


def two_sample_ks(a: np.ndarray, b: np.ndarray) -> float:
    """Two-sample KS statistic, with no binning and no range restriction.

    ``scipy`` is not required here: this is the asymptotic two-sample form,
    which is what the noise-floor comparison needs anyway.
    """
    a = np.sort(np.asarray(a, np.float64))
    b = np.sort(np.asarray(b, np.float64))
    if a.size == 0 or b.size == 0:
        return 0.0
    grid = np.concatenate([a, b])
    grid.sort(kind="stable")
    # Only compare at the sample points; ties handled by taking the right edge.
    fa = np.searchsorted(a, grid, side="right") / a.size
    fb = np.searchsorted(b, grid, side="right") / b.size
    return float(np.abs(fa - fb).max())


def noise_floor(n_a: int, n_b: int) -> float:
    """1.36 * sqrt(1/N_a + 1/N_b): the 5% two-sample KS quantile."""
    return 1.36 * np.sqrt(1.0 / n_a + 1.0 / n_b)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--before", type=Path, required=True,
                    help="n-tuple from the unmodified reference")
    ap.add_argument("--after", type=Path, required=True,
                    help="n-tuple from the modified reference")
    ap.add_argument("--expect-changed", nargs="*", default=[],
                    help="columns the edit was supposed to change")
    ap.add_argument("--include-run-state", action="store_true",
                    help=f"also compare {', '.join(RUN_STATE)} (not observables)")
    ap.add_argument("--min-moved", type=float, default=1.3,
                    help="a declared-changed column must reach this multiple of "
                         "the floor to count as having moved")
    ap.add_argument("--max-quiet", type=float, default=2.0,
                    help="an undeclared column above this multiple of the floor "
                         "counts as an unexpected move")
    args = ap.parse_args()

    a = load_fortran_ntuple(args.before)
    b = load_fortran_ntuple(args.after)
    n_a = len(a["es"])
    n_b = len(b["es"])
    floor = noise_floor(n_a, n_b)

    print(f"before : {args.before}  ({n_a} events)")
    print(f"after  : {args.after}  ({n_b} events)")
    print(f"noise floor (5% KS, two-sample) = {floor:.5f}\n")

    names = list(NTP_COLUMNS) + list(NTP_EXTRA)
    expected = set(args.expect_changed)
    unknown = expected - set(names)
    if unknown:
        raise SystemExit(f"unknown column(s) in --expect-changed: {sorted(unknown)}")
    if not args.include_run_state:
        skipped = [n for n in names if n in RUN_STATE]
        names = [n for n in names if n not in RUN_STATE]
        print(f"skipping run-state columns: {', '.join(skipped)}\n")

    print(f"{'column':<10s} {'KS':>9s} {'/floor':>7s}  {'mean(before)':>14s} "
          f"{'mean(after)':>14s}  verdict")
    print("-" * 78)

    problems = []
    for name in names:
        ks = two_sample_ks(a[name], b[name])
        ratio = ks / floor if floor > 0 else 0.0
        if name in expected:
            moved = ratio >= args.min_moved
            verdict = "CHANGED (expected)" if moved else "expected to change but did NOT"
            if not moved:
                problems.append(name)
        else:
            quiet = ratio < args.max_quiet
            verdict = "same" if quiet else "MOVED (unexpected)"
            if not quiet:
                problems.append(name)
        print(f"{name:<10s} {ks:9.5f} {ratio:7.2f}  {a[name].mean():14.6g} "
              f"{b[name].mean():14.6g}  {verdict}")

    print()
    if problems:
        print("FAILED: the following columns did not behave as declared:")
        for name in problems:
            print(f"  - {name}")
        return 1
    print(f"OK: every column outside {sorted(expected) or '[]'} stayed at the "
          "noise floor, and every declared column moved.")
    print("Note that a declared column need only clear the floor to count as "
          "moved: a stale value and a correct one can be drawn from "
          "overlapping pools, so the two distributions differ subtly rather "
          "than grossly. Confirm a fix by comparing against the port, not by "
          "how far the reference moved.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
