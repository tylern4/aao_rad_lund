#!/usr/bin/env python3
"""What does filling one node with 128 runs actually cost?

``bench_fortran.sbatch`` timed four runs on a 128-core node.  Four runs never
contend with each other, so what it measured is one idle core, and
``bench_report.py`` then divides by 128 to price a run.  The division is the
assumption; this measures it, on the same node type, with all 192 runs of the
grid in flight at one run per physical core.

Two things have to be right about the comparison, and both are reasons not to
use the projections ``bench_report`` prints.

**Not the projections.**  ``proj`` is ``wall * 20000 / events``, and events per
wall is trials-per-event divided by trials-per-second -- trials per event being
set by ``sigr_max``, which each run re-estimates from its own random points and
which varies about threefold between two runs of the *same* configuration.  A
projection ratio therefore carries the ceiling's noise as if it were packing;
cfg_000 projects to 710 s in the probe and 1,087 s in the reference while the
machine ran at 209,691 and 230,148 trials/s.  Trials per second is the quantity
that measures the machine.

**Not two medians over different sets.**  The probe times all 96
configurations and the reference the four ``bench_fortran`` chose; trials per
event span 400x across the grid, so medians taken over different sets are
different for reasons that have nothing to do with packing.  Rates are
compared configuration by configuration where both arms timed one, and overall
only through the fit below, which pools work rather than configurations.

The fit is ``wall = startup + trials / rate``: trials are known exactly from
``out.txt``, so the slope is throughput and the intercept is everything that
did not scale with work.  It matters because the probe's runs are ~17 s and the
reference's ~53 s, and a fixed ~1 s of process startup is therefore a larger
fraction of the short runs -- comparing raw medians would charge that to
packing.

Usage::

    packing_penalty.py <probe summary.json> <reference summary.json>
"""

from __future__ import annotations

import argparse
import json
import statistics
import sys
from pathlib import Path


def fit_rate(records: list[dict]) -> tuple[float | None, float | None, int]:
    """Least squares ``wall = startup + trials / rate`` over runs.

    Returns ``(rate, startup, n)``.  ``(None, None, n)`` when the fit is not
    identified -- fewer than two runs, or no spread in the trial counts, which
    makes the slope zero; a non-positive slope means throughput came out
    negative and is not a rate.

    ``startup`` is returned signed.  A negative intercept is not clamped away
    because it is evidence: it says the runs share work that does not scale
    with trials, and reporting it as zero would present a fit that does not
    describe the data as though it did.
    """
    pts = [
        (float(r["trials"]), float(r["wall"]))
        for r in records
        if r.get("trials") and r.get("wall")
    ]
    n = len(pts)
    if n < 2:
        return None, None, n
    sx = sum(t for t, _ in pts)
    sy = sum(w for _, w in pts)
    sxx = sum(t * t for t, _ in pts)
    sxy = sum(t * w for t, w in pts)
    denom = n * sxx - sx * sx
    if denom == 0:
        return None, None, n
    slope = (n * sxy - sx * sy) / denom
    if slope <= 0:
        return None, None, n
    return 1.0 / slope, (sy - slope * sx) / n, n


def rate_by_config(records: list[dict]) -> dict[str, float]:
    """Median trials/s per configuration -- the one figure both arms share."""
    out: dict[str, list[float]] = {}
    for rec in records:
        if rec.get("rate") and rec.get("complete"):
            cfg = rec["run_id"].rsplit("_s", 1)[0]
            out.setdefault(cfg, []).append(float(rec["rate"]))
    return {cfg: statistics.median(v) for cfg, v in out.items()}


def _line(tag: str, rate: float | None, startup: float | None, n: int) -> None:
    if rate is None:
        print(f"  {tag:<10}{n:>4} runs   fit not identified")
        return
    print(f"  {tag:<10}{n:>4} runs   {rate:>10,.0f} trials/s   "
          f"startup {startup:+.2f} s")


def compare(probe: dict, ref: dict) -> dict:
    """Print the penalty and the node-hours; return the numbers printed."""
    print("== packing penalty: 192 runs at one per core vs the 8-run bench")
    print("  throughput from wall = startup + trials / rate")

    p_rate, p_start, p_n = fit_rate(probe.get("runs", []))
    r_rate, r_start, r_n = fit_rate(ref.get("runs", []))
    _line("reference", r_rate, r_start, r_n)
    _line("probe", p_rate, p_start, p_n)
    ratio = (r_rate / p_rate) if (r_rate and p_rate) else None
    if ratio:
        print(f"  rate ratio {ratio:.2f} x   "
              f"(1.00 = packing a node costs nothing)")

    # Configuration by configuration where both timed one, as a check on the
    # fit: it uses runs the fit pools, but through raw trials/s rather than
    # through the regression.
    p_cfg, r_cfg = rate_by_config(probe.get("runs", [])), rate_by_config(ref.get("runs", []))
    shared = sorted(set(p_cfg) & set(r_cfg))
    per_cfg = [r_cfg[c] / p_cfg[c] for c in shared if p_cfg[c]]
    if per_cfg:
        print(f"  per-configuration ratio over the {len(shared)} both timed: "
              f"median {statistics.median(per_cfg):.2f} x   "
              f"range {min(per_cfg):.2f}-{max(per_cfg):.2f} x")

    print("\n== node-hours per run, 20,000 events, one run per physical core")
    out: dict = {"rate_ratio": ratio, "per_config_median": None}
    p_nh, r_nh = probe.get("node_hours"), ref.get("node_hours")
    if p_nh:
        print(f"  measured over {len(p_cfg)} configurations: {p_nh:.6f}")
        out["probe_node_hours"] = p_nh
    if r_nh:
        print(f"  uncontended, the {len(r_cfg)} the bench timed: "
              f"{r_nh:.6f}  (extrapolated to 128 in flight)")
        out["ref_node_hours"] = r_nh
    if per_cfg:
        out["per_config_median"] = statistics.median(per_cfg)
    return out


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("probe", type=Path, help="summary.json of the packed arm")
    p.add_argument("reference", type=Path, help="summary.json of the 8-run bench")
    args = p.parse_args(argv)

    missing = [str(f) for f in (args.probe, args.reference) if not f.is_file()]
    if missing:
        # The tables have already printed by here, so say what is absent rather
        # than raising -- but fail, because a missing summary means bench_report
        # did not finish and the numbers above have nothing behind them.
        print(f"\n== packing penalty: no summary ({', '.join(missing)})",
              file=sys.stderr)
        return 1

    probe = json.loads(args.probe.read_text())
    ref = json.loads(args.reference.read_text())
    compare(probe, ref)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
