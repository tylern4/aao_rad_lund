#!/usr/bin/env python3
"""Per-configuration event quotas for the 192-run packing probe.

``bench_fortran.sbatch`` times four runs that never contend with each other,
which is the rate the node-hours arithmetic divides by -- but production put
192 processes on the same 128-core node, and nothing has measured what that
costs.  Production's own history cannot answer it: the arm ran across three
cancelled jobs (14:27 + 5:04 + 10:44) that each restarted the runs the
previous one left unfinished, so its wall clocks mix contention with thrown-away
work and a queue gap.

This sizes every one of the 96 configurations to the same *uncontended* budget
so the whole grid can be run at production's exact layout inside one 30-minute
debug job.  Equal budgets matter: a uniform event quota would make the run at
delta 0.005 (median 58,566 trials per event) six times longer than the one at
delta 0.05 (9,118), so the packing penalty would be read off whichever configs
happened to be expensive.

Trials per event come from the production scan's own ``out.txt``, whose
counters are checked every ``n_events/25`` events.  ``ntries`` is ``integer*4``
(src/aao_rad.f90:154) and goes *negative* partway through for the 23
configurations whose trials exceed 2**31 -- after that the sample is no longer a
count, and taking the last one, as the verifier deliberately does for a
different question, would divide by a wrapped number.  Here only the leading run
of non-decreasing samples is read.

Usage::

    make_probe_quota.py --scan <production root> --out <override csv>
"""

from __future__ import annotations

import argparse
import csv
import re
import statistics
from pathlib import Path

import make_grid  # for N_EVENTS -- see event_quota

NTRIES_RE = re.compile(r"ntries, nevent, mcall_max:\s*(-?\d+)\s+(\d+)")

# Measured uncontended rate of one Fortran process on a PM-CPU node, from the
# eight runs of bench_fortran.sbatch (225,812-239,620 trials/s).
TRIALS_PER_SECOND = 234_000

# Seconds of uncontended work every run is sized for.  The probe's whole point
# is that it survives a large packing penalty: at 80x this is 1,600 s, still
# inside the 30-minute debug cap with room for the grid build.  It is also long
# enough that the ~1 s of process startup, table load and integration is a few
# percent rather than a rounding error.
TARGET_SECONDS = 20.0


def trials_per_event(text: str) -> float | None:
    """Trials per event from one ``out.txt``, or None if it says nothing usable.

    Cumulative trials over cumulative events, taken at the end of the longest
    prefix of strictly increasing ``ntries``.  The first sample that is negative
    or smaller than its predecessor ends the prefix, because from there on the
    numerator is a wrapped 32-bit counter; and a file with no positive sample at
    all yields None rather than a ratio against zero.
    """
    last: tuple[int, int] | None = None
    for match in NTRIES_RE.finditer(text):
        ntries, nevent = int(match.group(1)), int(match.group(2))
        if ntries < 0 or nevent <= 0:
            break
        if last is not None and ntries <= last[0]:
            break
        last = (ntries, nevent)
    if last is None:
        return None
    return last[0] / last[1]


def event_quota(tpe: float, target: float = TARGET_SECONDS, rate: int = TRIALS_PER_SECOND) -> int:
    """Events that buy *target* seconds of uncontended work at *tpe* trials/event.

    Clamped at both ends: a configuration expensive enough to need fewer than
    one event still has to produce one, and a cheap one -- trials per event near
    unity -- would otherwise ask for millions, which :func:`make_grid.load_overrides`
    rejects as *above* the standard quota and takes the whole probe down before
    it measures anything.  Sharing make_grid's constant rather than repeating
    20,000 keeps the two from drifting the way the run card and manifest did.
    """
    return min(make_grid.N_EVENTS, max(1, round(target * rate / tpe)))


def scan_rates(scan: Path) -> dict[str, float]:
    """Median trials per event per configuration across the production scan."""
    outs = sorted((scan / "fortran").glob("cfg_*_s*/out.txt"))
    if not outs:
        raise SystemExit(f"{scan}/fortran holds no run directories -- is that the scan root?")
    by_cfg: dict[str, list[float]] = {}
    for out in outs:
        tpe = trials_per_event(out.read_text(errors="replace"))
        if tpe is not None:
            by_cfg.setdefault(out.parent.name.rsplit("_s", 1)[0], []).append(tpe)
    return {cfg: statistics.median(vals) for cfg, vals in by_cfg.items()}


def write_quotas(
    scan: Path, out: Path, target: float = TARGET_SECONDS
) -> tuple[dict[str, int], dict[str, float]]:
    """Write ``cfg_id,n_events`` rows for every configuration in *scan*.

    Returns the quotas and the rates they were derived from, so a caller can
    report the serial work without reading 192 files twice.

    Every configuration must have a measurable rate.  One missing would silently
    fall back to the full 20,000-event quota, which is seven hours of serial
    work for cfg_034 -- the probe would then hit its wall limit and report
    nothing, having spent the job finding out that a config was absent.
    """
    rates = scan_rates(scan)
    quotas = {cfg: event_quota(tpe, target) for cfg, tpe in rates.items()}
    missing = [f"cfg_{i:03d}" for i in range(96) if f"cfg_{i:03d}" not in quotas]
    if missing:
        raise SystemExit(
            f"{scan}/fortran has no usable trials-per-event for "
            f"{len(missing)} configuration(s): {', '.join(missing)}"
        )
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["cfg_id", "n_events", "trials_per_event", "trials_per_second"])
        for cfg in sorted(quotas):
            writer.writerow([cfg, quotas[cfg], round(rates[cfg]), TRIALS_PER_SECOND])
    return quotas, rates


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--scan", type=Path, required=True, help="production scan root")
    ap.add_argument("--out", type=Path, required=True, help="override csv to write")
    ap.add_argument(
        "--target",
        type=float,
        default=TARGET_SECONDS,
        help=f"uncontended seconds of work per run (default {TARGET_SECONDS})",
    )
    args = ap.parse_args(argv)

    quotas, rates = write_quotas(args.scan, args.out, args.target)
    values = sorted(quotas.values())
    trials = 2 * sum(quotas[cfg] * rates[cfg] for cfg in quotas)  # both seeds
    print(
        f"{args.out}: {len(quotas)} configurations, quota "
        f"{values[0]}-{values[-1]} events (median {statistics.median(values):.0f}) "
        f"for {args.target:.0f}s of uncontended work each"
    )
    print(
        f"  both seeds: {sum(values) * 2} events, {trials:.3g} trials, "
        f"{trials / TRIALS_PER_SECOND / 3600:.2f} node-h of serial work -- "
        f"{trials / TRIALS_PER_SECOND / 128 / 60:.1f} min at 128-way"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
