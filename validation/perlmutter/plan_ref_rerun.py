#!/usr/bin/env python3
"""Size a rerun for configurations whose Fortran reference is unusable.

The reference is unusable for two different reasons, and they need different
treatments, which is the whole reason this is a separate step rather than a
blanket reduction of the event quota.

``integer*4 ntries`` (``src/aao_rad.f90:154``) counts every rejection-sampling
trial and is a divisor of the printed cross section.  Past 2**31 it wraps, and
the cross section comes back *negative*; if it wraps a second time it comes back
**positive and wrong** -- ``cfg_077_s0`` printed a positive 368,439,614 where
its own progress lines imply 4.66e9, inflating its cross section by 12.6x.
Those configurations need a smaller quota, sized so the wrapped counter never
happens: at 3e5 trials per event the original 20,000-event quota demands 6e9
trials.

A second group has a healthy cross section but a truncated n-tuple: those runs
were killed mid-generation, and because the n-tuple is written as fixed
800-byte records, the partial run overwrote only its own records and left a
stale tail from an earlier attempt.  The cross section still reads correctly --
it estimates a constant, and those runs agree with their own earlier estimates
to 0.1% -- but the *events* are two runs interleaved, so the n-tuple cannot be
used.  These need no reduction at all: they are the cheap configurations, at
~1,700 trials per event, and re-running them at the standard quota costs
nothing.

Prints a plan; ``--write`` saves it where :func:`make_grid.main` will pick it up
on its own.  That matters because every sbatch script re-runs ``make_grid.py``
with only ``--root``, so the quota has to live in the grid rather than in a
command-line flag that the next job would drop.

Usage::

    python plan_ref_rerun.py --root $SCRATCH/aao_rad_scan          # dry run
    python plan_ref_rerun.py --root $SCRATCH/aao_rad_scan --write
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))

from verify_statistics import default_root, defects  # noqa: E402

PROGRESS_RE = re.compile(r"ntries,\s+nevent,\s+mcall_max:\s*(-?\d+)\s+(\d+)")
SIGMA_RE = re.compile(
    r"Integrated cross section \(MC, numerical\) =\s*(-?[\d.Ee+-]+)\s*(-?[\d.Ee+-]+)"
)
NTYPE_MAX = 2**31 - 1

CAUSE_OVERFLOW = "counter_overflow"
CAUSE_SHORT = "short_n_tuple"
CAUSE_MISSING = "no_output"


def trials_per_event(text: str) -> float | None:
    """Trials per event from the last trustworthy progress line.

    ``ntries`` rises monotonically until it wraps, so the largest *positive*
    value is the last pre-wrap sample.  Using the first line instead would
    sample the run before its rejection ceiling has matured, when acceptance is
    still poor, and would understate the quota a configuration can safely take.
    """
    samples = [
        (int(n), int(e))
        for n, e in PROGRESS_RE.findall(text)
        if int(n) > 0 and int(e) > 0
    ]
    if not samples:
        return None
    n, e = max(samples, key=lambda t: t[0])
    return n / e


def run_causes(root: Path, run_id: str, n_events: int) -> tuple[list[str], str | None]:
    """(causes, text) for one expected run, naming a missing run rather than
    skipping it.

    Driven from the manifest rather than from ``collect()``, which skips runs
    that never printed a cross section or never produced an ``out.txt``.  A
    planner that skipped those would omit the configurations most in need of a
    rerun and then report that every reference was usable.
    """
    run_dir = root / "fortran" / run_id
    out_txt = run_dir / "out.txt"
    if not out_txt.is_file():
        return [CAUSE_MISSING], None
    text = out_txt.read_text(errors="replace")
    reasons = defects("fortran", run_dir, text, n_events)
    out: list[str] = []
    if "counter_overflow" in reasons or "non_positive_sigma" in reasons:
        out.append(CAUSE_OVERFLOW)
    if any(r.startswith("short_ntuple") for r in reasons):
        out.append(CAUSE_SHORT)
    return out, text


def plan(root: Path, margin: float) -> dict[str, dict]:
    """{cfg_id: {cause, trials_per_event, n_events, why}} for unusable configs."""
    manifest = root / "grid" / "manifest.csv"
    if not manifest.is_file():
        raise SystemExit(f"{manifest} not found -- generate the grid first")
    with open(manifest, newline="") as fh:
        rows = list(csv.DictReader(fh))

    cfg_info: dict[str, dict] = {}
    for row in rows:
        run_id, cfg, n_events = row["run_id"], row["cfg_id"], int(row["n_events"])
        found, text = run_causes(root, run_id, n_events)
        if not found:
            continue
        entry = cfg_info.setdefault(
            cfg, {"causes": set(), "trials": [], "n_events": n_events}
        )
        entry["causes"].update(found)
        if text:
            tpe = trials_per_event(text)
            if tpe:
                entry["trials"].append(tpe)

    out: dict[str, dict] = {}
    for cfg, entry in sorted(cfg_info.items()):
        n_events = entry["n_events"]
        tpe = max(entry["trials"]) if entry["trials"] else None
        quota = n_events
        why = ""
        if CAUSE_MISSING in entry["causes"]:
            why = "no output at all; rerun at the standard quota"
        elif CAUSE_OVERFLOW in entry["causes"] and tpe:
            safe = int(margin * NTYPE_MAX / tpe)
            if safe < n_events:
                quota = safe
                why = (
                    f"{tpe:,.0f} trials/event x {n_events:,} = "
                    f"{tpe * n_events:.3g} trials would pass 2**31"
                )
            else:
                why = "no reduction needed: the counter stays in range"
        elif CAUSE_SHORT in entry["causes"]:
            why = "killed mid-run; rerun at the standard quota"
        out[cfg] = {
            "causes": sorted(entry["causes"]),
            "trials_per_event": tpe,
            "n_events": quota,
            "original": n_events,
            "why": why,
        }
    return out


def main() -> int:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--root", type=Path, default=None)
    p.add_argument(
        "--margin", type=float, default=0.9,
        help="fraction of 2**31 the worst run may reach (default 0.9)",
    )
    p.add_argument(
        "--write", action="store_true",
        help=f"write <root>/grid/n_events_override.csv for make_grid.py to pick up",
    )
    args = p.parse_args()
    root = args.root or default_root()

    entries = plan(root, args.margin)
    if not entries:
        print("every configuration has a usable Fortran reference; nothing to do")
        return 0

    reduced = {c: e for c, e in entries.items() if e["n_events"] < e["original"]}
    short_only = {
        c: e for c, e in entries.items()
        if e["causes"] == [CAUSE_SHORT]
    }

    print(f"root        {root}")
    print(f"margin      {args.margin:.2f} of 2**31 ({args.margin * NTYPE_MAX:,.0f})")
    print(f"unusable configs: {len(entries)}")
    print(f"  need a reduced quota : {len(reduced)}")
    print(f"  rerun at standard    : {len(short_only)}")
    print()
    print(f"{'cfg':<10} {'cause':<9} {'trials/ev':>11} {'original':>9} {'rerun':>9}  note")
    for cfg, e in entries.items():
        cause = "overflow" if CAUSE_OVERFLOW in e["causes"] else "short"
        if CAUSE_OVERFLOW in e["causes"] and CAUSE_SHORT in e["causes"]:
            cause = "both"
        tpe = f"{e['trials_per_event']:,.0f}" if e["trials_per_event"] else "?"
        print(
            f"{cfg:<10} {cause:<9} {tpe:>11} {e['original']:>9,} "
            f"{e['n_events']:>9,}  {e['why']}"
        )

    if reduced:
        worst = max(
            (e["trials_per_event"] * e["n_events"] for e in reduced.values()),
            default=0.0,
        )
        print(
            f"\nworst projected ntries at these quotas: {worst:,.0f} "
            f"({100 * worst / NTYPE_MAX:.0f}% of 2**31)"
        )

    if not args.write:
        print("\n(dry run: add --write to save the override)")
        return 0

    out = root / "grid" / "n_events_override.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["cfg_id", "n_events"])
        for cfg, e in sorted(reduced.items()):
            w.writerow([cfg, e["n_events"]])
    print(f"\nwrote {len(reduced)} override(s) to {out}")
    print("re-run make_grid.py to apply them to the run cards and the manifest")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
