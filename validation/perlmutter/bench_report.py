#!/usr/bin/env python3
"""Wall clock, throughput and node-hours for the speedup benchmark.

The wall clock comes from ``timing.json``, which ``run_batch.execute()`` writes
after the child exits: ``start``, ``end``, ``wall`` and the return code.  That
is the driver's own ``time.time()`` bracket around the process, and it is the
number node-hours have to be computed from.

Reading the clock off the filesystem instead -- ``out.txt``'s birth time to its
mtime -- is kept only as a fallback, and the report flags when it was used,
because on Perlmutter it is wrong in a way the numbers do not reveal.  ``stat
-c %W`` reports the Lustre *server's* clock while mtime is stamped by the node,
and the two disagree by ~490 s; birth-to-mtime therefore measured a clock skew
and charged it to every run.  The first benchmark's pass 1 reported 490.0 s for
all eight runs -- identical, because they were all reading the same constant --
for runs that had each taken about 20 s, and the 490 s of "startup" derived
from it was that skew rather than a compilation.

Recording the timing in the driver also removes the conditions the filesystem
method had to impose.  A retried run overwrites its own ``timing.json``, so the
wall clock belongs to one attempt instead of spanning two attempts plus the
idle gap between them -- in the production scan 117 of 174 clean Fortran runs
had a median 19.5-hour gap of exactly that kind.  What still matters is:

*   **The run must have finished.**  A process stopped by the wall limit would
    otherwise report a wall clock equal to the limit and a rate made almost
    entirely of the time it did not finish in.  Completeness is therefore the
    same test the verifier uses -- the run printed its cross section -- rather
    than a file merely existing.

The port also reports its own loop time (``wall time`` in ``summary()``,
``generate.py:738``), measured from ``t_start`` inside ``generate.py`` after
configuration and table loading.  The difference between the two is everything
the port does before the sampling loop starts: interpreter startup, the JAX
import, loading the MAID tables, and -- when the XLA cache is cold -- compiling
the jitted kernels.  That difference is what ``bench_gpu.sbatch`` exists to
measure, so it is reported as ``startup`` rather than being folded into a rate.

Usage::

    bench_report.py [--label TEXT] [--concurrency N] [--target-events N] DIR
"""

from __future__ import annotations

import argparse
import json
import re
import statistics
import subprocess
import sys
from pathlib import Path

# The port's summary() block (generate.py:719).  All three are needed, not
# just the cross section: the counts give the trials per event, which is what
# scales a reduced quota back up to the production one, and the loop time is
# what separates compilation from sampling.
PORT_EVENTS_RE = re.compile(r"^events\s*:\s*([\d,]+)\s*$", re.M)
PORT_TRIALS_RE = re.compile(r"^trials\s*:\s*([\d,]+)\s*$", re.M)
PORT_LOOP_RE = re.compile(r"^wall time\s*:\s*([\d.]+)\s*s", re.M)
PORT_SUMMARY_RE = re.compile(r"^sigma \(MC\)\s*:", re.M)

# The Fortran's closing lines.  The counters are the last of their kind in the
# file, so the scan is a last-match search -- ``re.search`` takes the first, and
# an earlier progress line of a different shape made one attempt report a 25x
# speedup that did not exist.
FORT_NTRIES_RE = re.compile(r"ntries, nevent, mcall_max:\s*(\d+)\s+(\d+)")
FORT_SUMMARY_RE = re.compile(r"Integrated cross section \(MC, numerical\)")

DEFAULT_TARGET = 20_000


def _num(text: str) -> int:
    return int(text.replace(",", ""))


def _last(rx: re.Pattern, text: str, group: int = 1) -> str | None:
    found = None
    for m in rx.finditer(text):
        found = m
    return found.group(group) if found else None


def birth_time(path: Path) -> float | None:
    """Creation instant of *path* in epoch seconds, or None if unknowable.

    ``os.stat`` exposes ``st_birthtime`` on BSD and macOS but not on Linux, so
    the portable route there is coreutils' ``%W`` -- which is what the Perlmutter
    jobs use.  ``%W`` prints 0 when the filesystem cannot say, and 0 is not an
    epoch worth reporting, so it becomes None rather than a wall clock of 1970.

    On Lustre a non-zero value is not to be trusted either: it is the server's
    clock, and Pscratch's runs ~490 s behind the node's.  This is only ever used
    when ``timing.json`` is missing, and the report says so.
    """
    st = path.stat()
    bt = getattr(st, "st_birthtime", None)
    if bt:
        return float(bt)
    try:
        out = subprocess.run(
            ["stat", "-c", "%W", str(path)],
            capture_output=True,
            text=True,
            check=True,
        )
        value = int(out.stdout.strip())
    except (OSError, subprocess.CalledProcessError, ValueError):
        return None
    return float(value) or None


def driver_timing(path: Path) -> dict | None:
    """The ``timing.json`` ``run_batch.execute()`` wrote next to *path*.

    Absent for anything run before the driver recorded it, and unreadable
    rather than fatal if a run died before writing one -- a missing timing is a
    fallback, not an error.
    """
    try:
        data = json.loads((path.parent / "timing.json").read_text())
    except (OSError, ValueError):
        return None
    if not isinstance(data, dict):
        return None
    wall = data.get("wall")
    if not isinstance(wall, (int, float)) or wall <= 0:
        return None
    return data


def read_run(path: Path) -> dict:
    """Parse one run's ``out.txt`` into counts, completion and wall clock."""
    rec: dict = {
        "run_id": path.parent.name,
        "kind": None,
        "complete": False,
        "events": None,
        "trials": None,
        "loop": None,
        "wall": None,
        "wall_source": None,
    }
    try:
        text = path.read_text(errors="replace")
    except OSError as exc:
        rec["error"] = str(exc)
        return rec

    events = _last(PORT_EVENTS_RE, text)
    trials = _last(PORT_TRIALS_RE, text)
    loop = _last(PORT_LOOP_RE, text)
    if events is not None and trials is not None and loop is not None and PORT_SUMMARY_RE.search(text):
        # A negative ntries cannot come from these unsigned patterns, so a
        # wrapped counter lands in the "did not finish" branch below.
        rec.update(
            kind="port",
            events=_num(events),
            trials=_num(trials),
            loop=float(loop),
            complete=True,
        )
    else:
        last = None
        for m in FORT_NTRIES_RE.finditer(text):
            last = m
        if last is not None and FORT_SUMMARY_RE.search(text):
            trials_f, events_f = int(last.group(1)), int(last.group(2))
            if trials_f > 0 and events_f > 0:
                rec.update(kind="fortran", events=events_f, trials=trials_f, complete=True)
        if rec["kind"] is None:
            if FORT_NTRIES_RE.search(text) or FORT_SUMMARY_RE.search(text):
                rec["kind"] = "fortran"
            elif PORT_LOOP_RE.search(text) or "throughput" in text:
                rec["kind"] = "port"
            else:
                rec["kind"] = "?"

    # Driver timing wins: it is the driver's own clock around the child.
    # Birth-to-mtime is the fallback for anything run before the driver wrote
    # one, and it is labelled so the report can say how many runs it is guessing
    # at rather than presenting a skewed number as a measurement.
    timing = driver_timing(path)
    if timing is not None:
        rec["wall"] = float(timing["wall"])
        rec["wall_source"] = "driver"
        return rec
    born = birth_time(path)
    if born is not None:
        rec["wall"] = max(path.stat().st_mtime - born, 0.0)
        rec["wall_source"] = "filesystem"
    return rec


def derive(rec: dict, target: int = DEFAULT_TARGET) -> dict:
    """Fill in rate, startup and the projection back to a production-length run.

    The port's rate is trials per second *inside the loop*; charging compilation
    to it would make the sampling appear slower than it is, and the compilation
    is charged separately to ``startup``.  The Fortran has no loop timer, so its
    whole wall clock is the rate -- which is fair, because its startup is a
    serial table load and a few seconds of process start, not a compilation.

    ``proj`` is the wall clock the same run would have taken at ``target``
    events: the startup, which does not scale with the event count, plus the
    loop scaled by the run's own trials per event.  Sizing the projection from
    the run's own acceptance rather than from a nominal one matters because
    trials per event spans an order of magnitude across this grid and differs by
    a further factor of ~5 between the two implementations.
    """
    rec["rate"] = None
    rec["startup"] = None
    rec["proj"] = None
    if not rec["complete"] or not rec["events"] or not rec["trials"]:
        return rec

    wall, loop = rec["wall"], rec["loop"]
    if loop:
        rec["rate"] = rec["trials"] / loop
        rec["startup"] = max((wall - loop), 0.0) if wall is not None else None
    elif wall:
        rec["rate"] = rec["trials"] / wall
        rec["startup"] = 0.0

    if not rec["rate"]:
        return rec
    trials_per_event = rec["trials"] / rec["events"]
    rec["proj"] = (rec["startup"] or 0.0) + target * trials_per_event / rec["rate"]
    return rec


def _secs(value: float | None) -> str:
    return "n/a" if value is None else f"{value:.1f}s"


def _count(value: int | None) -> str:
    return "n/a" if value is None else f"{value:,}"


def report(label: str, recs: list[dict], concurrency: int = 1,
           grid: int = 192) -> dict:
    """Print one pass's table, medians and node-hours; return those medians.

    The summary is returned as well as printed because the only interesting
    check on it -- that node-hours scale as wall/concurrency -- is arithmetic on
    the medians, and reading them back out of the formatted table would be a
    test of the formatting.
    """
    print(f"\n== {label}")
    summary: dict = {"total": len(recs), "complete": 0, "concurrency": concurrency}
    if not recs:
        print("  no runs found")
        return summary

    print(
        f"  {'run':<14}{'kind':<9}{'events':>9}{'wall':>10}{'loop':>10}"
        f"{'startup':>10}{'trials/s':>14}{'node-h@N':>12}"
    )
    for rec in recs:
        mark = "" if rec["complete"] else "   <- INCOMPLETE, not timed"
        nodeh = f"{rec['proj'] / 3600 / concurrency:.6f}" if rec["proj"] else "n/a"
        print(
            f"  {rec['run_id']:<14}{str(rec['kind']):<9}{_count(rec['events']):>9}"
            f"{_secs(rec['wall']):>10}{_secs(rec['loop']):>10}"
            f"{_secs(rec['startup']):>10}"
            f"{_count(int(rec['rate'])) if rec['rate'] else 'n/a':>14}"
            f"{nodeh:>12}{mark}"
        )

    good = [r for r in recs if r["complete"] and r["proj"]]
    print(f"  {len(good)}/{len(recs)} complete")
    summary["complete"] = len(good)
    # Say how many of these are measurements.  A filesystem-timed run on Lustre
    # carries the server/node clock skew as a constant, which reads as a larger
    # wall and so a slower rate -- it biases the speedup *against* whichever arm
    # timed that way, quietly, with no number in the table looking wrong.
    guessed = [r for r in recs if r.get("wall_source") == "filesystem"]
    summary["filesystem_timed"] = len(guessed)
    if guessed:
        print(
            f"  WARNING: {len(guessed)}/{len(recs)} wall clock(s) read off the "
            f"filesystem because timing.json was missing; on Lustre those add the "
            f"~490 s server/node skew and are upper bounds, not measurements"
        )
    if not good:
        return summary

    def med(key: str) -> float | None:
        vals = [r[key] for r in good if r.get(key) is not None]
        return statistics.median(vals) if vals else None

    wall_m, loop_m = med("wall"), med("loop")
    start_m, rate_m, proj_m = med("startup"), med("rate"), med("proj")
    summary.update(
        median_wall=wall_m, median_loop=loop_m, median_startup=start_m,
        median_rate=rate_m, median_proj=proj_m,
        node_hours=(proj_m / 3600 / concurrency) if proj_m else None,
    )
    print(
        f"  median wall {_secs(wall_m)} | loop {_secs(loop_m)} | "
        f"startup {_secs(start_m)} | rate {_count(int(rate_m)) if rate_m else 'n/a'} trials/s"
    )
    nodeh = summary["node_hours"]
    if nodeh is not None:
        print(f"  node-h per run at concurrency {concurrency}: {nodeh:.6f}"
              f"  (measured wall: {wall_m / 3600 / concurrency:.6f})")
        if grid:
            print(f"  extrapolated {grid}-run grid at this run's own rate: "
                  f"{nodeh * grid:.2f} node-h")
    return summary


_JSON_FIELDS = (
    "run_id", "kind", "complete", "wall", "wall_source", "loop", "startup",
    "events", "trials", "rate", "proj",
)


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("dir", type=Path, help="directory of run directories, each holding an out.txt")
    p.add_argument("--label", default=None, help="heading for the table (default: the directory)")
    p.add_argument(
        "--concurrency",
        type=int,
        default=1,
        help="runs the production layout keeps in flight on one node; node-hours "
        "per run are wall/3600/this (4 for py_gpu, 128 for fortran) "
        "(default: %(default)s)",
    )
    p.add_argument(
        "--target-events",
        type=int,
        default=DEFAULT_TARGET,
        help="event count the projection scales back to (default: %(default)s)",
    )
    p.add_argument(
        "--grid",
        type=int,
        default=192,
        help="runs in the full grid, for the extrapolated total; 0 suppresses it "
        "(default: %(default)s)",
    )
    p.add_argument(
        "--json",
        type=Path,
        default=None,
        help="also write the summary (medians, node_hours, counts) as JSON here, "
        "for a script comparing two arms without scraping the table",
    )
    args = p.parse_args(argv)

    paths = sorted(args.dir.glob("*/out.txt"))
    if (args.dir / "out.txt").is_file():
        paths = [args.dir / "out.txt"]
    if not paths:
        print(f"no */out.txt under {args.dir}", file=sys.stderr)
        return 1

    recs = [derive(read_run(path), args.target_events) for path in paths]
    summary = report(args.label or str(args.dir), recs, args.concurrency, args.grid)
    if args.json is not None:
        # The per-run projections too, not just the medians.  Comparing two
        # arms of different breadth -- all 96 configurations against the four
        # bench_fortran timed -- has to be done configuration by configuration,
        # because trials per event span 400x across the grid and two medians
        # taken over different sets differ for that reason alone.
        summary["runs"] = [
            {k: rec.get(k) for k in _JSON_FIELDS} for rec in recs
        ]
        args.json.write_text(json.dumps(summary, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
