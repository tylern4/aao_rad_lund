#!/usr/bin/env python3
"""Generate k independent Fortran reference runs and combine them.

A single 8000-event reference run carries a ~0.3% statistical error on its
integrated cross section -- the same order as the port-vs-Fortran residual
being measured, so a single run cannot resolve it.  Independent runs combine
like any MC estimate: the trial-weighted mean of k runs is k times more
precise, and the run-to-run scatter is the honest error bar (it includes the
heavy tail of the weight distribution, which a naive ``1/sqrt(N)`` does not).

The same runs also concatenate into a larger reference n-tuple, which is what
the KS floors actually need: the two-sample 5% critical value
``1.36 * sqrt(1/N_py + 1/N_f77)`` is dominated by the reference sample, so a
``k x 8000``-event reference tightens every floor by roughly ``sqrt(k)``.

Seeding: ``myran`` is initialised from the unix clock at one-second
resolution (``src/unixtime.c:37``, read at ``aao_rad.f90:375``), so runs
launched inside the same second would be *identical*, not independent.
Launches are therefore staggered by ``--stagger`` seconds and every resulting
n-tuple is checksummed against the others; a collision is reported and the
offending runs are excluded from both the combination and the concatenation
rather than silently double-counted.

Each run needs its own directory: the MAID tables are opened relative to the
current directory (``src/maid_lee.f90:40``) and the n-tuple path is a
compile-time constant (``src/aao_rad.f90:199``), so per-run directories with
a ``spp_tbl`` symlink are the only way to run in parallel.

Usage
-----
    python validation/reference_fleet.py \
        --binary build/aao_rad \
        --run-card /tmp/frun_fixed/run_card.txt \
        --root /tmp/frun_fleet --n-runs 16

Resumable: a run directory whose ``out.txt`` already parses is kept, so an
interrupted fleet can be re-run with the same command.
"""

from __future__ import annotations

import argparse
import hashlib
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

SIG_RE = re.compile(
    r"Integrated cross section \(MC, numerical\) =\s*"
    r"(\S+)\s+(\S+)"
)
NTRIES_RE = re.compile(r"ntries, nevent, mcall_max:\s*(\d+)\s+(\d+)")


def parse_out_txt(path: Path) -> dict | None:
    """Pull the final ``sig_sum`` and trial count out of one run's stdout.

    ``sig_sum`` (the second number) is the importance-sampling estimate
    ``sum(w) / tries * phase_space`` -- the one the sigma tables quote.
    ``sig_int`` (the first) is the acceptance-fraction estimate and is kept
    only for the record.
    """
    text = path.read_text()
    sigs = SIG_RE.findall(text)
    trys = NTRIES_RE.findall(text)
    if not sigs or not trys:
        return None
    sig_int, sig_sum = float(sigs[-1][0]), float(sigs[-1][1])
    ntries, nevent = int(trys[-1][0]), int(trys[-1][1])
    return {
        "sig_int": sig_int,
        "sig_sum": sig_sum,
        "ntries": ntries,
        "nevent": nevent,
    }


def md5(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def launch_run(
    run_dir: Path,
    binary: Path,
    card: Path,
    spp_tbl: Path,
    env_extra: dict[str, str],
) -> subprocess.Popen:
    """Start one Fortran run in its own directory."""
    if not (run_dir / "spp_tbl").exists():
        (run_dir / "spp_tbl").symlink_to(spp_tbl)
    shutil.copy(card, run_dir / "run_card.txt")
    env = dict(os.environ, **env_extra)
    with (run_dir / "run_card.txt").open("rb") as stdin, (run_dir / "out.txt").open("wb") as out:
        return subprocess.Popen(
            [str(binary)],
            cwd=run_dir,
            stdin=stdin,
            stdout=out,
            stderr=subprocess.STDOUT,
            env=env,
        )


def combine(runs: list[dict]) -> dict:
    """Trial-weighted mean of ``sig_sum`` with the run-to-run scatter.

    Each run's ``sig_sum`` is ``sum(w) * PS / tries``, so the pooled estimate
    over all runs is the trial-weighted mean; the unweighted mean and the
    scatter are reported alongside because with equal event targets the two
    agree to much better than the scatter itself.
    """
    n = len(runs)
    total_tries = sum(r["ntries"] for r in runs)
    pooled = sum(r["sig_sum"] * r["ntries"] for r in runs) / total_tries
    sigs = [r["sig_sum"] for r in runs]
    mean = sum(sigs) / n
    if n > 1:
        var = sum((s - mean) ** 2 for s in sigs) / (n - 1)
        sd = var**0.5
    else:
        sd = 0.0
    return {
        "n": n,
        "pooled": pooled,
        "mean": mean,
        "sd": sd,
        "sd_rel": sd / mean if mean else 0.0,
        "se": sd / n**0.5 if n else 0.0,
        "total_tries": total_tries,
        "total_events": sum(r["nevent"] for r in runs),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--binary", type=Path, default=ROOT / "build" / "aao_rad")
    ap.add_argument("--run-card", type=Path, required=True)
    ap.add_argument("--root", type=Path, default=Path("/tmp/frun_fleet"))
    ap.add_argument("--n-runs", type=int, default=16)
    ap.add_argument("--stagger", type=float, default=2.0,
                    help="seconds between launches; the seed clock has 1 s "
                         "resolution, so this must stay above 1")
    ap.add_argument("--spp-tbl", type=Path, default=ROOT / "parms" / "spp_tbl")
    ap.add_argument("--n-python", type=int, default=200_000,
                    help="port sample size, used only to print the KS floor "
                         "the concatenated reference buys")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    card = args.run_card.resolve()
    binary = args.binary.resolve()
    spp_tbl = args.spp_tbl.resolve()
    if not binary.is_file():
        sys.exit(f"no Fortran binary at {binary}")
    if not card.is_file():
        sys.exit(f"no run card at {card}")

    env_extra = {"AAO_TRIALS": "0"}
    root = args.root
    root.mkdir(parents=True, exist_ok=True)

    # ---- launch, staggered ------------------------------------------------
    procs: list[tuple[int, subprocess.Popen]] = []
    t0 = time.perf_counter()
    for i in range(args.n_runs):
        run_dir = root / f"run_{i:02d}"
        run_dir.mkdir(exist_ok=True)
        done = parse_out_txt(run_dir / "out.txt") if (run_dir / "out.txt").exists() else None
        if done is not None and (run_dir / "aao_rad.ntuple").is_file():
            print(f"run_{i:02d}: already complete (sigma {done['sig_sum']:.6e}), skipping")
            continue
        if args.dry_run:
            print(f"run_{i:02d}: would launch in {run_dir}")
            continue
        procs.append((i, launch_run(run_dir, binary, card, spp_tbl, env_extra)))
        print(f"run_{i:02d}: launched in {run_dir} "
              f"(t+{time.perf_counter() - t0:.0f}s)")
        if i < args.n_runs - 1:
            time.sleep(args.stagger)
    if args.dry_run:
        return 0

    # ---- wait -------------------------------------------------------------
    failed = []
    for i, proc in procs:
        rc = proc.wait()
        if rc != 0:
            failed.append(i)
            print(f"run_{i:02d}: exit {rc}")

    # ---- collect ----------------------------------------------------------
    runs, collisions = [], []
    checksums: dict[str, list[int]] = {}
    for i in range(args.n_runs):
        run_dir = root / f"run_{i:02d}"
        ntuple = run_dir / "aao_rad.ntuple"
        parsed = parse_out_txt(run_dir / "out.txt") if (run_dir / "out.txt").exists() else None
        if parsed is None or not ntuple.is_file():
            if i not in failed:
                failed.append(i)
            print(f"run_{i:02d}: no parseable result, excluded")
            continue
        if parsed["nevent"] == 0:
            print(f"run_{i:02d}: zero events, excluded")
            continue
        checksums.setdefault(md5(ntuple), []).append(i)
        runs.append({"idx": i, **parsed})

    # identical checksums share a seed: keep the first of each group
    for idxs in checksums.values():
        if len(idxs) > 1:
            collisions.extend(idxs[1:])
            print(f"SEED COLLISION: runs {idxs} wrote identical n-tuples; "
                  f"keeping {idxs[0]}, excluding {idxs[1:]}")
    if collisions:
        runs = [r for r in runs if r["idx"] not in collisions]

    if not runs:
        sys.exit("no usable runs")

    # ---- combine ----------------------------------------------------------
    c = combine(runs)
    print()
    print(f"runs combined      : {c['n']} of {args.n_runs} "
          f"({c['total_events']} events, {c['total_tries']:.3g} trials)")
    for r in runs:
        print(f"  run_{r['idx']:02d}: sigma = {r['sig_sum']:.6e}  "
              f"sig_int = {r['sig_int']:.6e}  ntries = {r['ntries']}")
    print()
    print(f"pooled sigma       : {c['pooled']:.7e} micro-barn  "
          f"(trial-weighted)")
    print(f"mean of runs       : {c['mean']:.7e} micro-barn")
    print(f"run-to-run sd      : {c['sd']:.3e}  ({c['sd_rel'] * 100:.3f}% relative)")
    print(f"error on pooled    : {c['se']:.3e}  ({c['se'] / c['mean'] * 100:.3f}% relative)")

    # ---- concatenate the reference n-tuple --------------------------------
    if not collisions or len(runs) > 1:
        combined_ntuple = root / "aao_rad_combined.ntuple"
        with combined_ntuple.open("wb") as out:
            for r in runs:
                out.write((root / f"run_{r['idx']:02d}" / "aao_rad.ntuple").read_bytes())
        n_ref = c["total_events"]
        floor_old = 1.36 * (1.0 / args.n_python + 1.0 / 8000) ** 0.5
        floor_new = 1.36 * (1.0 / args.n_python + 1.0 / n_ref) ** 0.5
        print()
        print(f"combined n-tuple   : {combined_ntuple}  ({n_ref} events)")
        print(f"KS floor           : {floor_old:.4f} (single 8000) -> "
              f"{floor_new:.4f} ({n_ref} events), for N_py = {args.n_python}")

    summary = root / "combined_sigma.txt"
    summary.write_text(
        f"runs: {c['n']}\n"
        f"events: {c['total_events']}\n"
        f"trials: {c['total_tries']}\n"
        f"pooled sigma (micro-barn): {c['pooled']:.10e}\n"
        f"mean sigma (micro-barn): {c['mean']:.10e}\n"
        f"run-to-run sd: {c['sd']:.10e}\n"
        f"run-to-run sd relative: {c['sd_rel']:.6f}\n"
        f"error on pooled: {c['se']:.10e}\n"
        f"error on pooled relative: {c['se'] / c['mean']:.6f}\n"
    )
    print(f"summary            : {summary}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
