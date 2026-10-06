#!/usr/bin/env python3
"""The two numbers the benchmarks exist to produce: speedup and node-hours.

Both arms are priced as seconds for the same 20,000-event run, because they do
not do the same work per event -- the port accepts ~1.2-2.9x more often, so
trials per second alone would credit the Fortran's extra work to the machine and
seconds per event alone would credit the port's better acceptance to speed.

Each side's trials per event is taken where it is least noisy.  ``sigr_max``,
the rejection ceiling, is re-estimated from each run's own random points, so a
short run's cumulative trials per event is one draw of a quantity that moves
2-3x between two runs of the same configuration: the probe's cfg_000 runs read
7,453 and 9,592 against production's 16,181 and 11,434, and cfg_064 read 61,432
against 24,910.  The Fortran cost of a full run is therefore rebuilt as
``20000 * production_tpe / fitted_rate`` -- production's tpe comes from the full
scan's 20,000-event runs, whose checkpoints are stable within the run, and the
rate comes from fitting ``wall = startup + trials/rate`` over the uncontended
bench so the ~1 s of process startup is not counted as throughput.

Node-hours are quoted at the packing each arm actually runs at: 128 in flight
for the Fortran (one process per physical core, measured by ``probe_packing``),
4 for the GPU.  A 50x per-run speedup against a 32x concurrency deficit buys
about 1.6x, which is the whole economic question.

Usage::

    speedup_report.py            # uses $SCRATCH/aao_bench and $SCRATCH/aao_rad_scan
    speedup_report.py --root /path/to/aao_bench --scan /path/to/aao_rad_scan
"""

from __future__ import annotations

import argparse
import json
import os
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import make_probe_quota as q  # noqa: E402
import packing_penalty as pp  # noqa: E402

# The four configurations bench_fortran and bench_gpu both timed, spanning the
# delta axis that sets how expensive a configuration is.
CFGS = ("cfg_000", "cfg_034", "cfg_064", "cfg_082")

# Runs in flight, which is what node-hours per run is divided by.
FORT_PACK = 128
GPU_PACK = 4

N_EVENTS = 20_000


def load(path: Path) -> dict:
    """Summary json written by ``bench_report --json``, keyed by run id."""
    if not path.is_file():
        raise SystemExit(f"missing {path} -- run bench_report.py --json first")
    return json.loads(path.read_text())


def index(summary: dict) -> dict[str, dict]:
    return {r["run_id"]: r for r in summary.get("runs", [])}


def med_proj(runs: dict[str, dict], cfg: str) -> float | None:
    """Median projected 20,000-event wall over a configuration's seeds."""
    vals = [
        runs[f"{cfg}_{seed}"]["proj"]
        for seed in ("s0", "s1")
        if runs.get(f"{cfg}_{seed}", {}).get("proj")
    ]
    return statistics.median(vals) if vals else None


def section(title: str) -> None:
    print("\n" + "=" * 74)
    print(title)
    print("=" * 74)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--root", type=Path,
                    default=Path(os.environ.get("SCRATCH", "/tmp")) / "aao_bench")
    ap.add_argument("--scan", type=Path,
                    default=Path(os.environ.get("SCRATCH", "/tmp")) / "aao_rad_scan")
    args = ap.parse_args(argv)

    probe = load(args.root / "timing" / "pack128.json")
    ref = load(args.root / "timing" / "bench8.json")
    warm = load(args.root / "full" / "timing" / "pass_warm.json")
    cold = load(args.root / "full" / "timing" / "pass_cold.json")

    prod_tpe = q.scan_rates(args.scan)          # median over seeds, full runs
    warm_i, cold_i = index(warm), index(cold)

    # The Fortran rate is fitted here rather than remembered from the run that
    # produced it: it is the denominator of every Fortran second in the table,
    # so a stale value would move the whole speedup without anything looking
    # wrong.  bench_report --json writes trials and wall, so the fit has its
    # denominator; if it does not, the summary is from before that change.
    fort_rate, _, n_fit = pp.fit_rate(ref.get("runs", []))
    if not fort_rate:
        raise SystemExit(
            f"cannot fit wall = startup + trials/rate over {n_fit} run(s) with "
            "both -- regenerate the reference bench with bench_report.py --json")

    projs = [r["proj"] for r in probe["runs"] if r.get("proj")]
    grid_sum = sum(projs) / 3600 / FORT_PACK
    grid_fit = 2 * N_EVENTS * sum(prod_tpe.values()) / fort_rate / 3600 / FORT_PACK

    section("FORTRAN, one run per physical core (probe: 192 runs, 96 configs)")
    print(f"  20,000-event run, median config: {statistics.median(projs):6.0f} s of one"
          f" core  = {statistics.median(projs) / 3600 / FORT_PACK:.5f} node-h/run")
    print(f"  20,000-event run, grid mean:     {statistics.mean(projs):6.0f} s of one"
          f" core  = {statistics.mean(projs) / 3600 / FORT_PACK:.5f} node-h/run")
    print(f"  full 192-run grid (sum):          {grid_sum:6.2f} node-h")
    if ref.get("node_hours"):
        print(f"  uncontended, the 4 bench configs: {ref['node_hours']:.5f} node-h/run"
              f"   (cheaper than the grid average -- not the grid's price)")
    print(f"  cross-check, production trials/event x fitted rate: {grid_fit:.2f} node-h")

    section("PACKING PENALTY")
    # Fitted rather than ratioed, because the probe's runs are ~17 s against the
    # reference's ~53 s and a fixed ~1 s of startup is a larger share of the
    # short runs -- charging it to packing would overstate the penalty.  The fit
    # is packing_penalty's, so the numbers here cannot drift from that report.
    p_rate, p_start, p_n = pp.fit_rate(probe.get("runs", []))
    r_rate, r_start, r_n = pp.fit_rate(ref.get("runs", []))
    if p_rate and r_rate:
        print("  wall = startup + trials/rate, fitted over every run")
        print(f"    reference {r_rate:>10,.0f} trials/s  startup {r_start:+.2f} s"
              f"   ({r_n} runs, idle core)")
        print(f"    probe     {p_rate:>10,.0f} trials/s  startup {p_start:+.2f} s"
              f"   ({p_n} runs, one per core)")
        print(f"  rate ratio {r_rate / p_rate:.2f} x   "
              f"(1.00 = packing a node costs nothing)")
    p_cfg, r_cfg = pp.rate_by_config(probe.get("runs", [])), pp.rate_by_config(
        ref.get("runs", []))
    ratios = [r_cfg[c] / p_cfg[c] for c in sorted(set(p_cfg) & set(r_cfg)) if p_cfg[c]]
    if ratios:
        print(f"  per-configuration over the {len(ratios)} both timed: "
              f"median {statistics.median(ratios):.2f} x, "
              f"range {min(ratios):.2f}-{max(ratios):.2f} x")

    section("SPEEDUP AND NODE-HOURS PER RUN, THE 4 BENCHMARK CONFIGURATIONS")
    print(f"  {'config':<9} {'F s':>7} {'G warm':>7} {'G cold':>7} "
          f"{'warm x':>7} {'cold x':>7}   {'F nh':>8} {'Gw nh':>8} {'Gc nh':>8}")
    warm_sp, cold_sp, f_nh, gw_frac, gc_frac = [], [], [], [], []
    for cfg in CFGS:
        if cfg not in prod_tpe:
            print(f"  {cfg:<9} no production trials/event")
            continue
        gp, cp = med_proj(warm_i, cfg), med_proj(cold_i, cfg)
        if not (gp and cp):
            print(f"  {cfg:<9} incomplete in the GPU bench")
            continue
        f_sec = N_EVENTS * prod_tpe[cfg] / fort_rate
        f_n = f_sec / 3600 / FORT_PACK
        gw_n, gc_n = gp / 3600 / GPU_PACK, cp / 3600 / GPU_PACK
        print(f"  {cfg:<9} {f_sec:>7.0f} {gp:>7.0f} {cp:>7.0f} "
              f"{f_sec / gp:>7.1f} {f_sec / cp:>7.1f}   "
              f"{f_n:>8.5f} {gw_n:>8.5f} {gc_n:>8.5f}")
        warm_sp.append(f_sec / gp)
        cold_sp.append(f_sec / cp)
        f_nh.append(f_n)
        gw_frac.append(gw_n / f_n)
        gc_frac.append(gc_n / f_n)

    if warm_sp:
        print(f"\n  speedup      warm median {statistics.median(warm_sp):.1f} x "
              f"(range {min(warm_sp):.1f}-{max(warm_sp):.1f})   "
              f"cold median {statistics.median(cold_sp):.1f} x "
              f"(range {min(cold_sp):.1f}-{max(cold_sp):.1f})")
        print(f"  node-h/run   Fortran median {statistics.median(f_nh):.5f} at "
              f"{FORT_PACK} in flight "
              f"(grid median {statistics.median(projs) / 3600 / FORT_PACK:.5f})")
        print(f"               GPU warm {statistics.median(gw_frac):.2f} x Fortran's "
              f"(range {min(gw_frac):.2f}-{max(gw_frac):.2f}) -- lower is better")
        print(f"               GPU cold {statistics.median(gc_frac):.2f} x Fortran's "
              f"(range {min(gc_frac):.2f}-{max(gc_frac):.2f}) -- lower is better")
        print(f"\n  throughput   Fortran {fort_rate:,.0f} trials/s fitted, "
              f"GPU {warm['median_rate']:,.0f} warm / {cold['median_rate']:,.0f} cold"
              f" = {warm['median_rate'] / fort_rate:.1f} x / "
              f"{cold['median_rate'] / fort_rate:.1f} x")
        print("  the remainder is acceptance: the port rejects fewer points, so it")
        print("  needs 1.2-2.9x fewer trials for the same 20,000 events")

    paid = sum((14 * 3600 + 26 * 60 + 48, 5 * 3600 + 3 * 60 + 59,
                10 * 3600 + 44 * 60 + 4))
    section("WHAT PRODUCTION PAID")
    print(f"  3 cancelled jobs: 14:26:48 + 5:03:59 + 10:44:04 = {paid / 3600:.2f} node-h"
          f"   (1 node each, all CANCELLED)")
    print(f"  measured requirement {min(grid_sum, grid_fit):.2f}-"
          f"{max(grid_sum, grid_fit):.2f} node-h  ->  "
          f"{paid / 3600 / max(grid_sum, grid_fit):.1f}-"
          f"{paid / 3600 / min(grid_sum, grid_fit):.1f}x the compute")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
