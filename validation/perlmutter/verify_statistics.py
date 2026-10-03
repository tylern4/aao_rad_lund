#!/usr/bin/env python3
"""Statistical verification of the Perlmutter scan.

Three implementations sampled the same 96-configuration grid (beam energies
2-12 GeV crossed with channel, beam polarisation, photon-energy cut and
target thickness), two runs each: the Fortran generator, the JAX port on CPU
and the JAX port on GPU.  This script verifies them at two levels.

**Like-for-like** -- the two runs of one configuration *within* one
implementation.  This measures how much of any cross-code difference is just
Monte Carlo noise:

* run-to-run cross-section scatter, pooled across all configurations into one
  fractional sigma per implementation (the heavy photon tail makes the naive
  ``1/sqrt(n_trials)`` error optimistic, so the scatter is measured, not
  assumed);
* two-sample KS distances per observable, against the 95% noise floor
  ``c*sqrt(2/N)``.

**Cross-code** -- Fortran vs JAX-CPU, Fortran vs JAX-GPU and JAX-CPU vs
JAX-GPU, per configuration:

* cross-section ratio with an error bar built from each implementation's
  pooled run-to-run scatter, and the deviation in sigmas;
* KS distances per observable (the two runs of each side pooled) against
  ``c*sqrt(1/N_A + 1/N_B)``.

CPU and GPU share seeds per configuration, so any CPU-vs-GPU difference is
pure backend floating-point, not sampling; the Fortran seeds from ``unixtime``
and its two runs are simply independent samples.

**Overall agreement** aggregates everything: the fraction of
configuration-by-observable KS checks within the noise floor, the sigma
deviations, and the worst offenders, written to ``summary.txt``.

Both sigma estimators are the same one the sigma tables quote -- mean trial
weight over all trials times phase space (the Fortran's ``sig_sum`` from
``out.txt``, the port's ``sigma (MC)`` line) -- in micro-barns on both sides.

Usage (after the three batch jobs):
    python verify_statistics.py --root $SCRATCH/aao_rad_scan
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import sys
import warnings
from contextlib import contextmanager
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from compare_distributions import OBSERVABLES, ks_test, load_fortran_ntuple  # noqa: E402


@contextmanager
def quiet_scipy():
    """scipy warns when the exact KS p-value falls back to asymptotic; only
    the statistic is used here, so the warning is noise."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        yield


def ks_statistic(a: np.ndarray, b: np.ndarray) -> float:
    with quiet_scipy():
        stat, _ = ks_test(a, b)
    return stat


def fmt(values, spec: str) -> str:
    """Format the max of an array, tolerating all-NaN input."""
    finite = [v for v in np.asarray(values, dtype=float).ravel() if np.isfinite(v)]
    if not finite:
        return "n/a"
    return format(max(finite), spec)


def mean_or_na(values, spec: str) -> str:
    finite = [v for v in np.asarray(values, dtype=float).ravel() if np.isfinite(v)]
    if not finite:
        return "n/a"
    return format(float(np.mean(finite)), spec)


CODES = ("fortran", "py_cpu", "py_gpu")
PAIRS = (("fortran", "py_cpu"), ("fortran", "py_gpu"), ("py_cpu", "py_gpu"))

FORT_SIGMA_RE = re.compile(r"Integrated cross section \(MC, numerical\) =\s*(\S+)\s+(\S+)")
FORT_NTRIES_RE = re.compile(r"ntries, nevent, mcall_max:\s*(\d+)\s+(\d+)")
PORT_SIGMA_RE = re.compile(r"sigma \(MC\)\s*:\s*([0-9.eE+\-]+)")
PORT_TRIALS_RE = re.compile(r"\(([\d,]+) trials/s\)")


def default_root() -> Path:
    env = os.environ.get("AAO_SCAN_ROOT")
    if env:
        return Path(env)
    scratch = os.environ.get("SCRATCH", "/tmp")
    return Path(scratch) / "aao_rad_scan"


def parse_sigma(code: str, run_dir: Path) -> dict:
    """Pull the run's final cross section (micro-barn) and trial count."""
    text = (run_dir / "out.txt").read_text(errors="replace")
    out: dict = {"sigma": None, "trials": None}
    if code == "fortran":
        sigs = FORT_SIGMA_RE.findall(text)
        trys = FORT_NTRIES_RE.findall(text)
        if sigs:
            out["sigma"] = float(sigs[-1][1])  # sig_sum: sum(w)/tries * PS
        if trys:
            out["trials"] = int(trys[-1][0])
    else:
        sigs = PORT_SIGMA_RE.findall(text)
        trs = PORT_TRIALS_RE.findall(text)
        if sigs:
            out["sigma"] = float(sigs[-1])
        if trs:
            out["trials"] = int(trs[-1].replace(",", ""))
    return out


def load_events(code: str, run_dir: Path) -> dict[str, np.ndarray]:
    if code == "fortran":
        return load_fortran_ntuple(run_dir / "aao_rad.ntuple")
    with np.load(run_dir / "out.npz") as data:
        return {name: np.asarray(data[name], dtype=np.float64) for name in data.files}


def collect(root: Path) -> dict[str, dict[str, dict]]:
    """{code: {run_id: {sigma, trials, path, n_events}}} for every complete run."""
    manifest_path = root / "grid" / "manifest.csv"
    with open(manifest_path, newline="") as fh:
        runs = list(csv.DictReader(fh))

    results: dict[str, dict[str, dict]] = {code: {} for code in CODES}
    for code in CODES:
        for run in runs:
            run_dir = root / code / run["run_id"]
            if not (run_dir / "out.txt").is_file():
                continue
            info = parse_sigma(code, run_dir)
            if info["sigma"] is None:
                continue
            info["path"] = run_dir
            info["n_events"] = int(run["n_events"])
            info["cfg_id"] = run["cfg_id"]
            info["seed"] = run["seed"]
            results[code][run["run_id"]] = info
    return results


def pairs_by_config(results: dict[str, dict]) -> dict[str, dict[str, list[dict]]]:
    """{code: {cfg_id: [run, run]}} ordered by seed."""
    out: dict[str, dict[str, list[dict]]] = {}
    for code in CODES:
        by_cfg: dict[str, list[dict]] = {}
        for run in results[code].values():
            by_cfg.setdefault(run["cfg_id"], []).append(run)
        for cfg in by_cfg:
            by_cfg[cfg].sort(key=lambda r: r["seed"])
        out[code] = by_cfg
    return out


def rel_scatter(by_cfg: dict[str, list[dict]]) -> float:
    """Pooled fractional run-to-run scatter of one implementation.

    For two independent runs with fractional sd ``rel`` the difference has
    E[d^2] = 2 (rel sigma)^2, so rel^2 = mean(d^2 / (2 sigma_mean^2)).
    """
    terms = []
    for runs in by_cfg.values():
        if len(runs) < 2:
            continue
        s = np.array([r["sigma"] for r in runs])
        mean = s.mean()
        if mean <= 0:
            continue
        terms.append(((s.max() - s.min()) ** 2) / (2.0 * mean**2))
    if not terms:
        return float("nan")
    return float(np.sqrt(np.mean(terms)))


def ks_floor(n_a: int, n_b: int, c: float) -> float:
    return c * np.sqrt(1.0 / n_a + 1.0 / n_b)


def ks_row_set(events_a: dict, events_b: dict, floor: float) -> list[dict]:
    rows = []
    for title, column, *_ in OBSERVABLES:
        stat = ks_statistic(events_a[column], events_b[column])
        rows.append(
            {
                "observable": title,
                "ks": stat,
                "floor": floor,
                "ok": bool(np.isfinite(stat) and stat <= floor),
            }
        )
    return rows


def main() -> int:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--root", type=Path, default=None)
    p.add_argument(
        "--floor-c", type=float, default=1.36, help="KS noise-floor multiplier (95%% = 1.36)"
    )
    p.add_argument("--out-dir", type=Path, default=None)
    args = p.parse_args()

    root = args.root or default_root()
    out_dir = args.out_dir or root / "verify"
    out_dir.mkdir(parents=True, exist_ok=True)
    c = args.floor_c

    results = collect(root)
    by_cfg = pairs_by_config(results)

    print("== completeness")
    for code in CODES:
        n_cfg = sum(len(v) >= 2 for v in by_cfg[code].values())
        n_runs = len(results[code])
        total = 2 * sum(1 for _ in by_cfg[code])
        print(f"  {code}: {n_runs}/{total} runs complete, {n_cfg} configurations with both runs")

    # ---- like-for-like: run-to-run scatter and KS within each code
    scatter = {code: rel_scatter(by_cfg[code]) for code in CODES}
    print("== pooled run-to-run scatter (fractional, heavy-tail aware)")
    for code in CODES:
        print(f"  {code}: {scatter[code]:.4%}" if np.isfinite(scatter[code]) else f"  {code}: n/a")

    like_rows: list[dict] = []
    ks_rows: list[dict] = []
    for code in CODES:
        for cfg_id, runs in sorted(by_cfg[code].items()):
            if len(runs) < 2:
                continue
            s = [r["sigma"] for r in runs]
            mean = float(np.mean(s))
            z = float("nan")
            if np.isfinite(scatter[code]) and scatter[code] > 0:
                z = abs(s[1] - s[0]) / (np.sqrt(2.0) * scatter[code] * mean)
            events_a = load_events(code, runs[0]["path"])
            events_b = load_events(code, runs[1]["path"])
            n = min(len(events_a["es"]), len(events_b["es"]))
            kss = ks_row_set(events_a, events_b, ks_floor(n, n, c))
            worst = max(kss, key=lambda r: r["ks"])
            like_rows.append(
                {
                    "code": code,
                    "cfg_id": cfg_id,
                    "sigma_a": s[0],
                    "sigma_b": s[1],
                    "sigma_mean": mean,
                    "rel_diff": abs(s[1] - s[0]) / mean,
                    "z": z,
                    "ks_max": worst["ks"],
                    "ks_max_obs": worst["observable"],
                    "ks_over_floor": sum(0 if r["ok"] else 1 for r in kss),
                    "n_obs": len(kss),
                    "n_events": n,
                }
            )
            for row in kss:
                ks_rows.append({"level": "like", "pair": code, "cfg_id": cfg_id, **row})
            del events_a, events_b

    # ---- cross-code: sigma ratio and KS on the pooled runs
    cross_rows: list[dict] = []
    for code_a, code_b in PAIRS:
        common = sorted(set(by_cfg[code_a]) & set(by_cfg[code_b]))
        for cfg_id in common:
            ra, rb = by_cfg[code_a][cfg_id], by_cfg[code_b][cfg_id]
            if len(ra) < 2 or len(rb) < 2:
                continue
            sa = float(np.mean([r["sigma"] for r in ra]))
            sb = float(np.mean([r["sigma"] for r in rb]))
            ratio = sa / sb
            rels = [scatter[code_a], scatter[code_b]]
            z = float("nan")
            if all(np.isfinite(rels)) and min(rels) > 0:
                z = abs(ratio - 1.0) / np.sqrt(sum(r**2 for r in rels) / 2.0)
            events_a = {
                name: np.concatenate([load_events(code_a, r["path"])[name] for r in ra])
                for _, name, *_ in OBSERVABLES
            }
            events_b = {
                name: np.concatenate([load_events(code_b, r["path"])[name] for r in rb])
                for _, name, *_ in OBSERVABLES
            }
            n_a = len(events_a["es"])
            n_b = len(events_b["es"])
            kss = ks_row_set(events_a, events_b, ks_floor(n_a, n_b, c))
            worst = max(kss, key=lambda r: r["ks"])
            cross_rows.append(
                {
                    "pair": f"{code_a} vs {code_b}",
                    "cfg_id": cfg_id,
                    "sigma_a": sa,
                    "sigma_b": sb,
                    "ratio": ratio,
                    "z": z,
                    "ks_max": worst["ks"],
                    "ks_max_obs": worst["observable"],
                    "ks_over_floor": sum(0 if r["ok"] else 1 for r in kss),
                    "n_obs": len(kss),
                    "n_events_a": n_a,
                    "n_events_b": n_b,
                }
            )
            for row in kss:
                ks_rows.append(
                    {"level": "cross", "pair": f"{code_a} vs {code_b}", "cfg_id": cfg_id, **row}
                )
            del events_a, events_b

    # ---- write CSVs
    def write_csv(name: str, rows: list[dict]) -> Path:
        path = out_dir / name
        with open(path, "w", newline="") as fh:
            if rows:
                w = csv.DictWriter(fh, fieldnames=list(rows[0]))
                w.writeheader()
                w.writerows(rows)
        return path

    sigma_rows = [
        {
            "code": code,
            "run_id": run_id,
            "cfg_id": info["cfg_id"],
            "sigma": info["sigma"],
            "trials": info["trials"],
            "n_events": info["n_events"],
            "seed": info["seed"],
        }
        for code in CODES
        for run_id, info in sorted(results[code].items())
    ]
    write_csv("sigmas.csv", sigma_rows)
    write_csv("like_for_like.csv", like_rows)
    write_csv("cross_code.csv", cross_rows)
    write_csv("ks_full.csv", ks_rows)

    # ---- per-observable aggregation for the cross-code pairs
    per_obs: list[dict] = []
    for level in ("like", "cross"):
        codes = sorted({r["pair"] for r in ks_rows if r["level"] == level})
        for pair in codes:
            subset = [r for r in ks_rows if r["level"] == level and r["pair"] == pair]
            for title, *_ in OBSERVABLES:
                obs = [r for r in subset if r["observable"] == title]
                if not obs:
                    continue
                stats = np.array([r["ks"] for r in obs], dtype=float)
                floors = np.array([r["floor"] for r in obs], dtype=float)
                per_obs.append(
                    {
                        "level": level,
                        "pair": pair,
                        "observable": title,
                        "n_checks": len(obs),
                        "median_ks": float(np.nanmedian(stats)),
                        "max_ks": float(np.nanmax(stats)),
                        "frac_within_floor": float(np.mean(stats <= floors)),
                    }
                )
    write_csv("per_observable.csv", per_obs)

    # ---- overall summary
    lines: list[str] = []
    lines.append("Perlmutter scan: statistical verification")
    lines.append(
        f"grid: {len(by_cfg.get('fortran', {}))} configurations "
        f"(beam 2-12 GeV x channel x polarisation x delta-cut x target)"
    )
    lines.append("")
    lines.append("run-to-run scatter (pooled fractional):")
    for code in CODES:
        lines.append(
            f"  {code}: {scatter[code]:.4%}" if np.isfinite(scatter[code]) else f"  {code}: n/a"
        )
    lines.append("")
    lines.append("like-for-like (same code, same configuration, two runs):")
    for code in CODES:
        rows = [r for r in like_rows if r["code"] == code]
        if not rows:
            continue
        zs = [r["z"] for r in rows]
        ks = [r["ks_max"] for r in rows]
        over = sum(r["ks_over_floor"] for r in rows)
        total = sum(r["n_obs"] for r in rows)
        lines.append(
            f"  {code}: {len(rows)} cfgs, sigma |z| max {fmt(zs, '.2f')} "
            f"mean {mean_or_na(zs, '.2f')}; KS max {max(ks):.4f}; "
            f"{over}/{total} observable checks over the noise floor"
        )
    lines.append("")
    lines.append("cross-code agreement per pair:")
    verdict_lines = []
    for pair in [f"{a} vs {b}" for a, b in PAIRS]:
        rows = [r for r in cross_rows if r["pair"] == pair]
        if not rows:
            verdict_lines.append(f"  {pair}: no complete configurations")
            continue
        zs = [r["z"] for r in rows]
        ratios = np.array([r["ratio"] for r in rows], dtype=float)
        ratios = ratios[np.isfinite(ratios)]
        over = sum(r["ks_over_floor"] for r in rows)
        total = sum(r["n_obs"] for r in rows)
        worst = max(rows, key=lambda r: r["ks_max"])
        zvals = np.array(zs, dtype=float)
        zfin = zvals[np.isfinite(zvals)]
        under2 = float(np.mean(zfin < 2)) if zfin.size else float("nan")
        lines.append(
            f"  {pair}: {len(rows)} cfgs, "
            f"sigma ratio median {np.median(ratios):.5f} "
            f"(range {ratios.min():.5f}..{ratios.max():.5f}), "
            f"|z| max {fmt(zs, '.2f')}, "
            f"{under2:.0%} of cfgs under 2 sigma; "
            f"KS: {over}/{total} observable checks over the noise floor"
        )
        verdict_lines.append(
            f"    worst KS {worst['ks_max']:.4f} ({worst['ks_max_obs']}) in {worst['cfg_id']}"
        )
    lines.extend(verdict_lines)
    lines.append("")
    ok_frac = np.mean([r["ok"] for r in ks_rows if r["level"] == "cross"])
    lines.append(
        f"overall: {ok_frac:.1%} of cross-code observable checks within "
        f"the {args.floor_c}*sqrt(1/N1+1/N2) noise floor"
    )

    summary = "\n".join(lines)
    (out_dir / "summary.txt").write_text(summary + "\n")
    print()
    print(summary)
    print(f"\nwritten to {out_dir}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
