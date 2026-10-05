#!/usr/bin/env python3
"""Overlay the Fortran and JAX event distributions and show where they differ.

For a chosen set of configurations this draws, per observable, the two
normalised histograms on top of each other and -- because two overlaid curves
say *that* they differ but not *where* or *by how much* -- a second panel with
the signed bin-by-bin difference and its two-sample error band.

The same defect gate the statistical verification uses is applied here, so a
configuration whose Fortran reference overflowed its 32-bit counters is never
plotted: that would show the port disagreeing with a reference that could not
count its own trials.  Samples are also quantized to the precision the
reference reports (``es16.8``, 8 significant digits) so the plots agree with
``verify_statistics.py``; without that, ``E_s`` at 4.244 GeV is a float32-vs-
float64 artifact and looks like disjoint physics.

Usage::

    python make_plots.py --root $SCRATCH/aao_rad_scan --out-dir plots

Figures written:

``cfg_NNN.png``            per-configuration overlay + difference grid
``sigma_by_energy.png``    cross-section agreement by beam energy
``sigma_by_config.png``    per-configuration cross-section ratio
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))

from compare_distributions import OBSERVABLES, ks_test, load_fortran_ntuple  # noqa: E402
from verify_statistics import (  # noqa: E402
    CODES,
    collect,
    default_root,
    pairs_by_config,
    quantize_sig,
)

# Colours chosen to stay distinguishable in greyscale as well: the two arms
# differ in dash pattern as well as hue, so a black-and-white print still reads.
STYLE = {
    "fortran": {"color": "#1f4e9c", "label": "Fortran", "linestyle": "-", "hatch": None},
    "py_cpu": {"color": "#c2570a", "label": "JAX-CPU", "linestyle": "--", "hatch": "//"},
    "py_gpu": {"color": "#127a52", "label": "JAX-GPU", "linestyle": "-.", "hatch": None},
}


def load_events(code: str, run_dir: Path, observables: list[str]) -> dict[str, np.ndarray]:
    """Observable columns for one run, quantized to the reference's precision.

    Returns an empty mapping when the file is absent or unreadable.  A run counts
    as complete from its ``out.txt`` alone, so an arm can be 'complete' with a
    missing or half-written npz; that should cost that arm's curves, not kill
    the whole figure.
    """
    try:
        if code == "fortran":
            raw = load_fortran_ntuple(run_dir / "aao_rad.ntuple")
        else:
            with np.load(run_dir / "out.npz") as data:
                raw = {name: np.asarray(data[name]) for name in observables if name in data}
    except (OSError, ValueError, EOFError):
        return {}
    return {name: quantize_sig(raw[name]) for name in observables if name in raw}


def two_sample_bins(
    a: np.ndarray, b: np.ndarray, edges: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Normalised histograms and the per-bin difference with its error.

    Returns densities (probability per unit x) to plot, and the two-sample error
    on their difference.  The error is built from the per-bin *probabilities*
    p = n/N, not from the densities: sqrt(p(1-p)/N) is the variance of a
    multinomial bin.  Feeding densities in instead makes p exceed 1, p(1-p) go
    negative, and the band come back NaN.
    """
    na, _ = np.histogram(a, bins=edges)
    nb, _ = np.histogram(b, bins=edges)
    width = np.diff(edges)
    na_sum, nb_sum = max(int(na.sum()), 1), max(int(nb.sum()), 1)
    pa, pb = na / na_sum, nb / nb_sum
    dens_a, dens_b = pa / width, pb / width
    err = (
        np.sqrt(pa * (1.0 - pa) / na_sum) + np.sqrt(pb * (1.0 - pb) / nb_sum)
    ) / width
    centres = 0.5 * (edges[1:] + edges[:-1])
    return centres, dens_a, dens_b - dens_a, err


def bin_edges(a: np.ndarray, b: np.ndarray, target_bins: int = 60) -> np.ndarray:
    """Shared edges over the union of both samples, so neither is clipped."""
    lo = float(min(a.min(), b.min()))
    hi = float(max(a.max(), b.max()))
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        return np.linspace(lo - 0.5, lo + 0.5, target_bins + 1)
    return np.linspace(lo, hi, target_bins + 1)


def plot_config(
    cfg_id: str,
    events: dict[str, dict[str, np.ndarray]],
    observables: list[tuple],
    out_dir: Path,
    ks_floor_c: float,
    per_row: int = 2,
    dpi: int = 110,
    bins: int = 60,
) -> Path | None:
    """Overlay grid for one configuration: histogram on top, difference below.

    ``per_row`` observables are placed side by side.  Stacking all sixteen in one
    column produced a 4400-pixel-tall figure, which is technically a plot and
    practically unreadable.
    """
    codes = [c for c in CODES if c in events]
    per_row = max(1, per_row)
    n_cols = per_row * 2  # histogram | verdict, and difference | spare
    n_rows = 2 * math.ceil(len(observables) / per_row)
    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(5.5 * n_cols / 2, 2.6 * n_rows),
        squeeze=False,
        gridspec_kw={"height_ratios": [3, 1.15] * (n_rows // 2)},
    )

    for i, (title, name, lo, hi, logscale) in enumerate(observables):
        r = 2 * (i // per_row)
        c = 2 * (i % per_row)
        samples = {}
        for code in codes:
            v = events[code].get(name)
            if v is not None and v.size:
                samples[code] = np.asarray(v, float)
        if len(samples) < 2:
            for rr in (r, r + 1):
                for cc in (c, c + 1):
                    axes[rr][cc].axis("off")
            continue

        ref = codes[0]
        edges = bin_edges(samples[ref], samples[codes[-1]], bins)
        centres, fa, diff, err = two_sample_bins(samples[ref], samples[codes[-1]], edges)

        # ---- the distributions on top of each other
        ax = axes[r][c]
        for code, sample in samples.items():
            st = STYLE[code]
            counts, _ = np.histogram(sample, bins=edges)
            ax.step(
                centres,
                counts / counts.sum() / np.diff(edges),
                where="mid",
                color=st["color"],
                linestyle=st["linestyle"],
                linewidth=1.6,
                label=f"{st['label']} (n={sample.size})",
            )
        if logscale:
            ax.set_yscale("log")
        ax.set_ylabel("density")
        ax.legend(fontsize=7, frameon=False, loc="upper right")

        # ---- the statistic, so each panel states its own verdict
        axk = axes[r][c + 1]
        axk.axis("off")
        axk.set_title(title, fontsize=10)
        y = 0.95
        for code in codes[1:]:
            stat, _ = ks_test(samples[ref], samples[code])
            n_a, n_b = samples[ref].size, samples[code].size
            floor = ks_floor_c * np.sqrt(1.0 / n_a + 1.0 / n_b)
            verdict = "within floor" if stat <= floor else "OVER FLOOR"
            block = "\n".join(
                [
                    f"{STYLE[code]['label']} vs {STYLE[ref]['label']}",
                    f"  KS    {stat:.4f}",
                    f"  floor {floor:.4f}",
                    f"  {verdict}",
                ]
            )
            axk.text(
                0.02, y, block, family="monospace", fontsize=8, va="top",
                bbox=dict(boxstyle="round,pad=0.4", facecolor="#f4f4f4", edgecolor="#999"),
            )
            y -= 0.42

        # ---- where they differ, with the two-sample band
        axd = axes[r + 1][c]
        axd.bar(centres, diff, width=np.diff(edges), color="#444", alpha=0.75)
        axd.fill_between(
            centres, -err, err, color="#888", alpha=0.4, linewidth=0,
        )
        axd.axhline(0.0, color="k", linewidth=0.8)
        axd.set_yscale(
            "symlog", linthresh=max(float(np.max(np.abs(diff))) * 0.02, 1e-6)
        )
        axd.set_ylabel(f"{STYLE[codes[-1]]['label']} - {STYLE[ref]['label']}")
        axd.set_xlabel(name)
        axes[r + 1][c + 1].axis("off")

    # Blank off any unused cells in the last block row.  The grid holds
    # per_row * n_block_rows cells, not per_row * n_rows: each block is two rows.
    for i in range(len(observables), per_row * (n_rows // 2)):
        r = 2 * (i // per_row)
        c = 2 * (i % per_row)
        for rr in (r, r + 1):
            for cc in (c, c + 1):
                axes[rr][cc].axis("off")

    fig.suptitle(
        f"{cfg_id}: Fortran vs JAX event distributions "
        f"(both seeds pooled per arm; quantized to 8 significant digits)",
        fontsize=12,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.985))
    path = out_dir / f"{cfg_id}.png"
    fig.savefig(path, dpi=dpi)
    plt.close(fig)
    return path


def plot_sigma_by_energy(
    cross_rows: list[dict], ebeam_by_cfg: dict, out_dir: Path, dpi: int = 110
) -> Path:
    """Cross-section difference per beam energy, against the same-code noise."""
    pairs = sorted({r["pair"] for r in cross_rows})
    fig, ax = plt.subplots(figsize=(8.2, 5.0))
    for pair in pairs:
        rows = [r for r in cross_rows if r["pair"] == pair]
        by_e: dict[float, list[float]] = {}
        for r in rows:
            e = ebeam_by_cfg.get(r["cfg_id"])
            if e is None:
                continue
            by_e.setdefault(float(e), []).append(
                (r["sigma_a"] - r["sigma_b"]) / r["sigma_b"] * 100.0
            )
        if not by_e:
            continue
        es = sorted(by_e)
        med = [float(np.median(by_e[e])) for e in es]
        p90 = [float(np.percentile(np.abs(by_e[e]), 90)) for e in es]
        label = pair.replace(" vs ", " vs\n")
        ax.plot(es, med, marker="o", linewidth=1.6, label=f"{label}  median")
        ax.fill_between(
            es, np.array(med) - np.array(p90), np.array(med) + np.array(p90),
            alpha=0.18, linewidth=0,
        )
        ax.annotate(
            "shaded: median ± |d| p90",
            xy=(0.99, 0.02), xycoords="axes fraction", ha="right", fontsize=8,
            color="#555",
        )
    ax.axhline(0.0, color="k", linewidth=0.8)
    ax.set_xlabel("beam energy (GeV)")
    ax.set_ylabel("cross-section difference (%)")
    ax.set_title(
        "Cross-section agreement by beam energy\n"
        "valid configurations only; a pooled number hides this trend",
        fontsize=11,
    )
    ax.legend(fontsize=8, frameon=False)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    path = out_dir / "sigma_by_energy.png"
    fig.savefig(path, dpi=dpi)
    plt.close(fig)
    return path


def plot_sigma_by_config(
    cross_rows: list[dict], ebeam_by_cfg: dict, out_dir: Path, dpi: int = 110
) -> Path:
    """Every configuration's ratio, coloured by beam energy, so outliers are visible."""
    pairs = sorted({r["pair"] for r in cross_rows})
    fig, axes = plt.subplots(
        1, len(pairs), figsize=(5.0 * len(pairs), 4.4), squeeze=False
    )
    energies = sorted({float(v) for v in ebeam_by_cfg.values()})
    colours = plt.get_cmap("viridis")(
        np.linspace(0.08, 0.92, max(len(energies), 1))
    )
    for ax, pair in zip(axes[0], pairs):
        rows = [r for r in cross_rows if r["pair"] == pair]
        for colour, e in zip(colours, energies):
            sub = [
                r for r in rows if ebeam_by_cfg.get(r["cfg_id"]) and float(ebeam_by_cfg[r["cfg_id"]]) == e
            ]
            if not sub:
                continue
            ax.scatter(
                [r["cfg_id"] for r in sub],
                [(r["sigma_a"] - r["sigma_b"]) / r["sigma_b"] * 100.0 for r in sub],
                s=28, color=colour, label=f"{e:g} GeV", zorder=3,
            )
        ax.axhline(0.0, color="k", linewidth=0.8)
        ax.axhspan(-0.5, 0.5, color="#888", alpha=0.12, zorder=0)
        ax.set_title(pair, fontsize=10)
        ax.set_ylabel("difference (%)")
        ax.tick_params(axis="x", rotation=90, labelsize=6)
        ax.grid(alpha=0.25, axis="y")
    axes[0][0].legend(fontsize=7, frameon=False, ncol=2)
    fig.suptitle(
        "Per-configuration cross-section difference (shaded band: ±0.5%)",
        fontsize=11,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    path = out_dir / "sigma_by_config.png"
    fig.savefig(path, dpi=dpi)
    plt.close(fig)
    return path


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--root", type=Path, default=None)
    p.add_argument("--out-dir", type=Path, default=None)
    p.add_argument("--observables", default="", help="comma-separated subset of observable columns")
    p.add_argument("--per-energy", type=int, default=2, help="configurations to plot per beam energy")
    p.add_argument("--bins", type=int, default=60, help="bins per observable histogram")
    p.add_argument(
        "--per-row", type=int, default=4,
        help=(
            "observables placed side by side.  4 is the default because it is "
            "the smallest value that keeps a full 16-observable figure about "
            "as wide as it is tall; 2 leaves it 3.8x taller than wide, which is "
            "still a scrolling exercise"
        ),
    )
    p.add_argument("--dpi", type=int, default=110)
    p.add_argument("--floor-c", type=float, default=1.36)
    p.add_argument("--all-valid", action="store_true", help="plot every valid configuration")
    args = p.parse_args()

    root = args.root or default_root()
    out_dir = args.out_dir or root / "plots"
    out_dir.mkdir(parents=True, exist_ok=True)

    wanted = {s.strip() for s in args.observables.split(",") if s.strip()}
    observables = [
        o for o in OBSERVABLES if not wanted or o[1] in wanted or o[0] in wanted
    ]
    if not observables:
        print(f"no observables matched {sorted(wanted)}", file=sys.stderr)
        return 2

    ebeam_by_cfg: dict[str, str] = {}
    with open(root / "grid" / "manifest.csv", newline="") as fh:
        for row in csv.DictReader(fh):
            if row.get("ebeam"):
                ebeam_by_cfg[row["cfg_id"]] = row["ebeam"]

    results = collect(root)
    by_cfg = pairs_by_config(results)

    # A configuration is plottable when the *reference* is sound -- both seeds
    # clean -- and at least one port arm has both seeds to compare against.
    # Requiring every arm would throw away configurations that are perfectly
    # plottable Fortran-vs-GPU just because an unrelated arm is still running,
    # which is the normal state of a scan in progress.
    valid = [
        cfg
        for cfg in sorted(by_cfg["fortran"])
        if len(by_cfg["fortran"][cfg]) >= 2
        and any(len(by_cfg[c].get(cfg, [])) >= 2 for c in ("py_cpu", "py_gpu"))
    ]
    if not valid:
        print("no configuration has a clean two-seed reference to plot against", file=sys.stderr)
        return 1

    if args.all_valid:
        chosen = valid
    else:
        # Take the configurations with the most events at each energy, so the
        # plotted histograms are the best-resolved ones available.
        by_e: dict[float, list[str]] = {}
        for cfg in valid:
            by_e.setdefault(float(ebeam_by_cfg.get(cfg, 0.0)), []).append(cfg)
        chosen = []
        for e in sorted(by_e):
            ranked = sorted(by_e[e], key=lambda c: -by_cfg["fortran"][c][0]["n_events"])
            chosen.extend(ranked[: args.per_energy])

    print(f"root       {root}")
    print(f"out-dir    {out_dir}")
    # Count against the manifest, not the filtered set: deriving the denominator
    # from by_cfg would print "2 of 2" and hide that a configuration was dropped.
    n_total = len({cfg for cfg in ebeam_by_cfg})
    dropped = n_total - len(by_cfg["fortran"])
    print(
        f"valid cfgs {len(valid)} of {n_total} in the grid "
        f"({len(by_cfg['fortran'])} with a usable reference"
        + (f", {dropped} with no clean reference" if dropped else "")
        + f"); plotting {len(chosen)}"
    )
    for cfg in sorted({r["cfg_id"] for r in results["fortran"].values()}):
        if cfg in valid:
            continue
        defects = sorted(
            {
                d
                for r in results["fortran"].values()
                if r["cfg_id"] == cfg
                for d in (r.get("defects") or ("reference unusable",))
            }
        )
        short = [c for c in ("py_cpu", "py_gpu") if len(by_cfg[c].get(cfg, [])) < 2]
        parts = [f"fortran {', '.join(defects)}"] if defects else ["fortran not usable"]
        if short:
            parts.append(f"{'/'.join(short)} incomplete")
        print(f"  skipped {cfg}: {'; '.join(parts)}")

    made: list[Path] = []
    for i, cfg in enumerate(chosen, 1):
        arms = [c for c in CODES if len(by_cfg[c].get(cfg, [])) >= 2]
        events: dict[str, dict[str, np.ndarray]] = {}
        for code in arms:
            # Pool both seeds so each histogram has the full event count.
            pooled: dict[str, list[np.ndarray]] = {}
            for run in by_cfg[code][cfg][:2]:
                for name, v in load_events(
                    code, run["path"], [o[1] for o in observables]
                ).items():
                    pooled.setdefault(name, []).append(v)
            events[code] = {k: np.concatenate(v) for k, v in pooled.items()}
        path = plot_config(
            cfg, events, observables, out_dir, args.floor_c,
            per_row=args.per_row, dpi=args.dpi, bins=args.bins,
        )
        if path:
            made.append(path)
            print(f"  [{i}/{len(chosen)}] {path.name}  arms={','.join(arms)}")

    cross_rows = _cross_rows(results, by_cfg)
    made.append(plot_sigma_by_energy(cross_rows, ebeam_by_cfg, out_dir, args.dpi))
    made.append(plot_sigma_by_config(cross_rows, ebeam_by_cfg, out_dir, args.dpi))
    print(f"wrote {len(made)} figures to {out_dir}")
    return 0


def _cross_rows(results: dict, by_cfg: dict) -> list[dict]:
    """The same two-seed cross-section comparison verify_statistics reports."""
    rows = []
    for code_a in CODES:
        for code_b in CODES:
            if code_a == code_b:
                continue
            for cfg, ra in by_cfg[code_a].items():
                rb = by_cfg[code_b].get(cfg)
                if not ra or not rb or len(ra) < 2 or len(rb) < 2:
                    continue
                sa = float(np.mean([r["sigma"] for r in ra]))
                sb = float(np.mean([r["sigma"] for r in rb]))
                rows.append(
                    {
                        "pair": f"{code_a} vs {code_b}",
                        "cfg_id": cfg,
                        "sigma_a": sa,
                        "sigma_b": sb,
                    }
                )
    return rows


if __name__ == "__main__":
    raise SystemExit(main())
