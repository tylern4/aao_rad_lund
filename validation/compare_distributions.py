#!/usr/bin/env python3
"""Statistical comparison between the Fortran and Python event generators.

Two independent levels of evidence are produced.

**Pointwise (exact).**  ``dump_dsigma.f90`` / the instrumented ``aao_rad.f90``
write the internal CGLN and helicity amplitudes, so the whole amplitude chain
can be compared value by value.  See ``compare_dsigma_to_fortran.py`` and
``compare_sigma_points.py``; this script does not repeat them.

**Distributional (at scale).**  This script compares the *event* distributions.
Event-by-event equality is impossible -- the Fortran seeds ``myran`` from
``unixtime`` and there is no way to make the two generators draw the same
sequence -- so the criterion is that the two samples describe the same
distribution.  Each observable is summarised by

* a two-sample Kolmogorov-Smirnov statistic (distribution-free, sensitive to
  the whole shape, not just the mean), and
* a per-bin ratio histogram, plotted with Poisson error bars.

Inputs
------
``--fortran`` is the n-tuple written by this repository's instrumented
``aao_rad.f90`` (``aao_rad.ntuple``), *not* the LUND file.  The LUND writer in
the original only ever emits two tracks -- ``aao_rad.f90:236`` forces
``npart = 4`` down to 2, and the four-track blocks at lines 1091-1104 are
commented out -- so tracks 3 and 4, and hence the proton momentum,
``cos(theta*)``, ``phi*``, ``E_gamma`` and ``mm^2``, simply do not exist in it.
The n-tuple has all 32 columns and allows all of them to be checked.

Run card
--------
The run card is passed in explicitly (``--run-card``) rather than defaulted, so
the two generators are guaranteed to cover the same kinematic window.  A
mismatched window shows up immediately as a large KS statistic on ``Q^2`` and
``E'``, which is a useless diagnostic.

Usage
-----
    aao_rad_validate_run <run card>          # write a Fortran n-tuple
    python validation/compare_distributions.py \\
        --fortran /tmp/frun/aao_rad.ntuple --n-python 200000
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections.abc import Sequence
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
VALIDATION = ROOT / "validation"
PLOTS = VALIDATION / "plots"
PARMS = ROOT / "parms"
sys.path.insert(0, str(ROOT / "src"))

# ---------------------------------------------------------------------------
# n-tuple layout
# ---------------------------------------------------------------------------
#: The 32 n-tuple columns written by aao_rad.f90 (1-based ntp index -> name).
#: The port's EVENT_COLUMNS uses the same names in the same order, so the two
#: records can be compared by position.
NTP_COLUMNS: tuple[str, ...] = (
    "es", "ep", "theta", "w", "w_real",
    "ppx", "ppy", "ppz", "eprot",
    "ppix", "ppiy", "ppiz", "epi",
    "csthcm", "phicm", "mm2",
    "eg", "cstk", "phik",
    "qx", "qz", "q0", "csthe",
    "egx", "egy", "egz",
    "vx", "vy", "vz",
    "q2", "e_hel", "asym_p",
)

#: Extra columns appended by the validation instrumentation in aao_rad.f90.
NTP_EXTRA: tuple[str, ...] = (
    "weight", "sigr_max", "ntries", "mcall",
    "sigr_raw", "sigma_ek", "cstk_d", "csthcm_d", "phicm_d", "phik_d",
    "ehel_d", "es_d", "ep_d", "q2_d", "w_d",
)

#: (title, n-tuple name, lo, hi, log-y?)
#:
#: The ranges are deliberately generous -- the goal is to see the whole
#: distribution, not to crop it to where most events land.  ``vertex_*`` is
#: excluded: it is a uniform box and tests nothing.  ``asym_p`` and ``phik`` are
#: in degrees in the Fortran record and the port matches that.
#:
#: Each range must contain the data.  The reported ``outside`` column is what
#: checks that, and it is worth re-checking after changing a range: ``phi_k``
#: once read as ``[-30, 390]`` when the variable actually spans ``[-180, 180]``,
#: which quietly dropped 31% of the events from the histogram and the ratio
#: panel.  KS is unaffected (it runs on the raw samples) but the per-bin
#: comparison and the plot both were.
#: ``(title, n-tuple column, lo, hi, log y)`` for every observable compared.
#:
#: The ranges must cover the *whole* support of both samples.  Anything outside
#: is counted in the ``outside`` column and dropped from the histogram, so a
#: range that is too narrow does not merely clip a plot -- it quietly changes
#: the comparison.  ``phi_k`` was originally plotted over [-30, 360] when it
#: spans [-180, 180] and threw away 31.5% of the events.
#:
#: ``mm^2`` and ``q0`` are a little wider than their physical windows on
#: purpose.  The missing-mass cut is applied to the *pre-exit* electron energy
#: (aao_rad.f90:905-916) but the recorded column comes from the final state,
#: which uses the post-exit energy, so a handful of events land just outside
#: the cut window.  Two in 200k for ``mm^2``, one for ``q0``.
OBSERVABLES: tuple[tuple[str, str, float, float, bool], ...] = (
    ("E_s (GeV)", "es", 3.9, 4.25, False),
    ("E' (GeV)", "ep", 1.55, 2.95, False),
    ("Q^2 (GeV^2)", "q2", 0.0, 2.0, True),
    ("theta_e (deg)", "theta", 4.0, 60.0, True),
    ("W (GeV)", "w", 1.0, 2.4, True),
    ("W_real (GeV)", "w_real", 1.0, 2.5, False),
    ("E_pion (GeV)", "epi", 0.0, 3.6, True),
    ("cos(theta*)", "csthcm", -1.0, 1.0, True),
    ("phi* (deg)", "phicm", 0.0, 360.0, False),
    ("E_gamma (GeV)", "eg", 0.0, 0.6, True),
    ("mm^2 (GeV^2)", "mm2", 0.65, 1.15, True),
    ("cos(theta_k)", "cstk", -1.0, 1.0, True),
    ("phi_k (deg)", "phik", -180.0, 180.0, False),
    ("q0 (GeV)", "q0", 1.0, 2.75, True),
    ("cos(theta_e)", "csthe", 0.4, 1.0, True),
    ("asym_p", "asym_p", -0.9, 0.9, True),
)


def load_fortran_ntuple(path: Path, max_rows: int | None = None) -> dict[str, np.ndarray]:
    """Read the Fortran validation n-tuple into named columns."""
    raw = np.loadtxt(path, ndmin=2)
    n_expected = len(NTP_COLUMNS) + len(NTP_EXTRA)
    if raw.shape[1] < n_expected:
        raise SystemExit(
            f"{path} has {raw.shape[1]} columns, expected {n_expected} "
            f"({len(NTP_COLUMNS)} n-tuple + {len(NTP_EXTRA)} validation). "
            "Rebuild aao_rad.f90 with the instrumentation and re-run it."
        )
    if max_rows is not None:
        raw = raw[:max_rows]
    out = {name: raw[:, i].astype(np.float64) for i, name in enumerate(NTP_COLUMNS)}
    for j, name in enumerate(NTP_EXTRA):
        out[name] = raw[:, len(NTP_COLUMNS) + j].astype(np.float64)
    return out


def run_card_from_file(path: Path) -> dict[str, float]:
    """Read the numbers aao_rad.f90 reads, in the order it reads them.

    The reads are one value per line except the three two-value lines (the
    integration regions, the Q^2 window and the E' window).  gfortran's
    list-directed reads reject ``!`` comments, so the card is bare numbers.
    """
    values: list[float] = []
    for line in path.read_text().split("\n"):
        line = line.split("!")[0].strip()
        if not line:
            continue
        for tok in line.replace(",", " ").split():
            values.append(float(tok))

    it = iter(values)
    card: dict[str, float] = {}
    card["theory"] = next(it)
    card["flag_ehel"] = next(it)
    card["reg"] = [next(it) for _ in range(4)]
    card["npart"] = next(it)
    card["epirea"] = next(it)
    card["mm_cut"] = next(it)
    card["t_targ"] = next(it)
    card["r_targ"] = next(it)
    card["vertex_x"] = next(it)
    card["vertex_y"] = next(it)
    card["vz"] = next(it)
    card["beam_energy"] = next(it)
    card["q2_min"], card["q2_max"] = next(it), next(it)
    card["ep_min"], card["ep_max"] = next(it), next(it)
    card["delta"] = next(it)
    card["nmax"] = next(it)
    card["fmcall"] = next(it)
    if card["fmcall"] == 0.0:
        # fmcall == 0 makes aao_rad.f90 read an explicit sigr_max ... and then
        # multiply it by fmcall, i.e. zero it (aao_rad.f90:513).  A zero ceiling
        # makes sig_ratio infinite and the acceptance test meaningless, so the
        # reference run must use fmcall != 0.
        card["sigr_max_input"] = next(it)
        raise SystemExit(
            f"{path} sets fmcall = 0, which zeroes sigr_max at aao_rad.f90:513 "
            "and breaks the importance-sampling test.  Use fmcall = 1 to have "
            "the generator estimate the ceiling from its own 10k-point scan."
        )
    return card


def generate_python(card: dict[str, float], n_events: int, seed: int,
                    batch_size: int) -> tuple[dict[str, np.ndarray], object]:
    """Generate with the port using exactly the Fortran run card."""
    from aao_rad import EventGenerator, GeneratorConfig, build_grid

    is_pi0 = int(card["epirea"]) == 1
    cfg = GeneratorConfig(
        channel=1 if is_pi0 else 3,
        beam_energy=card["beam_energy"],
        q2_min=card["q2_min"],
        q2_max=card["q2_max"],
        ep_min=card["ep_min"],
        ep_max=card["ep_max"],
        min_photon_energy=card["delta"],
        missing_mass_cut=card["mm_cut"],
        k_exp=5.0,
        target_length_cm=card["t_targ"],
        target_radius_cm=card["r_targ"],
        beam_x_cm=card["vertex_x"],
        beam_y_cm=card["vertex_y"],
        beam_z_cm=card["vz"],
        # flag_ehel == 1 draws the helicity per trial with get_spin
        # (aao_rad.f90:521); the port's polarised_beam does the same.
        polarized_beam=bool(card["flag_ehel"]),
        regions=tuple(card["reg"]),
        n_events=n_events,
        batch_size=batch_size,
        seed=seed,
        # The Fortran saturates the table for out-of-range (W, Q^2) lookups
        # (multipole_amps.f90:16-28).  Reproduce that, or a large part of the
        # sampled window lands in different kinematics.
        w_max="clamp",
        ek_sampling="fortran",
    )
    grid = build_grid(cfg.channel, scheme=cfg.interp_scheme, parms_dir=str(PARMS))
    events, stats = EventGenerator(grid).generate(cfg)
    return {name: events[name].astype(np.float64) for name in events.dtype.names}, stats


def ks_test(a: np.ndarray, b: np.ndarray) -> tuple[float, float]:
    """Two-sample KS statistic and p-value, with a finite-value guard."""
    from scipy import stats

    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if a.size < 2 or b.size < 2:
        return float("nan"), float("nan")
    res = stats.ks_2samp(a, b)
    return float(res.statistic), float(res.pvalue)


def per_bin(py: np.ndarray, fr: np.ndarray, edges: np.ndarray) -> list[dict]:
    """Bin-by-bin comparison of two samples.

    A KS statistic is a worst-case over the whole distribution and says nothing
    about *where* two samples differ, and it is insensitive to a difference that
    cancels between bins.  So also compare bin by bin, which is what a
    histogram overlay is claiming to show:

    ``ratio``
        ``f_py / f_fr`` with the Poisson error a ratio of two counts carries,
        ``sigma(r)/r = sqrt(1/h_py + 1/h_fr)``.  Only defined where both bins
        are populated; the empty-reference bins are the informative ones and are
        reported with the count instead.
    ``pull``
        ``(f_py - f_fr) / sqrt(f_fr (1 - f_fr) (1/N_py + 1/N_fr))``.  This is the
        signed, per-bin version of "how many sigma apart are these two
        histograms", and unlike the KS it flags a difference that would not
        change the KS at all.  The reference term dominates ``1/N_fr`` whenever
        the reference sample is the smaller one, which is the usual case here.

    Fractions, not counts, so a bin is comparable across samples of different
    size.  Events outside ``edges`` are not counted here -- the caller reports
    them, because a plotting range that silently crops a third of the sample
    makes the comparison look better than the truth.  The ranges in
    :data:`OBSERVABLES` are checked against that column, which is how
    ``phi_k``'s ``[-30, 390]`` was caught hiding 31% of the events.
    """
    py = np.asarray(py, dtype=float)
    fr = np.asarray(fr, dtype=float)
    py = py[np.isfinite(py)]
    fr = fr[np.isfinite(fr)]
    n_py, n_fr = py.size, fr.size
    h_py, _ = np.histogram(py, bins=edges)
    h_fr, _ = np.histogram(fr, bins=edges)
    f_py = h_py / max(n_py, 1)
    f_fr = h_fr / max(n_fr, 1)
    # Poisson variance of the reference fraction, floored at zero for empty bins.
    var = np.where(f_fr > 0.0, f_fr * (1.0 - f_fr) * (1.0 / max(n_py, 1) + 1.0 / max(n_fr, 1)), 0.0)
    pull = np.where(var > 0.0, (f_py - f_fr) / np.sqrt(var), 0.0)

    rows = []
    for i in range(len(edges) - 1):
        both = h_py[i] > 0 and h_fr[i] > 0
        ratio = f_py[i] / f_fr[i] if both else float("nan")
        rel_err = (
            np.sqrt(1.0 / h_py[i] + 1.0 / h_fr[i]) if both else float("nan")
        )
        rows.append({
            "bin": i,
            "lo": float(edges[i]),
            "hi": float(edges[i + 1]),
            "count_py": int(h_py[i]),
            "count_fr": int(h_fr[i]),
            "frac_py": float(f_py[i]),
            "frac_fr": float(f_fr[i]),
            "ratio": float(ratio),
            "ratio_rel_err": float(rel_err),
            "pull": float(pull[i]),
        })
    return rows


def bin_summary(rows: list[dict], n_py: int, n_fr: int,
                outside: tuple[int, int]) -> dict:
    """Condense :func:`per_bin` output to the numbers worth printing."""
    pulled = np.array([r["pull"] for r in rows])
    counted = [r for r in rows if r["count_py"] and r["count_fr"]]
    ratios = np.array([r["ratio"] for r in counted])
    worst = int(np.argmax(np.abs(pulled))) if pulled.size else -1
    empty_ref = sum(1 for r in rows if r["count_fr"] == 0 and r["count_py"] > 0)
    return {
        "n_bins": len(rows),
        "n_bins_over_3sigma": int((np.abs(pulled) > 3.0).sum()),
        "max_abs_pull": float(np.abs(pulled).max()) if pulled.size else float("nan"),
        "worst_bin": worst,
        "median_ratio": float(np.median(ratios)) if ratios.size else float("nan"),
        "ratio_p05": float(np.percentile(ratios, 5)) if ratios.size else float("nan"),
        "ratio_p95": float(np.percentile(ratios, 95)) if ratios.size else float("nan"),
        "n_bins_ref_empty": empty_ref,
        "outside_py": outside[0],
        "outside_fr": outside[1],
        "n_py": n_py,
        "n_fr": n_fr,
    }


def make_figure(obs: Sequence[tuple], py: dict, fr: dict, summary: list[dict],
                bin_rows: dict, n_bins: int, title: str, out: Path) -> None:
    """One page per observable: overlaid histograms above, ratio below."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ncol = 4
    nrow = (len(obs) + ncol - 1) // ncol
    fig, axes = plt.subplots(nrow, 2 * ncol, figsize=(4.2 * ncol, 3.4 * nrow),
                             gridspec_kw={"width_ratios": [2, 1] * ncol})
    axes = axes.reshape(nrow, 2 * ncol)
    hist_axes = axes[:, 0::2].ravel()
    ratio_axes = axes[:, 1::2].ravel()

    by_title = {s["title"]: s for s in summary}
    for ax, rax, (t, _name, lo, hi, logy) in zip(hist_axes, ratio_axes, obs,
                                                 strict=False):
        edges = np.linspace(lo, hi, n_bins + 1)
        rows = bin_rows[t]
        ctr = 0.5 * (edges[:-1] + edges[1:])
        w = 0.4 * (edges[1] - edges[0])
        fp = np.array([r["frac_fr"] for r in rows])
        fq = np.array([r["frac_py"] for r in rows])

        ax.bar(ctr - w / 2, fp, width=w, color="0.55", label="Fortran")
        ax.bar(ctr + w / 2, fq, width=w, color="C0", alpha=0.7, label="Python")
        s = by_title.get(t)
        if s is not None:
            p_text = "<1e-4" if s["p"] < 1e-4 else f"{s['p']:.3g}"
            ax.text(0.03, 0.96, f"KS D={s['D']:.4f}\np={p_text}",
                    transform=ax.transAxes, va="top", ha="left", fontsize=8)
        ax.set_title(t, fontsize=10)
        if logy:
            ax.set_yscale("log")
        ax.set_xlim(lo, hi)
        ax.grid(alpha=0.25)
        ax.tick_params(labelsize=8)

        # Ratio panel: the per-bin question the KS cannot answer.
        both = (fq > 0) & (fp > 0)
        rax.axhline(1.0, color="0.4", lw=0.8, ls="--")
        if both.any():
            ratio = fq[both] / fp[both]
            err = ratio * np.sqrt(1.0 / fq[both] + 1.0 / fp[both])
            rax.errorbar(ctr[both], ratio, yerr=err, fmt="o", ms=2.2, lw=0.7,
                         color="C0", elinewidth=0.7, capsize=0)
        # Reference-empty bins: Python has events where Fortran has none.  Mark
        # them at the top of the axis rather than dropping them silently.
        empty = (fq > 0) & (fp == 0)
        if empty.any():
            rax.plot(ctr[empty], np.full(empty.sum(), 4.0), "^", ms=3,
                     color="C3", clip_on=False)
            rax.text(0.5, 0.92, f"{int(empty.sum())} ref-empty bin(s)",
                     transform=rax.transAxes, ha="center", fontsize=7,
                     color="C3")
        rax.set_yscale("log")
        rax.set_ylim(0.2, 8.0)
        rax.set_xlim(lo, hi)
        rax.grid(alpha=0.25)
        rax.tick_params(labelsize=8)

    for ax in list(hist_axes[len(obs):]) + list(ratio_axes[len(obs):]):
        ax.set_visible(False)

    handles = hist_axes[0].get_legend_handles_labels()[0]
    fig.legend(handles, ["Fortran", "Python"], loc="upper right", ncols=2)
    fig.suptitle(title, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.985))
    fig.savefig(out, dpi=140)
    plt.close(fig)


def print_bin_overview(rows: list[dict]) -> None:
    """One line per observable: how many bins disagree, and by how much."""
    print()
    print("per-bin agreement")
    print(f"{'observable':>18} {'>3sigma':>8} {'max|z|':>7} {'ratio p05':>10} "
          f"{'p50':>7} {'p95':>7} {'ref-empty':>10} {'outside':>9}")
    print("-" * 92)
    for s in rows:
        print(f"{s['title']:>18} {s['n_bins_over_3sigma']:4d}/{s['n_bins']:<4d} "
              f"{s['max_abs_pull']:7.2f} {s['ratio_p05']:10.3f} "
              f"{s['median_ratio']:7.3f} {s['ratio_p95']:7.3f} "
              f"{s['n_bins_ref_empty']:10d} "
              f"{s['outside_py']:4d}/{s['outside_fr']:<4d}")
    print("\n  >3sigma   bins whose pull exceeds 3, out of the bin count")
    print("  max|z|    largest absolute pull over all bins")
    print("  ratio     python fraction / fortran fraction, at the 5th/50th/95th")
    print("             percentile of the bins populated on both sides")
    print("  ref-empty bins with Python events and no Fortran event")
    print("  outside   events outside the plotted range, python/fortran")


def print_per_bin(obs: Sequence[tuple], bin_rows: dict,
                  summaries: list[dict]) -> None:
    """Print the full bin-by-bin table for every observable."""
    by_title = {s["title"]: s for s in summaries}
    for t, _name, _lo, _hi, _logy in obs:
        s = by_title[t]
        print()
        print("=" * 78)
        print(f"{t}   python {s['n_py']:,} events, fortran {s['n_fr']:,}")
        print(f"{'bin':>4} {'range':>18} {'n_py':>8} {'n_fr':>7} "
              f"{'frac_py':>11} {'frac_fr':>11} {'ratio':>8} {'+/-':>7} {'z':>7}")
        print("-" * 78)
        for r in bin_rows[t]:
            rng = f"{r['lo']:.4g}..{r['hi']:.4g}"
            ratio = "      -" if np.isnan(r["ratio"]) else f"{r['ratio']:8.3f}"
            err = "      -" if np.isnan(r["ratio_rel_err"]) else f"{r['ratio_rel_err']:7.3f}"
            z = f"{r['pull']:+7.2f}"
            flag = "  <<<" if abs(r["pull"]) > 3.0 else ""
            print(f"{r['bin']:4d} {rng:>18} {r['count_py']:8d} {r['count_fr']:7d} "
                  f"{r['frac_py']:11.3e} {r['frac_fr']:11.3e} {ratio} {err} {z}{flag}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fortran", type=Path, default=Path("/tmp/frun/aao_rad.ntuple"),
                    help="Fortran validation n-tuple (aao_rad.ntuple)")
    ap.add_argument("--run-card", type=Path, default=Path("/tmp/frun/run_card.txt"),
                    help="the run card aao_rad.f90 was given")
    ap.add_argument("--n-python", type=int, default=200_000)
    ap.add_argument("--max-fortran", type=int, default=None)
    ap.add_argument("--seed", type=int, default=1234)
    ap.add_argument("--batch-size", type=int, default=1 << 17)
    ap.add_argument("--bins", type=int, default=40)
    ap.add_argument("--per-bin", action="store_true",
                    help="print the full bin-by-bin table for every observable")
    ap.add_argument("--out", type=Path, default=PLOTS / "distributions.png")
    ap.add_argument("--csv", type=Path, default=VALIDATION / "distributions.csv")
    ap.add_argument("--bins-csv", type=Path,
                    default=VALIDATION / "distributions_bins.csv")
    ap.add_argument("--port-npz", type=Path, default=None,
                    help="cache the generated port sample here (.npz). On a "
                         "second run the sample is loaded instead of "
                         "regenerated -- the port draw is the expensive part "
                         "and, with a fixed --seed, identical every time. "
                         "The cache records the card, seed, event count and "
                         "batch size, and refuses to load a sample generated "
                         "under different ones.")
    ap.add_argument("--fortran-sigma", type=float, default=None,
                    help="the Fortran run's integrated cross section in "
                         "micro-barns -- the second number on the last "
                         "'Integrated cross section' line of its out.txt, or "
                         "the pooled value from reference_fleet.py. Prints "
                         "the port/Fortran sigma ratio after the tables.")
    args = ap.parse_args()

    import matplotlib

    matplotlib.use("Agg")

    if not args.fortran.exists():
        raise SystemExit(f"{args.fortran} not found -- run aao_rad with the "
                         "instrumented source first (see the Validation section of README.md)")
    if not args.run_card.exists():
        raise SystemExit(f"{args.run_card} not found -- the run card is needed so "
                         "both generators cover the same window")

    card = run_card_from_file(args.run_card)
    fr = load_fortran_ntuple(args.fortran, args.max_fortran)
    print(f"Fortran: {len(fr['es'])} events from {args.fortran}")
    print(f"  run card: beam {card['beam_energy']} GeV, Q^2 "
          f"[{card['q2_min']}, {card['q2_max']}], E' [{card['ep_min']}, {card['ep_max']}], "
          f"epirea={int(card['epirea'])}, delta={card['delta']}, mm_cut={card['mm_cut']}")

    print(f"\nGenerating {args.n_python} Python events on the same card ...")
    cache_meta = {
        "beam_energy": card["beam_energy"], "q2_min": card["q2_min"],
        "q2_max": card["q2_max"], "ep_min": card["ep_min"],
        "ep_max": card["ep_max"], "delta": card["delta"],
        "mm_cut": card["mm_cut"], "epirea": card["epirea"],
        "n_events": args.n_python, "seed": args.seed,
        "batch_size": args.batch_size,
    }
    py = stats = None
    if args.port_npz is not None and args.port_npz.exists():
        cached = np.load(args.port_npz, allow_pickle=False)
        if all(f"_{k}" in cached.files and float(cached[f"_{k}"]) == float(v)
               for k, v in cache_meta.items()):
            py = {k: cached[k] for k in cached.files if not k.startswith("_")}
            print(f"  loaded cached sample from {args.port_npz} "
                  f"(card/seed match)")
        else:
            raise SystemExit(
                f"{args.port_npz} holds a sample generated under a different "
                "card/seed/count -- delete it or pass a fresh --port-npz path")
    if py is None:
        py, stats = generate_python(card, args.n_python, args.seed, args.batch_size)
        if args.port_npz is not None:
            args.port_npz.parent.mkdir(parents=True, exist_ok=True)
            np.savez_compressed(args.port_npz, **py, **{f"_{k}": v for k, v in cache_meta.items()})
            print(f"  cached sample to {args.port_npz}")
    if stats is not None:
        print(f"  {stats.summary()}")
    else:
        print(f"  {len(next(iter(py.values())))} events, "
              f"sampler stats unavailable from cache")

    summary: list[dict] = []
    bin_rows: dict[str, list[dict]] = {}
    bin_summary_rows: list[dict] = []
    print()
    print(f"{'observable':>18} {'KS D':>9} {'p':>10} {'mean py':>11} {'mean for':>11}")
    print("-" * 64)
    for t, name, lo, hi, _logy in OBSERVABLES:
        a, b = py[name], fr[name]
        d, p = ks_test(a, b)
        finite_a, finite_b = a[np.isfinite(a)], b[np.isfinite(b)]
        summary.append({
            "title": t, "name": name, "D": d, "p": p,
            "mean_py": float(finite_a.mean()) if finite_a.size else float("nan"),
            "mean_fr": float(finite_b.mean()) if finite_b.size else float("nan"),
            "median_py": float(np.median(finite_a)) if finite_a.size else float("nan"),
            "median_fr": float(np.median(finite_b)) if finite_b.size else float("nan"),
            "n_py": int(finite_a.size), "n_fr": int(finite_b.size),
        })
        print(f"{t:>18} {d:9.4f} {p:10.2e} {summary[-1]['mean_py']:11.5g} "
              f"{summary[-1]['mean_fr']:11.5g}")

        rows = per_bin(a, b, np.linspace(lo, hi, args.bins + 1))
        bin_rows[t] = rows
        outside = (
            int((finite_a < lo).sum() + (finite_a >= hi).sum()),
            int((finite_b < lo).sum() + (finite_b >= hi).sum()),
        )
        bs = bin_summary(rows, finite_a.size, finite_b.size, outside)
        bs["title"] = t
        bs["name"] = name
        bin_summary_rows.append(bs)

    if args.per_bin:
        print_per_bin(OBSERVABLES, bin_rows, bin_summary_rows)
    print_bin_overview(bin_summary_rows)

    # sigma is the strongest single-number check available.
    sigma_f = args.fortran_sigma if args.fortran_sigma is not None else card.get("sigma_fortran")
    print()
    if stats is not None:
        print(f"port sigma (mean-weight) : {stats.sigma_mc:.6g} micro-barn")
        if sigma_f:
            print(f"fortran sigma            : {sigma_f:.6g} micro-barn")
            print(f"ratio                    : {stats.sigma_mc / sigma_f:.4f}")
    else:
        print("port sigma (mean-weight) : n/a (sample loaded from cache; "
              "run without --port-npz for the trial-level sigma)")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    make_figure(OBSERVABLES, py, fr, summary, bin_rows, args.bins,
                f"Fortran aao_rad vs JAX port   ({len(fr['es'])} vs {len(py['es'])} events)",
                args.out)
    print(f"\nwrote {args.out}")

    with args.csv.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(summary[0]))
        w.writeheader()
        w.writerows(summary)
    print(f"wrote {args.csv}")

    # The per-bin rows go in their own file so the summary stays one row per
    # observable; both are written for every run, not only --per-bin, because
    # the numbers are what anyone actually looks at afterwards.
    with args.bins_csv.open("w", newline="") as fh:
        # every observable produces the same per-bin columns, so the first
        # non-empty bin list defines the header for all of them
        first = next(r for _t, *_rest in OBSERVABLES for r in bin_rows[_t])
        w = csv.DictWriter(fh, fieldnames=["observable", *first])
        w.writeheader()
        for t, _name, _lo, _hi, _logy in OBSERVABLES:
            for r in bin_rows[t]:
                w.writerow({"observable": t, **r})
    print(f"wrote {args.bins_csv}")

    with args.csv.with_name("distributions_bins_summary.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(bin_summary_rows[0]))
        w.writeheader()
        w.writerows(bin_summary_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
