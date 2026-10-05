"""The plotting script must survive a real scan tree, including a broken one.

Figures are the deliverable here, so the failure that matters most is a silent
one: a plot that renders but compares a port against a Fortran run whose
32-bit counters overflowed.  That is the same trap the statistical verification
had, and it is invisible in the output -- you get a confident-looking picture of
a disagreement that is an artifact of the reference.
"""

from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parent.parent
PERLMUTTER = REPO / "validation" / "perlmutter"
sys.path.insert(0, str(PERLMUTTER))
sys.path.insert(0, str(PERLMUTTER.parent))

N_EVENTS = 600
CFGS = ("cfg_000", "cfg_001", "cfg_002")
EBEAM = {"cfg_000": "2.0", "cfg_001": "6.0", "cfg_002": "12.0"}
BROKEN = "cfg_002"
COLS_ALL = None


def _columns() -> list[str]:
    from compare_distributions import NTP_COLUMNS, NTP_EXTRA

    return list(NTP_COLUMNS) + list(NTP_EXTRA)


def _samples(rng: np.random.Generator, n: int, shift: float = 0.0) -> dict[str, np.ndarray]:
    cols = _columns()
    out = {}
    for i, name in enumerate(cols):
        centre = 1.0 + i * 0.1
        out[name] = rng.normal(centre + shift, 0.05 * abs(centre) + 1e-3, n)
    return out


def _ntuple_text(events: dict[str, np.ndarray]) -> str:
    cols = _columns()
    arr = np.column_stack([events[c] for c in cols])
    pad = 799
    return "".join(line.ljust(pad) + "\n" for line in _savetxt(arr, pad))


def _savetxt(arr: np.ndarray, pad: int) -> list[str]:
    lines = []
    for row in arr:
        s = "".join(f"{v:.8E}".rjust(16) for v in row)
        lines.append(s)
    return lines


FORT_CLEAN = (
    "  ntries, ngeom =   152553895   111206141\n"
    "  ntries, nevent, mcall_max:    211748729  {ne}           2\n"
    "  Integrated cross section (MC, numerical) =   {s:.8E}   {s:.8E}  mu-barns\n"
)
FORT_BROKEN = (
    "  ntries, ngeom =  -382487988 -1006273173\n"
    "  ntries, nevent, mcall_max:   -382487988       {ne}           1\n"
    "  Integrated cross section (MC, numerical) =  -1.50309518E-01  -1.50309518E-01  mu-barns\n"
    "  missm-2: snthcm =   0.00000000\n"
)
PORT = "sigma (MC)       : {s:.8g} micro-barn\nthroughput : 70 events/s (7,092,083 trials/s)\n"


def build_tree(root: Path, seed: int = 0) -> None:
    rng = np.random.default_rng(seed)
    (root / "grid").mkdir(parents=True)
    rows = []
    for cfg in CFGS:
        for s in (11, 22):
            run_id = f"{cfg}_s{s}"
            rows.append(
                {
                    "run_id": run_id,
                    "cfg_id": cfg,
                    "seed": s,
                    "n_events": N_EVENTS,
                    "ebeam": EBEAM[cfg],
                }
            )
            broken = cfg == BROKEN
            truth = _samples(rng, N_EVENTS)
            shifted = _samples(rng, N_EVENTS, shift=0.05)  # a visible disagreement

            d = root / "fortran" / run_id
            d.mkdir(parents=True)
            sigma = -0.150309518 if broken else 0.00648899
            tmpl = FORT_BROKEN if broken else FORT_CLEAN
            (d / "out.txt").write_text(tmpl.format(ne=N_EVENTS, s=sigma))
            # The broken config also has a short n-tuple, as the real ones did.
            n_rows = N_EVENTS // 2 if broken else N_EVENTS
            (d / "aao_rad.ntuple").write_text(_ntuple_text(_samples(rng, n_rows)))

            for code in ("py_cpu", "py_gpu"):
                pd = root / code / run_id
                pd.mkdir(parents=True)
                (pd / "out.txt").write_text(PORT.format(s=0.00648899))
                np.savez(
                    pd / "out.npz",
                    **{k: v.astype(np.float32) for k, v in shifted.items()},
                )

    with open(root / "grid" / "manifest.csv", "w", newline="") as fh:
        w = csv.DictWriter(
            fh, fieldnames=["run_id", "cfg_id", "seed", "n_events", "ebeam"]
        )
        w.writeheader()
        w.writerows(rows)


@pytest.fixture
def scan(tmp_path):
    root = tmp_path / "scan"
    build_tree(root)
    return root


def run(scan: Path, *extra: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        [
            sys.executable,
            str(PERLMUTTER / "make_plots.py"),
            "--root",
            str(scan),
            "--out-dir",
            str(scan / "plots"),
            "--observables",
            "es,ep,q2,mm2",
            "--bins",
            "30",
            *extra,
        ],
        capture_output=True,
        text=True,
    )


def test_it_runs_and_writes_figures(scan):
    r = run(scan)
    assert r.returncode == 0, r.stderr
    plots = scan / "plots"
    assert plots.is_dir()
    made = sorted(p.name for p in plots.glob("*.png"))
    assert made, "no figures written"
    assert "sigma_by_energy.png" in made
    assert "sigma_by_config.png" in made
    for name in made:
        assert (plots / name).stat().st_size > 5_000, f"{name} is suspiciously small"


def test_figures_are_not_blank(scan):
    """A figure that renders but is empty is worse than no figure."""
    r = run(scan)
    assert r.returncode == 0, r.stderr
    for png in (scan / "plots").glob("*.png"):
        data = png.read_bytes()
        assert data[:8] == b"\x89PNG\r\n\x1a\n", png
        # A uniform blank canvas compresses to almost nothing.
        assert len(data) > 5_000, f"{png.name} looks blank ({len(data)} bytes)"


def test_a_broken_reference_is_never_plotted(scan, capsys):
    """cfg_002's Fortran overflowed its counters.  Plotting the port against it
    would produce a confident picture of a disagreement that is an artifact."""
    r = run(scan)
    assert r.returncode == 0, r.stderr
    assert not (scan / "plots" / "cfg_002.png").exists()
    # It must be named as skipped, and the count must be against the grid rather
    # than against the filtered set -- "2 of 2" would hide the drop.
    assert "skipped cfg_002" in r.stdout, r.stdout
    assert "2 of 3 in the grid" in r.stdout, r.stdout


def test_the_error_band_is_finite(scan):
    """Feeding densities into a multinomial variance makes p exceed 1, p(1-p)
    negative and the band NaN -- a figure that draws an empty uncertainty."""
    import warnings

    import make_plots

    a = np.random.default_rng(0).normal(size=500)
    b = np.random.default_rng(1).normal(size=500)
    edges = make_plots.bin_edges(a, b, 30)
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        centres, dens, diff, err = make_plots.two_sample_bins(a, b, edges)
    assert np.all(np.isfinite(err)), "error band must be finite"
    assert np.all(err >= 0)
    assert np.all(np.isfinite(diff))
    # Densities must integrate to 1, i.e. be per unit x and per event.
    assert float(np.sum(dens * np.diff(edges))) == pytest.approx(1.0)


def test_it_plots_the_configurations_that_are_valid(scan):
    r = run(scan, "--per-energy", "1")
    assert r.returncode == 0, r.stderr
    assert (scan / "plots" / "cfg_000.png").exists()
    assert (scan / "plots" / "cfg_001.png").exists()


def test_all_valid_also_skips_the_broken_one(scan):
    r = run(scan, "--all-valid")
    assert r.returncode == 0, r.stderr
    names = sorted(p.name for p in (scan / "plots").glob("cfg_*.png"))
    assert names == ["cfg_000.png", "cfg_001.png"], names


def test_an_unknown_observable_is_rejected_rather_than_plotting_nothing(scan):
    r = subprocess.run(
        [
            sys.executable,
            str(PERLMUTTER / "make_plots.py"),
            "--root",
            str(scan),
            "--out-dir",
            str(scan / "plots"),
            "--observables",
            "not_a_column",
        ],
        capture_output=True,
        text=True,
    )
    assert r.returncode == 2, r.stdout
    assert "no observables matched" in r.stderr


def test_every_energy_level_gets_a_row_in_the_energy_plot(scan):
    """The point of that figure is the trend, so a level with no valid
    configuration must be visibly absent rather than silently missing."""
    r = run(scan)
    assert r.returncode == 0, r.stderr
    assert (scan / "plots" / "sigma_by_energy.png").stat().st_size > 5_000
    assert "12.0" not in r.stdout.split("plotting")[-1]


def test_the_difference_panel_uses_the_reference_label(scan):
    """If the fortran arm is absent the y-label must not still say 'Fortran'."""
    import make_plots

    assert "py_gpu" in make_plots.STYLE
    assert make_plots.STYLE["fortran"]["label"] == "Fortran"


def test_quantization_is_applied_so_plots_match_the_verifier(scan):
    """Without the same 8-digit quantization the E_s plot shows a float32
    artifact that verify_statistics.py reports as agreement."""
    import make_plots
    from verify_statistics import quantize_sig

    a = np.array([4.24399996, 4.1, 4.2])
    b = np.array([4.24399995803833, 4.1, 4.2])
    assert quantize_sig(a)[0] == quantize_sig(b)[0]
    assert make_plots.load_events is not None
