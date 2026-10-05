"""End-to-end check that the verifier's *report* is correct, not just its helpers.

Two real bugs got through unit tests of the helper functions and were caught
only by running the script:

* a ``NameError`` in the exclusion-summary condition (``d`` with no binding),
  which killed the job three seconds in with a traceback;
* a completeness denominator derived from the *filtered* set, so excluding 46
  broken Fortran runs printed ``192/16 runs complete``.

Both live in ``main``'s reporting, so these tests drive ``main`` itself over a
synthetic scan tree: two configurations, two seeds, three arms, with one
Fortran configuration broken in the way the real grid broke.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent.parent / "validation" / "perlmutter"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))

import verify_statistics as vs  # noqa: E402
from compare_distributions import NTP_COLUMNS, NTP_EXTRA, OBSERVABLES  # noqa: E402

N_EVENTS = 400
CFGS = ("cfg_000", "cfg_001")
BROKEN = "cfg_001"  # overflowed counters, like the real grid
COLS = list(NTP_COLUMNS) + list(NTP_EXTRA)


def sample_events(rng: np.random.Generator, n: int) -> dict[str, np.ndarray]:
    """Plausible values for every column the observables read."""
    out = {}
    for i, name in enumerate(COLS):
        centre = 1.0 + i * 0.1
        out[name] = rng.normal(centre, 0.05 * abs(centre) + 1e-3, n)
    return out


def write_npz(path: Path, events: dict[str, np.ndarray]) -> None:
    np.savez(path, **{k: v.astype(np.float32) for k, v in events.items()})


def write_ntuple(path: Path, events: dict[str, np.ndarray]) -> None:
    """The Fortran writes formatted text, 47 items per record, which loadtxt reads."""
    cols = np.column_stack([events[c] for c in COLS])
    np.savetxt(path, cols, fmt="%.8E")
    # savetxt does not reproduce the '(50(1x,es16.8))' field width, and the
    # record count is checked by counting lines, so pad each line out to the real
    # record width.  Without this the fixture's own n-tuple reads as short and
    # every fortran run gets excluded for the wrong reason.
    pad = vs.FORT_RECORD_BYTES - 1
    text = path.read_text()
    path.write_text("".join(line.ljust(pad) + "\n" for line in text.splitlines()))


FORT_CLEAN = """\
  ntries, ngeom =   152553895   111206141
  ntries, nevent, mcall_max:    211748729       {ne}           2
  Integrated cross section (MC, numerical) =   {s:.8E}   {s:.8E}  mu-barns
"""

FORT_BROKEN = """\
  ntries, ngeom =  -382487988 -1006273173
  ntries, nevent, mcall_max:   -382487988       {ne}           1
  Integrated cross section (MC, numerical) =  -1.50309518E-01  -1.50309518E-01  mu-barns
  missm-2: snthcm =   0.00000000
"""

PORT = """\
sigma (MC)       : {s:.8g} micro-barn
throughput       : 70 events/s (7,092,083 trials/s)
"""


def build_tree(root: Path, seed: int = 0) -> None:
    (root / "grid").mkdir(parents=True)
    rng = np.random.default_rng(seed)
    rows = []
    for cfg in CFGS:
        for s in (11, 22):
            run_id = f"{cfg}_s{s}"
            rows.append(
                {"run_id": run_id, "cfg_id": cfg, "seed": s, "n_events": N_EVENTS}
            )
            truth = sample_events(rng, N_EVENTS)
            broken = cfg == BROKEN

            d = root / "fortran" / run_id
            d.mkdir(parents=True)
            sigma = -0.150309518 if broken else 0.0064889900
            text = (FORT_BROKEN if broken else FORT_CLEAN).format(ne=N_EVENTS, s=sigma)
            (d / "out.txt").write_text(text)
            if broken:
                write_ntuple(d / "aao_rad.ntuple", sample_events(rng, N_EVENTS // 2))
            else:
                write_ntuple(d / "aao_rad.ntuple", truth)

            for code in ("py_cpu", "py_gpu"):
                pd = root / code / run_id
                pd.mkdir(parents=True)
                (pd / "out.txt").write_text(PORT.format(s=0.0064889900))
                write_npz(pd / "out.npz", truth)

    with open(root / "grid" / "manifest.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["run_id", "cfg_id", "seed", "n_events"])
        w.writeheader()
        w.writerows(rows)


@pytest.fixture
def tree(tmp_path):
    build_tree(tmp_path / "scan")
    return tmp_path / "scan"


def run_main(tree: Path, capsys) -> tuple[int, str]:
    rc = vs.main(["--root", str(tree), "--out-dir", str(tree / "verify")])
    return rc, capsys.readouterr().out


def test_main_runs_clean_over_a_tree_with_a_broken_reference(tree, capsys):
    """The regression: this raised NameError three seconds into the job."""
    rc, out = run_main(tree, capsys)
    assert rc == 0, out
    assert "Traceback" not in out


def test_completeness_denominator_is_the_manifest_not_the_filtered_set(tree, capsys):
    """The other regression: excluding the broken runs must not shrink the total."""
    _, out = run_main(tree, capsys)
    line = next(l for l in out.splitlines() if l.strip().startswith("fortran:"))
    assert "4/4 runs complete" in line, line
    assert "2 usable after defect checks" in line, line
    # The old bug produced "192/16"; the denominator must never exceed the count.
    reported = int(line.split()[1].split("/")[1])
    assert reported == 4, line


def test_broken_runs_are_reported_as_excluded(tree, capsys):
    _, out = run_main(tree, capsys)
    assert "excluded as unusable" in out
    assert "counter_overflow" in out
    assert "integer*4" in out, "the report must say why, not just drop runs"


def test_broken_configuration_does_not_reach_the_comparisons(tree, capsys):
    """cfg_001 has no usable fortran reference, so it cannot appear in a
    fortran-vs-port row -- and must not silently do so as an outlier."""
    _, out = run_main(tree, capsys)
    for line in out.splitlines():
        if "fortran vs" in line:
            assert "1 cfgs" in line, line
            assert "2 cfgs" not in line, line


def test_the_healthy_pairing_still_compares(tree, capsys):
    _, out = run_main(tree, capsys)
    assert "py_cpu vs py_gpu: 2 cfgs" in out


def test_identical_port_arms_give_a_ratio_of_unity(tree, capsys):
    """The port arms are fed the same events here, so their ratio must be 1."""
    _, out = run_main(tree, capsys)
    line = next(l for l in out.splitlines() if "py_cpu vs py_gpu" in l)
    assert "median 1.000" in line, line


def test_a_fully_clean_tree_excludes_nothing(tmp_path, capsys):
    """With no defects the exclusion section must not appear at all, and the
    denominator must still be right."""
    root = tmp_path / "clean"
    build_tree(root)
    # Rewrite the broken config's log as a healthy one.
    for s in (11, 22):
        d = root / "fortran" / f"cfg_001_s{s}"
        (d / "out.txt").write_text(FORT_CLEAN.format(ne=N_EVENTS, s=0.00648899))
        write_ntuple(d / "aao_rad.ntuple", sample_events(np.random.default_rng(1), N_EVENTS))

    rc, out = run_main(root, capsys)
    assert rc == 0, out
    assert "excluded as unusable" not in out
    assert "2 configurations with both runs" in out


def test_a_run_full_of_missm_warnings_is_still_compared(tmp_path, capsys):
    """The regression that cost the whole comparison.

    167 of the 192 real Fortran runs print 'missm-2' hundreds of times and are
    otherwise healthy -- full n-tuple, positive cross section.  It is a per-event
    clamp, not a failure.  With it counted as a defect the reference arm reported
    1 of 96 configurations usable.
    """
    root = tmp_path / "noisy"
    build_tree(root)
    for s in (11, 22):
        d = root / "fortran" / f"cfg_001_s{s}"
        healthy = (d / "out.txt").read_text()
        (d / "out.txt").write_text(healthy + "  missm-2: snthcm =   0.00000000\n" * 500)
        write_ntuple(d / "aao_rad.ntuple", sample_events(np.random.default_rng(1), N_EVENTS))

    rc, out = run_main(root, capsys)
    assert rc == 0, out
    # cfg_000 is still genuinely broken, cfg_001 is not.
    line = next(l for l in out.splitlines() if l.strip().startswith("fortran:"))
    assert "4/4 runs complete" in line, line
    assert "2 usable after defect checks" in line, line
    assert next(l for l in out.splitlines() if "fortran vs py_gpu" in l).count("1 cfgs")
    assert "'missm-2' in these logs is NOT a defect" in out


def test_observables_are_all_present_in_the_synthetic_ntuple():
    """Guards the fixture: if a column goes missing the KS rows silently shrink."""
    events = sample_events(np.random.default_rng(0), 8)
    missing = {c for _, c, *_ in OBSERVABLES} - set(events)
    assert not missing, missing
    assert set(COLS) >= {c for _, c, *_ in OBSERVABLES}
