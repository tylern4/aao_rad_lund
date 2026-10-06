"""The speedup benchmark's timing and throughput report.

``run_batch.py`` prints a run's elapsed time only every twentieth completion
(``run_batch.py:599``), so a four-run batch prints none -- the benchmark's
timings are read off the filesystem instead, from the ``out.txt`` that
``execute()`` opens with mode ``"wb"`` immediately before launching the
process.  These tests pin down what that reading means: which runs it is valid
for, how a startup is separated from a loop, and how a reduced quota is scaled
back to a production-length run so node-hours can be compared between two arms
that deliberately ran different amounts of work.
"""

from __future__ import annotations

import csv
import importlib.util
import os
import sys
import time
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
REPORT = ROOT / "validation" / "perlmutter" / "bench_report.py"
OVERRIDES = ROOT / "validation" / "perlmutter" / "bench_events.csv"
OVERRIDES_LOADER = ROOT / "validation" / "perlmutter" / "make_grid.py"


def _load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


br = _load(REPORT, "aao_bench_report")
mg = _load(OVERRIDES_LOADER, "aao_make_grid_for_bench")

PORT_SUMMARY = """\
events           : 1,000
trials           : 1,000,000
sampling accept. : 0.1000%
event yield      : 5.000%
weight max / mean: 1.2e+06 / 4.1e+03
above ceiling    : 0 trials, 0.00% of the cross section
sigma (MC)       : 3.41226 micro-barn
sigma (accepted) : 3.41101 micro-barn
throughput       : 100 events/s (100,000 trials/s)
wall time        : 10.0 s
"""

FORT_SUMMARY = """\
  12000  600  100000
ntries, nevent, mcall_max:   600000    1000    900000
Integrated cross section (MC, numerical) = 3.41226E+00 micro-barn
"""


def write(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return path


def test_port_summary_is_read_as_one_run(tmp_path: Path):
    out = write(tmp_path / "cfg_000_s0" / "out.txt", PORT_SUMMARY)
    rec = br.derive(br.read_run(out))
    assert rec["kind"] == "port"
    assert rec["complete"] is True
    assert rec["events"] == 1_000
    assert rec["trials"] == 1_000_000
    assert rec["loop"] == 10.0
    # 1,000,000 trials over a 10 s loop, not over whatever the file's whole
    # life happens to be: compilation is charged separately.
    assert rec["rate"] == pytest.approx(100_000)


def test_fortran_summary_has_no_loop_timer(tmp_path: Path):
    out = write(tmp_path / "cfg_000_s0" / "out.txt", FORT_SUMMARY)
    rec = br.derive(br.read_run(out))
    assert rec["kind"] == "fortran"
    assert rec["complete"] is True
    assert (rec["events"], rec["trials"]) == (1_000, 600_000)
    assert rec["loop"] is None
    # The serial binary reports nothing finer than its own wall clock, so the
    # rate is trials over the file's life.
    if rec["wall"]:
        assert rec["startup"] == 0.0


def test_last_counter_line_wins(tmp_path: Path):
    """re.search takes the first match; the closing counters are the last."""
    text = "ntries, nevent, mcall_max:            1       1\n" + FORT_SUMMARY
    out = write(tmp_path / "cfg_000_s0" / "out.txt", text)
    rec = br.read_run(out)
    assert rec["trials"] == 600_000


def test_wrapped_counter_is_not_a_finished_run(tmp_path: Path):
    """integer*4 ntries prints negative past 2**31 and the run is unusable."""
    text = (
        "ntries, nevent, mcall_max: -2147483648    1000    900000\n"
        "Integrated cross section (MC, numerical) = -1.0E+00 micro-barn\n"
    )
    out = write(tmp_path / "cfg_000_s0" / "out.txt", text)
    rec = br.derive(br.read_run(out))
    assert rec["complete"] is False
    assert rec["rate"] is None
    assert rec["proj"] is None


def test_a_run_stopped_short_is_not_timed(tmp_path: Path):
    """Killed at the wall limit it would report the limit as a real duration."""
    text = "events           : 12,000\ntrials           : 200,000\n"
    out = write(tmp_path / "cfg_000_s0" / "out.txt", text)
    rec = br.derive(br.read_run(out))
    assert rec["complete"] is False
    assert rec["proj"] is None


def test_port_projection_charges_startup_once(tmp_path: Path):
    out = write(tmp_path / "cfg_000_s0" / "out.txt", PORT_SUMMARY)
    rec = br.derive(br.read_run(out), target=20_000)
    if rec["wall"] is None:
        pytest.skip("this platform exposes no birth time")
    # startup does not scale with the event count; the loop does.
    expected = rec["startup"] + 20_000 * (rec["trials"] / rec["events"]) / rec["rate"]
    assert rec["proj"] == pytest.approx(expected)
    assert rec["proj"] >= 20_000 * (rec["trials"] / rec["events"]) / rec["rate"]


def test_fortran_projection_scales_linearly(tmp_path: Path):
    out = write(tmp_path / "cfg_000_s0" / "out.txt", FORT_SUMMARY)
    rec = br.derive(br.read_run(out), target=20_000)
    if rec["wall"] is None:
        pytest.skip("this platform exposes no birth time")
    assert rec["startup"] == 0.0
    assert rec["proj"] == pytest.approx(rec["wall"] * 20_000 / rec["events"])


def test_wall_clock_runs_from_birth_to_last_write(tmp_path: Path):
    out = write(tmp_path / "cfg_000_s0" / "out.txt", PORT_SUMMARY)
    born = br.birth_time(out)
    if born is None:
        pytest.skip("this platform exposes no birth time")
    assert born > 0
    # mtime is what the child's last write set; move it to make the arithmetic
    # checkable without sleeping.
    future = time.time() + 100
    os.utime(out, (future, future))
    rec = br.read_run(out)
    assert 95 <= rec["wall"] <= 110


def test_report_renders_medians_and_node_hours(tmp_path: Path, capsys):
    recs = []
    for name in ("cfg_000_s0", "cfg_034_s0"):
        out = write(tmp_path / name / "out.txt", PORT_SUMMARY)
        recs.append(br.derive(br.read_run(out)))
    if any(r["wall"] is None for r in recs):
        pytest.skip("this platform exposes no birth time")
    br.report("pass 1: cold XLA cache", recs, concurrency=4, grid=192)
    text = capsys.readouterr().out
    assert "pass 1: cold XLA cache" in text
    assert "2/2 complete" in text
    assert "median wall" in text
    assert "node-h per run at concurrency 4" in text
    assert "192-run grid" in text


def test_report_counts_an_unfinished_run_as_incomplete(tmp_path: Path, capsys):
    out = write(tmp_path / "cfg_000_s0" / "out.txt", "events           : 3\n")
    br.report("partial", [br.derive(br.read_run(out))], concurrency=4)
    text = capsys.readouterr().out
    assert "0/1 complete" in text
    assert "INCOMPLETE" in text


def test_main_finds_runs_and_reports_a_missing_directory(tmp_path: Path, capsys):
    write(tmp_path / "runs" / "cfg_000_s0" / "out.txt", PORT_SUMMARY)
    assert br.main([str(tmp_path / "runs"), "--label", "warm"]) == 0
    assert "warm" in capsys.readouterr().out

    empty = tmp_path / "empty"
    empty.mkdir()
    assert br.main([str(empty)]) == 1


def test_concurrency_scales_the_node_hours(tmp_path: Path):
    out = write(tmp_path / "cfg_000_s0" / "out.txt", PORT_SUMMARY)
    rec = br.derive(br.read_run(out))
    if rec["proj"] is None:
        pytest.skip("this platform exposes no birth time")
    alone = br.report("run", [rec], concurrency=1)
    packed = br.report("run", [rec], concurrency=128)
    assert alone["node_hours"] == pytest.approx(packed["node_hours"] * 128, rel=1e-9)
    assert alone["node_hours"] == pytest.approx(rec["proj"] / 3600, rel=1e-9)


def test_bench_overrides_are_reductions_within_the_standard_quota():
    """The quota exists to keep integer*4 ntries inside 2**31, so a raise is
    never useful and make_grid rejects one; every bench row has to be a
    reduction, and has to buy the ~600 s of serial work the job was sized for."""
    rows = list(csv.DictReader(open(OVERRIDES)))
    assert rows, "no bench quotas"
    assert mg.load_overrides(OVERRIDES) == {
        r["cfg_id"]: int(r["n_events"]) for r in rows
    }
    for row in rows:
        n = int(row["n_events"])
        assert 0 < n <= mg.N_EVENTS
        tpe = int(row["trials_per_event"])
        rate = int(row["trials_per_second"])
        seconds = n * tpe / rate
        assert 590 <= seconds <= 600 + 300, (
            f"{row['cfg_id']}: {n} events x {tpe} trials/event at {rate} "
            f"trials/s is {seconds:.0f}s, not the ~600s the debug job is sized for"
        )


def test_bench_configurations_are_four_energies_and_both_deltas():
    """delta, not beam energy, is what decides how expensive a configuration
    is -- median trials per event 58,566 at 0.005 against 9,118 at 0.05 -- so a
    benchmark covering only one of them measures a fifth of the grid."""
    rows = list(csv.DictReader(open(OVERRIDES)))
    cfgs = {r["cfg_id"] for r in rows}
    assert cfgs == {"cfg_000", "cfg_034", "cfg_064", "cfg_082"}
