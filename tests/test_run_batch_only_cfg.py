"""Selecting a subset of configurations to re-run.

Every run in this scan has printed a cross section -- including the ones whose
32-bit counter wrapped and returned a negative one -- so ``run_complete`` is
true for all 192 of them.  A filtered batch that honoured completeness would
select 70 runs, run none, and exit 0.  That silent no-op is the thing these
tests exist to prevent.
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

ROOT = Path(__file__).resolve().parent.parent
DRIVER = ROOT / "validation" / "perlmutter" / "run_batch.py"


def _load_driver():
    spec = importlib.util.spec_from_file_location("aao_run_batch_only_cfg", DRIVER)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


rb = _load_driver()


def write_manifest(root: Path, rows: list[tuple[str, str]]) -> list[dict]:
    (root / "grid").mkdir(parents=True, exist_ok=True)
    with open(root / "grid" / "manifest.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["run_id", "cfg_id", "seed", "n_events"])
        for run_id, cfg in rows:
            w.writerow([run_id, cfg, 0, 20_000])
    return rb.load_manifest(root)


def mark_complete(root: Path, code: str, run_id: str, cfg_id: str) -> None:
    d = root / code / run_id
    d.mkdir(parents=True, exist_ok=True)
    marker = rb.SIGMA_MARKERS[code]
    (d / "out.txt").write_text(f" ... {marker} ... \n")
    # run_complete wants the data file alongside the marker.
    (d / ("aao_rad.ntuple" if code == "fortran" else "out.npz")).write_bytes(b"")
    (root / "grid").mkdir(parents=True, exist_ok=True)
    (root / "grid" / f"{cfg_id}.txt").write_text("7\n")  # card, for build_job


def args_for(**kw) -> SimpleNamespace:
    base = dict(code="fortran", only_cfg=None, limit=0)
    base.update(kw)
    return argparse.Namespace(**base)


# --------------------------------------------------------------- parse_only_cfg


def test_an_unset_selector_selects_nothing_in_particular():
    assert rb.parse_only_cfg(None) is None
    assert rb.parse_only_cfg("") is None


def test_a_comma_list_is_split_and_trimmed():
    assert rb.parse_only_cfg("cfg_001,cfg_002 , cfg_003") == {
        "cfg_001",
        "cfg_002",
        "cfg_003",
    }


def test_a_bare_file_lists_one_id_per_line(tmp_path):
    f = tmp_path / "ids.txt"
    f.write_text("# the 35 broken ones\ncfg_040\n\ncfg_041\n")
    assert rb.parse_only_cfg(f"@{f}") == {"cfg_040", "cfg_041"}


def test_a_csv_override_file_is_read_through_its_header(tmp_path):
    """``plan_ref_rerun.py --write`` produces exactly this shape, and the point
    of ``@file`` is that the plan passes through untouched."""
    f = tmp_path / "n_events_override.csv"
    f.write_text("cfg_id,n_events\ncfg_040,5412\ncfg_041,15080\n")
    assert rb.parse_only_cfg(f"@{f}") == {"cfg_040", "cfg_041"}


def test_a_missing_file_fails_loudly(tmp_path):
    with pytest.raises(OSError):
        rb.parse_only_cfg(f"@{tmp_path / 'nope.csv'}")


# ------------------------------------------------------------------ select_runs


def test_selecting_nothing_keeps_every_run():
    runs = [{"cfg_id": "cfg_001"}, {"cfg_id": "cfg_002"}]
    assert rb.select_runs(runs, None) == runs


def test_select_runs_keeps_both_seeds_of_the_named_config():
    runs = [
        {"cfg_id": "cfg_001"},
        {"cfg_id": "cfg_001"},
        {"cfg_id": "cfg_002"},
    ]
    kept = rb.select_runs(runs, "cfg_001")
    assert len(kept) == 2


def test_an_unknown_id_aborts_instead_of_running_fewer_runs(tmp_path):
    """Silently dropping an id would shrink the batch to 34 and still exit 0."""
    runs = [{"cfg_id": "cfg_001"}]
    with pytest.raises(SystemExit) as exc:
        rb.select_runs(runs, "cfg_001,cfg_999")
    assert "cfg_999" in str(exc.value)


# --------------------------------------------------------------- select_pending


def test_a_complete_run_is_skipped_when_nothing_was_selected(tmp_path):
    root = tmp_path / "scan"
    mark_complete(root, "fortran", "cfg_001_s0", "cfg_001")
    runs = [{"run_id": "cfg_001_s0", "cfg_id": "cfg_001"}]
    pending, rerun = rb.select_pending(args_for(), runs, root)
    assert pending == [] and rerun is False


def test_selecting_configurations_forces_them_to_run_again(tmp_path):
    """The trap: all 192 runs are complete by this test's own definition."""
    root = tmp_path / "scan"
    mark_complete(root, "fortran", "cfg_001_s0", "cfg_001")
    runs = [{"run_id": "cfg_001_s0", "cfg_id": "cfg_001"}]
    pending, rerun = rb.select_pending(args_for(only_cfg="cfg_001"), runs, root)
    assert rerun is True
    assert [r["run_id"] for r in pending] == ["cfg_001_s0"]


def test_an_incomplete_run_is_pending_without_a_selector(tmp_path):
    root = tmp_path / "scan"
    runs = [{"run_id": "cfg_001_s0", "cfg_id": "cfg_001"}]
    pending, _ = rb.select_pending(args_for(), runs, root)
    assert [r["run_id"] for r in pending] == ["cfg_001_s0"]


def test_limit_still_applies_after_a_selection(tmp_path):
    root = tmp_path / "scan"
    for i in range(4):
        mark_complete(root, "fortran", f"cfg_{i:03d}_s0", f"cfg_{i:03d}")
    runs = [{"run_id": f"cfg_{i:03d}_s0", "cfg_id": f"cfg_{i:03d}"} for i in range(4)]
    pending, _ = rb.select_pending(
        args_for(only_cfg="cfg_000,cfg_001,cfg_002,cfg_003", limit=2), runs, root
    )
    assert len(pending) == 2


# ------------------------------------------------------------------- through main


def test_a_filtered_dry_run_reports_a_rerun_not_a_skip(tmp_path):
    root = tmp_path / "scan"
    mark_complete(root, "py_cpu", "cfg_001_s0", "cfg_001")
    mark_complete(root, "py_cpu", "cfg_001_s1", "cfg_001")
    mark_complete(root, "py_cpu", "cfg_002_s0", "cfg_002")
    write_manifest(
        root, [("cfg_001_s0", "cfg_001"), ("cfg_001_s1", "cfg_001"),
               ("cfg_002_s0", "cfg_002")]
    )

    r = subprocess.run(
        [sys.executable, str(DRIVER), "--code", "py_cpu", "--root", str(root),
         "--only-cfg", "cfg_001", "--dry-run"],
        capture_output=True, text=True,
    )
    assert r.returncode == 0, r.stderr or r.stdout
    assert "selected run(s) to re-run" in r.stdout, r.stdout
    assert "complete already" not in r.stdout, r.stdout
    # Both seeds of the selected configuration, and neither of the other's.
    assert r.stdout.count("would run:") == 2, r.stdout


def test_without_a_selector_the_dry_run_reports_the_skip(tmp_path):
    root = tmp_path / "scan"
    mark_complete(root, "py_cpu", "cfg_001_s0", "cfg_001")
    write_manifest(root, [("cfg_001_s0", "cfg_001")])

    r = subprocess.run(
        [sys.executable, str(DRIVER), "--code", "py_cpu", "--root", str(root),
         "--dry-run"],
        capture_output=True, text=True,
    )
    assert r.returncode == 0, r.stderr or r.stdout
    assert "1/1 complete already" in r.stdout, r.stdout
    assert "0 to run" in r.stdout, r.stdout


def test_the_override_file_can_be_passed_straight_through(tmp_path):
    """No transcribing 35 ids by hand."""
    root = tmp_path / "scan"
    for cfg in ("cfg_001", "cfg_002"):
        for s in (0, 1):
            mark_complete(root, "py_cpu", f"{cfg}_s{s}", cfg)
    write_manifest(
        root,
        [(f"{cfg}_s{s}", cfg) for cfg in ("cfg_001", "cfg_002") for s in (0, 1)],
    )
    plan = root / "grid" / "n_events_override.csv"
    plan.write_text("cfg_id,n_events\ncfg_001,5412\n")

    r = subprocess.run(
        [sys.executable, str(DRIVER), "--code", "py_cpu", "--root", str(root),
         "--only-cfg", f"@{plan}", "--dry-run"],
        capture_output=True, text=True,
    )
    assert r.returncode == 0, r.stderr or r.stdout
    assert r.stdout.count("would run:") == 2, r.stdout


# ------------------------------------------------------- stale outputs on a rerun


def test_the_previous_n_tuple_is_cleared_before_a_rerun(tmp_path):
    """The Fortran opens its n-tuple without status='replace', so a re-run
    rewinds and overwrites from the start, leaving the old tail behind.  A
    reduced quota over a longer old file would hand the verifier records from a
    run that no longer exists -- and it would pass, because the file is longer
    than the quota."""
    root = tmp_path / "scan"
    run = {"run_id": "cfg_001_s0", "cfg_id": "cfg_001"}
    d = root / "fortran" / "cfg_001_s0"
    d.mkdir(parents=True)
    (d / "aao_rad.ntuple").write_bytes(b"x" * 800 * 20_000)
    (d / "out.txt").write_text("previous run's cross section\n")
    (d / "unrelated").write_text("keep me")

    removed = rb.clear_stale_outputs([run], "fortran", root)
    assert removed == 2
    assert not (d / "aao_rad.ntuple").exists()
    assert not (d / "out.txt").exists()
    assert (d / "unrelated").read_text() == "keep me"


def test_only_the_runs_being_run_are_cleared(tmp_path):
    root = tmp_path / "scan"
    for rid in ("cfg_001_s0", "cfg_002_s0"):
        d = root / "fortran" / rid
        d.mkdir(parents=True)
        (d / "aao_rad.ntuple").write_bytes(b"x")
    rb.clear_stale_outputs([{"run_id": "cfg_001_s0", "cfg_id": "cfg_001"}], "fortran", root)
    assert not (root / "fortran" / "cfg_001_s0" / "aao_rad.ntuple").exists()
    assert (root / "fortran" / "cfg_002_s0" / "aao_rad.ntuple").exists()


def test_nothing_left_to_clear_is_not_an_error(tmp_path):
    assert rb.clear_stale_outputs([], "fortran", tmp_path) == 0


# ------------------------------------------------------------- the run card


def test_a_reduced_quota_reaches_the_run_card(tmp_path):
    """Without this the Fortran would keep generating the original quota, wrap
    its counter exactly as before, and print a cross section that looks no
    different from the rejected one."""
    root = tmp_path / "scan"
    mark_complete(root, "py_cpu", "cfg_001_s0", "cfg_001")
    write_manifest(root, [("cfg_001_s0", "cfg_001")])
    (root / "grid" / "cfg_001.txt").write_text("7\n5412\n")
    stale = root / "py_cpu" / "cfg_001_s0" / "run_card.txt"
    stale.write_text("7\n20000\n")

    r = subprocess.run(
        [sys.executable, str(DRIVER), "--code", "py_cpu", "--root", str(root),
         "--only-cfg", "cfg_001", "--dry-run"],
        capture_output=True, text=True,
    )
    assert r.returncode == 0, r.stderr or r.stdout
    assert stale.read_text() == "7\n5412\n", "the stale card survived"
    assert "run card refreshed" in r.stdout, r.stdout


def test_an_unchanged_card_is_left_quiet(tmp_path):
    root = tmp_path / "scan"
    mark_complete(root, "py_cpu", "cfg_001_s0", "cfg_001")
    write_manifest(root, [("cfg_001_s0", "cfg_001")])
    card = root / "grid" / "cfg_001.txt"
    (root / "py_cpu" / "cfg_001_s0" / "run_card.txt").write_text(card.read_text())

    r = subprocess.run(
        [sys.executable, str(DRIVER), "--code", "py_cpu", "--root", str(root),
         "--only-cfg", "cfg_001", "--dry-run"],
        capture_output=True, text=True,
    )
    assert r.returncode == 0, r.stderr or r.stdout
    assert "run card refreshed" not in r.stdout, r.stdout
