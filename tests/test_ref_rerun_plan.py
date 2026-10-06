"""Reducing an event quota is only safe if every consumer sees the same number.

The run card and the manifest are written by the same script but from different
expressions -- the card used the ``N_EVENTS`` constant while the manifest used
``cfg["n_events"]``.  An override written to one and not the other would have
the Fortran generate 20,000 events while the verifier expected 5,412, reporting
every run as short, or worse, leave the original quota in place and let the
32-bit counter wrap again exactly as it did before.
"""

from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
PERLMUTTER = REPO / "validation" / "perlmutter"
sys.path.insert(0, str(PERLMUTTER))

import make_grid  # noqa: E402
import plan_ref_rerun as plan_mod  # noqa: E402


# ---------------------------------------------------------------- make_grid


def test_card_reads_the_quota_from_the_config_not_the_constant():
    """The manifest already used cfg["n_events"]; the card did not."""
    cfg = make_grid.build_configs()[0]
    cfg["n_events"] = 7306
    card = make_grid.card_text(cfg).splitlines()
    assert card[15] == "7306", f"n_events is card line 16, got {card[15]!r}"


def test_card_defaults_to_the_standard_quota():
    cfg = make_grid.build_configs()[0]
    assert cfg["n_events"] == make_grid.N_EVENTS
    assert make_grid.card_text(cfg).splitlines()[15] == str(make_grid.N_EVENTS)


def test_an_override_touches_only_the_named_configuration():
    cfgs = make_grid.build_configs({"cfg_007": 15_000})
    by_id = {c["cfg_id"]: c["n_events"] for c in cfgs}
    assert by_id["cfg_007"] == 15_000
    assert by_id["cfg_006"] == make_grid.N_EVENTS
    assert by_id["cfg_008"] == make_grid.N_EVENTS
    assert len(by_id) == 96


def test_overrides_are_validated(tmp_path):
    """A quota above the standard would reintroduce the overflow the override
    exists to prevent, and zero would make ntries a divisor of nothing."""
    good = tmp_path / "ok.csv"
    good.write_text("cfg_id,n_events\ncfg_000,5412\n")
    assert make_grid.load_overrides(good) == {"cfg_000": 5412}

    too_big = tmp_path / "big.csv"
    too_big.write_text(f"cfg_id,n_events\ncfg_000,{make_grid.N_EVENTS + 1}\n")
    with pytest.raises(ValueError, match="more than the standard quota"):
        make_grid.load_overrides(too_big)

    zero = tmp_path / "zero.csv"
    zero.write_text("cfg_id,n_events\ncfg_000,0\n")
    with pytest.raises(ValueError, match="non-positive"):
        make_grid.load_overrides(zero)

    assert make_grid.load_overrides(None) == {}


def test_card_and_manifest_agree_after_a_regenerated_grid(tmp_path):
    """The point of the whole exercise: both consumers see one quota."""
    root = tmp_path / "scan"
    (root / "grid").mkdir(parents=True)
    (root / "grid" / "n_events_override.csv").write_text(
        "cfg_id,n_events\ncfg_000,5412\n"
    )

    r = subprocess.run(
        [sys.executable, str(PERLMUTTER / "make_grid.py"), "--root", str(root)],
        capture_output=True, text=True,
    )
    assert r.returncode == 0, r.stderr
    assert "overridden for 1 configuration" in r.stdout, r.stdout

    card = (root / "grid" / "cfg_000.txt").read_text().splitlines()
    assert card[15] == "5412", card[15]
    # An untouched configuration must not inherit it.
    assert (root / "grid" / "cfg_001.txt").read_text().splitlines()[15] == str(
        make_grid.N_EVENTS
    )

    rows = list(csv.DictReader(open(root / "grid" / "manifest.csv")))
    by_cfg = {r["cfg_id"]: int(r["n_events"]) for r in rows}
    assert by_cfg["cfg_000"] == 5412
    assert by_cfg["cfg_001"] == make_grid.N_EVENTS
    assert len(rows) == 192
    # Every run of the reduced configuration carries the reduced quota, since
    # the verifier checks each run against it.
    assert {int(r["n_events"]) for r in rows if r["cfg_id"] == "cfg_000"} == {5412}


def test_regenerating_without_a_flag_still_honours_the_saved_override(tmp_path):
    """Every sbatch script calls make_grid with only --root.  An override that
    survived only on the command line would be dropped by the next job to start,
    and the quota would revert silently between the Fortran arm and the verifier."""
    root = tmp_path / "scan"
    (root / "grid").mkdir(parents=True)
    (root / "grid" / "n_events_override.csv").write_text(
        "cfg_id,n_events\ncfg_003,8269\n"
    )
    r = subprocess.run(
        [sys.executable, str(PERLMUTTER / "make_grid.py"), "--root", str(root)],
        capture_output=True, text=True,
    )
    assert r.returncode == 0, r.stderr
    rows = list(csv.DictReader(open(root / "grid" / "manifest.csv")))
    assert {int(x["n_events"]) for x in rows if x["cfg_id"] == "cfg_003"} == {8269}


# ------------------------------------------------------------------ planner

OK_TUPLE = b"".join(b"x" * 799 + b"\n" for _ in range(20))


def write_run(
    root: Path,
    run_id: str,
    n_events: int,
    *,
    progress: list[tuple[int, int]] | None = None,
    final: tuple[int, int] = (50_000, 40_000),
    sigma: str = "1.00000000E-03   1.00000000E-03",
    records: int | None = 20,
    with_out: bool = True,
) -> None:
    d = root / "fortran" / run_id
    d.mkdir(parents=True, exist_ok=True)
    if not with_out:
        return
    lines = []
    for n, e in progress or [(1_000, 800), (2_000, 1_600)]:
        lines.append(f" ntries, nevent, mcall_max: {n:>10} {e:>8}           1")
        lines.append(
            f" Integrated cross section (MC, numerical) = {sigma}  mu-barns"
        )
    lines.append(f" ntries, ngeom = {final[0]:>10} {final[1]:>10}")
    (d / "out.txt").write_text("\n".join(lines) + "\n")
    if records is not None:
        (d / "aao_rad.ntuple").write_bytes(b"".join(b"x" * 799 + b"\n" for _ in range(records)))


def write_manifest(root: Path, rows: list[tuple[str, str, int]]) -> None:
    (root / "grid").mkdir(parents=True, exist_ok=True)
    with open(root / "grid" / "manifest.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["run_id", "cfg_id", "seed", "n_events"])
        for run_id, cfg, n_events in rows:
            w.writerow([run_id, cfg, 0, n_events])


def test_trials_per_event_uses_the_last_positive_sample_not_the_first():
    """The first line samples the run before its rejection ceiling matures, when
    acceptance is still poor; it overstates trials per event and would cut the
    quota further than necessary."""
    text = "\n".join(
        f" ntries, nevent, mcall_max: {n:>10} {e:>8}           1"
        for n, e in [(80_000, 800), (400_000, 4_000), (1_000_000, 10_000)]
    )
    assert plan_mod.trials_per_event(text) == pytest.approx(100.0)


def test_trials_per_event_ignores_wrapped_lines():
    """A negative counter is wrapped, so it is not a sample of anything."""
    text = (
        " ntries, nevent, mcall_max:      100000        800           1\n"
        " ntries, nevent, mcall_max:   -90000000       1600           1\n"
    )
    assert plan_mod.trials_per_event(text) == pytest.approx(125.0)


def test_a_wrapped_reference_gets_a_smaller_quota(tmp_path):
    root = tmp_path / "scan"
    # 3e5 trials/event at the standard quota would need 6e9 trials.
    write_run(root, "cfg_040_s0", 20_000, progress=[(30_000_000, 100), (300_000_000, 1_000)],
              final=(-1_000, -1_000), records=20)
    write_run(root, "cfg_040_s1", 20_000, progress=[(30_000_000, 100), (300_000_000, 1_000)],
              final=(-1_000, -1_000), records=20)
    write_manifest(root, [("cfg_040_s0", "cfg_040", 20_000),
                          ("cfg_040_s1", "cfg_040", 20_000)])
    out = plan_mod.plan(root, 0.9)
    assert out["cfg_040"]["n_events"] == int(0.9 * 2**31 / 300_000.0)
    assert out["cfg_040"]["n_events"] < 20_000
    assert "counter_overflow" in out["cfg_040"]["causes"]


def test_a_truncated_ntuple_costs_no_reduction(tmp_path):
    """Those are the cheap configurations.  Their cross section is fine -- it
    estimates a constant and agrees with their own earlier estimates -- so
    cutting their quota would only throw away statistics to fix nothing."""
    root = tmp_path / "scan"
    write_run(root, "cfg_002_s0", 20_000, progress=[(3_000, 800)], records=10)
    write_run(root, "cfg_002_s1", 20_000, progress=[(3_000, 800)], records=10)
    write_manifest(root, [("cfg_002_s0", "cfg_002", 20_000),
                          ("cfg_002_s1", "cfg_002", 20_000)])
    out = plan_mod.plan(root, 0.9)
    assert out["cfg_002"]["causes"] == ["short_n_tuple"]
    assert out["cfg_002"]["n_events"] == 20_000


def test_a_healthy_configuration_is_not_in_the_plan(tmp_path):
    """A full n-tuple, a positive counter and a cross section: nothing to rerun.

    The quota and the record count have to match, or the run would be short and
    land in the plan for the wrong reason.
    """
    root = tmp_path / "scan"
    for s in (0, 1):
        write_run(root, f"cfg_001_s{s}", 20, progress=[(3_000, 800), (9_000, 2_400)],
                  final=(50_000, 40_000), records=20)
    write_manifest(root, [("cfg_001_s0", "cfg_001", 20),
                          ("cfg_001_s1", "cfg_001", 20)])
    assert plan_mod.plan(root, 0.9) == {}


def test_a_run_with_no_output_is_reported_not_skipped(tmp_path):
    """collect() drops runs with no output; a planner built on it would omit the
    configurations most in need of a rerun and then report success."""
    root = tmp_path / "scan"
    write_run(root, "cfg_050_s0", 20_000, with_out=False)
    write_manifest(root, [("cfg_050_s0", "cfg_050", 20_000)])
    out = plan_mod.plan(root, 0.9)
    assert out["cfg_050"]["causes"] == ["no_output"]
    assert out["cfg_050"]["n_events"] == 20_000


def test_the_worst_projected_counter_stays_inside_the_margin(tmp_path):
    root = tmp_path / "scan"
    rows = []
    for i, tpe in enumerate([128_000, 310_000, 357_000, 1_000_000]):
        cfg = f"cfg_{i:03d}"
        for s in (0, 1):
            run = f"{cfg}_s{s}"
            write_run(root, run, 20_000, progress=[(tpe * 100, 100)],
                      final=(-1, -1), records=20)
            rows.append((run, cfg, 20_000))
    write_manifest(root, rows)

    out = plan_mod.plan(root, 0.9)
    for cfg, e in out.items():
        projected = e["trials_per_event"] * e["n_events"]
        assert projected <= 0.9 * 2**31 + 1, f"{cfg} projects {projected:,.0f}"


def test_the_cli_writes_where_make_grid_reads(tmp_path):
    root = tmp_path / "scan"
    write_run(root, "cfg_040_s0", 20_000, progress=[(30_000_000, 100)], final=(-1, -1),
              records=20)
    write_manifest(root, [("cfg_040_s0", "cfg_040", 20_000)])
    r = subprocess.run(
        [sys.executable, str(PERLMUTTER / "plan_ref_rerun.py"),
         "--root", str(root), "--write"],
        capture_output=True, text=True,
    )
    assert r.returncode == 0, r.stderr
    written = root / "grid" / "n_events_override.csv"
    assert written.is_file(), "make_grid would never see the quota"
    rows = list(csv.DictReader(open(written)))
    assert rows == [{"cfg_id": "cfg_040", "n_events": str(int(0.9 * 2**31 / 300_000.0))}]


def test_a_dry_run_writes_nothing(tmp_path):
    root = tmp_path / "scan"
    write_run(root, "cfg_040_s0", 20_000, progress=[(30_000_000, 100)], final=(-1, -1),
              records=20)
    write_manifest(root, [("cfg_040_s0", "cfg_040", 20_000)])
    r = subprocess.run(
        [sys.executable, str(PERLMUTTER / "plan_ref_rerun.py"), "--root", str(root)],
        capture_output=True, text=True,
    )
    assert r.returncode == 0, r.stderr
    assert not (root / "grid" / "n_events_override.csv").exists()
    assert "add --write" in r.stdout


def test_a_clean_scan_says_there_is_nothing_to_do(tmp_path):
    root = tmp_path / "scan"
    for s in (0, 1):
        write_run(root, f"cfg_001_s{s}", 20,
                  progress=[(3_000, 800), (9_000, 2_400)], final=(50_000, 40_000),
                  records=20)
    write_manifest(root, [("cfg_001_s0", "cfg_001", 20),
                          ("cfg_001_s1", "cfg_001", 20)])
    r = subprocess.run(
        [sys.executable, str(PERLMUTTER / "plan_ref_rerun.py"), "--root", str(root)],
        capture_output=True, text=True,
    )
    assert r.returncode == 0, r.stderr
    assert "nothing to do" in r.stdout


def test_the_write_selects_every_unusable_configuration_not_just_the_reduced(tmp_path):
    """--only-cfg has to pick up the eleven killed mid-run as well; they keep the
    standard quota, so they appear in no override file and would otherwise be
    silently left out of the batch that is meant to rebuild them."""
    root = tmp_path / "scan"
    # wrapped counter: needs a smaller quota
    write_run(root, "cfg_040_s0", 20_000, progress=[(30_000_000, 100)],
              final=(-1, -1), records=20)
    # truncated n-tuple only: needs a rerun, but no reduction
    write_run(root, "cfg_002_s0", 20_000, progress=[(3_000, 800)], records=10)
    write_manifest(root, [("cfg_040_s0", "cfg_040", 20_000),
                          ("cfg_002_s0", "cfg_002", 20_000)])

    r = subprocess.run(
        [sys.executable, str(PERLMUTTER / "plan_ref_rerun.py"),
         "--root", str(root), "--write"],
        capture_output=True, text=True,
    )
    assert r.returncode == 0, r.stderr

    override = list(csv.DictReader(open(root / "grid" / "n_events_override.csv")))
    assert [row["cfg_id"] for row in override] == ["cfg_040"], "only the quota changes"

    selected = (root / "grid" / "rerun_configs.txt").read_text().split()
    assert sorted(selected) == ["cfg_002", "cfg_040"], selected


def test_the_ids_the_plan_writes_are_the_ids_the_batch_accepts(tmp_path):
    """The two halves of the workflow have to agree on format."""
    root = tmp_path / "scan"
    write_run(root, "cfg_040_s0", 20_000, progress=[(30_000_000, 100)],
              final=(-1, -1), records=20)
    write_manifest(root, [("cfg_040_s0", "cfg_040", 20_000)])
    subprocess.run(
        [sys.executable, str(PERLMUTTER / "plan_ref_rerun.py"),
         "--root", str(root), "--write"],
        capture_output=True, text=True, check=True,
    )

    sys.path.insert(0, str(PERLMUTTER))
    import run_batch  # noqa: E402
    ids = run_batch.parse_only_cfg(f"@{root / 'grid' / 'rerun_configs.txt'}")
    assert ids == {"cfg_040"}
    assert run_batch.parse_only_cfg(
        f"@{root / 'grid' / 'n_events_override.csv'}"
    ) == {"cfg_040"}
