"""Speedup and node-hours: the arithmetic the two benchmarks exist to produce.

The fixture is deliberately hand-checkable -- one trials-per-event for every
configuration, one wall for every arm -- so a regression in how the two sides
are priced shows up as a number rather than as a reshuffling of output.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
PERLMUTTER = REPO / "validation" / "perlmutter"
sys.path.insert(0, str(PERLMUTTER))

import speedup_report as sr  # noqa: E402

TPE = 10_000.0
# The reference bench's fitted rate: the report must read it out of bench8.json
# rather than keep its own copy, so the fixture sets it and the assertions use
# the same value -- if the two were pinned independently, a hard-coded rate in
# the report would still pass here.
REF_RATE = 240_004.0
# 20,000 events at 10,000 trials/event on the reference's fitted 240,004 trials/s.
F_SEC = 20_000 * TPE / REF_RATE
F_NH = F_SEC / 3600 / sr.FORT_PACK


def write_scan(root: Path) -> Path:
    """Production scan: 25 checkpoints per run, stable trials per event."""
    for cfg in sr.CFGS:
        for seed in ("s0", "s1"):
            d = root / "fortran" / f"{cfg}_{seed}"
            d.mkdir(parents=True)
            lines = [
                f" ntries, nevent, mcall_max: {int(TPE * 800 * k)} {800 * k}"
                for k in range(1, 26)
            ]
            (d / "out.txt").write_text("\n".join(lines) + "\n")
    return root


def summary(root: Path, rel: str, rate: float, projs: dict[str, float],
            startup: float = 0.0, node_hours: float | None = None) -> dict:
    """A bench_report summary whose runs obey wall = startup + trials/rate.

    Trial counts vary across runs or the fit's denominator is zero -- every
    combination of startup and rate would then pass through the points.
    """
    runs = []
    for i, (rid, proj) in enumerate(projs.items()):
        trials = 1e6 * (i + 1)
        runs.append({
            "run_id": rid, "proj": proj, "rate": rate, "complete": True,
            "trials": trials, "wall": startup + trials / rate,
        })
    payload = {"median_rate": rate, "runs": runs}
    if node_hours is not None:
        payload["node_hours"] = node_hours
    path = root / rel
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload))
    return payload


def build(tmp_path: Path) -> Path:
    """A bench root whose four arms all parse, plus the scan they are priced on."""
    root = tmp_path / "bench"
    every = {f"{cfg}_{s}" for cfg in sr.CFGS for s in ("s0", "s1")}
    summary(root, "timing/pack128.json", 220_582.0,
            {rid: 3600.0 for rid in sorted(every)}, startup=1.05)
    summary(root, "timing/bench8.json", REF_RATE,
            {rid: 914.0 for rid in sorted(every)}, startup=1.27,
            node_hours=914.0 / 3600 / sr.FORT_PACK)
    # 100 s warm, 200 s cold for every configuration: no per-config noise to
    # disentangle, so the printed speedup is F_SEC / 100 and F_SEC / 200.
    summary(root, "full/timing/pass_warm.json", 7_385_400.0,
            {rid: 100.0 for rid in sorted(every)})
    summary(root, "full/timing/pass_cold.json", 4_237_787.0,
            {rid: 200.0 for rid in sorted(every)})
    return root


def run(tmp_path: Path, capsys) -> str:
    root = build(tmp_path)
    scan = write_scan(tmp_path / "scan")
    assert sr.main(["--root", str(root), "--scan", str(scan)]) == 0
    return capsys.readouterr().out


# ------------------------------------------------------------------- pricing


def test_speedup_is_fortran_seconds_over_gpu_seconds_for_the_same_work(tmp_path, capsys):
    text = run(tmp_path, capsys)
    # Every configuration is identical in the fixture, so the median is exact.
    expected_warm = F_SEC / 100.0
    expected_cold = F_SEC / 200.0
    assert f"warm median {expected_warm:.1f} x" in text
    assert f"cold median {expected_cold:.1f} x" in text


def test_node_hours_divide_by_the_packing_each_arm_runs_at(tmp_path, capsys):
    """A 50x faster run has to beat a 32x packing deficit; the fixture's GPU is
    slower per run, so its node-hours must come out *worse*, not better."""
    text = run(tmp_path, capsys)
    gw_nh = 100.0 / 3600 / sr.GPU_PACK
    assert f"{F_NH:.5f}" in text
    assert f"GPU warm {gw_nh / F_NH:.2f} x Fortran's" in text
    assert gw_nh / F_NH > 1, "fixture should show the GPU using more node-hours"


def test_the_grid_is_priced_from_production_trials_per_event(tmp_path, capsys):
    """The probe's own short runs draw a different sigr_max and read 0.67x the
    work of production's full runs, so the grid total must not use them."""
    text = run(tmp_path, capsys)
    grid_fit = 2 * 20_000 * len(sr.CFGS) * TPE / REF_RATE / 3600 / sr.FORT_PACK
    assert f"production trials/event x fitted rate: {grid_fit:.2f} node-h" in text
    # The probe's own projections are still reported, alongside it.
    probe_sum = 8 * 3600.0
    assert "full 192-run grid (sum):" in text
    assert f"{probe_sum / 3600 / sr.FORT_PACK:.2f} node-h" in text


def test_the_four_bench_configs_are_priced_against_the_grid_median(tmp_path, capsys):
    text = run(tmp_path, capsys)
    assert "uncontended, the 4 bench configs" in text
    assert "not the grid's price" in text
    grid_median = 3600.0 / 3600 / sr.FORT_PACK
    assert f"grid median {grid_median:.5f}" in text
    assert "0.00345" not in text, "the grid median must be measured, not typed in"


def test_the_production_bill_is_compared_against_the_measurement(tmp_path, capsys):
    text = run(tmp_path, capsys)
    assert "30.25 node-h" in text
    assert "x the compute" in text


def test_the_packing_penalty_is_fitted_from_the_runs_rather_than_typed_in(tmp_path, capsys):
    """Startup is ~1 s on runs of ~17 s (probe) and ~53 s (reference), so a raw
    rate ratio charges the probe's startup share to packing.  The section must
    report the fit packing_penalty computes, not a constant copied from it."""
    text = run(tmp_path, capsys)
    assert "wall = startup + trials/rate, fitted over every run" in text
    assert f"reference {REF_RATE:>10,.0f} trials/s  startup +1.27 s" in text
    assert f"probe     {220_582:>10,.0f} trials/s  startup +1.05 s" in text
    assert f"rate ratio {REF_RATE / 220_582:.2f} x" in text
    # Every configuration carries the same rate in the fixture, so the
    # per-configuration spread is a single number rather than a range.
    assert f"over the {len(sr.CFGS)} both timed: median {REF_RATE / 220_582:.2f} x" in text


def test_a_summary_the_fit_cannot_identify_is_refused(tmp_path, capsys):
    """bench_report wrote trials and wall once it grew --json, so a summary
    without them is stale, and pricing the table off a remembered rate instead
    would silently produce numbers nobody measured."""
    root = build(tmp_path)
    path = root / "timing" / "bench8.json"
    stale = json.loads(path.read_text())
    for rec in stale["runs"]:
        rec.pop("trials", None)
        rec.pop("wall", None)
    path.write_text(json.dumps(stale))
    with pytest.raises(SystemExit) as exc:
        sr.main(["--root", str(root), "--scan", str(write_scan(tmp_path / "scan"))])
    assert "bench_report.py --json" in str(exc.value)
    assert "startup + trials/rate" in str(exc.value)


def test_startup_must_differ_between_the_arms_or_the_fit_cannot_see_it():
    """Guards the fixture itself: if both arms were given the same startup the
    test above would still pass while no longer pinning anything down."""
    import packing_penalty as pp
    runs = [{"wall": 1.05 + t / 220_582.0, "trials": t}
            for t in (1e6, 2e6, 3e6, 4e6)]
    rate, startup, _ = pp.fit_rate(runs)
    assert rate == pytest.approx(220_582.0, rel=1e-9)
    assert startup == pytest.approx(1.05, abs=1e-9)


# ------------------------------------------------------------------- helpers


def test_index_and_median_ignore_runs_that_have_no_projection():
    summary = {"runs": [
        {"run_id": "cfg_000_s0", "proj": 10.0},
        {"run_id": "cfg_000_s1", "proj": None},
        {"run_id": "cfg_034_s0", "proj": 40.0},
        {"run_id": "cfg_034_s1", "proj": 20.0},
    ]}
    idx = sr.index(summary)
    assert sr.med_proj(idx, "cfg_000") == 10.0
    assert sr.med_proj(idx, "cfg_034") == 30.0
    assert sr.med_proj(idx, "cfg_999") is None


def test_a_missing_summary_is_named_rather_than_raised(tmp_path):
    with pytest.raises(SystemExit) as exc:
        sr.load(tmp_path / "nope.json")
    assert "nope.json" in str(exc.value)
    assert "bench_report" in str(exc.value)
