"""Measuring the packing penalty, and the two ways it can be measured wrong.

The penalty is the ratio of one machine's throughput to another's, so the
estimator has to ignore everything that is not throughput.  Two candidates
fail: the projections bench_report prints, because they scale by trials per
event and ``sigr_max`` re-estimates that per run; and a ratio of medians over
the two arms, because one times 96 configurations and the other four, and
trials per event span 400x across the grid.  What is left is trials per second,
with the fixed process startup fitted out -- the probe's runs are ~17 s against
the reference's ~53 s, so ~1 s of startup is a larger share of the short runs
and would otherwise be charged to packing.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
PERLMUTTER = REPO / "validation" / "perlmutter"
sys.path.insert(0, str(PERLMUTTER))

import packing_penalty as pp  # noqa: E402

RATE, STARTUP = 230_000.0, 1.25


def runs(trials: list[float], rate: float = RATE, startup: float = STARTUP) -> list[dict]:
    """Runs that obey wall = startup + trials / rate exactly."""
    return [
        {"run_id": f"cfg_{i:03d}_s0", "complete": True, "trials": t,
         "wall": startup + t / rate, "rate": t / (startup + t / rate)}
        for i, t in enumerate(trials)
    ]


# ------------------------------------------------------------------- the fit


def test_the_fit_recovers_the_rate_and_the_startup_it_was_built_with():
    rate, startup, n = pp.fit_rate(runs([1e6, 2e6, 5e6, 1.1e7]))
    assert n == 4
    assert rate == pytest.approx(RATE, rel=1e-9)
    assert startup == pytest.approx(STARTUP, rel=1e-9)


def test_startup_is_not_clamped_to_zero():
    """A negative intercept is information -- the runs share work that does not
    scale with trials.  Reporting it as zero would present a fit that does not
    describe the data as though it did."""
    pts = [{"trials": t, "wall": -0.5 + t / RATE} for t in (1e6, 4e6, 9e6)]
    _, startup, _ = pp.fit_rate(pts)
    assert startup == pytest.approx(-0.5, rel=1e-9)


def test_a_single_run_does_not_identify_a_rate():
    assert pp.fit_rate(runs([3e6])[:1]) == (None, None, 1)
    assert pp.fit_rate([]) == (None, None, 0)


def test_trial_counts_with_no_spread_do_not_identify_a_rate():
    """Every run doing the same work leaves the slope denominator zero: any
    startup and any rate would pass through the points."""
    assert pp.fit_rate(runs([4e6, 4e6, 4e6]))[0] is None


def test_runs_without_wall_or_trials_are_skipped_rather_than_crashing():
    broken = [{"run_id": "cfg_000_s0", "complete": False}, {"wall": None, "trials": 1e6}]
    rate, _, n = pp.fit_rate(runs([1e6, 3e6]) + broken)
    assert n == 2
    assert rate == pytest.approx(RATE, rel=1e-9)


# ------------------------------------------------------- configuration pairing


def test_rates_pair_by_configuration_not_by_position():
    """cfg_000_s0 in one arm must meet cfg_000_s1 in the other, or the ratio is
    between different physics."""
    probe = [
        {"run_id": "cfg_000_s0", "complete": True, "rate": 200_000.0},
        {"run_id": "cfg_000_s1", "complete": True, "rate": 240_000.0},
        {"run_id": "cfg_034_s0", "complete": True, "rate": 190_000.0},
    ]
    got = pp.rate_by_config(probe)
    assert got == {"cfg_000": 220_000.0, "cfg_034": 190_000.0}


def test_an_incomplete_run_contributes_nothing():
    assert pp.rate_by_config([{"run_id": "cfg_000_s0", "complete": False, "rate": 1.0}]) == {}


# ------------------------------------------------------------------ comparing


def test_identical_arms_report_no_penalty(capsys):
    summary = {"runs": runs([1e6, 4e6, 9e6]), "node_hours": 0.002}
    out = pp.compare(summary, json.loads(json.dumps(summary)))
    text = capsys.readouterr().out

    assert out["rate_ratio"] == pytest.approx(1.0, rel=1e-9)
    assert out["per_config_median"] == pytest.approx(1.0, rel=1e-9)
    assert "rate ratio 1.00 x" in text


def test_the_penalty_is_the_reference_rate_over_the_probe_rate(capsys):
    """The reference is one idle core; a slower probe means packing cost
    something, so the ratio has to be reference/probe and not its inverse."""
    slow = [{"run_id": f"cfg_{i:03d}_s0", "complete": True, "trials": t,
             "wall": 0.0 + t / (RATE / 1.5), "rate": RATE / 1.5}
            for i, t in enumerate([1e6, 4e6, 9e6])]
    out = pp.compare({"runs": slow, "node_hours": 0.003}, {"runs": runs([1e6, 4e6, 9e6]),
                                                           "node_hours": 0.002})
    assert out["rate_ratio"] == pytest.approx(1.5, rel=1e-6)
    assert "rate ratio 1.50 x" in capsys.readouterr().out


def test_the_node_hours_are_labelled_as_measured_and_as_extrapolated(capsys):
    """The probe prices all 96 configurations at 128 in flight; the reference
    only ever ran four at once, so its figure is an assumption."""
    pp.compare({"runs": runs([1e6]), "node_hours": 0.00345},
               {"runs": runs([1e6]), "node_hours": 0.001984})
    text = capsys.readouterr().out
    assert "measured over 1 configurations: 0.003450" in text
    assert "0.001984  (extrapolated to 128 in flight)" in text


# ----------------------------------------------------------------------- main


def test_a_missing_summary_fails_instead_of_printing_a_comparison(tmp_path: Path, capsys):
    probe = tmp_path / "probe.json"
    probe.write_text(json.dumps({"runs": []}))
    assert pp.main([str(probe), str(tmp_path / "absent.json")]) == 1
    assert "no summary" in capsys.readouterr().err


def test_main_reads_two_summaries(tmp_path: Path, capsys):
    payload = {"runs": runs([1e6, 4e6, 9e6]), "node_hours": 0.002}
    a, b = tmp_path / "a.json", tmp_path / "b.json"
    a.write_text(json.dumps(payload))
    b.write_text(json.dumps(payload))
    assert pp.main([str(a), str(b)]) == 0
    assert "rate ratio 1.00 x" in capsys.readouterr().out
