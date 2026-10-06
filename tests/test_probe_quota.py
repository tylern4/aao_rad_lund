"""Sizing the packing probe: what each configuration must be asked for.

``bench_fortran.sbatch`` measures four runs that never contend with each other,
and the node-hours arithmetic divides by that rate.  Production ran 192
processes on the same 128-core node, so the probe exists to measure the
difference -- which means its quotas have to be right before it starts.

Two properties are load-bearing.  The quota must be a *reduction* everywhere,
because ``make_grid.load_overrides`` refuses anything above the standard quota
and a probe run would die before measuring anything.  And trials per event must
be read from the leading run of increasing ``ntries`` samples: 23 configurations
wrap ``integer*4`` partway through, and the last sample -- the one the verifier
uses for a different question -- is a negative number there, so dividing by it
would size those configurations as if they cost nothing and strand them at the
wall limit for the full 30 minutes.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
PERLMUTTER = REPO / "validation" / "perlmutter"
sys.path.insert(0, str(PERLMUTTER))

import make_grid  # noqa: E402
import make_probe_quota as q  # noqa: E402

# The production scan's own final samples: cfg_000 never wraps, cfg_040 goes
# negative at event 16,800 while still counting (its true total is 2.25e9).
CLEAN = """\
ntries, nevent, mcall_max:    16343056         800           1
ntries, nevent, mcall_max:    32362611        1600           1
ntries, nevent, mcall_max:   323626117       20000           1
"""

WRAPPED = """\
ntries, nevent, mcall_max:    105478181         800           1
ntries, nevent, mcall_max:   1951302916       15200           2
ntries, nevent, mcall_max:   2050608146       16000           2
ntries, nevent, mcall_max:  -2138719157       16800           2
ntries, nevent, mcall_max:  -1732787474       20000           2
"""


def fake_scan(root: Path, cfgs: range | list[int], text: str = CLEAN) -> Path:
    """A scan-shaped tree with one seeded run per configuration."""
    for i in cfgs:
        d = root / "fortran" / f"cfg_{i:03d}_s0"
        d.mkdir(parents=True)
        (d / "out.txt").write_text(text)
    return root


# ------------------------------------------------------------ trials per event


def test_trials_per_event_is_the_cumulative_ratio_at_the_end():
    assert q.trials_per_event(CLEAN) == pytest.approx(323_626_117 / 20_000)


def test_trials_per_event_stops_before_the_wrapped_sample():
    """cfg_040's last sample is a negative int*4.  The verifier takes the last
    match because it wants the run's closing count; sizing wants a ratio."""
    assert q.trials_per_event(WRAPPED) == pytest.approx(2_050_608_146 / 16_000)


def test_trials_per_event_stops_when_the_counter_goes_backwards():
    """A wrap need not cross zero to stop being a count."""
    text = "ntries, nevent, mcall_max:  100  10  1\nntries, nevent, mcall_max:  60  20  1\n"
    assert q.trials_per_event(text) == pytest.approx(10.0)


def test_a_file_with_no_usable_sample_is_none_not_zero():
    assert q.trials_per_event("") is None
    assert q.trials_per_event("no counters here\n") is None
    assert q.trials_per_event("ntries, nevent, mcall_max:  -5  10  1\n") is None


# ------------------------------------------------------------------- the quota


def test_the_quota_buys_the_target_seconds_of_work():
    # 20 s x 234,000 trials/s at 16,181 trials per event = 289 events.
    assert q.event_quota(323_626_117 / 20_000) == 289


def test_the_quota_is_always_a_reduction_of_the_standard_quota():
    for tpe in (1.0, 100.0, 16_181.0, 3e5, 1e9):
        assert 1 <= q.event_quota(tpe) <= make_grid.N_EVENTS


def test_an_expensive_configuration_gets_fewer_events_than_a_cheap_one():
    """The whole point of per-configuration sizing: one shared quota would
    either finish delta 0.05 before the clock has resolution or strand delta
    0.005 at the wall limit."""
    cheap, dear = q.event_quota(9_118), q.event_quota(58_566)
    assert cheap > dear


# ------------------------------------------------------------- writing the file


def test_write_quotas_covers_every_configuration(tmp_path: Path):
    scan = fake_scan(tmp_path / "scan", range(96))
    out = tmp_path / "grid" / "n_events_override.csv"
    quotas, rates = q.write_quotas(scan, out)

    assert sorted(quotas) == [f"cfg_{i:03d}" for i in range(96)]
    assert rates["cfg_000"] == pytest.approx(323_626_117 / 20_000)

    rows = list(csv.DictReader(open(out, newline="")))
    assert len(rows) == 96
    assert all(1 <= int(r["n_events"]) <= make_grid.N_EVENTS for r in rows)
    # What make_grid will actually load from this file.
    assert make_grid.load_overrides(out)["cfg_000"] == 289


def test_a_configuration_without_a_rate_aborts_the_probe(tmp_path: Path):
    """Silently defaulting one configuration to 20,000 events would spend the
    whole job waiting for it and report nothing."""
    scan = fake_scan(tmp_path / "scan", [i for i in range(96) if i != 40])
    with pytest.raises(SystemExit) as exc:
        q.write_quotas(scan, tmp_path / "out.csv")
    assert "cfg_040" in str(exc.value)


def test_a_missing_scan_tree_aborts(tmp_path: Path):
    with pytest.raises(SystemExit):
        q.write_quotas(tmp_path / "nope", tmp_path / "out.csv")
