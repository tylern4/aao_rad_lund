"""The trial budget must be a knob, not a constant.

Acceptance is ``mean(weight) / ceiling``, and with a heavy-tailed integrand that
lands far below the original's few tenths of a percent.  On the validation grid
the ``pi0``-channel 8-12 GeV configurations needed 207,000-231,000 trials per
event, so the original hardcoded ``200_000 * n_events`` guard cut them off a few
percent short of their quota and they raised instead of finishing.

These tests pin the two halves of that: the guard still fires on a genuinely
empty phase space (so a broken config cannot silently emit a short run), and
raising the budget lets a slow-but-legitimate config through.
"""

from __future__ import annotations

import numpy as np
import pytest
from conftest import requires_tables

from aao_rad.config import GeneratorConfig
from aao_rad.generate import EventGenerator, build_grid


@pytest.fixture(scope="module")
def generator(parms_dir):
    return EventGenerator(build_grid(3, parms_dir=parms_dir))


def test_max_trials_per_event_defaults_to_the_original_value():
    assert GeneratorConfig(channel=3).max_trials_per_event == 200_000


def test_max_trials_per_event_must_be_positive():
    with pytest.raises(ValueError, match="max_trials_per_event"):
        GeneratorConfig(channel=3, max_trials_per_event=0)


def test_cli_exposes_the_budget():
    """The driver passes it explicitly, so a run is reproducible from its cmd."""
    from aao_rad import cli

    parser = cli.build_parser()
    args = parser.parse_args(["--channel", "3", "--max-trials-per-event", "600000"])
    assert args.max_trials_per_event == 600_000
    # And it is left alone when not given, rather than pinned to None.
    assert parser.parse_args(["--channel", "3"]).max_trials_per_event is None


def test_an_unreachable_budget_still_raises_rather_than_emitting_a_short_run(generator):
    """The guard's real job: refuse to pass off a truncated run as a complete one."""
    cfg = GeneratorConfig(
        channel=3,
        beam_energy=4.244,
        q2_min=0.2,
        q2_max=1.9,
        ep_min=1.6,
        ep_max=2.9,
        n_events=100_000,
        seed=1,
        batch_size=4096,
        max_trials_per_event=2,
    )
    with pytest.raises(RuntimeError, match="acceptance collapsed"):
        list(generator.stream(cfg, chunk_events=4096, n_events=100_000))


@requires_tables
def test_a_raised_budget_completes_what_the_default_would_abort(generator, small_config):
    """Same run, same seed: the budget decides whether it completes."""
    import dataclasses

    tight = dataclasses.replace(
        small_config, n_events=4000, batch_size=4096, max_trials_per_event=2
    )
    with pytest.raises(RuntimeError, match="acceptance collapsed"):
        list(generator.stream(tight, chunk_events=4096, n_events=4000))

    roomy = dataclasses.replace(tight, max_trials_per_event=10_000_000)
    events = sum(
        b.size for b, _ in generator.stream(roomy, chunk_events=4096, n_events=4000)
    )
    assert events == 4000


@requires_tables
def test_the_error_names_the_likely_cause_and_the_way_out(generator, small_config):
    """An operator hitting this needs to know both readings: empty cuts, or a
    budget that is merely too tight for a slow configuration."""
    import dataclasses

    cfg = dataclasses.replace(
        small_config, n_events=4000, batch_size=4096, max_trials_per_event=2
    )
    with pytest.raises(RuntimeError) as exc:
        list(generator.stream(cfg, chunk_events=4096, n_events=4000))
    msg = str(exc.value)
    assert "empty" in msg
    assert "max-trials-per-event" in msg
    # The progress so far is reported, so it is clear this was not a hard zero.
    assert "events from" in msg and "trials per event" in msg
    assert np.isfinite(float(msg.split("trials")[0].split()[-1]))
