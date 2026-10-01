"""Tests for the configuration object: validation, round-tripping and presets.

Validation happens in ``__post_init__``, so *constructing* a config is the
validation -- every one of these tests is really "this configuration is
accepted" or "this one is refused, and says why".
"""

from __future__ import annotations

import json
import math

import numpy as np
import pytest

from aao_rad.config import EXPERIMENTS, GeneratorConfig, ep_window_for_w
from aao_rad.constants import M_N, TABLE_Q2_MAX, TABLE_W_MAX, TABLE_W_MIN


def _replace(cfg: GeneratorConfig, **kw) -> GeneratorConfig:
    return GeneratorConfig(**{**cfg.to_dict(), **kw})


class TestValidation:
    def test_defaults_are_usable(self):
        GeneratorConfig()  # must not raise

    @pytest.mark.parametrize(
        ("kw", "match"),
        [
            ({"min_photon_energy": 0.0}, "delta"),
            ({"min_photon_energy": -0.1}, "delta"),
            ({"q2_min": 0.0}, "q2_min"),
            ({"q2_max": 0.0}, "q2_min"),
            ({"q2_min": 2.0, "q2_max": 1.0}, "q2_min"),
            ({"ep_min": 0.0}, "ep_min"),
            ({"ep_max": 0.0}, "ep_min"),
            ({"ep_min": 2.0, "ep_max": 1.0}, "ep_min"),
            ({"batch_size": 512}, "batch_size"),
            ({"regions": (0.3, 0.3, 0.3, 0.3)}, "region sizes"),
            ({"n_tracks": 3}, "n_tracks"),
            ({"theory": 5}, "theory=5"),
            ({"channel": 4}, "channel=4"),
            ({"channel": 2}, "channel=2"),
            ({"w_max": 2.5}, "outside the MAID07 table"),
        ],
    )
    def test_rejects_nonsense(self, kw, match):
        with pytest.raises(ValueError, match=match):
            _replace(GeneratorConfig(), **kw)

    def test_accepts_the_table_edges(self):
        for w in (TABLE_W_MIN, TABLE_W_MAX):
            _replace(GeneratorConfig(), w_max=w)

    def test_clamp_disables_the_cut(self):
        _replace(GeneratorConfig(), w_max="clamp")

    def test_frozen(self):
        cfg = GeneratorConfig()
        with pytest.raises(Exception):
            cfg.beam_energy = 9.0  # type: ignore[misc]

    def test_ek_sampling_choices(self):
        assert _replace(GeneratorConfig(), ek_sampling="fortran").ek_sampling == "fortran"
        with pytest.raises(ValueError):
            _replace(GeneratorConfig(), ek_sampling="nonsense")

    def test_interp_scheme_choices(self):
        for s in ("linear", "spline"):
            assert _replace(GeneratorConfig(), interp_scheme=s).interp_scheme == s
        with pytest.raises(ValueError):
            _replace(GeneratorConfig(), interp_scheme="cubic")


class TestRoundTrip:
    def test_to_dict_from_dict_is_the_identity(self):
        cfg = GeneratorConfig(seed=1234, n_events=99, channel=1)
        assert GeneratorConfig.from_dict(cfg.to_dict()) == cfg

    def test_survives_a_json_round_trip(self, tmp_path):
        cfg = GeneratorConfig(seed=7, regions=(0.1, 0.15, 0.2, 0.25))
        path = tmp_path / "cfg.json"
        path.write_text(json.dumps(cfg.to_dict()))
        again = GeneratorConfig.from_dict(json.loads(path.read_text()))
        assert again == cfg
        # The regions must come back as a tuple of floats, not a list of ints.
        assert again.regions == (0.1, 0.15, 0.2, 0.25)
        assert all(isinstance(x, float) for x in again.regions)

    def test_from_file_reads_json(self, tmp_path):
        cfg = GeneratorConfig(seed=11, n_events=5)
        path = tmp_path / "cfg.json"
        path.write_text(json.dumps(cfg.to_dict()))
        assert GeneratorConfig.from_file(path) == cfg

    def test_unknown_keys_are_rejected(self):
        with pytest.raises(ValueError, match="unknown configuration keys"):
            GeneratorConfig.from_dict({"beam_energy": 4.0, "nonsense": 1})

    def test_clam_mode_survives_the_round_trip(self):
        cfg = _replace(GeneratorConfig(), w_max="clamp")
        assert GeneratorConfig.from_dict(cfg.to_dict()).w_max == "clamp"

    def test_merge_ignores_none(self):
        cfg = GeneratorConfig(seed=1)
        assert cfg.merge(seed=None, n_events=5).n_events == 5
        assert cfg.merge(seed=None).seed == 1


class TestPresets:
    @pytest.mark.parametrize("name", sorted(EXPERIMENTS))
    def test_each_preset_is_constructible(self, name):
        GeneratorConfig(**EXPERIMENTS[name])

    @pytest.mark.parametrize("name", sorted(EXPERIMENTS))
    def test_presets_keep_the_sample_inside_the_table(self, name):
        """The shipped presets must not push most of the sample past the edge
        of the MAID07 grid, which is what the old driver's windows did."""
        import sys

        import jax
        import numpy as np

        import aao_rad.generate  # noqa: F401  (module shadows a function)
        gen = sys.modules["aao_rad.generate"]

        cfg = GeneratorConfig(**EXPERIMENTS[name])
        kin = gen.build_kinematics(cfg)
        kv = gen.draw_kinematics(jax.random.PRNGKey(0), kin, 200_000)
        ok = np.asarray(kv["ok"])
        assert ok.mean() > 0.95, f"{name}: only {ok.mean():.1%} of trials are legal"
        # epw is the true (pre-exit) W, which is what the w_cut applies to.
        w = np.asarray(kv["epw"])[ok]
        q2 = np.asarray(kv["q2"])[ok]
        assert np.mean((w >= TABLE_W_MIN) & (w <= TABLE_W_MAX)) > 0.95, name
        assert np.mean((q2 >= 0.0) & (q2 <= TABLE_Q2_MAX)) > 0.99, name

    @pytest.mark.parametrize("name", sorted(EXPERIMENTS))
    def test_presets_use_supported_channels(self, name):
        assert EXPERIMENTS[name]["channel"] in (1, 3)

    def test_region_sizes_sum_below_the_limit(self):
        for name, preset in EXPERIMENTS.items():
            assert sum(preset["regions"]) < 1.0, name

    @pytest.mark.parametrize("name", sorted(EXPERIMENTS))
    def test_presets_reach_their_own_q2_window(self, name):
        """A preset that cannot produce its requested Q^2 window is useless."""
        import sys

        import jax
        import numpy as np

        import aao_rad.generate  # noqa: F401
        gen = sys.modules["aao_rad.generate"]

        cfg = GeneratorConfig(**EXPERIMENTS[name])
        kin = gen.build_kinematics(cfg)
        kv = gen.draw_kinematics(jax.random.PRNGKey(1), kin, 200_000)
        q2 = np.asarray(kv["q2"])[np.asarray(kv["ok"])]
        assert np.mean((q2 >= cfg.q2_min) & (q2 <= cfg.q2_max)) > 0.9, name


class TestLegacyInputCards:
    """The old driver read a positional run card; that path must keep working."""

    CARD = """7
1
.20 .12 .20 .20
4
3
.2
5.0
0.486
0.3
0.03
0.0
4.244
0.2 1.9
1.6 2.9
.005
100
0.5
"""

    def test_parses(self):
        cfg = GeneratorConfig.from_legacy_input(self.CARD)
        assert cfg.theory == 7
        assert cfg.polarized_beam is True
        assert cfg.channel == 3
        assert cfg.regions == (0.20, 0.12, 0.20, 0.20)
        assert cfg.n_tracks == 4
        assert cfg.missing_mass_cut == pytest.approx(0.2)
        assert cfg.target_length_cm == pytest.approx(5.0)
        assert cfg.target_radius_cm == pytest.approx(0.486)
        assert cfg.beam_x_cm == pytest.approx(0.3)
        assert cfg.beam_y_cm == pytest.approx(0.03)
        assert cfg.beam_z_cm == pytest.approx(0.0)
        assert cfg.beam_energy == pytest.approx(4.244)
        assert cfg.q2_min == pytest.approx(0.2)
        assert cfg.q2_max == pytest.approx(1.9)
        assert cfg.ep_min == pytest.approx(1.6)
        assert cfg.ep_max == pytest.approx(2.9)
        assert cfg.min_photon_energy == pytest.approx(0.005)
        assert cfg.n_events == 100

    def test_comments_and_blank_lines_are_ignored(self):
        card = "! leading comment\n\n" + self.CARD.replace("4.244", "4.244  ! beam")
        assert GeneratorConfig.from_legacy_input(card) == GeneratorConfig.from_legacy_input(
            self.CARD
        )

    def test_ntracks_3_maps_to_2(self):
        cfg = GeneratorConfig.from_legacy_input(self.CARD.replace("\n4\n", "\n3\n"))
        assert cfg.n_tracks == 2

    def test_a_truncated_card_is_rejected(self):
        with pytest.raises(ValueError, match="at least 21 values"):
            GeneratorConfig.from_legacy_input("7\n1\n.2 .1 .2 .2\n")

    def test_shipped_run_card_is_refused_with_a_clear_reason(self):
        """``test.inp`` asks for MAID2000, which this port does not implement.

        The original could not run it either -- ``maid.F``'s MAID2000 path is
        dead in the shipped tree -- so refusing loudly beats producing nothing.
        """
        from pathlib import Path

        path = Path(__file__).resolve().parent.parent / "test.inp"
        if not path.is_file():
            pytest.skip("test.inp not present")
        with pytest.raises(ValueError, match="theory=5 is not available"):
            GeneratorConfig.from_legacy_input(path.read_text())

    def test_the_same_card_works_with_theory_7(self, tmp_path):
        from pathlib import Path

        path = Path(__file__).resolve().parent.parent / "test.inp"
        if not path.is_file():
            pytest.skip("test.inp not present")
        text = path.read_text().replace("5  ", "7  ", 1)
        cfg = GeneratorConfig.from_legacy_input(text)
        assert cfg.theory == 7
        assert cfg.channel == 1
        assert cfg.n_events == 500
        assert math.isfinite(cfg.beam_energy)
        # ...and the W window has to be pulled inside the table for it to run.
        cfg = cfg.merge(ep_min=1.7, ep_max=3.6)
        assert cfg.ep_max < cfg.beam_energy

class TestEpWindowForW:
    """The ``(Q^2, E')`` rectangle maps to a diagonal band in ``W``."""

    def test_reproduces_the_banana(self):
        """The returned window must keep every corner of the rectangle inside."""
        rng = np.random.default_rng(0)
        for beam, q0, q1 in [(4.244, 0.2, 1.9), (4.8, 0.9, 2.5), (5.0, 0.5, 1.2),
                             (6.0, 1.0, 3.0)]:
            lo, hi = ep_window_for_w(beam, q0, q1)
            q2 = rng.uniform(q0, q1, 500)
            ep = rng.uniform(lo, hi, 500)
            w_sq = M_N**2 + 2 * M_N * (beam - ep) - q2
            w = np.sqrt(np.maximum(w_sq, 0.0))
            assert w.min() >= TABLE_W_MIN - 1e-6, (beam, q0, q1, w.min())
            assert w.max() <= TABLE_W_MAX + 1e-6, (beam, q0, q1, w.max())

    def test_the_edges_are_tight(self):
        """Widening the window by any amount at either end must push W out."""
        beam, q0, q1 = 4.244, 0.2, 1.9
        lo, hi = ep_window_for_w(beam, q0, q1)
        w_of = lambda ep, q2: np.sqrt(M_N**2 + 2 * M_N * (beam - ep) - q2)  # noqa: E731
        assert w_of(lo - 1e-3, q0) > TABLE_W_MAX
        assert w_of(hi + 1e-3, q1) < TABLE_W_MIN

    def test_a_too_wide_q2_window_is_reported(self):
        """A Q^2 window wide enough leaves no E' window at all."""
        with pytest.raises(ValueError, match="narrow the Q"):
            ep_window_for_w(4.8, 0.9, 4.5)

    def test_a_wide_but_feasible_window_is_merely_narrow(self):
        """Q^2 up to 3.5 still works for a 4.8 GeV beam, but only just."""
        lo, hi = ep_window_for_w(4.8, 0.9, 3.5)
        assert 0.0 < hi - lo < 0.15
        # This is why the "default" preset caps Q^2 at 2.5: at 3.5 the sample
        # would be squeezed into a 0.1 GeV E' window.
        wide_lo, wide_hi = ep_window_for_w(4.8, 0.9, 2.5)
        assert wide_hi - wide_lo > hi - lo

    def test_an_inverted_window_is_rejected(self):
        with pytest.raises(ValueError):
            ep_window_for_w(4.244, 0.2, 1.9, w_min=1.9, w_max=1.95)

    @pytest.mark.parametrize("name", sorted(EXPERIMENTS))
    def test_presets_use_a_window_from_this_helper(self, name):
        """Each preset must at least sit inside the in-table rectangle."""
        p = EXPERIMENTS[name]
        lo, hi = ep_window_for_w(p["beam_energy"], p["q2_min"], p["q2_max"])
        assert lo - 1e-6 <= p["ep_min"]
        assert p["ep_max"] <= hi + 1e-6
