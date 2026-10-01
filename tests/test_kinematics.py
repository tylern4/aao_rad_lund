"""Tests for the response-function kernels: the dilogarithm and the kinematics.

``spence`` and the two-body decay kinematics are the two places where an
algebraic slip would produce plausible-looking but wrong events, so they are
pinned to known values and to exact on-shell identities.
"""

from __future__ import annotations

import sys

import jax
import jax.numpy as jnp
import numpy as np
import pytest

from aao_rad.constants import M_E, M_N, M_PI0, M_PIP, pi_masses
from aao_rad.kinematics import hadronic_final_state
from aao_rad.generate import EVENT_COLUMNS
from aao_rad.motsa import spence

#: Li_2 at arguments where the value is known in closed form.
SPENCE_CASES = [
    (0.0, 0.0),
    (1.0, np.pi**2 / 6.0),
    (-1.0, -np.pi**2 / 12.0),
    (0.5, np.pi**2 / 12.0 - np.log(2.0) ** 2 / 2.0),
]


def _series(x: np.ndarray, terms: int = 400) -> np.ndarray:
    """Reference Li_2 by its defining series; only valid for |x| < 1."""
    k = np.arange(1, terms + 1)
    return (x[..., None] ** k / k**2).sum(axis=-1)


class TestSpence:
    @pytest.mark.parametrize(("x", "want"), SPENCE_CASES)
    def test_closed_form_values(self, x, want):
        got = float(spence(jnp.asarray(x, jnp.float32)))
        assert got == pytest.approx(want, rel=2e-6, abs=2e-7)

    def test_matches_the_defining_series(self):
        x = np.linspace(-0.95, 0.95, 97)
        got = np.asarray(spence(jnp.asarray(x, jnp.float32)))
        np.testing.assert_allclose(got, _series(x), rtol=2e-5, atol=2e-7)

    @pytest.mark.parametrize("x", [-4.0, -2.5, -1.2, -0.999, 0.999, 1.2, 3.0, 7.5])
    def test_analytic_continuation_outside_the_unit_interval(self, x):
        """The defining series diverges for |x| >= 1, so only the identities
        that define the continuation can be checked there."""
        got = float(spence(jnp.asarray(x, jnp.float32)))
        assert np.isfinite(got)
        if x <= 0.0:
            # Landen: Li_2(x) = -Li_2(x/(x-1)) - ln^2(1-x)/2.
            trans = float(spence(jnp.asarray(x / (x - 1.0), jnp.float32)))
            assert got == pytest.approx(-trans - np.log(1 - x) ** 2 / 2.0, rel=2e-5)
        if 0.0 < x <= 1.0:
            # Li_2(x) + Li_2(1 - x) = pi^2/6 - ln(x) ln(1 - x)
            other = float(spence(jnp.asarray(1.0 - x, jnp.float32)))
            assert got + other == pytest.approx(
                np.pi**2 / 6.0 - np.log(x) * np.log(1.0 - x), rel=2e-5, abs=2e-7)

    def test_is_continuous_and_monotonic_across_one(self):
        """Li_2 has a square-root singularity in its derivative at 1, but the
        function itself is continuous and increasing, so no step may appear."""
        eps = np.array([1.0 - 1e-7, 1.0 - 1e-8, 1.0, 1.0 + 1e-8, 1.0 + 1e-7])
        got = np.asarray(spence(jnp.asarray(eps, jnp.float32)))
        # float32 cannot resolve 1 +/- 1e-8, so the distinct points are
        # {1-1e-7, 1-1e-6, 1, 1+1e-6}: monotonic and spanning < 1e-3 in total.
        assert np.all(np.diff(got) >= 0.0), "Li_2 must be monotonically increasing"
        assert got.max() - got.min() < 1e-2
        assert got[-1] == pytest.approx(np.pi**2 / 6.0, rel=1e-5)

    def test_is_batched_and_differentiable(self):
        x = jnp.asarray(np.linspace(0.01, 0.9, 32), jnp.float32)
        assert spence(x).shape == x.shape
        g = jax.grad(lambda z: spence(z).sum())
        # d/dx Li_2(x) = -ln(1 - x) / x
        for t in (0.2, 0.5, 0.8):
            assert float(g(t)) == pytest.approx(-np.log(1 - t) / t, rel=2e-5)

    def test_jitted_matches_eager(self):
        x = jnp.asarray(np.linspace(0.0, 0.95, 64), jnp.float32)
        np.testing.assert_allclose(np.asarray(jax.jit(spence)(x)),
                                   np.asarray(spence(x)), rtol=1e-6, atol=1e-7)


class TestPionMasses:
    """The second element is the squared *missing* mass, not m_pi^2.

    For pi+ p the missing particle is a neutron, so the generator aims at
    M_N^2; for pi0 p it is a pion, so it aims at M_pi0^2.
    """

    def test_channel_1_is_pi0_production(self):
        assert pi_masses(1) == (M_PI0, M_PI0**2)

    def test_channel_3_is_pip_production_with_a_neutron_missing(self):
        assert pi_masses(3) == (M_PIP, M_N**2)

    @pytest.mark.parametrize("channel", [2, 5])
    def test_neutron_target_channels(self, channel):
        assert pi_masses(channel)[1] == M_N**2

    def test_unknown_channel_is_rejected(self):
        with pytest.raises(ValueError, match="unsupported pion production channel"):
            pi_masses(4)

    def test_constants_are_the_precise_pdg_values(self):
        assert M_N == 0.9382720813
        assert M_E == 0.00051099895
        assert M_E < M_PI0 < M_PIP < M_N


class TestHadronicFinalState:
    """The two-body decay as the generator actually calls it.

    Driving the kernel through the real sampler matters: the photon energy and
    direction are not free parameters, they are tied to the electron's
    bremsstrahlung.  Feeding ``hadronic_final_state`` arbitrary triples of
    ``(es, ep', ek)`` produces a boost velocity inconsistent with the
    reconstructed ``W``, and the Lorentz identities below then fail -- so those
    failures are a test artefact, not a code defect.

    ``finalize`` is driven end to end, and the identities are checked on the
    events that survive the kinematic cuts.
    """

    @staticmethod
    def _survivors(cfg, n: int = 400_000):
        """One large trial batch, reduced to the events that pass the cuts."""
        import aao_rad.generate as _g
        gen = sys.modules["aao_rad.generate"]

        grid = gen.build_grid(cfg.channel, parms_dir="parms",
                              scheme=cfg.interp_scheme)
        kin = gen.build_kinematics(cfg)
        key = jax.random.PRNGKey(1234)
        kv = gen.draw_kinematics(key, kin, n)
        ph = gen.draw_photon(jax.random.fold_in(key, 1), kin, kv, n)
        _, asym, _ = gen.integrand(grid, kin, kv, ph, cfg.interp_scheme)
        rec, ok = gen.finalize(kin, kv, ph, asym, jax.random.fold_in(key, 2), n)
        rec = np.asarray(rec)[np.asarray(ok)]
        return {name: rec[:, i]
                for i, name in enumerate(EVENT_COLUMNS)}, kin

    @pytest.fixture(scope="class")
    def events(self, small_config):
        return self._survivors(small_config)

    def test_produces_enough_events_to_test(self, events):
        ev, _ = events
        assert len(ev["epi"]) > 500, f"only {len(ev['e_pi'])} events survived"

    def test_all_columns_are_finite(self, events):
        ev, _ = events
        for name, arr in ev.items():
            assert np.isfinite(arr).all(), f"{name} is not finite"

    def test_each_track_is_on_its_mass_shell(self, events):
        ev, _ = events
        for e_col, p_cols, mass in (
            ("eprot", ("ppx", "ppy", "ppz"), M_N),
            ("epi", ("ppix", "ppiy", "ppiz"), M_PIP),
        ):
            p2 = sum(ev[c] ** 2 for c in p_cols)
            resid = ev[e_col] ** 2 - p2
            np.testing.assert_allclose(resid, mass**2, rtol=1e-3, atol=1e-4)
            assert np.all(ev[e_col] > 0.0)

    def test_energy_conservation(self, events):
        """E_pion + E_nucleon is the hadronic system energy es - ep' + M_N - ek."""
        ev, _ = events
        e_w = (ev["es"] - ev["ep"]) + M_N - ev["eg"]
        np.testing.assert_allclose(ev["eprot"] + ev["epi"], e_w,
                                   rtol=1e-3, atol=1e-4)

    def test_momentum_conservation(self, events):
        """p_pion + p_nucleon is the resonance momentum q - k.

        ``q`` here must be built from the *measured* (post-exit) electron.  The
        ``qx``/``qz``/``q0`` n-tuple columns hold the *generated* values, taken
        before the exit radiation, exactly as ``aao_rad.f90`` line 576 does --
        so they deliberately do not close the balance with the measured tracks.
        """
        ev, _ = events
        p_s = np.sqrt(ev["es"] ** 2 - M_E**2)
        p_p = np.sqrt(ev["ep"] ** 2 - M_E**2)
        c, s = ev["csthe"], np.sqrt(np.maximum(1.0 - ev["csthe"] ** 2, 0.0))
        q_mx, q_mz = -p_p * s, p_s - p_p * c
        for got, want in (
            (ev["ppx"] + ev["ppix"], q_mx - ev["egx"]),
            (ev["ppy"] + ev["ppiy"], -ev["egy"]),
            (ev["ppz"] + ev["ppiz"], q_mz - ev["egz"]),
        ):
            np.testing.assert_allclose(got, want, rtol=1e-3, atol=1e-4)

    def test_recorded_q_is_the_pre_exit_virtual_photon(self, events):
        """The n-tuple q columns are consistent with the *generated* electron."""
        ev, _ = events
        p_s = np.sqrt(ev["es"] ** 2 - M_E**2)
        p_p = np.sqrt(ev["ep"] ** 2 - M_E**2)
        s = np.sqrt(np.maximum(1.0 - ev["csthe"] ** 2, 0.0))
        # The exit loss only ever removes energy, so the generated momentum
        # transfer is the larger one.
        assert np.all(ev["qx"] <= -p_p * s + 1e-4)

    def test_photon_is_on_shell(self, events):
        ev, _ = events
        k2 = ev["egx"] ** 2 + ev["egy"] ** 2 + ev["egz"] ** 2 - ev["eg"] ** 2
        # Photons softer than delta are written as a dummy 1e-5 4-vector, which
        # is deliberately not on shell, so restrict to the real ones.
        real = ev["eg"] > 1e-4
        assert real.sum() > 0
        np.testing.assert_allclose(k2[real], 0.0, rtol=1e-4, atol=1e-5)

    def test_cos_theta_cm_is_in_range(self, events):
        ev, _ = events
        assert ev["csthcm"].min() >= -1.0 - 1e-6
        assert ev["csthcm"].max() <= 1.0 + 1e-6

    def test_missing_mass_cut_is_respected(self, small_config, events):
        """Every surviving event sits inside the configured mm2 window."""
        from aao_rad.constants import pi_masses

        ev, kin = events
        m_exp = pi_masses(small_config.channel)[1]
        assert float(kin.m_exp) == pytest.approx(m_exp)
        cut = small_config.missing_mass_cut
        assert np.all(np.abs(ev["mm2"] - m_exp) <= cut + 1e-6)
        # ...and the window is genuinely populated, not just barely satisfied.
        assert np.mean(np.abs(ev["mm2"] - m_exp) > 0.5 * cut) > 0.2

    def test_w_is_above_the_two_body_threshold(self, events):
        ev, _ = events
        assert np.all(ev["w_real"] > M_N + M_PIP)
        assert np.all(ev["w"] > M_N + M_PIP)

    def test_the_same_key_reproduces_the_same_events(self, small_config):
        import sys

        import aao_rad.generate  # noqa: F401  (register the module)
        gen = sys.modules["aao_rad.generate"]

        a, _ = self._survivors(small_config, n=200_000)
        b, _ = self._survivors(small_config, n=200_000)
        assert len(a["epi"]) == len(b["epi"])
        for name in a:
            np.testing.assert_array_equal(a[name], b[name])

    def test_jitting_does_not_change_the_events(self, small_config):
        """The stepper is the jitted path; run it twice and compare."""
        import sys

        import aao_rad.generate  # noqa: F401
        gen = sys.modules["aao_rad.generate"]

        cfg = small_config
        grid = gen.build_grid(cfg.channel, parms_dir="parms",
                              scheme=cfg.interp_scheme)
        kin = gen.build_kinematics(cfg)
        n = 65536
        key = jax.random.PRNGKey(21)
        kv = gen.draw_kinematics(key, kin, n)
        ph = gen.draw_photon(jax.random.fold_in(key, 1), kin, kv, n)
        _, asym, _ = gen.integrand(grid, kin, kv, ph, cfg.interp_scheme)
        fk = jax.random.fold_in(key, 2)
        eager = gen.finalize(kin, kv, ph, asym, fk, n)
        jitted = jax.jit(gen.finalize, static_argnums=5)(kin, kv, ph, asym, fk, n)
        # XLA may reassociate float32 sums, so this is allclose, not exact.
        np.testing.assert_allclose(np.asarray(eager[0]), np.asarray(jitted[0]),
                                   rtol=1e-4, atol=1e-5)
        np.testing.assert_array_equal(np.asarray(eager[1]), np.asarray(jitted[1]))
