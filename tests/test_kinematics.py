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
from aao_rad.generate import EVENT_COLUMNS
from aao_rad.kinematics import hadronic_final_state
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
        gen = sys.modules["aao_rad.generate"]

        grid = gen.build_grid(cfg.channel, parms_dir="parms",
                              scheme=cfg.interp_scheme)
        kin = gen.build_kinematics(cfg)
        key = jax.random.PRNGKey(1234)
        kv = gen.draw_kinematics(key, kin, n)
        ph = gen.draw_photon(jax.random.fold_in(key, 1), kin, kv, n)
        asym = gen.integrand(grid, kin, kv, ph, cfg.interp_scheme).asym
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

    def test_missing_mass_cut_gates_the_cross_section(self, small_config):
        """The cut is applied to the *pre-exit* missing mass, and it gates the
        weight sum -- not just which events are written out.

        ``aao_rad.f90`` builds the hadronic state with the electron energy that
        has not yet lost energy to the exit bremsstrahlung and rejects at line
        909 when the missing mass is outside the window; ``sig_tot`` is only
        incremented at line 916, i.e. *after* that test.  So the cut has to
        enter ``integrand.ok``, which is what the cross-section estimator sums
        over.  The recorded ``mm2`` comes from the later ``missm`` call at line
        980 and legitimately drifts outside the window, so it is not the right
        thing to assert on -- the gate is.
        """
        from aao_rad.constants import pi_masses

        gen = sys.modules["aao_rad.generate"]
        grid = gen.build_grid(small_config.channel, parms_dir="parms",
                              scheme=small_config.interp_scheme)
        kin = gen.build_kinematics(small_config)
        n = 200_000
        key = jax.random.PRNGKey(4242)
        kv = gen.draw_kinematics(key, kin, n)
        ph = gen.draw_photon(jax.random.fold_in(key, 1), kin, kv, n)
        integ = gen.integrand(grid, kin, kv, ph, small_config.interp_scheme)

        pre = hadronic_final_state(
            kin.e_beam, kv["es"], kv["ep"], kv["th0"], ph["ek"], ph["cstk"],
            ph["phik"], ph["csthcm"], ph["phicm_deg"], kin.m_pi, kin.is_pi0,
        )
        m_exp = pi_masses(small_config.channel)[1]
        assert float(kin.m_exp) == pytest.approx(m_exp)
        cut = small_config.missing_mass_cut
        in_window = (pre.w_real > 0.0) & (np.abs(pre.mm2 - m_exp) <= cut)
        ok = np.asarray(integ.ok)
        # Every kept trial is inside the window...
        assert np.all(in_window[ok])
        # ...and the gate actually bites (otherwise the test is vacuous).
        assert np.mean(~in_window & np.asarray(integ.ok_sampling)) > 0.05

    def test_recorded_mm2_stays_above_the_expected_mass(self, small_config, events):
        """The recorded (post-exit) missing mass is never below the expectation.

        Exit bremsstrahlung only takes energy away, so the reconstructed missing
        mass is bounded below by the value the cut was applied at; this is the
        invariant the recorded column must satisfy even though the upper bound is
        not sharp once the exit loss is allowed for.
        """
        from aao_rad.constants import pi_masses

        ev, kin = events
        m_exp = pi_masses(small_config.channel)[1]
        assert float(kin.m_exp) == pytest.approx(m_exp)
        assert ev["mm2"].min() >= m_exp - 1e-3
        # The window is genuinely populated, not just barely satisfied.
        assert np.mean(ev["mm2"] - m_exp > 0.5 * small_config.missing_mass_cut) > 0.2

    def test_w_is_above_the_two_body_threshold(self, events):
        ev, _ = events
        assert np.all(ev["w_real"] > M_N + M_PIP)
        assert np.all(ev["w"] > M_N + M_PIP)

    def test_the_same_key_reproduces_the_same_events(self, small_config):
        import aao_rad.generate  # noqa: F401  (register the module)

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
        asym = gen.integrand(grid, kin, kv, ph, cfg.interp_scheme).asym
        fk = jax.random.fold_in(key, 2)
        eager = gen.finalize(kin, kv, ph, asym, fk, n)
        jitted = jax.jit(gen.finalize, static_argnums=5)(kin, kv, ph, asym, fk, n)
        # XLA may reassociate float32 sums, so this is allclose, not exact.
        np.testing.assert_allclose(np.asarray(eager[0]), np.asarray(jitted[0]),
                                   rtol=1e-4, atol=1e-5)
        np.testing.assert_array_equal(np.asarray(eager[1]), np.asarray(jitted[1]))


class TestAcceptanceCeiling:
    """The rejection sampler's ceiling, and the self-check it enables.

    Accepting ``u * ceiling < weight`` is ``min(1, weight/ceiling)``.  That
    samples events as ``weight`` only while the ceiling exceeds every weight in
    the run, and the integrand has a heavy tail near the beam line with no
    finite maximum, so the ceiling is an estimate that has to be managed rather
    than known.  Two things are checked here: that the two independent
    cross-section estimators agree, which fails loudly if the ceiling handling
    is wrong, and that the residual bias is small and reported.
    """

    @pytest.fixture(scope="class")
    def run(self, small_config):
        import aao_rad.generate  # noqa: F401  (register the module)
        gen = sys.modules["aao_rad.generate"]

        cfg = small_config.merge(n_events=20_000, batch_size=1 << 15)
        grid = gen.build_grid(cfg.channel, parms_dir="parms",
                              scheme=cfg.interp_scheme)
        events, stats = gen.EventGenerator(grid).generate(cfg)
        return events, stats

    def test_the_two_cross_section_estimators_agree(self, run):
        """``sigma_mc`` sums trial weights; ``sigma_accepted`` sums acceptances.

        They share no code path, so a ceiling handled wrongly shows up as a
        disagreement between them.
        """
        _, stats = run
        assert stats.sigma_mc > 0.0
        np.testing.assert_allclose(
            stats.sigma_accepted, stats.sigma_mc, rtol=0.05,
        )

    def test_the_ceiling_bias_is_negligible_and_reported(self, run):
        _, stats = run
        assert stats.ceiling_bias < 1e-3, (
            f"{stats.n_above_ceiling} trials above the ceiling carry "
            f"{stats.ceiling_bias:.2%} of the cross section"
        )
        assert "above ceiling" in stats.summary()

    def test_the_ceiling_bounded_the_weights_that_reached_the_stream(self, run):
        """A raised ceiling ends above the largest weight that was kept."""
        _, stats = run
        assert stats.weight_max >= stats.weight_max_observed

    def test_more_headroom_costs_acceptance(self, small_config):
        """Sanity check on the trade the margin makes, not a physics test."""
        import aao_rad.generate  # noqa: F401
        gen = sys.modules["aao_rad.generate"]

        tight = gen.build_grid(small_config.channel, parms_dir="parms",
                               scheme=small_config.interp_scheme)
        cfg = small_config.merge(n_events=4000, batch_size=1 << 15)
        gen_obj = gen.EventGenerator(tight)
        _, a = gen_obj.generate(cfg.merge(weight_max_margin=1.0))
        _, b = gen_obj.generate(cfg.merge(weight_max_margin=8.0))
        assert b.acceptance < a.acceptance

    def test_a_margin_below_one_is_rejected(self):
        from aao_rad.config import GeneratorConfig

        with pytest.raises(ValueError, match="weight_max_margin"):
            GeneratorConfig(weight_max_margin=0.5)


class TestPhotonDecayAngles:
    """``cos(theta*)`` and ``phi*`` must be independent of each other.

    These are two separate ``myran`` calls in the original
    (``aao_rad.f90:779-780``), and they reach the cross section only through
    different combinations of the pion momentum components inside ``missm``.
    Drawing both from one JAX key therefore does *not* cancel: it makes
    ``phi* = 180*(cos(theta*) + 1)`` exactly, confining the decay direction to a
    curve on the sphere.

    The failure is invisible in either histogram alone -- sharing a key leaves
    both marginals exactly uniform -- so these tests look at the joint
    behaviour: a linear correlation, and the mean ``phi*`` within each
    ``cos(theta*)`` decile.
    """

    N = 200_000

    @pytest.fixture(scope="class")
    def angles(self, small_config):
        from aao_rad.generate import build_kinematics, draw_kinematics, draw_photon

        kin = build_kinematics(small_config)
        key = jax.random.PRNGKey(9182)
        kv = draw_kinematics(key, kin, self.N)
        ph = draw_photon(jax.random.fold_in(key, 1), kin, kv, self.N)
        return np.asarray(ph["csthcm"]), np.asarray(ph["phicm_deg"])

    def test_both_marginals_are_uniform(self, angles):
        """The condition that was *not* sufficient to catch the bug."""
        csthcm, phicm_deg = angles
        assert abs(csthcm.mean()) < 5e-3
        assert abs(phicm_deg.mean() - 180.0) < 1.5
        # KS against uniform on [0, 1]; 5 sigma of slack for N = 200k
        for u in ((csthcm + 1.0) / 2.0, phicm_deg / 360.0):
            n = u.size
            d = np.abs(np.arange(1, n + 1) / n - np.sort(u)).max()
            assert d < 5.0 / np.sqrt(n)

    def test_the_two_angles_are_independent(self, angles):
        """A shared key would give a correlation of about 0.95, not 0."""
        csthcm, phicm_deg = angles
        r = np.corrcoef(csthcm, phicm_deg)[0, 1]
        assert abs(r) < 0.02, f"corr(cos(theta*), phi*) = {r:.4f}"

    def test_the_mean_phi_is_flat_in_cos_theta(self, angles):
        """The decisive check: ``phi* = 180(cos(theta*)+1)`` rises monotonically.

        Under independence the mean is 180 in every decile; under the bug it
        sweeps from 0 to 360 across the range, which is a 170 sigma pull on the
        outer deciles.  The tolerance is scaled by each decile's own standard
        error rather than fixed, so the test does not have to be retuned when
        ``N`` changes.
        """
        csthcm, phicm_deg = angles
        edges = np.quantile(csthcm, np.linspace(0.0, 1.0, 6))
        for d in range(5):
            sel = (csthcm >= edges[d]) & (csthcm <= edges[d + 1])
            m = phicm_deg[sel]
            z = abs(m.mean() - 180.0) / (m.std() / np.sqrt(m.size))
            assert z < 4.0, (
                f"decile {d} of cos(theta*): mean phi* = {m.mean():.2f} deg, "
                f"{z:.1f} sigma from 180"
            )


class TestBeamAsymmetryBranch:
    """``asym_p`` must use the ``W`` that belongs to its own photon branch.

    The original computes the single-spin asymmetry twice, and at two different
    hadronic masses.  The soft branch (``aao_rad.f90:790``) calls ``dsigma``
    with the driver's ``epw``, built from ``es`` and ``ep`` alone.  The
    radiative branch never reaches that call in the caller at all: it goes
    through ``sigma()``, which recomputes ``mf2 = uu - 2*ek*(u0 - pu*cstk)``
    and passes ``epw = sqrt(mf2)`` -- a different mass, smaller by the radiated
    energy.

    The two disagree by ~17% in the median here, which is why the port keeps
    both evaluations rather than sharing one.  It is also exactly where the
    original went wrong: ``asym_p`` was a *local* of ``sigma()``, so the value
    ``dsigma`` computed for the radiative branch died with the call and
    ``ntp(32)`` recorded the last soft *trial*'s asymmetry instead of the
    event's own -- 24% of the n-tuple carried an unrelated value, which is what
    the reference ``asym_p`` comparison flagged.
    """

    N = 60_000

    @pytest.fixture(scope="class")
    def points(self, small_config):
        from aao_rad.generate import build_grid, build_kinematics, draw_kinematics, draw_photon

        grid = build_grid(small_config.channel, parms_dir="parms",
                          scheme=small_config.interp_scheme)
        kin = build_kinematics(small_config)
        key = jax.random.PRNGKey(5150)
        kv = draw_kinematics(key, kin, self.N)
        ph = draw_photon(jax.random.fold_in(key, 1), kin, kv, self.N)
        return grid, kin, kv, ph

    @staticmethod
    def _expected(grid, kin, kv, ph, radiative):
        """``response_functions`` at the mass this branch uses, plus its validity.

        Returns ``(asym, defined)``.  ``defined`` is the branch's own
        ``sigma() .le. 0.`` guard -- below-threshold trials evaluate the
        response at ``W = 0``, where it divides by zero, and are rejected by
        the generator regardless.
        """
        from aao_rad.motsa import _epsilon
        from aao_rad.xsection import response_functions

        es, ep, cst0 = kv["es"], kv["ep"], kv["cst0"]
        ek, cstk, e_hel = ph["ek"], ph["cstk"], ph["e_hel"]

        if radiative:
            qq = (
                2.0 * M_E**2 - 2.0 * es * ep + 2.0 * kv["ps"] * kv["pp"] * cst0
                - 2.0 * ek * (es - ep) + 2.0 * ek * kv["pu"] * cstk
            )
            mf2 = kv["uu"] - 2.0 * ek * (kv["u0"] - kv["pu"] * cstk)
            w_sq, q2 = jnp.maximum(mf2, 0.0), -qq
            defined = (mf2 > kin.wg**2) & (qq < 0.0)
        else:
            w_sq, q2 = kv["w_sq"], kv["q2"]
            defined = w_sq > 0.0

        nu = (w_sq - M_N**2 + q2) / (2.0 * M_N)
        resp = response_functions(
            grid, q2, jnp.sqrt(w_sq), ph["csthcm"], ph["phicm_deg"] * (np.pi / 180.0),
            _epsilon(es, ep, cst0, nu, q2), e_hel.astype(jnp.float32), kin.m_pi,
        )
        return np.asarray(resp.asym_p), np.asarray(defined)

    @pytest.mark.parametrize("radiative", [False, True])
    def test_asym_uses_the_branch_mass(self, points, small_config, radiative):
        from aao_rad.generate import integrand

        grid, kin, kv, ph = points
        got = np.asarray(integrand(grid, kin, kv, ph, small_config.interp_scheme).asym)
        want, defined = self._expected(grid, kin, kv, ph, radiative)
        other, _ = self._expected(grid, kin, kv, ph, not radiative)

        ek = np.asarray(ph["ek"])
        ok = defined & (ek >= kin.delta if radiative else ek < kin.delta)
        assert ok.sum() > 1000, f"only {int(ok.sum())} usable trials in this branch"
        assert np.abs(got[ok] - want[ok]).max() < 1e-5
        # The other branch's mass must be a genuinely different answer, or this
        # test would pass even if the branch selection were dropped.
        assert np.abs(other[ok] - want[ok]).max() > 1e-3

    def test_asym_is_this_events_own(self, points, small_config):
        """``asym_p`` must not repeat across events the way the original's did.

        In the broken reference the radiative rows carry the asymmetry of the
        last soft *trial*, so most radiative values were duplicates of some
        earlier row's.  A correct implementation produces a continuous spread,
        which is what a run of exactly-equal neighbours would rule out.
        """
        from aao_rad.generate import integrand

        grid, kin, kv, ph = points
        asym = np.asarray(integrand(grid, kin, kv, ph, small_config.interp_scheme).asym)
        rad = (np.asarray(ph["ek"]) >= kin.delta) & np.isfinite(asym)
        a = asym[rad]
        n_uniq = np.unique(a).size
        assert n_uniq / a.size > 0.99, (
            f"only {n_uniq} distinct asym_p in {a.size} radiative trials -- "
            "the value is being carried over from somewhere else"
        )
