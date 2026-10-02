"""Tests for the interpolation kernels.

The interpolator is the numerical core of the response functions, so the two
properties that matter are checked here: it must reproduce the tabulated values
exactly at grid nodes, and it must stay finite and continuous everywhere in
between and at the boundaries.
"""

from __future__ import annotations

import jax
import jax.numpy as jnp
import numpy as np
import pytest
from conftest import requires_tables

from aao_rad.interpolate import (
    InterpolationGrid,
    _splint,
    bilinear,
    natural_spline_2nd_derivs,
    spline_tensor,
)

SCHEMES = ["linear", "spline"]


def _toy_grid(n1: int = 7, n2: int = 5, k: int = 3) -> InterpolationGrid:
    """A small, exactly-representable grid so node values are exact."""
    a1 = jnp.linspace(0.0, 1.0, n1)
    a2 = jnp.linspace(2.0, 4.0, n2)
    # A separable polynomial in the axis indices, so we know the answer.
    i = jnp.arange(n1, dtype=jnp.float32)[:, None]
    j = jnp.arange(n2, dtype=jnp.float32)[None, :]
    values = jnp.stack(
        [1.0 + i + 2.0 * j + 0.5 * i * j,
         2.0 - i + j,
         0.25 * i * i - 0.75 * j * j],
        axis=0,
    )
    return InterpolationGrid(
        values=values, axis1=a1, axis2=a2, uniform1=True, uniform2=True,
        d2_1=None, d2_2=None,
    )


class TestBilinear:
    def test_reproduces_nodes_exactly(self):
        g = _toy_grid()
        # Sample every node of the grid.
        x1 = jnp.repeat(g.axis1, len(g.axis2))
        x2 = jnp.tile(g.axis2, len(g.axis1))
        got = bilinear(g.values, g.axis1, g.axis2, x1, x2)
        want = jnp.moveaxis(g.values, 0, -1).reshape(len(x1), -1)
        np.testing.assert_allclose(np.asarray(got), np.asarray(want), rtol=1e-6, atol=1e-6)

    def test_shape_is_batch_by_amplitude(self):
        g = _toy_grid()
        x = jnp.linspace(0.1, 0.9, 17)
        out = bilinear(g.values, g.axis1, g.axis2, x, x)
        assert out.shape == (17, 3)

    def test_is_affine_in_a_bilinear_function(self):
        """The interpolant of an affine function is that function."""
        n1 = n2 = 9
        a1 = np.linspace(-1.0, 1.0, n1)
        a2 = np.linspace(3.0, 5.0, n2)
        f = lambda x, y: 2.0 * x - 3.0 * y + 7.0  # noqa: E731
        # Tabulate the function in axis coordinates, not index coordinates.
        values = f(a1[:, None], a2[None, :])[None, :, :]

        x = jnp.linspace(-1.0, 1.0, 23)
        y = jnp.linspace(3.0, 5.0, 23)
        got = bilinear(values, jnp.asarray(a1), jnp.asarray(a2), x, y)[:, 0]
        np.testing.assert_allclose(
            np.asarray(got), np.asarray(f(x, y)), rtol=1e-5, atol=1e-5
        )

    def test_clamps_outside_the_grid(self):
        """Outside the grid the value is held at the edge, not extrapolated.

        This is the behaviour that keeps a lookup at W > 2 GeV from running off
        the end of the table (the original does the same, by hand, in
        ``multipole_amps.f90``).
        """
        g = _toy_grid()
        # Query well outside axis 1 but at the same axis-2 node, so only the
        # axis-1 clamping is under test.
        far = bilinear(g.values, g.axis1, g.axis2,
                       jnp.array([-5.0, 10.0]), g.axis2[1])
        edge_lo = bilinear(g.values, g.axis1, g.axis2, g.axis1[0], g.axis2[1])
        edge_hi = bilinear(g.values, g.axis1, g.axis2, g.axis1[-1], g.axis2[1])
        np.testing.assert_allclose(np.asarray(far[0]), np.asarray(edge_lo),
                                   rtol=1e-5, atol=1e-6)
        np.testing.assert_allclose(np.asarray(far[1]), np.asarray(edge_hi),
                                   rtol=1e-5, atol=1e-6)

    def test_is_finite_everywhere_on_a_dense_grid(self):
        g = _toy_grid()
        x = jnp.linspace(-1.0, 2.0, 64)
        x1, x2 = jnp.meshgrid(x, x, indexing="ij")
        out = bilinear(g.values, g.axis1, g.axis2, x1.ravel(), x2.ravel())
        assert np.isfinite(np.asarray(out)).all()

    def test_is_jittable(self):
        """The jitted result must match the eager one (XLA may reassociate)."""
        g = _toy_grid()
        f = jax.jit(lambda x: bilinear(g.values, g.axis1, g.axis2, x, x))
        x = jnp.linspace(0.0, 1.0, 11)
        np.testing.assert_allclose(np.asarray(f(x)),
                                   np.asarray(bilinear(g.values, g.axis1, g.axis2, x, x)),
                                   rtol=1e-6, atol=1e-6)


class TestSpline:
    def test_second_derivatives_of_a_quadratic_with_matching_slopes(self):
        """A spline through a quadratic whose end slopes are pinned *is* that
        quadratic, so its second derivatives are the constant 2a everywhere --
        unlike the natural spline, which forces them to zero at the ends.
        """
        n = 9
        x = np.linspace(0.0, 1.0, n)
        a, b, c = 3.0, -2.0, 1.0
        y = a * x**2 + b * x + c
        d2 = np.asarray(
            natural_spline_2nd_derivs(y, x, yp1=2 * a * x[0] + b,
                                      ypn=2 * a * x[-1] + b)
        )
        np.testing.assert_allclose(d2, 2 * a, rtol=1e-5, atol=1e-5)

    def test_natural_spline_of_a_quadratic_is_exact_at_the_nodes(self):
        """Natural boundary conditions do not spoil node reproduction."""
        n = 9
        x = np.linspace(0.0, 1.0, n)
        y = 3.0 * x**2 - 2.0 * x + 1.0
        d2 = np.asarray(natural_spline_2nd_derivs(y, x))
        # The curvature is right in the interior; the ends are pinned to zero.
        assert d2[0] == 0.0 and d2[-1] == 0.0
        np.testing.assert_allclose(d2[1:-1], 6.0, rtol=0.35)
        out = np.asarray(
            _splint(jnp.asarray(x, jnp.float32), jnp.asarray(y, jnp.float32),
                    jnp.asarray(d2, jnp.float32), jnp.asarray(x, jnp.float32))
        )
        np.testing.assert_allclose(out, y, rtol=1e-6, atol=1e-6)

    def test_second_derivatives_of_a_line_vanish(self):
        n = 11
        x = np.linspace(-2.0, 3.0, n)
        d2 = np.asarray(natural_spline_2nd_derivs(4.0 * x + 1.0, x))
        np.testing.assert_allclose(d2, 0.0, atol=1e-4)

    def test_second_derivatives_vanish_at_the_ends(self):
        # Natural spline: the boundary second derivatives are zero by definition.
        n = 12
        x = np.linspace(0.0, 2.0, n)
        y = np.sin(x) + 0.3 * x**3
        d2 = np.asarray(natural_spline_2nd_derivs(y, x))
        assert abs(d2[0]) < 1e-6
        assert abs(d2[-1]) < 1e-6

    def test_reproduces_nodes(self):
        g = _toy_grid()
        d2_1 = np.moveaxis(
            natural_spline_2nd_derivs(np.moveaxis(np.asarray(g.values), 1, 0),
                                      np.asarray(g.axis1)),
            0, 1,
        )
        d2_2 = np.moveaxis(
            natural_spline_2nd_derivs(np.moveaxis(np.asarray(g.values), 2, 0),
                                      np.asarray(g.axis2)),
            0, 2,
        )
        assert d2_1.shape == g.values.shape
        assert d2_2.shape == g.values.shape

        x1 = jnp.repeat(g.axis1, len(g.axis2))
        x2 = jnp.tile(g.axis2, len(g.axis1))
        got = spline_tensor(g.values, g.axis1, g.axis2,
                            jnp.asarray(d2_1), jnp.asarray(d2_2), x1, x2)
        want = jnp.moveaxis(g.values, 0, -1).reshape(len(x1), -1)
        np.testing.assert_allclose(np.asarray(got), np.asarray(want), rtol=1e-5, atol=1e-5)

    def test_second_derivatives_have_the_table_shape(self, pnpi_table, parms_dir):
        """The loaded derivative grids must line up with the amplitude grid.

        ``natural_spline_2nd_derivs`` splines along its *leading* axis, so
        getting the transpose wrong silently yields an array of the right
        size in the wrong order -- which is exactly the kind of bug that
        produces plausible-looking but wrong events.
        """
        from aao_rad.generate import build_grid

        g = build_grid(3, parms_dir=parms_dir, scheme="spline")
        assert g.d2_1 is not None and g.d2_2 is not None
        assert g.d2_1.shape == g.values.shape
        assert g.d2_2.shape == g.values.shape
        assert g.d2_1.shape[0] == 62

    def test_spline_interpolates_the_physical_grid(self, parms_dir):
        """A real table: the spline and bilinear schemes must broadly agree.

        They are different interpolants of the same tabulated data, so they may
        differ -- but only slightly, and never by a large factor.
        """
        from aao_rad.generate import build_grid

        gl = build_grid(3, parms_dir=parms_dir, scheme="linear")
        gs = build_grid(3, parms_dir=parms_dir, scheme="spline")
        n = 4000
        rng = np.random.default_rng(7)
        q1 = jnp.asarray(rng.uniform(0.05, 4.95, n), jnp.float32)
        q2 = jnp.asarray(rng.uniform(1.12, 1.98, n), jnp.float32)

        b = np.asarray(bilinear(gl.values, gl.axis1, gl.axis2, q1, q2))
        s = np.asarray(spline_tensor(gs.values, gs.axis1, gs.axis2,
                                     gs.d2_1, gs.d2_2, q1, q2))
        assert np.isfinite(s).all()
        assert b.shape == s.shape
        # Compare where the tabulated amplitudes are not numerically zero.
        mask = np.abs(b) > 1e-6 * np.nanmax(np.abs(b))
        rel = np.abs(b[mask] - s[mask]) / np.abs(b[mask])
        assert np.median(rel) < 0.02, f"median relative difference {np.median(rel):.3f}"


class TestSchemesAgree:
    """Both schemes must track each other to within spline accuracy.

    The original hard-selected bilinear (``method_spline = 2``), so ``linear``
    is the default; ``spline`` is offered as a smoother alternative and is only
    trustworthy if it does not wander far from it.
    """

    @staticmethod
    def _points(n: int = 2000, seed: int = 11):
        """Random kinematics spread across the whole tabulated grid."""
        k1, k2, k3, k4, k5 = jax.random.split(jax.random.PRNGKey(seed), 5)
        return {
            "q2": jax.random.uniform(k1, (n,), minval=0.05, maxval=4.95,
                                     dtype=jnp.float32),
            "w": jax.random.uniform(k2, (n,), minval=1.12, maxval=1.98,
                                    dtype=jnp.float32),
            "cos_theta_cm": jax.random.uniform(k3, (n,), minval=-0.95, maxval=0.95,
                                               dtype=jnp.float32),
            "phi_cm_rad": jax.random.uniform(k4, (n,), minval=0.0,
                                             maxval=2 * np.pi, dtype=jnp.float32),
            # Polarised beam, as the default configuration uses.
            "epsilon": jax.random.uniform(k5, (n,), minval=0.3, maxval=0.99,
                                          dtype=jnp.float32),
            "e_hel": jnp.where(jax.random.uniform(k5, (n,)) < 0.5, -1.0, 1.0),
        }

    @requires_tables
    @pytest.mark.parametrize("scheme", SCHEMES)
    def test_response_is_finite_over_the_whole_table(self, parms_dir, scheme):
        from aao_rad.constants import pi_masses
        from aao_rad.generate import build_grid
        from aao_rad.xsection import Response, response_functions

        grid = build_grid(3, parms_dir=parms_dir, scheme=scheme)
        m_pi, m_exp = pi_masses(3)
        p = self._points()

        resp = response_functions(
            grid, p["q2"], p["w"], p["cos_theta_cm"], p["phi_cm_rad"],
            p["epsilon"], p["e_hel"], m_pi, scheme=scheme,
        )
        for name in Response._fields:
            val = np.asarray(getattr(resp, name))
            assert np.isfinite(val).all(), f"{name} has non-finite values"
            assert np.nanmax(np.abs(val)) < 1e6, f"{name} is wildly large"

    @requires_tables
    def test_unpolarised_beam_drops_the_single_spin_term(self, parms_dir):
        """With no beam polarisation only the unpolarised response is left."""
        from aao_rad.constants import pi_masses
        from aao_rad.generate import build_grid
        from aao_rad.xsection import response_functions

        grid = build_grid(3, parms_dir=parms_dir)
        m_pi, _ = pi_masses(3)
        p = self._points(500)

        one = jnp.ones_like(p["e_hel"])
        args = (p["q2"], p["w"], p["cos_theta_cm"], p["phi_cm_rad"], p["epsilon"])
        unpolarised = np.asarray(
            response_functions(grid, *args, jnp.zeros_like(one), m_pi).sigma0)
        plus = np.asarray(response_functions(grid, *args, one, m_pi).sigma0)
        minus = np.asarray(response_functions(grid, *args, -one, m_pi).sigma0)
        # sigma0(+1) + sigma0(-1) == 2 * sigma0(0) by construction.
        np.testing.assert_allclose(plus + minus, 2.0 * unpolarised,
                                   rtol=1e-5, atol=1e-9)
        # ...and the asymmetry is a real antisymmetric signal, not noise.
        asym = np.asarray(response_functions(grid, *args, one, m_pi).asym_p)
        assert np.abs(asym).max() > 0.0

    @requires_tables
    def test_spline_stays_close_to_linear(self, parms_dir):
        from aao_rad.constants import pi_masses
        from aao_rad.generate import build_grid
        from aao_rad.xsection import response_functions

        gl = build_grid(3, parms_dir=parms_dir, scheme="linear")
        gs = build_grid(3, parms_dir=parms_dir, scheme="spline")
        m_pi, _ = pi_masses(3)
        p = self._points(1500, seed=3)

        a = np.asarray(response_functions(
            gl, p["q2"], p["w"], p["cos_theta_cm"], p["phi_cm_rad"],
            p["epsilon"], p["e_hel"], m_pi, scheme="linear").sigma0)
        b = np.asarray(response_functions(
            gs, p["q2"], p["w"], p["cos_theta_cm"], p["phi_cm_rad"],
            p["epsilon"], p["e_hel"], m_pi, scheme="spline").sigma0)
        assert np.isfinite(a).all() and np.isfinite(b).all()
        # Compare where the response is not numerically zero.
        mask = np.abs(a) > 1e-3 * np.max(np.abs(a))
        assert mask.sum() > 100
        rel = np.abs(a[mask] - b[mask]) / np.abs(a[mask])
        assert np.median(rel) < 0.05, f"median relative difference {np.median(rel):.3f}"
