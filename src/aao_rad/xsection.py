"""Exclusive single-pion electroproduction response functions.

Batched port of ``maid_lee.f90`` + ``xsection.f90`` + ``dsigma.f90`` restricted
to the MAID07 model (``theory_opt = 7``), which is the only model the original
driver could actually reach.

The seven quantities returned are the *unpolarised* invariant cross section and
the six structure/response functions the radiative generator needs::

    sigma_0   dsigma / dOmega_h  (unpolarised)
    sigma_t   (|H1|^2 + |H2|^2 + |H3|^2 + |H4|^2) / 2
    sigma_l   |H5|^2 + |H6|^2
    sigma_tt  Re(H3 H2* - H4 H1*)
    sigma_lt  sqrt(2) Re(H5* (H1 - H4) + H6* (H2 + H3))
    sigma_ltp sqrt(2) Im(H5* (H4 - H1) - H6* (H2 + H3))
    asym_p    beam-helicity asymmetry sigma_p / sigma_u
"""

from __future__ import annotations

from typing import NamedTuple

import jax.numpy as jnp

from .amplitudes import (
    cgln_amplitudes,
    helicity_amplitudes,
    legendre_polynomials,
    multipole_amplitudes,
)
from .constants import M_N
from .interpolate import InterpolationGrid

__all__ = ["Response", "response_functions", "ResponseModel"]


class Response(NamedTuple):
    """The seven response functions at a batch of kinematic points."""

    sigma0: jnp.ndarray
    sigma_t: jnp.ndarray
    sigma_l: jnp.ndarray
    sigma_tt: jnp.ndarray
    sigma_lt: jnp.ndarray
    sigma_ltp: jnp.ndarray
    asym_p: jnp.ndarray


def _center_of_mass(w: jnp.ndarray, q2: jnp.ndarray, m_pi: float) -> tuple[jnp.ndarray, ...]:
    """Centre-of-mass kinematics of the pion (``maid_lee.f90`` lines 100-107).

    ``w`` and ``q2`` must be the *unclamped* kinematics.  ``maid_lee`` evaluates
    these before it ever reaches the table, and ``multipole_amps.f90`` later
    clamps only the ``interp`` argument, so the original really does combine a
    clamped lookup with unclamped ``nu_cm``/``qv_mag_cm``.  See
    :func:`response_functions`.
    """
    m_p2 = M_N * M_N
    w2 = w * w
    e_pi_cm = 0.5 * (w2 + m_pi * m_pi - m_p2) / w
    p_pi_cm = jnp.sqrt(jnp.maximum(e_pi_cm * e_pi_cm - m_pi * m_pi, 0.0))
    qv_cm = jnp.sqrt(jnp.maximum(((w2 + q2 + m_p2) / (2.0 * w)) ** 2 - m_p2, 0.0))
    nu_cm = (w2 - m_p2 - q2) / (2.0 * w)
    return p_pi_cm, qv_cm, nu_cm


def _fkt(w: jnp.ndarray, p_pi_cm: jnp.ndarray) -> jnp.ndarray:
    """Pion phase-space factor ``fkt`` (``xsection.f90`` line 28).

    ``xsection.f90`` recomputes ``fkt`` from ``W`` *after* calling
    ``multipole_amps``, and ``multipole_amps.f90`` line 16 assigns ``w = 1.1``
    into the shared COMMON when ``w <= 1.1``.  So the denominator's ``W`` carries
    that floor while ``ppi_mag_cm`` and ``nu_cm`` do not.  Below threshold the
    original therefore evaluates a different (and dimensionally odd) expression
    from the one above threshold; that asymmetry is reproduced here rather than
    smoothed over.
    """
    w_floor = jnp.maximum(w, 1.1)
    return 2.0 * w_floor * p_pi_cm / jnp.maximum(w_floor * w_floor - M_N * M_N, 1e-12)


def response_functions(
    grid: InterpolationGrid,
    q2: jnp.ndarray,
    w: jnp.ndarray,
    cos_theta_cm: jnp.ndarray,
    phi_cm_rad: jnp.ndarray,
    epsilon: jnp.ndarray,
    e_hel: jnp.ndarray,
    m_pi: float,
    *,
    scheme: str = "linear",
) -> Response:
    """Evaluate the MAID07 response functions at a batch of points.

    Parameters
    ----------
    grid
        Interpolation grid built from a :class:`~aao_rad.table.MaidTable`.
    q2, w
        Momentum transfer [GeV^2] and hadronic invariant mass [GeV].
    cos_theta_cm
        ``cos(theta*)`` of the pion in the hadronic CM frame.
    phi_cm_rad
        ``phi*`` **in radians** (the Fortran interface took degrees; radians
        avoids a multiply/divide pair inside the hot loop).
    epsilon
        Virtual photon polarisation, ``1 / (1 + 2(1 + nu^2/Q2) tan^2(theta_e/2))``.
    e_hel
        Beam helicity, ``+1``/``-1``/``0`` (0 selects the unpolarised average).
    m_pi
        Pion mass in GeV.
    scheme
        Interpolation scheme, ``"linear"`` (default, matches the original) or
        ``"spline"``.
    """
    # The original's clamping is *not* a single up-front clip, and getting that
    # wrong is visible.  multipole_amps.f90:16-28 raises W to 1.1 in the shared
    # COMMON and then picks one of four `interp` corners -- (Q2 or 5) x (W or
    # 2.0) -- so the *lookup* saturates in both directions while maid_lee.f90
    # has already computed nu_cm, qv_mag_cm and ppi_mag_cm from the unclamped
    # kinematics.  xsection.f90:27-28 then re-reads W (now floored at 1.1, but
    # not capped above) for ekin and fkt.  Clipping W before the kinematics
    # instead, as an earlier revision did, inflates sigma_l by (ekin_clamped /
    # ekin_unclipped)^2: 1.35x at W = 2.2, Q^2 = 0.1, rising to 5x by
    # W = 3, because ekin = sqrt(Q^2)/nu_cm grows as nu_cm shrinks.
    w_c = jnp.clip(w, 1.1, 2.0)
    q2_c = jnp.minimum(q2, 5.0)

    p_pi_cm, qv_cm, nu_cm = _center_of_mass(w, q2, m_pi)
    fkt = _fkt(w, p_pi_cm)

    amps = grid(q2_c, w_c, scheme=scheme)  # (batch, 62)
    sp, sm, ep, em, mp, mm = multipole_amplitudes(amps, nu_cm, qv_cm)

    pol = legendre_polynomials(cos_theta_cm)
    ff = cgln_amplitudes(pol, sp, sm, ep, em, mp, mm)
    hh = helicity_amplitudes(*ff, cos_theta_cm)

    hh1, hh2, hh3, hh4, hh5, hh6 = hh
    sqrt2 = jnp.sqrt(2.0)

    sigma_t = (jnp.abs(hh1) ** 2 + jnp.abs(hh2) ** 2 + jnp.abs(hh3) ** 2 + jnp.abs(hh4) ** 2) / 2.0
    sigma_l = jnp.abs(hh5) ** 2 + jnp.abs(hh6) ** 2
    sigma_tt = (hh3 * jnp.conj(hh2) - hh4 * jnp.conj(hh1)).real
    sigma_lt = sqrt2 * (jnp.conj(hh5) * (hh1 - hh4) + jnp.conj(hh6) * (hh2 + hh3)).real
    sigma_ltp = sqrt2 * (jnp.conj(hh5) * (hh4 - hh1) - jnp.conj(hh6) * (hh2 + hh3)).imag

    # Longitudinal terms are suppressed by ekin = |Q| / nu_cm, with both taken
    # unclamped (xsection.f90 line 27).
    ekin = jnp.sqrt(jnp.maximum(q2, 0.0)) / nu_cm
    sigma_l = sigma_l * ekin**2
    sigma_lt = sigma_lt * ekin
    sigma_ltp = sigma_ltp * ekin

    # vlt / vltp are a factor 2 smaller than in spp_int_e1 in order to match the
    # AO formalism (see the comment in xsection.f90).  epsilon is at most 1 by
    # construction, but rounding at forward scattering (tan(theta_e/2) -> 0)
    # can push it a few ulps past 1, and sqrt of a negative would poison the
    # whole batch with NaN, so clip before taking the root.
    eps_c = jnp.clip(epsilon, 0.0, 1.0)
    vlt = jnp.sqrt(eps_c * (1.0 + eps_c) / 2.0)
    vltp = jnp.sqrt(eps_c * (1.0 - eps_c) / 2.0)

    cos2 = jnp.cos(2.0 * phi_cm_rad)
    cos1 = jnp.cos(phi_cm_rad)
    sin1 = jnp.sin(phi_cm_rad)

    sigma_u = fkt * (
        sigma_t + epsilon * sigma_l + epsilon * sigma_tt * cos2 + vlt * sigma_lt * cos1
    )
    sigma_p = fkt * vltp * sigma_ltp * sin1

    sigma_0 = sigma_u + e_hel * sigma_p
    # An unpolarised beam cannot observe the single-spin asymmetry.
    safe_u = jnp.where(sigma_u == 0.0, 1.0, sigma_u)
    asym = jnp.where(jnp.abs(e_hel) > 0.5, sigma_p / safe_u, 0.0)

    return Response(sigma_0, sigma_t, sigma_l, sigma_tt, sigma_lt, sigma_ltp, asym)


class ResponseModel:
    """Callable wrapper binding a table, scheme and pion mass to the physics."""

    __slots__ = ("grid", "m_pi", "scheme", "_call")

    def __init__(self, grid: InterpolationGrid, m_pi: float, scheme: str = "linear") -> None:
        self.grid = grid
        self.m_pi = m_pi
        self.scheme = scheme
        import jax

        self._call = jax.jit(
            lambda q2, w, cth, phi, eps, ehel: response_functions(
                self.grid, q2, w, cth, phi, eps, ehel, self.m_pi, scheme=self.scheme
            ),
            static_argnums=(),
        )

    def __call__(
        self, q2, w, cos_theta_cm, phi_cm_rad, epsilon, e_hel
    ) -> Response:
        return self._call(q2, w, cos_theta_cm, phi_cm_rad, epsilon, e_hel)
