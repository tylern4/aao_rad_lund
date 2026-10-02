"""Mo & Tsai radiative cross section with the AO hadronic response.

Batched port of the ``sigma()`` function at the bottom of ``aao_rad.f90``:
the bremsstrahlung part is the exact Mo & Tsai expression, with the hadronic
part obtained by converting the MAID response functions into the ``f``/``g``
amplitudes the formula expects (the original author's "fake up" prescription,
preserved verbatim).

The whole function is elementwise in the kinematic variables, so a batch of
points costs exactly one pass.
"""

from __future__ import annotations

import jax.numpy as jnp

from .constants import ALPHA, M_E, M_N, PI
from .xsection import Response, response_functions

__all__ = ["motsa_sigma", "non_radiative_sigma", "spence", "polarization_modulation"]

_TINY = 1.0e-31
"""Sentinel the Fortran returned when a kinematic check failed."""

_SERIES_K = jnp.arange(1, 41, dtype=jnp.float32)


def _safe(x: jnp.ndarray) -> jnp.ndarray:
    """Replace exact zeros by 1.0 so divisions stay finite (results are masked)."""
    return jnp.where(x == 0.0, 1.0, x)


# --------------------------------------------------------------------------
# Dilogarithm
# --------------------------------------------------------------------------
def _li2_series(v: jnp.ndarray) -> jnp.ndarray:
    """``sum_k v**k / k**2`` -- accurate for ``|v| <= 1/2``."""
    return jnp.sum(v[..., None] ** _SERIES_K / _SERIES_K**2, axis=-1)


def _li2_unit(u: jnp.ndarray) -> jnp.ndarray:
    """``Li2(u)`` for ``0 <= u < 1`` via the direct series or Euler reflection."""
    v = 1.0 - u
    direct = jnp.where(u <= 0.5, _li2_series(u), 0.0)
    reflected = PI**2 / 6.0 - jnp.log(jnp.maximum(u, 1e-30)) * jnp.log(jnp.maximum(v, 1e-30)) - (
        _li2_series(v)
    )
    return jnp.where(u <= 0.5, direct, reflected)


def spence(x: jnp.ndarray) -> jnp.ndarray:
    r"""Real dilogarithm :math:`\mathrm{Li}_2(x)` for ``x <= 1``.

    The Fortran used a 100-step trapezoid rule with hard-coded special cases
    around ``x = +/-1``; this uses the convergent series plus the standard
    reflection/continuation identities, which is exact, has no branch
    discontinuities and is differentiable (so ``jit`` and ``grad`` work).
    """
    x = jnp.asarray(x, jnp.float32)
    positive = x >= 0.0
    xp = jnp.where(positive, x, 1.0)
    xn = jnp.where(positive, 0.0, x)

    # 0 <= x <= 1
    val_small = _li2_unit(jnp.minimum(xp, 1.0))
    # x > 1 :  Li2(x) = pi^2/3 - ln(x)^2/2 - Li2(1/x)
    val_big = PI**2 / 3.0 - 0.5 * jnp.log(xp) ** 2 - _li2_unit(1.0 / xp)
    val_pos = jnp.where(xp > 1.0, val_big, val_small)

    # x <= 0 :  Li2(x) = -Li2(x/(x-1)) - ln(1-x)^2/2     with x/(x-1) in [0, 1)
    val_neg = -_li2_unit(xn / _safe(xn - 1.0)) - 0.5 * jnp.log1p(-xn) ** 2

    return jnp.where(positive, val_pos, val_neg)


# --------------------------------------------------------------------------
# Shared helpers
# --------------------------------------------------------------------------
def polarization_modulation(
    resp: Response, epeps: jnp.ndarray, phi: jnp.ndarray, e_hel: jnp.ndarray
) -> jnp.ndarray:
    """Relative transverse/longitudinal interference term.

    Reproduces the ``1 + (...)`` factor that modulates both the radiative and
    the non-radiative cross section::

        (eps*sig_tt*cos(2 phi) + sqrt(eps(1+eps)/2)*sig_lt*cos(phi)
             + e_hel*sqrt(eps(1-eps)/2)*sig_ltp*sin(phi))
        ---------------------------------------------------------
                          sig_t + eps*sig_l
    """
    denom = resp.sigma_t + epeps * resp.sigma_l
    return (
        epeps * resp.sigma_tt * jnp.cos(2.0 * phi)
        + jnp.sqrt(epeps * (1.0 + epeps) / 2.0) * resp.sigma_lt * jnp.cos(phi)
        + e_hel * jnp.sqrt(epeps * (1.0 - epeps) / 2.0) * resp.sigma_ltp * jnp.sin(phi)
    ) / _safe(denom)


def _epsilon(es, ep, cst0, nu, q2):
    """Virtual-photon polarisation ``1 / (1 + 2(1 + nu^2/Q2) tan^2(theta_e/2))``."""
    s2 = (1.0 - cst0) / 2.0
    return 1.0 / (1.0 + 2.0 * (1.0 + nu**2 / _safe(q2)) * s2 / _safe(1.0 - s2))


def _pion_kinematic_factor(w_sq, epw, m_pi):
    """``2 W |q_pi*| / (W^2 - m_N^2)`` evaluated from the two-body kinematics."""
    fkt_raw = (w_sq - M_N**2 + m_pi**2) / (2.0 * _safe(epw))
    return (
        jnp.sqrt(jnp.maximum(fkt_raw**2 - m_pi**2, 0.0))
        * 2.0
        * _safe(epw)
        / _safe(w_sq - M_N**2)
    )


def _fg_amplitudes(resp: Response, w_sq, q2, m_pi):
    """Convert the response functions into Mo & Tsai's ``f`` and ``g``."""
    epw = jnp.sqrt(jnp.maximum(w_sq, 0.0))
    nu = (w_sq + q2 - M_N**2) / (2.0 * M_N)
    kfac = (w_sq - M_N**2) / (2.0 * M_N)
    f = (
        1.0
        / (2.0 * PI**2 * ALPHA * M_N)
        * (kfac / (1.0 + nu**2 / _safe(q2)))
        * (resp.sigma_t + resp.sigma_l)
    )
    g = M_N / (2.0 * PI**2 * ALPHA) * kfac * resp.sigma_t
    fkt = _pion_kinematic_factor(w_sq, epw, m_pi)
    return f * fkt, g * fkt, nu


# --------------------------------------------------------------------------
# Radiative (ek > delta) part
# --------------------------------------------------------------------------
def motsa_sigma(
    grid,
    *,
    es,
    ep,
    th0,
    cst0,
    ps,
    pp,
    sdotk,
    pdotk,
    u0,
    pu,
    uu,
    cstk,
    ek,
    csthcm,
    phi,
    e_hel,
    wg,
    m_pi,
    scheme: str = "linear",
) -> tuple[jnp.ndarray, jnp.ndarray]:
    """Mo & Tsai differential cross section for one batch of phase-space points.

    All array arguments share the same leading batch axis.  ``phi`` is the
    hadronic-frame azimuth in **radians**.

    Returns
    -------
    (sigma, asym_p)
        Cross section in the units the MAID tables are normalised to
        (microbarn), and the single-spin beam asymmetry.  ``sigma`` collapses
        to a negligible sentinel wherever a kinematic check rejects the point.
    """
    qq = 2.0 * M_E**2 - 2.0 * es * ep + 2.0 * ps * pp * cst0 - 2.0 * ek * (es - ep) + 2.0 * ek * pu * cstk
    mf2 = uu - 2.0 * ek * (u0 - pu * cstk)

    ok = (mf2 > wg**2) & (qq < 0.0)
    w_sq = jnp.where(ok, mf2, 0.0)

    sp = es * ep - ps * pp * cst0
    qsq = -qq

    sdotk_s, pdotk_s = _safe(sdotk), _safe(pdotk)

    ffac = (
        -(M_E / pdotk_s) ** 2 * (2.0 * es * (ep + ek) + qq / 2.0)  # ffac1
        - (M_E / sdotk_s) ** 2 * (2.0 * ep * (es - ek) + qq / 2.0)  # ffac2
        - 2.0  # ffac3
        + 2.0 / sdotk_s / pdotk_s
        * (M_E**2 * (sp - ek**2) + sp * (2.0 * es * ep - sp + ek * (es - ep)))  # ffac4
        + (2.0 * (es * ep + es * ek + ep**2) + qq / 2.0 - sp - M_E**2) / pdotk_s  # ffac5
        - (2.0 * (es * ep - ep * ek + es**2) + qq / 2.0 - sp - M_E**2) / sdotk_s  # ffac6
    )
    gfac = (
        M_E**2 * (2.0 * M_E**2 + qq) * (1.0 / pdotk_s**2 + 1.0 / sdotk_s**2)  # gfac1
        + 4.0  # gfac2
        + 4.0 * sp * (sp - 2.0 * M_E**2) / pdotk_s / sdotk_s  # gfac3
        + (2.0 * sp + 2.0 * M_E**2 - qq) * (1.0 / pdotk_s - 1.0 / sdotk_s)  # gfac4
    )

    epw = jnp.sqrt(jnp.maximum(w_sq, 0.0))
    eps = _epsilon(es, ep, cst0, (w_sq - M_N**2 + qsq) / (2.0 * M_N), qsq)
    resp = response_functions(grid, qsq, epw, csthcm, phi, eps, e_hel, m_pi, scheme=scheme)

    epeps = _epsilon(es, ep, cst0, es - ep, qsq)
    f, g, nu = _fg_amplitudes(resp, w_sq, qsq, m_pi)

    sig_r = ((ALPHA**3 / (2.0 * PI * qq) ** 2) / M_N) * (ep / _safe(es)) * ek
    sig_f = M_N**2 * f * ffac + g * gfac
    denom = resp.sigma_t + epeps * resp.sigma_l
    sig_f = sig_f * (1.0 + polarization_modulation(resp, epeps, phi, e_hel))

    good = (
        ok
        & (ffac > 0.0)
        & (gfac > 0.0)
        & (resp.sigma0 > 0.0)
        & (epeps > 0.0)
        & (epeps <= 1.0)
        & (nu > 0.0)
        & (denom > 0.0)
        & (sig_f > 0.0)
    )
    return jnp.where(good, sig_r * sig_f, _TINY), jnp.where(ok, resp.asym_p, 0.0)


# --------------------------------------------------------------------------
# Non-radiative (ek < delta) part
# --------------------------------------------------------------------------
def non_radiative_sigma(
    grid,
    *,
    es,
    ep,
    th0,
    cst0,
    ps,
    pp,
    q2,
    w_sq,
    csthcm,
    phi,
    e_hel,
    delta: float,
    m_pi: float,
    scheme: str = "linear",
) -> tuple[jnp.ndarray, jnp.ndarray]:
    """Non-radiative cross section folded with the soft-photon correction.

    Covers the ``ek < delta`` region of the integrand, where the radiated
    photon is unresolvable: the event is instead weighted by the
    vertex-corrected cross section times the standard soft factor
    ``exp(del_inf)`` (Mo & Tsai), averaged over ``0 < ek < delta`` and the full
    hadronic solid angle.

    Returns
    -------
    (sigma, asym_p)
        Cross section in the same units as :func:`motsa_sigma`, and the
        single-spin beam asymmetry.
    """
    sp = es * ep - ps * pp * cst0
    q0 = es - ep
    nu = (w_sq + q2 - M_N**2) / (2.0 * M_N)
    eps = _epsilon(es, ep, cst0, nu, q2)

    epw = jnp.sqrt(jnp.maximum(w_sq, 0.0))
    resp = response_functions(grid, q2, epw, csthcm, phi, eps, e_hel, m_pi, scheme=scheme)
    epeps = _epsilon(es, ep, cst0, q0, q2)

    f, g, _ = _fg_amplitudes(resp, w_sq, q2, m_pi)

    signr = 2.0 * (ALPHA * ep / _safe(q2)) ** 2 * (
        f * M_N * jnp.cos(th0 / 2.0) ** 2 + 2.0 * g * jnp.sin(th0 / 2.0) ** 2 / M_N
    )
    signr = signr * (1.0 + polarization_modulation(resp, epeps, phi, e_hel))

    # Soft-gluon and vertex corrections (Mo & Tsai, Appendix C).
    log_2sp = jnp.log(2.0 * sp / M_E**2)
    deltar = -(ALPHA / PI) * (
        28.0 / 9.0 - (13.0 / 6.0) * log_2sp - spence((ep - es) / ep) - spence((es - ep) / es)
    )
    delinf = -(ALPHA / PI) * jnp.log(es * ep / delta**2) * (log_2sp - 1.0)

    sigr1 = signr * (1.0 + deltar) * jnp.exp(jnp.clip(delinf, -80.0, 80.0))
    # Average over the soft-photon energy range and the full hadronic solid angle.
    return sigr1 / (delta * 4.0 * PI), resp.asym_p
