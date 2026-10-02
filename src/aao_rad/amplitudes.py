"""Multipole -> CGLN -> helicity amplitude chain.

This is a batched transcription of ``multipole_amps.f90`` + ``cgln_amps.f90`` +
``helicity_amps.f90`` + ``legendre.f90``.  The original called
``multipole_amps`` twice per event (once directly, once from ``cgln_amps``);
here it is called once.

All arrays carry a leading batch axis, so every routine works on an arbitrary
number of kinematic points at once.
"""

from __future__ import annotations

import jax.numpy as jnp

__all__ = [
    "N_WAVE",
    "MAX_L",
    "cgln_amplitudes",
    "helicity_amplitudes",
    "legendre_polynomials",
    "multipole_amplitudes",
]

N_WAVE = 6
"""Partial waves L = 0..5 present in the MAID07 tables."""

MAX_L = 5
"""Highest L summed (``wave_L`` in the Fortran code)."""

#: ``nu_cm * S = |q| * L`` conversion.  MAID tabulates S-multipoles in units of
#: 1e-3/m_pi+; Tiator's (and hence the AO) formalism wants L-multipoles scaled
#: to sqrt(microbarn).
_MAID_TO_SQRT_MUB = 0.141383

_ROOT2 = jnp.sqrt(2.0)
_LF = jnp.arange(N_WAVE, dtype=jnp.float32)
_L = jnp.arange(N_WAVE)
#: Fortran loop bounds, as masks.  ``cgln_amps.f90`` starts ff3/ff6 at ``l = 1``
#: and ff4 at ``l = 2``; the skipped terms are *not* zero (P_1^2 is the constant
#: 3 and sp(0) is tabulated), so they have to be masked off explicitly.
_L_GE_1 = (_L >= 1)[None, :]
_L_GE_2 = (_L >= 2)[None, :]


def _shift(a: jnp.ndarray, amount: int) -> jnp.ndarray:
    """``out[:, l] = a[:, l - amount]`` (Fortran ``pol(l-amount, k)``).

    The ``l < amount`` entries hold junk from negative Fortran indices; callers
    mask them out with a precomputed constant.
    """
    if amount == 1:
        return jnp.concatenate([a[:, :1], a[:, :-1]], axis=1)
    if amount == 2:
        return jnp.concatenate([a[:, :2], a[:, :-2]], axis=1)
    raise ValueError("only shifts of 1 and 2 occur")


def legendre_polynomials(cos_theta: jnp.ndarray) -> jnp.ndarray:
    """Associated Legendre polynomials used by the AO helicity formalism.

    Returns shape ``(batch, N_WAVE + 2, 2)`` where ``[..., 0]`` is
    ``P_l^1(cos theta)`` and ``[..., 1]`` is ``P_l^2(cos theta)`` for
    ``l = 0 .. N_WAVE + 1``.
    """
    x = jnp.clip(cos_theta, -1.0, 1.0)
    x2 = x * x

    p1 = [
        jnp.zeros_like(x),
        jnp.ones_like(x),
        3.0 * x,
        (15.0 * x2 - 3.0) / 2.0,
        (35.0 * x2 * x - 15.0 * x) / 2.0,
        (315.0 * x2 * x2 - 210.0 * x2 + 15.0) / 8.0,
        (693.0 * x2 * x2 * x - 630.0 * x2 * x + 105.0 * x) / 8.0,
        (3003.0 * x2 * x2 * x2 - 3465.0 * x2 * x2 + 945.0 * x2 - 35.0) / 16.0,
    ]
    p2 = [
        jnp.zeros_like(x),
        jnp.zeros_like(x),
        jnp.full_like(x, 3.0),
        15.0 * x,
        (105.0 * x2 - 15.0) / 2.0,
        (315.0 * x2 * x - 105.0 * x) / 2.0,
        (3465.0 * x2 * x2 - 1890.0 * x2 + 105.0) / 8.0,
        (9009.0 * x2 * x2 * x - 6930.0 * x2 * x + 945.0 * x) / 8.0,
    ]
    n = N_WAVE + 2
    return jnp.stack(
        [jnp.stack(p1[:n], axis=-1), jnp.stack(p2[:n], axis=-1)], axis=-1
    )


def multipole_amplitudes(
    amps: jnp.ndarray, nu_cm: jnp.ndarray, qv_mag_cm: jnp.ndarray
) -> tuple[jnp.ndarray, ...]:
    """Split the 62 tabulated components into the six multipole families.

    Parameters
    ----------
    amps
        Interpolated table values, shape ``(batch, 62)``.
    nu_cm, qv_mag_cm
        Centre-of-mass energy transfer and three-momentum transfer, shape
        ``(batch,)``.  Multipoles are converted from Sato-Lee's S convention
        to Tiator's L convention via ``nu_cm / |q|``.

    Returns
    -------
    (sp, sm, ep, em, mp, mm)
        Complex arrays of shape ``(batch, N_WAVE)``.
    """
    z = jnp.asarray(amps[:, 0::2], jnp.complex64) + 1j * jnp.asarray(amps[:, 1::2], jnp.complex64)

    factor = (_MAID_TO_SQRT_MUB * (nu_cm / qv_mag_cm))[:, None]
    efactor = _MAID_TO_SQRT_MUB

    # Index layout transcribed from multipole_amps.f90 lines 41-84:
    #   z = [S0+..S5+, S1-..S5-, E0+..E5+, E2-..E5-, M1+..M5+, M1-..M5-]
    sp, sm = z[:, 0:6], z[:, 6:11]
    el, em = z[:, 11:17], z[:, 17:21]
    mp, mm = z[:, 21:26], z[:, 26:31]

    # The Fortran hard-sets the waves the table does not contain to zero:
    # S0-, E0-, E1-, M0+ and M0- are absent (see mpintp.inc), so each of these
    # families gets leading zeros padded on to reach N_WAVE entries.
    sm = jnp.concatenate([jnp.zeros_like(sm[:, :1]), sm], axis=1)
    em = jnp.concatenate([jnp.zeros_like(em[:, :2]), em], axis=1)
    mp = jnp.concatenate([jnp.zeros_like(mp[:, :1]), mp], axis=1)
    mm = jnp.concatenate([jnp.zeros_like(mm[:, :1]), mm], axis=1)

    return (sp * factor, sm * factor, el * efactor, em * efactor, mp * efactor, mm * efactor)


def cgln_amplitudes(
    pol: jnp.ndarray, sp: jnp.ndarray, sm: jnp.ndarray, ep: jnp.ndarray,
    em: jnp.ndarray, mp: jnp.ndarray, mm: jnp.ndarray
) -> tuple[jnp.ndarray, ...]:
    """CGLN amplitudes ``ff1..ff6`` in Tiator's helicity formalism.

    ``pol`` is the output of :func:`legendre_polynomials`.  The Fortran's
    ``if (l < 2) ... else ...`` branches are merged by summing both
    contributions and masking the out-of-range one, which is exactly what the
    original intended (the masked terms would read uninitialised ``pol(-1, k)``).
    """
    p1 = pol[..., 0]  # P_l^1
    p2 = pol[..., 1]  # P_l^2

    same1 = p1[:, :N_WAVE]  # P_l^1
    same2 = p2[:, :N_WAVE]  # P_l^2
    up1 = p1[:, 1 : N_WAVE + 1]  # P_{l+1}^1
    up2 = p2[:, 1 : N_WAVE + 1]  # P_{l+1}^2
    # Down-shifted: index l reads P_{l-1}, valid only for l >= 2.  Truncating to
    # N_WAVE leaves the leading entries holding P_0/P_1, which the masks drop.
    dn1 = _shift(p1, 1)[:, :N_WAVE]  # P_{l-1}^1
    dn2 = _shift(p2, 1)[:, :N_WAVE]  # P_{l-1}^2

    ge1 = _L_GE_1.astype(up1.dtype)
    ge2 = _L_GE_2.astype(up1.dtype)

    # cgln_amps.f90 lines 23-57, loop by loop.
    ff1 = jnp.sum(
        (_LF * mp + ep) * up1 + ge2 * (((_LF + 1.0) * mm + em) * dn1), axis=1
    )
    ff2 = jnp.sum(((_LF + 1.0) * mp + _LF * mm) * same1, axis=1)
    ff3 = jnp.sum(ge1 * (ep - mp) * up2 + ge2 * ((em + mm) * dn2), axis=1)
    ff4 = jnp.sum(ge2 * (mp - ep - mm - em) * same2, axis=1)
    ff5 = jnp.sum(((_LF + 1.0) * sp) * up1 + ge2 * (-_LF * sm * dn1), axis=1)
    ff6 = jnp.sum(ge1 * (_LF * sm - (_LF + 1.0) * sp) * same1, axis=1)

    return ff1, ff2, ff3, ff4, ff5, ff6


def helicity_amplitudes(
    ff1: jnp.ndarray, ff2: jnp.ndarray, ff3: jnp.ndarray, ff4: jnp.ndarray,
    ff5: jnp.ndarray, ff6: jnp.ndarray, cos_theta_cm: jnp.ndarray
) -> tuple[jnp.ndarray, ...]:
    """Helicity amplitudes ``hh1..hh6`` from the CGLN amplitudes."""
    c = jnp.clip(cos_theta_cm, -1.0, 1.0)
    theta = jnp.arccos(c)
    s = jnp.sin(theta)
    s2 = jnp.sin(theta / 2.0)
    c2 = jnp.cos(theta / 2.0)

    hh1 = -s * c2 * (ff3 + ff4) / _ROOT2
    hh2 = c2 * ((ff2 - ff1) + 0.5 * (1.0 - c) * (ff3 - ff4)) * _ROOT2
    hh3 = s * s2 * (ff3 - ff4) / _ROOT2
    hh4 = s2 * ((ff2 + ff1) + 0.5 * (1.0 + c) * (ff3 + ff4)) * _ROOT2
    hh5 = c2 * (ff5 + ff6)
    hh6 = s2 * (ff6 - ff5)
    return hh1, hh2, hh3, hh4, hh5, hh6
