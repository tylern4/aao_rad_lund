"""Event kinematics: hadronic final state, missing mass, photon momentum.

Batched port of ``missm()`` from ``aao_rad.f90``.  Everything here is pure
elementwise arithmetic, so it maps straight onto ``jax.vmap``.
"""

from __future__ import annotations

from typing import NamedTuple

import jax.numpy as jnp

from .constants import M_E, M_N, PI

__all__ = ["HadronicFinalState", "hadronic_final_state", "photon_momentum"]


class HadronicFinalState(NamedTuple):
    """Lab-frame final state for one event."""

    p_p: jnp.ndarray  # proton (or the other nucleon) momentum, (B, 3)
    p_pi: jnp.ndarray  # pion momentum, (B, 3)
    e_prot: jnp.ndarray
    e_pi: jnp.ndarray
    w_real: jnp.ndarray
    mm2: jnp.ndarray


def photon_momentum(
    ep_s: jnp.ndarray,
    ep_p: jnp.ndarray,
    theta_e: jnp.ndarray,
    e_gamma: jnp.ndarray,
    cos_theta_k: jnp.ndarray,
    phi_k: jnp.ndarray,
) -> tuple[jnp.ndarray, jnp.ndarray]:
    """Photon 3-momentum and the hadronic-invariant-mass ingredients.

    Returns
    -------
    (k_vec, w2)
        ``k_vec`` has shape ``(B, 3)``; ``w2`` is the *non-radiative* hadronic
        invariant mass squared ``W^2``.
    """
    p_s = jnp.sqrt(jnp.maximum(ep_s**2 - M_E**2, 0.0))
    p_p = jnp.sqrt(jnp.maximum(ep_p**2 - M_E**2, 0.0))
    c_the = jnp.cos(theta_e)
    s_the = jnp.sin(theta_e)

    qx = -p_p * s_the
    qz = p_s - p_p * c_the
    q_vec = jnp.sqrt(qx**2 + qz**2)

    q2 = 2.0 * ep_s * ep_p - 2.0 * p_s * p_p * c_the - 2.0 * M_E**2
    w2 = M_N**2 - q2 + 2.0 * M_N * (ep_s - ep_p)

    cstk = jnp.clip(cos_theta_k, -1.0, 1.0)
    sntk = jnp.sqrt(jnp.maximum(1.0 - cstk**2, 0.0))
    cstq = jnp.where(q_vec > 0.0, qz / jnp.where(q_vec > 0.0, q_vec, 1.0), 0.0)
    sntq = jnp.sqrt(jnp.maximum(1.0 - cstq**2, 0.0))
    csphk, snphk = jnp.cos(phi_k), jnp.sin(phi_k)

    kx = e_gamma * (sntk * csphk * cstq - cstk * sntq)
    ky = e_gamma * sntk * snphk
    kz = e_gamma * (cstk * cstq + sntk * csphk * sntq)
    return jnp.stack([kx, ky, kz], axis=-1), w2


def hadronic_final_state(
    e_beam: jnp.ndarray,
    ep_s: jnp.ndarray,
    ep_p: jnp.ndarray,
    theta_e: jnp.ndarray,
    e_gamma: jnp.ndarray,
    cos_theta_k: jnp.ndarray,
    phi_k: jnp.ndarray,
    cos_theta_cm: jnp.ndarray,
    phi_cm_deg: jnp.ndarray,
    m_pi: float,
    is_pi0: jnp.ndarray | float = 1.0,
) -> HadronicFinalState:
    """Two-body decay in the hadronic frame, boosted to the lab.

    Mirrors ``missm()``: the pion is generated isotropically in the hadronic
    CM frame with ``cos theta*`` / ``phi*``, then boosted along the resonance
    momentum (the vector momentum of the hadronic system) and finally rotated
    into the beam frame.

    Notes
    -----
    The Fortran version emitted ``mm2 = 0`` and bailed out when the hadronic
    mass fell below threshold.  Here the same events simply get ``w_real = 0``
    and ``mm2 = 0`` and are rejected downstream by the missing-mass cut, which
    keeps the kernel branch-free.
    """
    p_s = jnp.sqrt(jnp.maximum(ep_s**2 - M_E**2, 0.0))
    p_p = jnp.sqrt(jnp.maximum(ep_p**2 - M_E**2, 0.0))
    c_the = jnp.cos(theta_e)
    s_the = jnp.sin(theta_e)
    qx = -p_p * s_the
    qz = p_s - p_p * c_the

    k_vec, w2 = photon_momentum(ep_s, ep_p, theta_e, e_gamma, cos_theta_k, phi_k)

    # True hadronic invariant mass once the photon is accounted for.
    k_dot_q = k_vec[..., 0] * qx + k_vec[..., 2] * qz
    w_real_sq = w2 - 2.0 * e_gamma * ((ep_s - ep_p) + M_N) + 2.0 * k_dot_q
    w_min = M_N + m_pi
    w_real = jnp.where(
        w_real_sq > w_min**2, jnp.sqrt(jnp.maximum(w_real_sq, 0.0)), 0.0
    )
    e_w = (ep_s - ep_p) + M_N - e_gamma  # hadronic system lab energy

    # Lab frame of the resonance momentum.
    pwr_x = qx - k_vec[..., 0]
    pwr_y = -k_vec[..., 1]
    pwr_z = qz - k_vec[..., 2]
    pwr = jnp.sqrt(pwr_x**2 + pwr_y**2 + pwr_z**2)

    beta = jnp.where(e_w != 0.0, pwr / jnp.where(e_w != 0.0, e_w, 1.0), 0.0)
    gamma = jnp.where(w_real > 0.0, e_w / jnp.where(w_real > 0.0, w_real, 1.0), 0.0)

    # Rotation basis with the resonance momentum along the local z axis.
    pfac = jnp.sqrt(pwr_y**2 + pwr_z**2)
    safe_pfac = jnp.where(pfac > 0.0, pfac, 1.0)
    safe_pwr = jnp.where(pwr > 0.0, pwr, 1.0)
    cxx = pfac / safe_pwr
    cxy = -pwr_x * pwr_y / safe_pfac / safe_pwr
    cxz = -pwr_x * pwr_z / safe_pfac / safe_pwr
    cyx = jnp.zeros_like(pwr_x)
    cyy = pwr_z / safe_pfac
    cyz = -pwr_y / safe_pfac
    czx = pwr_x / safe_pwr
    czy = pwr_y / safe_pwr
    czz = pwr_z / safe_pwr

    # Two-body decay momentum in the hadronic CM frame.
    pstar_sq = (
        (w_real**2 - M_N**2 - m_pi**2) ** 2 / 4.0 - (M_N * m_pi) ** 2
    ) / jnp.maximum(w_real**2, 1e-12)
    pstar = jnp.sqrt(jnp.maximum(pstar_sq, 0.0))
    e_pcm = jnp.sqrt(pstar**2 + M_N**2)
    e_picm = jnp.sqrt(pstar**2 + m_pi**2)

    snthcm = jnp.sqrt(jnp.maximum(1.0 - jnp.clip(cos_theta_cm, -1.0, 1.0) ** 2, 0.0))
    phi_cm = phi_cm_deg * (PI / 180.0)
    csphi, snphi = jnp.cos(phi_cm), jnp.sin(phi_cm)

    cthcm = jnp.clip(cos_theta_cm, -1.0, 1.0)
    ppi_wx = pstar * snthcm * csphi
    ppi_wy = pstar * snthcm * snphi
    ppi_wz = gamma * (pstar * cthcm + beta * e_picm)
    e_pi = gamma * (e_picm + beta * pstar * cthcm)

    # Boosted pion, then rotated into the beam frame.
    ppix = ppi_wx * cxx + ppi_wy * cyx + ppi_wz * czx
    ppiy = ppi_wx * cxy + ppi_wy * cyy + ppi_wz * czy
    ppiz = ppi_wx * cxz + ppi_wy * cyz + ppi_wz * czz

    pp_wx, pp_wy = -ppi_wx, -ppi_wy
    pp_wz = gamma * beta * w_real - ppi_wz
    e_prot = gamma * w_real - e_pi

    ppx = pp_wx * cxx + pp_wy * cyx + pp_wz * czx
    ppy = pp_wx * cxy + pp_wy * cyy + pp_wz * czy
    ppz = pp_wx * cxz + pp_wy * cyz + pp_wz * czz

    # Experimental ("reconstructed") missing mass from the charged hadron.
    q2_exp = 2.0 * e_beam * ep_p * (1.0 - c_the)
    qz_exp = e_beam - p_p * c_the
    nu = e_beam - ep_p

    # Which track is charged depends only on the production channel, but inside
    # jit the pion mass is a traced array, so this cannot be a Python `if`.
    # Evaluate both reconstructions and select with jnp.where.
    q_dot_prot = qx * ppx + qz_exp * ppz
    mm2_pi0 = (
        -q2_exp + 2.0 * M_N**2 + 2.0 * M_N * (nu - e_prot) - 2.0 * nu * e_prot
        + 2.0 * q_dot_prot
    )
    q_dot_pi = qx * ppix + qz_exp * ppiz
    mm2_pip = (
        -q2_exp + M_N**2 + m_pi**2 + 2.0 * M_N * (nu - e_pi) - 2.0 * nu * e_pi
        + 2.0 * q_dot_pi
    )
    mm2 = jnp.where(is_pi0 > 0.5, mm2_pi0, mm2_pip)
    mm2 = jnp.where(w_real > 0.0, mm2, 0.0)

    return HadronicFinalState(
        p_p=jnp.stack([ppx, ppy, ppz], axis=-1),
        p_pi=jnp.stack([ppix, ppiy, ppiz], axis=-1),
        e_prot=e_prot,
        e_pi=e_pi,
        w_real=w_real,
        mm2=mm2,
    )
