"""Batched, GPU-resident Monte-Carlo event generator.

The Fortran version was a scalar importance-sampling loop: draw one point,
evaluate the integrand, accept/reject, repeat.  Almost all of the cost lives in
the response-function evaluation -- a table interpolation followed by small sums
over six partial waves -- which is perfectly parallel and perfectly
embarrassing for a CPU.

This module keeps the *same* sampling and weighting scheme (so the generated
physics is unchanged) but evaluates a whole batch of trial points at once:

* every kinematic quantity becomes an elementwise array operation;
* the kinematic rejections of the original ``go to`` chain become boolean masks;
* acceptance becomes a single ``uniform < weight / weight_max`` comparison;
* accepted events are compacted and appended to an output buffer.

Two further changes remove real inefficiencies in the original:

* the weight maximum is estimated on the device from a large batch of trial
  points instead of being read from a file, which removes a whole class of
  "the run stalled because sigr_max was too small" failures;
* the photon energy is drawn from the *truncated* exponential (see
  :attr:`~aao_rad.config.GeneratorConfig.ek_sampling`), so no trial is thrown
  away merely for landing above ``ek_max``.
"""

from __future__ import annotations

import logging
import os
import time
from collections.abc import Iterator
from dataclasses import dataclass
from typing import Any, NamedTuple

import jax
import jax.numpy as jnp
import numpy as np

from .config import GeneratorConfig
from .constants import (
    BFAC,
    HYDROGEN_RAD,
    M_E,
    M_N,
    PI,
    TABLE_Q2_MAX,
    TABLE_W_MAX,
    pi_masses,
)
from .interpolate import InterpolationGrid
from .kinematics import hadronic_final_state
from .motsa import motsa_sigma, non_radiative_sigma
from .table import MaidTable, load_table

__all__ = [
    "EVENT_COLUMNS",
    "EVENT_DTYPE",
    "EventGenerator",
    "GenerationStats",
    "Integrand",
    "Kinematics",
    "build_grid",
    "build_kinematics",
    "draw_kinematics",
    "draw_photon",
    "finalize",
    "integrand",
]

log = logging.getLogger(__name__)

EVENT_COLUMNS: tuple[str, ...] = (
    "es", "ep", "theta", "w", "w_real",
    "ppx", "ppy", "ppz", "eprot",
    "ppix", "ppiy", "ppiz", "epi",
    "csthcm", "phicm", "mm2",
    "eg", "cstk", "phik",
    "qx", "qz", "q0", "csthe",
    "egx", "egy", "egz",
    "vx", "vy", "vz",
    "q2", "e_hel", "asym_p",
)
"""The 32 n-tuple variables of the original generator, in the same order."""

EVENT_DTYPE = np.dtype([(name, np.float32) for name in EVENT_COLUMNS])

#: Below this the integrand is numerically zero and the trial is discarded.
#: (The Fortran returned 1e-31 in those branches and let the point through.)
_NEGLIGIBLE = 1.0e-30


# ---------------------------------------------------------------------------
# Static parameters
# ---------------------------------------------------------------------------
class Kinematics(NamedTuple):
    """Constants derived once from the configuration; a JAX pytree."""

    e_beam: jnp.ndarray
    m_pi: jnp.ndarray
    is_pi0: jnp.ndarray
    m_exp: jnp.ndarray
    wg: jnp.ndarray
    ep_min: jnp.ndarray
    ep_max: jnp.ndarray
    ep_range: jnp.ndarray
    uq2_min: jnp.ndarray
    uq2_range: jnp.ndarray
    q2_hard_max: jnp.ndarray
    t_len_cm: jnp.ndarray
    t_len_rl: jnp.ndarray
    r_targ: jnp.ndarray
    vertex_x: jnp.ndarray
    vertex_y: jnp.ndarray
    vertex_z0: jnp.ndarray
    reg: jnp.ndarray  # cumulative region edges, (4,)
    cos_rng: jnp.ndarray
    k_exp: jnp.ndarray
    delta: jnp.ndarray
    w_cut: jnp.ndarray
    """Reject trials with ``W`` above this; <= 0 disables the cut (``clamp``)."""
    q2_cut: jnp.ndarray
    """Reject trials with ``Q^2`` above this; <= 0 disables the cut (``clamp``)."""
    mm_cut: jnp.ndarray
    polarized: jnp.ndarray
    truncated_ek: jnp.ndarray


def build_grid(
    channel: int,
    *,
    parms_dir: str | None = None,
    table_dir: str | None = None,
    scheme: str = "linear",
    dtype: Any = jnp.float32,
    use_cache: bool = True,
    table: MaidTable | None = None,
) -> InterpolationGrid:
    """Load a MAID table and package it for jitted, batched interpolation."""
    from .interpolate import natural_spline_2nd_derivs

    if table is None:
        table = load_table(
            channel, parms_dir=parms_dir, table_dir=table_dir, use_cache=use_cache
        )

    values = jnp.asarray(table.amps, dtype=dtype)
    axis1 = jnp.asarray(table.q2, dtype=dtype)
    axis2 = jnp.asarray(table.w, dtype=dtype)

    d2_1 = d2_2 = None
    if scheme == "spline":
        log.info("precomputing spline second derivatives (slowest part of start-up)")
        # natural_spline_2nd_derivs works on a 2-D array whose *first* axis is
        # the one being splined, so the target axis has to lead.  Move it back
        # afterwards to give (62, nq2, nw), matching `values`.
        d2_2_np = np.moveaxis(
            natural_spline_2nd_derivs(np.moveaxis(table.amps, 2, 0), table.w), 0, 2
        )
        d2_1_np = np.moveaxis(
            natural_spline_2nd_derivs(np.moveaxis(table.amps, 1, 0), table.q2), 0, 1
        )
        assert d2_1_np.shape == table.amps.shape, d2_1_np.shape
        assert d2_2_np.shape == table.amps.shape, d2_2_np.shape
        d2_1 = jnp.asarray(d2_1_np, dtype=dtype)
        d2_2 = jnp.asarray(d2_2_np, dtype=dtype)

    return InterpolationGrid(
        values=values,
        axis1=axis1,
        axis2=axis2,
        uniform1=table.uniform,
        uniform2=table.uniform,
        d2_1=d2_1,
        d2_2=d2_2,
    )


def build_kinematics(cfg: GeneratorConfig, dtype: Any = jnp.float32) -> Kinematics:
    """Derive the pre-computed kinematic constants from a configuration."""
    m_pi, m_exp = pi_masses(cfg.channel)
    wg = M_N + m_pi + 0.0005
    e_beam = float(cfg.beam_energy)

    # Cut off q2 at the value for 90 degree elastic scattering.
    q2_hard_max = 4.0 * e_beam**2 * 0.5 / (1.0 + 2.0 * e_beam * 0.5 / M_N)
    q2_max = min(cfg.q2_max, q2_hard_max)

    # Scattered electron window, also limited by the pion threshold.
    ep_max = min(cfg.ep_max, e_beam - (wg**2 + cfg.q2_min - M_N**2) / (2.0 * M_N))
    if ep_max <= cfg.ep_min:
        raise ValueError(
            f"the requested scattered-electron window ({cfg.ep_min}, {cfg.ep_max}] GeV) "
            f"is empty for a beam energy of {e_beam} GeV"
        )

    j = lambda v: jnp.asarray(v, dtype)  # noqa: E731
    return Kinematics(
        e_beam=j(e_beam),
        m_pi=j(m_pi),
        is_pi0=j(1.0 if cfg.channel in (1, 5) else 0.0),
        m_exp=j(m_exp),
        wg=j(wg),
        ep_min=j(cfg.ep_min),
        ep_max=j(ep_max),
        ep_range=j(ep_max - cfg.ep_min),
        uq2_min=j(1.0 / q2_max),
        uq2_range=j(1.0 / cfg.q2_min - 1.0 / q2_max),
        q2_hard_max=j(q2_hard_max),
        t_len_cm=j(cfg.target_length_cm),
        # The bremsstrahlung exponent needs the path length in radiation
        # lengths, not cm: the original converts with
        #   targs = bfac * targs[cm] / hydrogen_rad
        # before evaluating eloss = xs**(1/targs).
        t_len_rl=j(cfg.target_length_cm * BFAC / HYDROGEN_RAD),
        r_targ=j(cfg.target_radius_cm),
        vertex_x=j(cfg.beam_x_cm),
        vertex_y=j(cfg.beam_y_cm),
        vertex_z0=j(cfg.beam_z_cm),
        reg=j(np.cumsum(np.asarray(cfg.regions, dtype=np.float64))),
        cos_rng=j(cfg.cos_step),
        k_exp=j(cfg.k_exp),
        delta=j(cfg.min_photon_energy),
        # Table edges.  The original saturates outside them (see multipole_amps
        # .f90); we reject instead so the sample stays physical.  w_cut <= 0
        # disables the W cut, which is how "clamp" is requested.
        w_cut=j(0.0 if cfg.w_max == "clamp" else float(cfg.w_max or TABLE_W_MAX)),
        q2_cut=j(0.0 if cfg.w_max == "clamp" else TABLE_Q2_MAX),
        mm_cut=j(cfg.missing_mass_cut),
        polarized=j(1.0 if cfg.polarized_beam else 0.0),
        truncated_ek=j(1.0 if cfg.ek_sampling == "truncated" else 0.0),
    )


# ---------------------------------------------------------------------------
# Trial sampling (all operations are elementwise over the batch)
# ---------------------------------------------------------------------------
def draw_kinematics(key, kin: Kinematics, n: int) -> dict[str, jnp.ndarray]:
    """Electron-side kinematics plus the in-target energy loss.

    Mirrors statements 20-14 of the original main loop.  Rejected points are
    flagged in ``ok`` rather than re-drawn, which keeps the kernel branch-free;
    wasted slots simply contribute zero weight.
    """
    k_pos, k_x, k_test, k_uq2, k_ep, _ = jax.random.split(key, 6)
    dt = kin.e_beam.dtype
    u = lambda k: jax.random.uniform(k, (n,), dtype=dt)  # noqa: E731

    # --- interaction point and in-target bremsstrahlung loss ----------------
    # The original picks targs uniformly on [0, L] and then loops
    #     xs = rand;  eloss = xs**(1/targs);  if (rand > 1-eloss) retry
    # at fixed targs.  That retry has acceptance probability 1/(1+targs), so
    # the *accepted* density of targs is proportional to 1/(1+targs) on [0, L],
    # whose inverse CDF is the closed form below.  Sampling it directly is
    # exact, branch-free (so it survives jit) and needs no retry loop.
    #
    # Given an accepted targs, the survival uniforms give
    #     v = 1 - xtest ~ U(0,1),  xs ~ U(0, v**targs)
    #     => eloss = xs**(1/targs) = v * u3**(1/targs)
    # which reproduces the original's joint (targs, eloss) distribution while
    # guaranteeing acceptance by construction.
    len_rl = kin.t_len_rl
    targs = jnp.expm1(u(k_pos) * jnp.log1p(len_rl))
    # Back to cm for the vertex position (the original samples in cm and only
    # converts to radiation lengths afterwards).
    cm_per_rl = jnp.where(len_rl > 0.0, kin.t_len_cm / jnp.maximum(len_rl, 1e-30), 0.0)
    vertex_z = targs * cm_per_rl + kin.vertex_z0 - kin.t_len_cm / 2.0

    inv_targs = 1.0 / jnp.maximum(targs, 1e-12)
    eloss = (1.0 - u(k_test)) * jnp.power(u(k_x), inv_targs)
    eloss = jnp.where(eloss <= 1.0, eloss, 1.0)

    es = kin.e_beam * (1.0 - eloss)
    ok = es >= kin.e_beam / 4.0

    # --- Q^2 and scattered electron energy ---------------------------------
    uq2 = kin.uq2_min + kin.uq2_range * u(k_uq2)
    q2 = 1.0 / uq2
    ep = kin.ep_max - kin.ep_range * u(k_ep)
    ok = ok & (ep >= kin.ep_min) & (q2 > 0.0) & (q2 <= kin.q2_hard_max)

    s_kin = q2 / (4.0 * es * ep)
    ok = ok & (s_kin <= 0.5)  # cut off scattering at 90 degrees

    th0 = 2.0 * jnp.arcsin(jnp.sqrt(jnp.clip(s_kin, 0.0, 0.5)))
    cst0 = jnp.cos(th0)
    snt0 = jnp.sin(th0)

    ps = jnp.sqrt(jnp.maximum(es**2 - M_E**2, 0.0))
    pp = jnp.sqrt(jnp.maximum(ep**2 - M_E**2, 0.0))

    # Scattered electron above the pion threshold for this angle?
    ep_test = (M_N**2 + 2.0 * M_N * es - kin.wg**2) / (2.0 * (M_N + 2.0 * es * s_kin))
    ok = ok & (ep <= ep_test)

    w_sq = M_N**2 + 2.0 * M_N * (es - ep) - q2
    ok = ok & (w_sq >= M_N**2)
    epw = jnp.sqrt(jnp.maximum(w_sq, 0.0))
    ok = ok & (epw >= kin.wg + 0.002)
    # Stay inside the tabulated grid.  The original instead saturates the
    # lookup at the W = 2 / Q^2 = 5 edges, which returns a frozen response for
    # most of its own sample (see GeneratorConfig.w_max).
    ok = ok & ((kin.w_cut <= 0.0) | (epw <= kin.w_cut))
    ok = ok & ((kin.q2_cut <= 0.0) | (q2 <= kin.q2_cut))

    # --- Mo & Tsai intermediate quantities ---------------------------------
    u0 = es - ep + M_N
    pu_sq = ps**2 + pp**2 - 2.0 * ps * pp * cst0
    ok = ok & (pu_sq > 0.0)
    pu = jnp.sqrt(jnp.maximum(pu_sq, 0.0))
    uu = u0**2 - pu**2
    safe_pu = jnp.where(pu > 0.0, pu, 1.0)
    csths = (ps - pp * cst0) / safe_pu
    csthp = (ps * cst0 - pp) / safe_pu
    snths = jnp.sqrt(jnp.maximum(1.0 - csths**2, 0.0))
    snthp = jnp.sqrt(jnp.maximum(1.0 - csthp**2, 0.0))
    ok = ok & (snths > 0.0) & (snthp > 0.0)

    return dict(
        es=es, ep=ep, th0=th0, theta=th0 * (180.0 / PI), cst0=cst0, snt0=snt0,
        ps=ps, pp=pp, q2=q2, q0=es - ep, w_sq=w_sq, epw=epw,
        u0=u0, pu=pu, uu=uu, csths=csths, csthp=csthp, snths=snths, snthp=snthp,
        qvecx=-pp * snt0, qvecz=ps - pp * cst0,
        vertex_z=vertex_z, targs=targs, targs_cm=targs * cm_per_rl,
        s_kin=s_kin, ok=ok,
    )


def draw_photon(key, kin: Kinematics, kv: dict, n: int) -> dict[str, jnp.ndarray]:
    """Importance-sample the photon direction and energy (statements 10-14)."""
    # One key per independent variate.  ``k_cm`` and ``k_cm_phi`` must stay
    # separate: cos(theta*) and phi* are two independent uniforms in the
    # Fortran (aao_rad.f90:779-780) and they enter mm2 through different
    # combinations of the pion momentum components, so sharing one key does not
    # cancel -- it confines the decay direction to a curve instead of the
    # sphere.  Both marginals stay uniform under such a constraint, which is
    # why the error is invisible in either histogram alone.
    keys = jax.random.split(key, 10)
    dt = kin.e_beam.dtype
    u = lambda k: jax.random.uniform(k, (n,), dtype=dt)  # noqa: E731
    k_reg, k_sign, k_pos, k_pos2, k_phi, k_ek, k_cm, k_cm_phi, k_spin, k_redraw = keys

    cstk1 = jnp.maximum(kv["csths"], kv["csthp"])
    cstk2 = jnp.minimum(kv["csths"], kv["csthp"])

    cs_range = jnp.minimum(jnp.minimum(kin.cos_rng, 1.0 - cstk1), 0.5 * (cstk1 - cstk2))
    cs_range = jnp.maximum(cs_range, 0.0)
    cs_rngb = jnp.maximum(jnp.minimum(kin.cos_rng / 40.0, cs_range / 5.0), 1e-9)

    r1, r2, r3, r4 = kin.reg[0], kin.reg[1], kin.reg[2], kin.reg[3]

    csran = u(k_reg)
    rn1 = jnp.where(u(k_sign) > 0.5, 1.0, -1.0)
    rn2 = u(k_pos)
    rn3 = u(k_pos2)

    band = (2.0 * rn2 - 1.0) * cs_rngb
    spread = rn1 * (cs_rngb + rn3 * (cs_range - cs_rngb))

    cstk = jnp.where(
        csran < r1, cstk1 + band,
        jnp.where(csran < r2, cstk2 + band,
        jnp.where(csran < r3, cstk1 + spread,
        jnp.where(csran < r4, cstk2 + spread, 2.0 * rn3 - 1.0))),
    )
    # Region 5 gives back everything the narrow regions did not claim, less the
    # overlap with the bands it must avoid: (1 - csrnge*delphi/pi)/(1 - reg4).
    mcfac = jnp.where(
        csran < r1, cs_rngb / r1,
        jnp.where(csran < r2, cs_rngb / (r2 - r1),
        jnp.where(csran < r3, (cs_range - cs_rngb) / (r3 - r2),
        jnp.where(csran < r4, (cs_range - cs_rngb) / (r4 - r3),
                  (1.0 - cs_range / 9.0) / (1.0 - r4)))),
    )

    del_phi = PI / 9.0
    narrow = csran < r4
    ur = u(k_phi)
    phik = jnp.where(narrow, (ur - 0.5) * del_phi, 2.0 * PI * (ur - 0.5))
    mpfac = jnp.where(narrow, del_phi / (2.0 * PI), 1.0)

    # Region 5 must not double-sample the narrow bands around the beam lines.
    # aao_rad.f90:721 tests the two bands with ``.or.``: being inside *either*
    # one is enough to collide.
    collide = (
        (jnp.abs(cstk - cstk1) < cs_range)
        | (jnp.abs(cstk - cstk2) < cs_range)
    ) & (jnp.abs(phik) < del_phi / 2.0)
    # Branch-free redraw: no Python-level `if` on a traced array.  Each pass
    # needs its own key, otherwise every iteration draws the same uniform and
    # the retry does nothing.  Points that still collide after the last pass
    # are rejected below, which is conservative and probability-negligible.
    for redraw_key in jax.random.split(k_redraw, 4):
        redraw = (~narrow) & collide
        cstk = jnp.where(
            redraw, 2.0 * jax.random.uniform(redraw_key, (n,), dtype=dt) - 1.0, cstk
        )
        collide = (
            (~narrow)
            & ((jnp.abs(cstk - cstk1) < cs_range) | (jnp.abs(cstk - cstk2) < cs_range))
            & (jnp.abs(phik) < del_phi / 2.0)
        )

    ok = kv["ok"] & ~((~narrow) & collide)

    cstk = jnp.clip(cstk, -1.0, 1.0)
    tk = jnp.arccos(cstk)
    sntk = jnp.sin(tk)

    # --- photon energy -----------------------------------------------------
    ek_max = 0.5 * (kv["uu"] - kin.wg**2) / (kv["u0"] - kv["pu"] * cstk)
    # aao_rad.f90:747-750 only caps the top end; a non-positive ek_max makes
    # every draw fail ``ek > ekmax`` at line 763, so it has to be a rejection
    # rather than a clipped floor.
    ok = ok & (ek_max > 0.0)
    ek_max_cap = jnp.clip(ek_max, 1e-6, kin.e_beam)

    uek = u(k_ek)
    tail = jnp.exp(-kin.k_exp * ek_max_cap)
    # Truncated-exponential inverse CDF: no trial is discarded for ek > ek_max.
    # The *unclipped* draw is what gets tested against ``ek_max``: the original
    # draws ek = -log(uek)/kexp and sends the trial back to label 20 when
    # ``ek .gt. ekmax`` (aao_rad.f90:757-763).  Clipping first and then
    # comparing would make the test vacuous, silently replacing a rejection with
    # a spike of full-weight points piled up on ek_max -- which is worth about
    # 1% of the cross section and shifts <ek> by 2%.
    ek_drawn = jnp.where(
        kin.truncated_ek > 0.5,
        -jnp.log1p(-uek * (1.0 - tail)) / kin.k_exp,
        -jnp.log(jnp.maximum(uek, 1e-30)) / kin.k_exp,
    )
    ok = ok & (ek_drawn > 0.0) & (ek_drawn <= ek_max)
    # Only round-off can put the truncated draw outside the window; clamp it so
    # the recorded point stays physical, which cannot affect the mask above.
    ek = jnp.minimum(jnp.maximum(ek_drawn, 0.0), ek_max_cap)
    jac_ek = jnp.where(
        kin.truncated_ek > 0.5,
        jnp.exp(kin.k_exp * ek) / (kin.k_exp * jnp.maximum(1.0 - jnp.exp(-kin.k_exp * ek_max_cap), 1e-12)),
        jnp.exp(kin.k_exp * ek) / kin.k_exp,
    )

    csthcm = 2.0 * u(k_cm) - 1.0
    phicm_deg = 360.0 * u(k_cm_phi)
    # flag_ehel = 1 draws the helicity per trial (aao_rad.f90:521, get_spin);
    # the unpolarised case averages over it.
    e_hel = jnp.where(kin.polarized > 0.5, jnp.where(u(k_spin) < 0.5, -1.0, 1.0), 0.0)

    return dict(
        cstk=cstk, tk=tk, sntk=sntk, phik=phik, ek=ek, ek_max=ek_max_cap,
        csthcm=csthcm, phicm_deg=phicm_deg, e_hel=e_hel,
        mcfac=mcfac, mpfac=mpfac, jac_ek=jac_ek, ok=ok,
        # The region-selection uniform, kept so that validation can recover
        # ``intreg`` from aao_rad.f90:657-731 without reconstructing the
        # regions from the sampled point (see validation/compare_trials.py).
        csran=csran, cs_rngb=cs_rngb, cs_range=cs_range,
        cstk1=cstk1, cstk2=cstk2,
    )


class Integrand(NamedTuple):
    """Result of :func:`integrand` for a batch of trial points.

    Attributes
    ----------
    weight:
        Trial weight masked by :attr:`ok_sampling` only, so its maximum is the
        unconditional maximum of the integrand -- the ceiling the original
        scanned for at ``aao_rad.f90:475-513``.  Callers that want the
        cross-section estimate must mask it with :attr:`ok`.
    asym:
        Single-spin beam asymmetry, broadcast to the batch.
    ok:
        Trials that survive *every* selection, including the missing-mass cut.
        These, and only these, contribute to the cross-section estimate.
    ok_sampling:
        Trials for which the integrand itself is defined and positive.
    sigr:
        The cross section alone, before the region, multipole and Jacobian
        factors are applied.  ``weight`` is ``sigr * mcfac * mpfac * jacob``, and
        the three factors after ``sigr`` are pure geometry that the port
        reproduces exactly, so this is what validation has to compare in order
        to attribute a :attr:`weight` mismatch to the physics rather than to the
        sampling.
    """

    weight: Any
    asym: Any
    ok: Any
    ok_sampling: Any
    sigr: Any


def integrand(grid, kin: Kinematics, kv, ph, scheme: str) -> Integrand:
    """Integrand weight and beam asymmetry for every trial point."""
    es, ep, th0, cst0 = kv["es"], kv["ep"], kv["th0"], kv["cst0"]
    ek, cstk, sntk, phik = ph["ek"], ph["cstk"], ph["sntk"], ph["phik"]
    cos_phik = jnp.cos(phik)
    phi_cm = ph["phicm_deg"] * (PI / 180.0)
    e_hel = ph["e_hel"]

    sdotk = (
        es * ek
        - kv["ps"] * ek * cstk * kv["csths"]
        - kv["ps"] * ek * sntk * kv["snths"] * cos_phik
    )
    pdotk = (
        ep * ek
        - kv["pp"] * ek * cstk * kv["csthp"]
        - kv["pp"] * ek * sntk * kv["snthp"] * cos_phik
    )

    common = dict(
        es=es, ep=ep, th0=th0, cst0=cst0, ps=kv["ps"], pp=kv["pp"],
        csthcm=ph["csthcm"], phi=phi_cm, e_hel=e_hel, scheme=scheme,
    )
    sigr_rad, asym_rad = motsa_sigma(
        grid, sdotk=sdotk, pdotk=pdotk, u0=kv["u0"], pu=kv["pu"], uu=kv["uu"],
        cstk=cstk, ek=ek, wg=kin.wg, m_pi=kin.m_pi, **common,
    )
    sigr_nr, asym_nr = non_radiative_sigma(
        grid, q2=kv["q2"], w_sq=kv["w_sq"], delta=kin.delta, m_pi=kin.m_pi, **common,
    )
    low = ek < kin.delta
    sigr = jnp.where(low, sigr_nr, sigr_rad)
    asym = jnp.where(low, asym_nr, asym_rad)

    jacob = ph["jac_ek"] / (2.0 * es * ep) * kv["q2"] ** 2
    weight = sigr * ph["mcfac"] * ph["mpfac"] * jacob
    ok_sampling = ph["ok"] & (sigr > _NEGLIGIBLE) & jnp.isfinite(weight)

    # aao_rad.f90:905-909 builds the hadronic final state with the *pre-exit*
    # electron energy and throws the trial away when the missing mass falls
    # outside the cut -- and that test sits *before* ``sig_tot = sig_tot + sigr``
    # at line 916.  The cut therefore has to gate the trial weight, not just
    # the event record, or the cross-section estimate picks up the rejected
    # non-resonant weight.
    fs = hadronic_final_state(
        kin.e_beam, kv["es"], kv["ep"], kv["th0"], ph["ek"], ph["cstk"], ph["phik"],
        ph["csthcm"], ph["phicm_deg"], kin.m_pi, kin.is_pi0,
    )
    ok = ok_sampling & (fs.w_real > 0.0) & (jnp.abs(fs.mm2 - kin.m_exp) <= kin.mm_cut)
    return Integrand(
        jnp.where(ok_sampling, weight, 0.0), asym, ok, ok_sampling,
        jnp.where(ok_sampling, sigr, 0.0),
    )


def _photon_vector(kin: Kinematics, kv, ep_out, ph):
    """Photon momentum using the pre-exit beam energy and the exit electron."""
    p_s = jnp.sqrt(jnp.maximum(kv["es"] ** 2 - M_E**2, 0.0))
    p_p = jnp.sqrt(jnp.maximum(ep_out**2 - M_E**2, 0.0))
    c_the, s_the = kv["cst0"], kv["snt0"]
    qx = -p_p * s_the
    qz = p_s - p_p * c_the
    q_vec = jnp.sqrt(qx**2 + qz**2)
    cstk, sntk = ph["cstk"], ph["sntk"]
    safe_q = jnp.where(q_vec > 1e-9, q_vec, 1.0)
    cstq = jnp.where(q_vec > 1e-9, qz / safe_q, 0.0)
    sntq = jnp.sqrt(jnp.maximum(1.0 - cstq**2, 0.0))
    csphk, snphk = jnp.cos(ph["phik"]), jnp.sin(ph["phik"])
    e = ph["ek"]
    return jnp.stack(
        [
            e * (sntk * csphk * cstq - cstk * sntq),
            e * sntk * snphk,
            e * (cstk * cstq + sntk * csphk * sntq),
        ],
        axis=-1,
    )


def finalize(kin: Kinematics, kv, ph, asym, key, n: int):
    """Exit radiation loss, missing-mass test, and the event record.

    Returns ``(record, ok)`` where ``record`` has shape ``(n, 32)`` in the
    :data:`EVENT_COLUMNS` order.
    """
    k_x, k_test, k_fall = jax.random.split(key, 3)
    dt = kin.e_beam.dtype

    # Path the scattered electron still has to travel inside the target (cm).
    c_the, s_the = kv["cst0"], kv["snt0"]
    safe_s = jnp.where(s_the > 1e-6, s_the, 1e-6)
    safe_c = jnp.where(jnp.abs(c_the) > 1e-6, c_the, jnp.sign(c_the) * 1e-6 + 1e-12)
    t_geom = kin.r_targ / safe_s
    t_remain = (kin.t_len_cm - kv["targs_cm"]) / safe_c
    targp_cm = jnp.where(t_geom > 0.0, jnp.minimum(t_geom, t_remain), t_remain)
    targp_cm = jnp.where(jnp.isfinite(targp_cm) & (targp_cm > 1e-6), targp_cm, 1e-6)
    # The original converts the exit path back to radiation lengths before
    # forming the bremsstrahlung exponent.
    targp = targp_cm * (BFAC / HYDROGEN_RAD)

    # Same exact construction as the in-target loss above: the original loops
    # (label 222) rejecting until the electron survives the exit, which at fixed
    # targp accepts with probability 1/(1+targp).  Inverting that gives
    # eloss = (1-u2)*u3**(1/targp), with no retry loop and no tracer bool.
    eloss = (1.0 - jax.random.uniform(k_test, (n,), dtype=dt)) * jnp.power(
        jax.random.uniform(k_x, (n,), dtype=dt), 1.0 / jnp.maximum(targp, 1e-12)
    )
    eloss = jnp.where(eloss <= 1.0, eloss, 1.0)

    ep_out = kv["ep"] * (1.0 - eloss)
    ok = ep_out >= kin.ep_min

    fs = hadronic_final_state(
        kin.e_beam, kv["es"], ep_out, kv["th0"], ph["ek"], ph["cstk"], ph["phik"],
        ph["csthcm"], ph["phicm_deg"], kin.m_pi, kin.is_pi0,
    )

    # Reconstructed W from the measured electron, using exactly the same
    # convention as the missing-mass reconstruction.
    nu = kin.e_beam - ep_out
    q2_meas = 2.0 * kin.e_beam * ep_out * (1.0 - c_the) - 2.0 * M_E**2
    w_meas = jnp.sqrt(jnp.maximum(M_N**2 + 2.0 * M_N * nu - q2_meas, 0.0))

    # aao_rad.f90:984 rejects only ``mm2 == 0`` (below pion threshold) and
    # ``ep < ep_min`` here; the missing-mass *cut* was already applied, with the
    # pre-exit energy, before the trial was accepted at all.  Re-applying it
    # would drop a second, uncounted set of events.
    ok = ok & (fs.w_real > 0.0)

    k_vec = _photon_vector(kin, kv, ep_out, ph)
    # Column order must match EVENT_COLUMNS.
    cols = [
        kv["es"], ep_out, kv["theta"], w_meas, fs.w_real,
        fs.p_p[..., 0], fs.p_p[..., 1], fs.p_p[..., 2], fs.e_prot,
        fs.p_pi[..., 0], fs.p_pi[..., 1], fs.p_pi[..., 2], fs.e_pi,
        ph["csthcm"], ph["phicm_deg"], fs.mm2,
        ph["ek"], ph["cstk"], ph["phik"] * (180.0 / PI),
        kv["qvecx"], kv["qvecz"], kv["q0"], c_the,
        k_vec[..., 0], k_vec[..., 1], k_vec[..., 2],
        kin.vertex_x, kin.vertex_y, kv["vertex_z"],
        kv["q2"], ph["e_hel"], asym,
    ]
    # vertex_x/vertex_y are per-run scalars; broadcast so every column has the
    # batch length.
    rec = jnp.stack(
        [jnp.broadcast_to(c, (n,)) if jnp.ndim(c) == 0 else c for c in cols], axis=1
    )
    return rec, ok


# ---------------------------------------------------------------------------
# Public driver
# ---------------------------------------------------------------------------
@dataclass
class GenerationStats:
    """Run diagnostics, reported at the end of a generation."""

    n_events: int = 0
    n_trials: int = 0
    weight_max: float = 0.0
    """The acceptance ceiling, estimated on the device before the run."""
    weight_max_observed: float = 0.0
    """Largest legal trial weight actually seen during the run.

    If this exceeds :attr:`weight_max` the ceiling was too low: those trials
    could never be accepted.  :attr:`n_above_ceiling` counts them.
    """
    mean_weight: float = 0.0
    weight_sum: float = 0.0
    n_sampling_accepted: int = 0
    """Trials the importance sampler accepted, before any kinematic cut."""
    n_above_ceiling: int = 0
    """Legal trials whose weight exceeded :attr:`weight_max`.

    Such a trial is accepted with probability 1 rather than ``w / weight_max``,
    so it is over-represented in the event stream.  :attr:`ceiling_bias` is the
    size of that effect on the sampled cross section.
    """
    ceiling_bias: float = 0.0
    """Fraction of the cross section carried by :attr:`n_above_ceiling`.

    :attr:`sigma_mc` sums the trial weight over *every* trial, so it is
    unaffected by the ceiling; only the event-by-event distributions are.  This
    bounds that error, and it should be small.
    """
    acceptance: float = 0.0
    """Importance-sampler efficiency, ``n_sampling_accepted / n_trials``."""
    event_yield: float = 0.0
    """Recorded events per trial, i.e. ``n_events / n_trials``.

    Lower than :attr:`acceptance` by whatever the kinematic cuts reject.
    """
    sigma_mc: float = 0.0
    """Integrated cross section from the mean trial weight [microbarn]."""
    sigma_accepted: float = 0.0
    """The same, from the *sampling* acceptance against ``weight_max``.

    ``E[weight | accepted] = weight_max * P(u < weight / weight_max) =
    mean(weight)``, so summing ``accepted * weight_max`` over the trials gives an
    independent unbiased estimator of :attr:`sigma_mc` that shares none of its
    code path -- which makes the agreement between the two a useful self-check.
    Accumulated per batch, so it stays exact when the ceiling is raised
    mid-run.
    """
    phase_space: float = 0.0
    """(4 pi)^2 * 2 pi * d(1/Q^2) * dE': converts mean weight to a cross section."""
    seconds: float = 0.0
    events_per_second: float = 0.0
    trials_per_second: float = 0.0

    def summary(self) -> str:
        clipped = (
            f" / {self.weight_max_observed:.6g} seen"
            if self.weight_max_observed > self.weight_max
            else ""
        )
        return (
            f"events           : {self.n_events:,}\n"
            f"trials           : {self.n_trials:,}\n"
            f"sampling accept. : {100.0 * self.acceptance:.3f}%\n"
            f"event yield      : {100.0 * self.event_yield:.3f}%\n"
            f"weight max / mean: {self.weight_max:.6g}{clipped} / "
            f"{self.mean_weight:.6g}\n"
            f"above ceiling    : {self.n_above_ceiling:,} trials, "
            f"{self.ceiling_bias:.2%} of the cross section\n"
            f"sigma (MC)       : {self.sigma_mc:.6g} micro-barn\n"
            f"sigma (accepted) : {self.sigma_accepted:.6g} micro-barn\n"
            f"throughput       : {self.events_per_second:,.0f} events/s "
            f"({self.trials_per_second:,.0f} trials/s)\n"
            f"wall time        : {self.seconds:.1f} s"
        )


def _make_sampler(grid, kin: Kinematics, scheme: str, n: int):
    """Build the jitted kernel that evaluates one batch of trial weights.

    ``n`` is baked in as a Python int, so a given batch size compiles once and
    is then reused for every step of the production loop.
    """

    def sample(key):
        kv = draw_kinematics(key, kin, n)
        ph = draw_photon(jax.random.fold_in(key, 1), kin, kv, n)
        weight = integrand(grid, kin, kv, ph, scheme)[0]
        return weight

    return jax.jit(sample)


def _make_stepper(grid, kin: Kinematics, scheme: str, n: int, capacity: int):
    """Build the jitted kernel that accepts, compacts and buffers a batch.

    ``capacity`` is the number of rows the device buffer can hold; accepted
    events beyond it are dropped, which is what bounds host and device memory
    during a streaming run.
    """

    def step(key, buffer, write_pos, weight_max):
        kv = draw_kinematics(key, kin, n)
        ph = draw_photon(jax.random.fold_in(key, 1), kin, kv, n)
        integ = integrand(grid, kin, kv, ph, scheme)
        asym, ok = integ.asym, integ.ok
        # The cross-section estimator sums the trial weight over *every* trial
        # that survives the missing-mass cut (aao_rad.f90:916), zero elsewhere.
        weight = jnp.where(ok, integ.weight, 0.0)
        raw_weight = integ.weight

        u_acc = jax.random.uniform(jax.random.fold_in(key, 2), (n,), dtype=weight.dtype)
        # The sampling acceptance is the pure importance-sampler efficiency: it
        # must be counted *before* the kinematic cuts are applied downstream, or
        # the cross-section estimate built from it inherits their rejection.
        sampled = ok & (u_acc * weight_max < weight)
        # A weight above the ceiling can never be accepted (weight/weight_max
        # > 1), so the trial is simply lost.  Counting these says whether the
        # ceiling estimated on the device was good enough.
        above = ok & (raw_weight > weight_max)

        rec, final_ok = finalize(kin, kv, ph, asym, jax.random.fold_in(key, 3), n)
        accept = sampled & final_ok

        # Compaction: a stable argsort puts accepted points first, preserving
        # their relative order, so the buffer stays dense without a scan.
        order = jnp.argsort(~accept, stable=True)
        n_acc = jnp.minimum(jnp.sum(accept), capacity - write_pos)
        rows = jnp.arange(n)
        idx = write_pos + rows
        valid = (rows < n_acc)[:, None]
        payload = jnp.take(rec, order, axis=0)
        buffer = buffer.at[idx].set(
            jnp.where(valid, payload, buffer[jnp.clip(idx, 0, capacity - 1)]),
            mode="drop",
        )
        counters = jnp.array(
            [
                jnp.sum(sampled),
                jnp.sum(above),
                # Maximum over the trials that actually reach the event stream.
                # A trial the missing-mass cut rejects contributes nothing to the
                # accepted events, so its weight must not push the ceiling up.
                jnp.max(jnp.where(ok, raw_weight, 0.0)),
                # Weight, not count: this is the share of the cross section the
                # ceiling cannot sample correctly.  Accepting u * W < w is
                # exactly min(1, w/W), so a trial below the ceiling is sampled
                # as w while one above it is over-represented by W/w; the total
                # weight above the ceiling is therefore the size of the error.
                jnp.sum(jnp.where(above, jnp.where(ok, weight, 0.0), 0.0)),
            ],
            dtype=jnp.float32,
        )
        return buffer, write_pos + n_acc, weight, counters

    return jax.jit(step)


class EventGenerator:
    """Vectorised radiative pion electroproduction generator.

    Examples
    --------
    >>> grid = build_grid(3, parms_dir="parms")            # doctest: +SKIP
    >>> events, stats = EventGenerator(grid).generate(cfg)  # doctest: +SKIP
    """

    def __init__(self, grid: InterpolationGrid) -> None:
        self.grid = grid

    def estimate_weight_max(
        self, cfg: GeneratorConfig, kin: Kinematics | None = None, n: int | None = None
    ) -> tuple[float, float]:
        """Return ``(weight_max, mean_weight)`` from a large batch of trials.

        The maximum is what the rejection sampler needs as its ceiling; the
        original read it from ``stdin`` after a serial 10,000-point scan.
        """
        kin = kin if kin is not None else build_kinematics(cfg)
        n = n or max(int(cfg.batch_size), 1 << 16)
        sample = _make_sampler(self.grid, kin, cfg.interp_scheme, n)
        weights = sample(jax.random.PRNGKey(cfg.seed if cfg.seed is not None else 0))
        return float(jnp.max(weights)), float(jnp.mean(weights))

    def stream(
        self,
        cfg: GeneratorConfig,
        *,
        chunk_events: int = 500_000,
        n_events: int | None = None,
    ) -> Iterator[tuple[np.ndarray, GenerationStats]]:
        """Yield blocks of accepted events, and the final stats on the last block.

        The device buffer only ever holds ``chunk_events`` rows, so a billion
        event run uses the same memory as a million event one.  Blocks come out
        in *production* order, not sorted, because each block is filled in the
        order its accepted trials appeared.
        """
        t_start = time.perf_counter()
        kin = build_kinematics(cfg)
        total = int(n_events if n_events is not None else cfg.n_events)
        if total <= 0:
            raise ValueError("n_events must be positive")

        seed = cfg.seed if cfg.seed is not None else _random_seed()
        key = jax.random.PRNGKey(seed)
        batch = int(cfg.batch_size)
        chunk = max(min(int(chunk_events), total), batch)

        weight_max, _ = self.estimate_weight_max(cfg, kin, max(batch, 1 << 16))
        if not np.isfinite(weight_max) or weight_max <= 0.0:
            raise RuntimeError(
                "could not find a positive integrand maximum -- check that the MAID "
                "table matches the requested channel and that the kinematic cuts "
                "are not empty"
            )
        # A pilot batch only bounds the maximum of the points it actually
        # contains.  ``accept iff u * W < w`` is exactly ``min(1, w/W)``, so as
        # long as W exceeds the largest weight in the run the event stream is
        # distributed as w; a weight that slips over the ceiling instead gets
        # accepted with probability 1, which biases the sample towards the
        # heavy tail.  Scaling the pilot estimate gives room for the tail
        # beyond it, and the running maximum below closes the gap for good.
        weight_max *= cfg.weight_max_margin
        log.info("weight maximum estimated as %.6g (seed %d)", weight_max, seed)

        step = _make_stepper(self.grid, kin, cfg.interp_scheme, batch, chunk)
        buffer = jnp.zeros((chunk, len(EVENT_COLUMNS)), dtype=jnp.float32)
        write_pos = jnp.zeros((), dtype=jnp.int32)

        n_trials = 0
        weight_sum = 0.0
        n_sampled = 0
        n_above = 0
        weight_above = 0.0
        # sum(acceptances) * ceiling, accumulated per batch with the ceiling
        # that batch actually used, so it stays exact when the ceiling moves.
        sampled_weight = 0.0
        weight_max_seen = 0.0
        n_done = 0
        next_report = cfg.verbose_every or 0
        stats: GenerationStats | None = None

        while n_done < total:
            key, subkey = jax.random.split(key)
            buffer, write_pos, weights, counters = step(
                subkey, buffer, write_pos, weight_max
            )
            weight_sum += float(jnp.sum(weights))
            n_trials += batch
            n_sampled += int(counters[0])
            n_above += int(counters[1])
            weight_max_seen = max(weight_max_seen, float(counters[2]))
            weight_above += float(counters[3])
            sampled_weight += int(counters[0]) * weight_max
            if weight_max_seen > weight_max:
                # The pilot estimate was short.  Raising the ceiling fixes every
                # later batch; the batch that found the outlier is already
                # over-weighted, which is what ``ceiling_bias`` measures.
                log.warning(
                    "weight ceiling %.6g was exceeded (%.6g seen); raising it. "
                    "Increase --weight-max-margin to avoid this.",
                    weight_max, weight_max_seen,
                )
                weight_max = weight_max_seen
            produced = int(write_pos)
            n_done += produced

            if cfg.verbose_every and n_done >= next_report:
                log.info(
                    "%d/%d events after %d trials (%.3f%% sampled, "
                    "%.3f%% recorded)",
                    min(n_done, total), total, n_trials,
                    100.0 * n_sampled / max(n_trials, 1),
                    100.0 * n_done / max(n_trials, 1),
                )
                next_report += cfg.verbose_every

            # The original ran at a few tenths of a percent acceptance, so the
            # budget has to be generous; this guard is for a genuinely empty
            # or near-empty phase space, not for the normal case.  Real
            # configurations reach 231,000 trials per event, so this is a
            # config field rather than a constant -- see max_trials_per_event.
            if n_trials > cfg.max_trials_per_event * max(total, 1):
                raise RuntimeError(
                    f"acceptance collapsed: {n_done}/{total} events from {n_trials} "
                    f"trials ({n_trials / max(n_done, 1):,.0f} trials per event). "
                    "Either the kinematic cuts are empty for this beam energy, or "
                    "this configuration needs a larger --max-trials-per-event."
                )

            # Trim the last block so we never emit more events than requested.
            keep = min(produced, total - (n_done - produced))
            # Copy out only the rows that were actually written, and take an
            # owning copy of them.  ``buffer`` is (chunk_events, 32) float32 --
            # 64 MB at the default chunk -- and a run needs thousands of
            # iterations to reach 20k events at the original's sub-percent
            # acceptance.  Transferring the whole buffer and then slicing it
            # hands the caller a *view*, so every yielded block pinned all
            # 64 MB for as long as the caller kept it, and a caller that
            # accumulates blocks (the CLI does, to write one npz) grew by
            # ~192 GB per run.  That was the 103 OOM kills on the 64-worker
            # arm and the ~100 GB per process on the GPU arm.  Slicing before
            # device_get also keeps the transfer proportional to the events
            # produced rather than to the buffer.
            block = np.array(jax.device_get(buffer[:keep]), copy=True)

            is_last = n_done >= total
            if is_last:
                stats = self._stats(
                    cfg, kin, weight_max, weight_sum, n_trials, n_done, t_start,
                    n_sampled, n_above, weight_max_seen, weight_above, sampled_weight,
                )
            yield block.view(EVENT_DTYPE).reshape(-1), stats

            if is_last:
                return
            # Reset the buffer for the next block; the tail is re-filled.
            write_pos = jnp.zeros((), dtype=jnp.int32)

    def generate(
        self,
        cfg: GeneratorConfig,
        *,
        chunk_events: int | None = None,
    ) -> tuple[np.ndarray, GenerationStats]:
        """Generate ``cfg.n_events`` accepted events and return them all at once.

        For large samples prefer :meth:`stream` and write each block as it
        arrives, so host memory stays bounded.

        Returns
        -------
        (events, stats)
            ``events`` is a structured array with the fields in
            :data:`EVENT_COLUMNS`; ``stats`` carries the run diagnostics.
        """
        blocks: list[np.ndarray] = []
        stats = None
        for block, blk_stats in self.stream(cfg, chunk_events=chunk_events or 10**9):
            blocks.append(block)
            stats = blk_stats or stats
        assert stats is not None  # stream always yields at least one block
        events = (
            np.concatenate(blocks) if len(blocks) > 1 else blocks[0]
        )
        return events, stats

    @staticmethod
    def _stats(
        cfg: GeneratorConfig,
        kin: Kinematics,
        weight_max: float,
        weight_sum: float,
        n_trials: int,
        n_events: int,
        t_start: float,
        n_sampled: int = 0,
        n_above: int = 0,
        weight_max_seen: float = 0.0,
        weight_above: float = 0.0,
        sampled_weight: float = 0.0,
    ) -> GenerationStats:
        seconds = time.perf_counter() - t_start
        phase_space = (4.0 * PI) ** 2 * 2.0 * PI * float(kin.uq2_range) * float(kin.ep_range)
        mean_w = weight_sum / max(n_trials, 1)
        acceptance = n_sampled / n_trials if n_trials else 0.0
        return GenerationStats(
            n_events=n_events,
            n_trials=n_trials,
            weight_max=float(weight_max),
            weight_max_observed=float(weight_max_seen),
            mean_weight=mean_w,
            weight_sum=weight_sum,
            n_sampling_accepted=int(n_sampled),
            n_above_ceiling=int(n_above),
            ceiling_bias=(weight_above / weight_sum) if weight_sum > 0.0 else 0.0,
            acceptance=acceptance,
            event_yield=n_events / n_trials if n_trials else 0.0,
            sigma_mc=mean_w * phase_space,
            sigma_accepted=(sampled_weight / n_trials * phase_space) if n_trials else 0.0,
            phase_space=phase_space,
            seconds=seconds,
            events_per_second=n_events / seconds if seconds > 0 else 0.0,
            trials_per_second=n_trials / seconds if seconds > 0 else 0.0,
        )


def _random_seed() -> int:
    return int.from_bytes(os.urandom(4), "little")
