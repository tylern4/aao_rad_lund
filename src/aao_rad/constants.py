"""Physical constants.

The Fortran original carried several *different* sets of constants in
different subroutines (``spp.inc`` used ``m_p = .9382799`` while
``missm`` in ``aao_rad.f90`` hard-coded ``mp = .938`` and ``MEL = 0.511e-3``,
and the fine-structure constant was ``1/137`` in the generator but
``0.729735e-02`` in the response-function code).  Those inconsistencies are a
few times ``1e-4`` in relative terms, i.e. far below the accuracy of the
underlying phenomenological fit, so this port uses a single consistent set of
PDG values throughout.  See ``README.md`` for the full list of differences.
"""

from __future__ import annotations

import jax.numpy as jnp

__all__ = [
    "ALPHA",
    "BFAC",
    "HYDROGEN_RAD",
    "M_E",
    "M_N",
    "M_PI0",
    "M_PIP",
    "PI",
    "TABLE_Q2_MAX",
    "TABLE_W_MAX",
    "TABLE_W_MIN",
    "pi_masses",
]

PI = 3.141592653589793
"""Mathematical pi (the Fortran code used the truncated ``3.1415926``)."""

#: Nucleon (proton) mass [GeV]
M_N = 0.9382720813
#: Charged pion mass [GeV]
M_PIP = 0.13957039
#: Neutral pion mass [GeV]
M_PI0 = 0.1349766
#: Electron mass [GeV]
M_E = 0.51099895000e-3
#: Fine-structure constant (CODATA 2018).  The Fortran used 1/137 and 0.729735e-2.
ALPHA = 7.2973525693e-3

#: Target lengthening factor used to convert a length in cm into radiation
#: lengths, ``t_rl = bfac * t_cm / hydrogen_rad`` (``aao_rad.f90``).
BFAC = 4.0 / 3.0
#: Hydrogen radiation length [cm], exactly as the original hard-coded it.
HYDROGEN_RAD = 865.0

#: Upper edges of the shipped MAID07 table grid (``mpintp.inc``: NVAR1=101
#: from 0 to 5 GeV^2, NVAR2=93 from 1.08 to 2 GeV).  Lookups outside these
#: bounds have no tabulated value; see :attr:`GeneratorConfig.w_max`.
TABLE_Q2_MAX = 5.0
TABLE_W_MAX = 2.0
#: Lower edge of the W grid.  ``multipole_amps.f90`` raises it to 1.1 anyway.
TABLE_W_MIN = 1.1

#: jax-typed aliases for use inside traced functions.
PI_J = jnp.asarray(PI)


def pi_masses(channel: int) -> tuple[float, float]:
    """Return ``(m_pi, m_expected)`` for a production channel.

    ``channel`` follows the original ``epirea`` convention:

    ======  =============================
    value   reaction
    ======  =============================
    ``1``   ``p -> pi0 p``   (m_exp = m_pi0^2)
    ``2``   ``n -> pi- p``
    ``3``   ``p -> pi+ n``   (m_exp = m_p^2)
    ``5``   ``n -> pi0 n``
    ======  =============================

    Returns
    -------
    (m_pi, m_exp_sq)
        Pion mass in GeV and the squared expected *missing* mass in GeV^2
        that the generator aims at (the missing mass cut is applied against
        this value).
    """
    if channel == 1:
        return M_PI0, M_PI0**2
    if channel == 2:
        return M_PIP, M_N**2
    if channel == 3:
        return M_PIP, M_N**2
    if channel == 5:
        return M_PI0, M_N**2
    raise ValueError(
        f"unsupported pion production channel {channel!r}; expected one of 1, 2, 3, 5"
    )
