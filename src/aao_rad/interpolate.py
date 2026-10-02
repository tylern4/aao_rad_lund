"""Batched grid interpolation for the MAID multipole tables.

Everything here is written so that a single call evaluates an *arbitrary
number* of ``(Q2, W)`` points at once, which is what makes the GPU port
possible: the table is uploaded once and the interpolation becomes four
gathers plus a weighted sum.

Two schemes are provided:

``linear``
    Multilinear (bilinear) interpolation.  This is what the Fortran code
    actually used (``read_sf_file`` hard-sets ``method_spline = 2``) and is the
    default here so that results match the original.

``spline``
    Tensor-product natural cubic spline, reproducing Numerical Recipes'
    ``SPLIN2``.  The second-derivative tables along the *W* axis are built once
    at load time; the ones along *Q2* are, unlike in the original, cached too
    instead of being rebuilt on every query.
"""

from __future__ import annotations

from typing import NamedTuple

import jax.numpy as jnp
import numpy as np

__all__ = [
    "InterpolationGrid",
    "bilinear",
    "natural_spline_2nd_derivs",
    "spline_tensor",
]

_NATURAL = 1.0e30  # sentinel: >= this means "use natural boundary condition"


def _bracket(
    axis: jnp.ndarray, x: jnp.ndarray, *, uniform: bool
) -> tuple[jnp.ndarray, jnp.ndarray]:
    """Return ``(i, t)`` with ``x`` between ``axis[i]`` and ``axis[i+1]``.

    ``t`` is the linear position inside that interval.  Out-of-range values are
    clamped onto the grid instead of aborting the run (the original Fortran
    called ``STOP``; this generator already clamps ``Q2`` and ``W`` to the
    table bounds, so clamping here is a safety net rather than a behaviour
    change).
    """
    n = axis.size
    xc = jnp.clip(x, axis[0], axis[-1])
    if uniform:
        # Closed form; avoids a binary search and is a single fused expression.
        dx = axis[1] - axis[0]
        i = jnp.clip(jnp.floor((xc - axis[0]) / dx).astype(jnp.int32), 0, n - 2)
        t = (xc - axis[i]) / (axis[i + 1] - axis[i])
    else:
        i = jnp.clip(jnp.searchsorted(axis, xc, side="right") - 1, 0, n - 2)
        t = (xc - axis[i]) / (axis[i + 1] - axis[i])
    return i, t


def bilinear(
    values: jnp.ndarray, axis1: jnp.ndarray, axis2: jnp.ndarray, x1: jnp.ndarray,
    x2: jnp.ndarray, *, uniform1: bool = True, uniform2: bool = True
) -> jnp.ndarray:
    """Bilinear interpolation of ``values[k, i, j]`` at points ``(x1, x2)``.

    Parameters
    ----------
    values
        Array of shape ``(k, n1, n2)``.
    x1, x2
        Query coordinates, any broadcastable shape (e.g. ``(batch,)``).

    Returns
    -------
    jnp.ndarray
        Shape ``broadcast(x1, x2).shape + (k,)``.
    """
    x1, x2 = jnp.broadcast_arrays(jnp.asarray(x1), jnp.asarray(x2))
    shape = x1.shape
    flat1, flat2 = x1.ravel(), x2.ravel()

    i, t = _bracket(axis1, flat1, uniform=uniform1)
    j, u = _bracket(axis2, flat2, uniform=uniform2)

    w00 = (1.0 - t) * (1.0 - u)
    w01 = (1.0 - t) * u
    w10 = t * (1.0 - u)
    w11 = t * u

    # Gather as (k, n): the amplitude axis is the *leading* one in `values` and
    # the weights run along the flattened batch.  Transpose back at the end.
    out = (
        w00[None, :] * values[:, i, j]
        + w01[None, :] * values[:, i, j + 1]
        + w10[None, :] * values[:, i + 1, j]
        + w11[None, :] * values[:, i + 1, j + 1]
    )
    return jnp.moveaxis(out, 0, -1).reshape(*shape, values.shape[0])


def natural_spline_2nd_derivs(
    y: np.ndarray, x: np.ndarray, yp1: float = _NATURAL, ypn: float = _NATURAL
) -> np.ndarray:
    """Second derivatives of the natural cubic spline through ``(x, y)``.

    Vectorised transcription of Numerical Recipes' ``SPLINE``: the *same*
    tridiagonal system applies to every column of ``y``, so all leading
    coefficients are computed at once and only the back-substitution loop is
    run in Python.
    """
    n = x.size
    if n < 4:
        raise ValueError("natural cubic spline needs at least 4 points")

    h = np.diff(x)  # (n-1,)
    y2 = np.zeros_like(y)
    u = np.zeros_like(y)

    natural_lo = yp1 >= 0.99e30
    natural_hi = ypn >= 0.99e30

    if natural_lo:
        y2[0] = 0.0
        u[0] = 0.0
    else:
        y2[0] = -0.5
        u[0] = (3.0 / h[0]) * ((y[1] - y[0]) / h[0] - yp1)

    for i in range(1, n - 1):  # tridiagonal decomposition
        sig = h[i - 1] / (h[i - 1] + h[i])
        p = sig * y2[i - 1] + 2.0
        y2[i] = (sig - 1.0) / p
        u[i] = (
            6.0
            * ((y[i + 1] - y[i]) / h[i] - (y[i] - y[i - 1]) / h[i - 1])
            / (h[i - 1] + h[i])
            - sig * u[i - 1]
        ) / p

    if natural_hi:
        qn = 0.0
        un = 0.0
    else:
        qn = 0.5
        un = (3.0 / h[-1]) * (ypn - (y[-1] - y[-2]) / h[-1])
    y2[-1] = (un - qn * u[-2]) / (qn * y2[-2] + 1.0)

    for j in range(n - 2, -1, -1):  # back substitution
        y2[j] = y2[j] * y2[j + 1] + u[j]
    return y2


def _splint(xa: jnp.ndarray, ya: jnp.ndarray, y2a: jnp.ndarray, x: jnp.ndarray) -> jnp.ndarray:
    """Evaluate a pre-conditioned cubic spline at a batch of points."""
    klo, u = _bracket(xa, x, uniform=True)
    khi = klo + 1
    h = xa[khi] - xa[klo]
    a = (xa[khi] - x) / h
    b = (x - xa[klo]) / h
    return (
        a * ya[klo]
        + b * ya[khi]
        + ((a * a * a - a) * y2a[klo] + (b * b * b - b) * y2a[khi]) * (h * h) / 6.0
    )


def _spline_coeffs(
    xa: jnp.ndarray, x: jnp.ndarray, lo: jnp.ndarray, hi: jnp.ndarray
) -> tuple[jnp.ndarray, jnp.ndarray, jnp.ndarray]:
    """Weights of the cubic-spline basis on ``[lo, hi]`` for each query point.

    Returns ``(a, b, h2/6)`` where the interpolated value is
    ``a*y[lo] + b*y[hi] + (a^3-a)*y2[lo]*h2/6 + (b^3-b)*y2[hi]*h2/6``.
    """
    h = xa[hi] - xa[lo]
    a = (xa[hi] - x) / h
    b = (x - xa[lo]) / h
    return a, b, (h * h) / 6.0


def spline_tensor(
    values: jnp.ndarray,
    axis1: jnp.ndarray,
    axis2: jnp.ndarray,
    d2_1: jnp.ndarray,
    d2_2: jnp.ndarray,
    x1: jnp.ndarray,
    x2: jnp.ndarray,
) -> jnp.ndarray:
    """Tensor-product natural cubic spline evaluation.

    ``values``/``d2_1``/``d2_2`` all have shape ``(k, n1, n2)``; ``d2_2`` holds
    second derivatives along axis 2 (built once at load time) and ``d2_1``
    those along axis 1.

    Both passes are plain fancy-index gathers rather than ``vmap``: the bracket
    index differs per query point, so a vmapped ``_splint`` would have to
    diagonal-index anyway, and doing it directly keeps the shapes obvious.
    """
    x1, x2 = jnp.broadcast_arrays(jnp.asarray(x1), jnp.asarray(x2))
    shape = x1.shape
    flat1, flat2 = x1.ravel(), x2.ravel()
    k = values.shape[0]

    # Pass 1: spline along axis 2, for every (k, i) at once.  Indexing the
    # last axis with a (nbatch,) index array yields (k, n1, nbatch).
    j, _ = _bracket(axis2, flat2, uniform=True)
    j = j.astype(jnp.int32)
    a2, b2, s2 = _spline_coeffs(axis2, flat2, j, j + 1)
    y2 = (
        a2 * values[:, :, j]
        + b2 * values[:, :, j + 1]
        + ((a2**3 - a2) * d2_2[:, :, j] + (b2**3 - b2) * d2_2[:, :, j + 1]) * s2
    )  # (k, n1, nbatch)

    # Pass 2: spline along axis 1.  d2_1 depends on (i, j), so the axis-2
    # bracket j is reused rather than re-derived.  y2 is (k, n1, nbatch) and
    # the gather has to be diagonal -- row i[n] of column n -- which needs an
    # explicit arange alongside i; d2_1 needs no such trick because its two
    # advanced indices are adjacent.
    i, _ = _bracket(axis1, flat1, uniform=True)
    i = i.astype(jnp.int32)
    idx = jnp.arange(flat1.size)
    a1, b1, s1 = _spline_coeffs(axis1, flat1, i, i + 1)
    out = (
        a1 * y2[:, i, idx]
        + b1 * y2[:, i + 1, idx]
        + ((a1**3 - a1) * d2_1[:, i, j] + (b1**3 - b1) * d2_1[:, i + 1, j]) * s1
    )  # (k, nbatch)
    return jnp.moveaxis(out, 0, -1).reshape(*shape, k)


class InterpolationGrid(NamedTuple):
    """Immutable, jittable bundle of a table plus everything needed to read it."""

    values: jnp.ndarray
    axis1: jnp.ndarray
    axis2: jnp.ndarray
    uniform1: bool
    uniform2: bool
    d2_1: jnp.ndarray | None
    d2_2: jnp.ndarray | None

    def __call__(
        self, x1: jnp.ndarray, x2: jnp.ndarray, *, scheme: str = "linear"
    ) -> jnp.ndarray:
        """Interpolate at ``(x1, x2)``; result has shape ``x1.shape + (k,)``."""
        if scheme == "linear" or self.d2_1 is None or self.d2_2 is None:
            return bilinear(
                self.values,
                self.axis1,
                self.axis2,
                x1,
                x2,
                uniform1=self.uniform1,
                uniform2=self.uniform2,
            )
        if scheme == "spline":
            return spline_tensor(self.values, self.axis1, self.axis2, self.d2_1, self.d2_2, x1, x2)
        raise ValueError(f"unknown interpolation scheme {scheme!r}")
