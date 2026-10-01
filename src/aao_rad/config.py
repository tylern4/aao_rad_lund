"""Run configuration for the radiative event generator.

Everything the Fortran driver read from ``stdin`` is a field here.  The
original prompt sequence is preserved by :meth:`GeneratorConfig.from_legacy_input`
so old ``.inp`` files can still be used, but the normal entry points are
keyword arguments on the CLI, a JSON/YAML file, or a Python API call.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field, fields, replace
from pathlib import Path
from typing import Any, Literal

from .constants import M_N, TABLE_W_MAX, TABLE_W_MIN

__all__ = ["EXPERIMENTS", "GeneratorConfig", "ep_window_for_w"]


def ep_window_for_w(
    beam_energy: float,
    q2_min: float,
    q2_max: float,
    w_min: float = TABLE_W_MIN,
    w_max: float = TABLE_W_MAX,
) -> tuple[float, float]:
    """Scattered-electron window that keeps ``W`` inside the tabulated range.

    The generator samples ``Q^2`` and ``E'`` independently inside their windows,
    so ``W^2 = M_N^2 + 2 M_N (E_beam - E') - Q^2`` varies over the whole
    rectangle.  ``W`` falls as ``E'`` rises and as ``Q^2`` rises, so the corners
    ``(Q^2_min, E'_max)`` and ``(Q^2_max, E'_min)`` bracket the extremes and the
    rectangle maps to a diagonal band in ``W``.  Inverting at those two corners
    gives the window that keeps the whole rectangle inside ``[w_min, w_max]``::

        E'_min = E_beam - (w_max^2 - M_N^2 + Q^2_min) / (2 M_N)
        E'_max = E_beam - (w_min^2 - M_N^2 + Q^2_max) / (2 M_N)

    Ignoring the electron mass (as the generator does) shifts these by
    ``m_e^2 / (2 M_N) ~ 1.4e-7`` GeV, far below the sampling resolution, so
    the two agree to well within float precision.

    Parameters
    ----------
    beam_energy, q2_min, q2_max
        Beam energy and momentum-transfer window [GeV, GeV^2].
    w_min, w_max
        Hadronic invariant-mass band to stay within [GeV].  Defaults to the
        MAID07 table's own edges.

    Returns
    -------
    (ep_min, ep_max)
        The usable window.  It may come back empty or inverted when the
        ``Q^2`` window is too wide for ``W`` to stay in range; that is reported
        by raising ``ValueError`` rather than silently returning nonsense.
    """
    nu_hi = (w_max**2 - M_N**2 + q2_min) / (2.0 * M_N)
    nu_lo = (w_min**2 - M_N**2 + q2_max) / (2.0 * M_N)
    ep_min = beam_energy - nu_hi
    ep_max = beam_energy - nu_lo
    if not 0.0 < ep_min < ep_max < beam_energy:
        raise ValueError(
            f"no scattered-electron window keeps W in [{w_min:g}, {w_max:g}] GeV "
            f"for a {beam_energy:g} GeV beam and Q^2 in [{q2_min:g}, {q2_max:g}] "
            f"GeV^2 (the rectangle maps to [{ep_max:.3f}, {ep_min:.3f}] GeV, which is "
            "empty); narrow the Q^2 window first"
        )
    return ep_min, ep_max


@dataclass(frozen=True, slots=True)
class GeneratorConfig:
    """All knobs of the generator.

    The defaults reproduce the ``rgb`` preset of the original ``aao_rad``
    driver script.
    """

    # --- physics / channel -------------------------------------------------
    theory: int = 7
    """Response-function model.  Only ``7`` (MAID07) is implemented."""

    channel: int = 3
    """Pion production channel: ``1`` = pi0 p, ``3`` = pi+ n."""

    polarized_beam: bool = True
    """Flip the beam helicity randomly per event."""

    # --- target geometry ---------------------------------------------------
    target_length_cm: float = 5.0
    target_radius_cm: float = 2.5
    beam_x_cm: float = 0.0
    beam_y_cm: float = 0.0
    beam_z_cm: float = -0.4

    # --- beam --------------------------------------------------------------
    beam_energy: float = 4.244
    """Incident electron energy [GeV]."""

    # --- acceptance --------------------------------------------------------
    q2_min: float = 0.2
    q2_max: float = 1.9
    ep_min: float = 0.3
    ep_max: float = 1.8
    """Scattered electron energy window [GeV]."""

    w_max: float | str | None = None
    """Reject trials whose hadronic mass exceeds this [GeV].  Defaults to the
    table's own upper edge (2.0); pass ``"clamp"`` for the original's
    saturating lookup.

    Reject trials whose hadronic mass exceeds this [GeV].

    The shipped MAID07 tables only span ``W in [1.08, 2.0]`` and
    ``Q^2 in [0, 5]``.  The original saturates out-of-range lookups
    (``multipole_amps.f90`` replaces ``W > 2`` with ``W = 2`` and ``Q^2 > 5``
    with ``Q^2 = 5``), so any kinematics reaching past the table edge returns
    the boundary row -- a frozen, unphysical response that also wrecks the
    importance sampling.  The shipped ``test.inp`` / ``clas12_test.inp`` run
    cards sit almost entirely outside the table for exactly this reason.

    By default this port rejects such trials instead, which keeps the sampled
    cross section inside the region the response actually describes.  Set to a
    number to choose the edge yourself, or to ``"clamp"`` to reproduce the
    original's saturating behaviour exactly.
    """

    # --- radiative sampling -------------------------------------------------
    min_photon_energy: float = 0.005
    """``delta`` [GeV]: photons softer than this are folded in analytically."""

    regions: tuple[float, float, float, float] = (0.20, 0.12, 0.20, 0.20)
    """Relative sizes of the four importance-sampling regions in ``cos(theta_k)``.

    Their sum must stay below 1; the remainder is sampled uniformly.
    """

    cos_step: float = 0.04
    """``csrng``: width of the narrow bands around the beam directions."""

    k_exp: float = 5.0
    """Slope of the logarithmic photon-energy variable ``uek = exp(-k_exp * ek)``."""

    ek_sampling: Literal["truncated", "fortran"] = "truncated"
    """How the photon energy is sampled.

    ``truncated``
        Draw ``ek`` from the exponential distribution *truncated* to
        ``[0, ek_max]``.  No trial is wasted on ``ek > ek_max``; the weight
        carries the corresponding Jacobian, so the sampled distribution is
        unchanged.  This is the default and is typically ~2x more efficient
        than the original.

    ``fortran``
        Reproduce the original's ``uek = exp(-k_exp * ek)`` transform driven by
        a 1e-9-resolution uniform RNG (which silently capped ``ek`` at
        ``ln(1e9)/k_exp``).
    """

    # --- output / selection -------------------------------------------------
    n_events: int = 1_000_000
    missing_mass_cut: float = 0.2
    """Half-width of the ``(mm^2 - mm_exp^2)`` window [GeV^2]."""

    n_tracks: int = 4
    """Particles declared in the LUND header (2 or 4)."""

    # --- run-time ----------------------------------------------------------
    seed: int | None = None
    batch_size: int = 65_536
    """Number of trial points evaluated per vectorised step."""

    interp_scheme: Literal["linear", "spline"] = "linear"
    """Table interpolation.  ``linear`` matches the original (which set
    ``method_spline = 2``)."""

    write_tracks: bool = True
    """Write all ``n_tracks`` particle lines.  The Fortran wrote only two
    track lines while declaring up to four in the header, which produces a
    file that most readers cannot load; set to ``False`` to reproduce it."""

    verbose_every: int = 0
    """Print progress every N accepted events (0 disables progress reporting)."""

    extra: dict[str, Any] = field(default_factory=dict)

    # ------------------------------------------------------------------
    def __post_init__(self) -> None:
        if self.theory != 7:
            raise ValueError(
                f"theory={self.theory} is not available; this port implements "
                "MAID07 (theory=7), the only model the original driver could reach"
            )
        if self.channel not in (1, 3):
            raise ValueError(
                f"channel={self.channel} is not available; use 1 (pi0 p) or 3 (pi+ n)"
            )
        if sum(self.regions) >= 1.0:
            raise ValueError(
                f"the region sizes must sum to less than 1, got {sum(self.regions):g}"
            )
        if self.n_tracks not in (2, 4):
            raise ValueError("n_tracks must be 2 or 4")
        if self.ek_sampling not in ("truncated", "fortran"):
            raise ValueError(
                f"ek_sampling={self.ek_sampling!r} is not available; "
                "use 'truncated' (default) or 'fortran'"
            )
        if self.interp_scheme not in ("linear", "spline"):
            raise ValueError(
                f"interp_scheme={self.interp_scheme!r} is not available; "
                "use 'linear' (matches the original) or 'spline'"
            )
        if self.min_photon_energy <= 0.0:
            raise ValueError("min_photon_energy (delta) must be positive")
        if not 0.0 < self.q2_min < self.q2_max:
            raise ValueError("require 0 < q2_min < q2_max")
        if not 0.0 < self.ep_min < self.ep_max:
            raise ValueError("require 0 < ep_min < ep_max")
        if self.batch_size < 1024:
            raise ValueError("batch_size must be at least 1024 to fill the GPU")
        if self.w_max is not None and self.w_max != "clamp":
            if not 1.1 <= float(self.w_max) <= 2.0:
                raise ValueError(
                    f"w_max={self.w_max} is outside the MAID07 table, which spans "
                    "W in [1.08, 2.0]; use 'clamp' to reproduce the original's "
                    "saturating lookup instead"
                )

    # ------------------------------------------------------------------
    def to_dict(self) -> dict[str, Any]:
        d = asdict(self)
        d["regions"] = list(self.regions)
        return d

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> GeneratorConfig:
        known = {f.name for f in fields(cls)}
        kwargs = {k: v for k, v in data.items() if k in known}
        if "regions" in kwargs:
            kwargs["regions"] = tuple(float(x) for x in kwargs["regions"])
        unknown = set(data) - known
        if unknown:
            raise ValueError(f"unknown configuration keys: {sorted(unknown)}")
        return cls(**kwargs)

    @classmethod
    def from_file(cls, path: str | Path) -> GeneratorConfig:
        """Load from JSON, or from YAML/TOML if the optional parser is present."""
        path = Path(path)
        text = path.read_text()
        suffix = path.suffix.lower()
        if suffix in (".yaml", ".yml"):
            try:
                import yaml  # type: ignore[import-untyped]
            except ModuleNotFoundError as exc:  # pragma: no cover
                raise ImportError("PyYAML is required to read YAML configs") from exc
            data = yaml.safe_load(text)
        elif suffix == ".toml":
            import tomllib

            data = tomllib.loads(text)
        else:
            data = json.loads(text)
        return cls.from_dict(data)

    def merge(self, **overrides: Any) -> GeneratorConfig:
        """Return a copy with ``overrides`` applied (``None`` values ignored)."""
        clean = {k: v for k, v in overrides.items() if v is not None}
        return replace(self, **clean)

    # ------------------------------------------------------------------
    @classmethod
    def from_legacy_input(cls, text: str) -> GeneratorConfig:
        """Parse the interactive prompt sequence of the original ``aao_rad``.

        The values, in order, are::

            theory, polarized, reg1..reg4, n_tracks, channel, mm_cut,
            t_target, r_target, x_beam, y_beam, z_beam, e_beam,
            q2_min, q2_max, ep_min, ep_max, delta, n_events, sigr_factor, sigr_max

        Blank lines and ``!`` comments are ignored.  ``sigr_factor``/``sigr_max``
        are dropped (the maximum is now estimated on the device).
        """
        tokens: list[float] = []
        for raw in text.splitlines():
            line = raw.split("!")[0].strip()
            if line:
                tokens.extend(float(tok) for tok in line.replace(",", " ").split())
        if len(tokens) < 21:
            raise ValueError(
                f"legacy input needs at least 21 values, found {len(tokens)}"
            )
        theory, pol = tokens[0], tokens[1]
        (
            r1, r2, r3, r4, n_tracks, channel, mm_cut, t_len, t_rad,
            bx, by, bz, e_beam, q2min, q2max, epmin, epmax, delta,
        ) = tokens[2:20]
        return cls(
            theory=int(theory),
            channel=int(channel),
            polarized_beam=bool(int(pol)),
            regions=(r1, r2, r3, r4),
            n_tracks=4 if int(n_tracks) == 4 else 2,
            missing_mass_cut=mm_cut,
            target_length_cm=t_len,
            target_radius_cm=t_rad,
            beam_x_cm=bx,
            beam_y_cm=by,
            beam_z_cm=bz,
            beam_energy=e_beam,
            q2_min=q2min,
            q2_max=q2max,
            ep_min=epmin,
            ep_max=epmax,
            min_photon_energy=delta,
            n_events=int(tokens[20]),
        )


#: Named presets matching the ``--experiment`` values of the old driver script.
#:
#: The ``(q2_min, q2_max)`` and ``(ep_min, ep_max)`` pairs are not independent.
#: ``Q^2`` and ``E'`` are sampled independently, so the pair of windows maps to a
#: diagonal band in ``W``; :func:`ep_window_for_w` gives the ``E'`` window that
#: keeps the whole rectangle inside the MAID07 table (``W in [1.08, 2.0]``).  The
#: values below come from that, rounded inwards.  The old driver script's
#: windows ignored this, which is why most of *its* sample landed in the frozen
#: edge row of the response tables.
EXPERIMENTS: dict[str, dict[str, Any]] = {
    "default": {
        "polarized_beam": True,
        "regions": (0.20, 0.12, 0.20, 0.20),
        "n_tracks": 4,
        "channel": 3,
        "missing_mass_cut": 0.2,
        "target_length_cm": 5.0,
        "target_radius_cm": 0.486,
        "beam_x_cm": 0.3,
        "beam_y_cm": 0.03,
        "beam_z_cm": 0.0,
        "beam_energy": 4.8,
        # The old driver asked for Q^2 up to 3.5 GeV^2 with E' in 0.1-4.25 GeV.
        # For a 4.8 GeV beam that rectangle reaches W ~ 3 GeV, and no E' window
        # at all can keep W <= 2 GeV across it (see ep_window_for_w), so the Q^2
        # window has to narrow first.
        "q2_min": 0.9,
        "q2_max": 2.5,
        "ep_min": 2.66,
        "ep_max": 3.29,
        "min_photon_energy": 0.005,
    },
    "rgb": {
        "polarized_beam": True,
        "regions": (0.20, 0.12, 0.20, 0.20),
        "n_tracks": 4,
        "channel": 3,
        "missing_mass_cut": 0.2,
        "target_length_cm": 5.0,
        "target_radius_cm": 2.5,
        "beam_x_cm": 0.0,
        "beam_y_cm": 0.0,
        "beam_z_cm": -0.4,
        "beam_energy": 4.244,
        "q2_min": 0.2,
        "q2_max": 1.9,
        # The old driver used 0.3-1.8 GeV here, which puts essentially the whole
        # sample at W > 2 GeV -- outside the MAID07 table and into its frozen
        # edge row.  ep_window_for_w gives 2.475-3.056 GeV for this Q^2 window.
        "ep_min": 2.48,
        "ep_max": 3.05,
        "min_photon_energy": 0.005,
    },
}
