"""Writers for the CLAS12-style LUND generator file and the raw event record.

The original wrote three lines per event::

    <npart> <version> <1> <1> <flag_ehel> <11> <ebeam> <1> <1> <0.0>
    1  0.0  <parent> <pid>  <status> <0>  px py pz E  <mass> vx vy vz
    2  0.0  <parent> <pid>  <status> <0>  px py pz E  <mass> vx vy vz

with the momenta rotated by a random azimuth ``phir`` about the beam axis and
photons softer than ``delta`` replaced by a dummy ``1e-5`` 4-vector so the slot
always exists.

Deviations from the original, all deliberate and all affecting only the *file*,
never the sampled physics:

* **All ``npart`` tracks are written.**  The Fortran declared ``npart = 4`` in
  the header but emitted only two track lines, so files it produced were
  self-inconsistent.  The second nucleon and the photon momenta are known, so
  the missing lines are now emitted.  Pass ``write_tracks=False`` to reproduce
  the original's truncated file exactly.
* **The mass column holds the particle mass in GeV.**  The original wrote the
  PDG code there for the hadron line.
* **Fixed-width numeric formatting** instead of Fortran list-directed output.
  gfortran's list-directed layout is not portable and cannot be reproduced
  faithfully across compilers; the fields are whitespace delimited, so any
  reader sees the same numbers.

The writer is streaming: it takes ``(n, 32)`` blocks, so a billion-event run uses
the same memory as a million-event one.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterator, TextIO

import numpy as np

from .config import GeneratorConfig
from .constants import M_E, M_N, pi_masses
from .generate import EVENT_COLUMNS

__all__ = ["LundWriter", "PDG_CODES", "PHOTON_PDG", "write_npz"]

#: (pion PDG, nucleon PDG) for each production channel.
PDG_CODES: dict[int, tuple[int, int]] = {
    1: (111, 2212),  # p -> pi0 p
    2: (211, 2212),  # n -> pi- p
    3: (211, 2112),  # p -> pi+ n
    5: (111, 2112),  # n -> pi0 n
}

#: PDG code for a real photon.
PHOTON_PDG = 22

#: Energy written for a photon too soft to be resolved.
_DUMMY_ENERGY = 1.0e-5

_REAL_FMT = "{:16.8E}"
_INT_FMT = "{:>6d}"


def _reals(*values: float) -> str:
    """Format a run of reals in aligned fixed-width columns."""
    return "".join(_REAL_FMT.format(v) for v in values)


def _ints(*values: int) -> str:
    return "".join(_INT_FMT.format(v) for v in values)


class LundWriter:
    """Stream accepted events into a LUND generator file.

    Parameters
    ----------
    path
        Destination file; ``"-"`` writes to stdout.
    cfg
        Run configuration (supplies the channel, track count, beam energy and
        the soft-photon threshold).
    write_tracks
        Override :attr:`GeneratorConfig.write_tracks`.

    Examples
    --------
    >>> with LundWriter("out.lund", cfg) as w:          # doctest: +SKIP
    ...     for block, _ in generator.stream(cfg):
    ...         w.write(block)
    """

    def __init__(
        self, path: str | Path, cfg: GeneratorConfig, *, write_tracks: bool | None = None
    ) -> None:
        self.path = path
        self.cfg = cfg
        self.write_tracks = cfg.write_tracks if write_tracks is None else write_tracks
        self.n_written = 0
        self._fh: TextIO | None = None
        self._closeable = True
        self._rng: np.random.Generator | None = None

    # -- context management ---------------------------------------------
    def __enter__(self) -> LundWriter:
        if str(self.path) in ("-", ":-"):
            import sys

            self._fh = sys.stdout
            self._closeable = False
        else:
            self._fh = Path(self.path).open("w")
            self._closeable = True
        return self

    def __exit__(self, *exc: object) -> None:
        self.close()

    def close(self) -> None:
        if self._fh is not None and self._closeable:
            self._fh.close()
        self._fh = None

    # -- writing ---------------------------------------------------------
    def write(self, events: np.ndarray) -> int:
        """Append a block of events; returns the number written."""
        if self._fh is None:
            raise RuntimeError("use LundWriter as a context manager")
        n = len(events)
        if n == 0:
            return 0
        self._fh.write(self.format(events))
        self.n_written += n
        return n

    def format(self, events: np.ndarray) -> str:
        """Render a block of events as LUND text."""
        n = len(events)
        if n == 0:
            return ""

        col = {
            name: np.asarray(events[name], dtype=np.float64).reshape(-1)
            for name in EVENT_COLUMNS
        }
        cfg = self.cfg
        m_pi, _ = pi_masses(cfg.channel)
        pion_pdg, nucleon_pdg = PDG_CODES[cfg.channel]
        # Header: npart 1 1 1 flag_ehel 11 ebeam 1 1 0.0
        header = (
            _ints(cfg.n_tracks, 1, 1, 1, int(cfg.polarized_beam), 11)
            + _reals(float(cfg.beam_energy))
            + _ints(1, 1)
            + _reals(0.0)
            + "\n"
        )

        soft = col["eg"] <= cfg.min_photon_energy
        theta = col["theta"] * (np.pi / 180.0)
        # Rotate every momentum by the same random azimuth about the beam axis,
        # exactly as the original did with its per-event `phir`.
        phir = self._draw_rotation(len(events))
        rotc, rots = np.cos(phir), np.sin(phir)

        def rot(px: np.ndarray, py: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
            return px * rotc + py * rots, py * rotc - px * rots

        ep = col["ep"]
        epx, epy = rot(ep * np.sin(theta), np.zeros_like(ep))
        epz = ep * np.cos(theta)

        # The scattered electron always occupies track slot 1.
        trk = (epx, epy, epz, ep, M_E, 11)
        # `had` is the charged hadron, which is the proton for pi0 and the
        # pion for pi+; `oth` is the neutral partner, so the two always come
        # out in the slot order the header declares.
        if cfg.channel == 1:  # pi0: the proton is the charged track
            had_x, had_y, had_z, had_e, had_m, had_pid = (
                *rot(col["ppx"], col["ppy"]), col["ppz"], col["eprot"], M_N, nucleon_pdg
            )
            oth_x, oth_y, oth_z, oth_e, oth_m, oth_pid = (
                *rot(col["ppix"], col["ppiy"]), col["ppiz"], col["epi"], m_pi, pion_pdg
            )
        else:  # pi+: the pion is the charged track
            had_x, had_y, had_z, had_e, had_m, had_pid = (
                *rot(col["ppix"], col["ppiy"]), col["ppiz"], col["epi"], m_pi, pion_pdg
            )
            oth_x, oth_y, oth_z, oth_e, oth_m, oth_pid = (
                *rot(col["ppx"], col["ppy"]), col["ppz"], col["eprot"], M_N, nucleon_pdg
            )

        gam_x, gam_y = rot(np.where(soft, 0.0, col["egx"]), np.where(soft, 0.0, col["egy"]))
        gam_z = np.where(soft, _DUMMY_ENERGY, col["egz"])
        gam_e = np.where(soft, _DUMMY_ENERGY, col["eg"])

        vx, vy, vz = col["vx"], col["vy"], col["vz"]

        tracks = [
            (1, trk[5], trk[0], trk[1], trk[2], trk[3], trk[4]),
            (2, had_pid, had_x, had_y, had_z, had_e, had_m),
        ]
        if self.write_tracks and cfg.n_tracks == 4:
            tracks.append((3, oth_pid, oth_x, oth_y, oth_z, oth_e, oth_m))
            tracks.append((4, PHOTON_PDG, gam_x, gam_y, gam_z, gam_e, 0.0))

        # Track masses may be scalars (the electron and pion/proton masses are
        # constants for the run); broadcast them to per-event arrays.
        def column(v: np.ndarray | float) -> np.ndarray:
            return np.broadcast_to(np.asarray(v, dtype=np.float64), (n,))

        out: list[str] = []
        for i in range(n):
            out.append(header)
            for track_no, pid, x, y, z, e, mass in tracks:
                out.append(
                    _ints(track_no, 0)
                    + _reals(0.0)
                    + _ints(1, pid, 0, 0)
                    + _reals(
                        column(x)[i], column(y)[i], column(z)[i], column(e)[i],
                        column(mass)[i], vx[i], vy[i], vz[i],
                    )
                    + "\n"
                )
        return "".join(out)

    # -- helpers ---------------------------------------------------------
    def _draw_rotation(self, n: int) -> np.ndarray:
        """Per-event azimuth used to rotate momenta about the beam axis.

        The original drew this from the same stream as the event, so it is
        uniform on ``[0, 2 pi)``.  A dedicated generator seeded from the run
        seed keeps the rotation reproducible without perturbing the physics
        stream.
        """
        if self._rng is None:
            self._rng = np.random.default_rng(
                0 if self.cfg.seed is None else int(self.cfg.seed) + 1
            )
        return self._rng.uniform(0.0, 2.0 * np.pi, size=n)


def write_npz(path: str | Path, blocks: Iterator[np.ndarray], *, extra: dict | None = None) -> int:
    """Accumulate event blocks into a compressed ``.npz``.

    Returns the number of events written.  This is the lossless format: every
    n-tuple variable the original tracked is preserved.
    """
    chunks = [np.asarray(b) for b in blocks if len(b)]
    if not chunks:
        chunks = [np.zeros(0, dtype=_event_dtype())]
    data = np.concatenate(chunks)
    payload: dict[str, np.ndarray] = {name: data[name] for name in EVENT_COLUMNS}
    if extra:
        payload.update(extra)
    np.savez_compressed(path, **payload)
    return len(data)


def _event_dtype() -> np.dtype:
    from .generate import EVENT_DTYPE

    return EVENT_DTYPE
