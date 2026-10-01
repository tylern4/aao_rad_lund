"""A reader for CLAS12-style LUND generator files.

Used by the test suite and by ``validation/`` to compare this port against the
Fortran original's output.

Two quirks of the original have to be handled:

* It declares ``npart = 4`` in the header but emits only two track lines, so its
  files have three lines per event.  This port writes all four tracks.
* It writes the four-momenta as ``px py pz E`` (``aao_rad.f90`` line 1092; its
  own comment there calls the order unresolved), whereas the CLAS12 convention
  -- and this port -- is ``E px py pz``.  The order is auto-detected from the
  mass shells, so callers always get ``energy, px, py, pz`` in that order.

Header layout (10 fields)::

    npart version subversion track_recl_version polarized beam_pid
    beam_energy target_pid target_mass beam_vertex_z

Track layout (14 fields)::

    trackid parentid parentpid pid status status2  f7 f8 f9 f10  mass  vx vy vz

where ``(f7..f10)`` is the four-momentum in one of the two orders above.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

HEADER_FIELDS = (
    "npart", "version", "subversion", "track_recl_version", "polarized",
    "beam_pid", "beam_energy", "target_pid", "target_mass", "beam_vertex_z",
)

#: Track fields, before the four-momentum order has been resolved.
_TRACK_FIELDS = (
    "track", "parentid", "parentpid", "pid", "status", "status2",
    "f7", "f8", "f9", "f10", "mass", "vx", "vy", "vz",
)

_N_TRACK_FIELDS = len(_TRACK_FIELDS)
_INT_KEYS = ("track", "parentid", "parentpid", "pid", "status", "status2")

#: Field order written by this port (``LundWriter``).
E_FIRST = "E_pxpy_pz"
#: Field order written by the Fortran original.
MOMENTA_FIRST = "pxpy_pz_E"


@dataclass
class LundFile:
    """A parsed LUND file: the run header plus one row per track.

    The four-momentum columns are exposed in the canonical CLAS12 order
    (:data:`E_FIRST`) whatever the file's own order was.
    """

    header: dict[str, float]
    tracks: dict[str, np.ndarray]
    #: The order the file itself used: :data:`E_FIRST` or :data:`MOMENTA_FIRST`.
    order: str = E_FIRST

    def __len__(self) -> int:
        return int(len(self.tracks["px"])) if self.tracks else 0

    @property
    def n_events(self) -> int:
        """Number of events, i.e. the number of track-1 lines."""
        if not len(self.tracks["px"]):
            return 0
        return int(np.count_nonzero(self.tracks["track"] == 1))

    @property
    def n_tracks_written(self) -> int:
        """Highest track number actually present in the file."""
        return int(self.tracks["track"].max()) if len(self.tracks["track"]) else 0

    def events(self) -> dict[str, np.ndarray]:
        """Collapse the per-track rows into per-event arrays.

        Slots the file never wrote come back as NaN, which is how the original's
        missing third and fourth tracks show up.
        """
        n = self.n_events
        # A row belongs to the event that opened on the preceding track-1 line.
        event_of_row = np.cumsum(self.tracks["track"] == 1) - 1
        out: dict[str, np.ndarray] = {
            "npart": np.full(n, self.header["npart"]),
            "beam_energy": np.full(n, self.header["beam_energy"]),
        }
        for slot in range(1, int(self.header["npart"]) + 1):
            here = self.tracks["track"] == slot
            target = event_of_row[here]
            for field in ("pid", "energy", "px", "py", "pz", "mass"):
                col = np.full(n, np.nan)
                col[target] = self.tracks[field][here]
                out[f"trk{slot}_{field}"] = col
        return out


def read_lund(
    path: str | Path,
    *,
    max_events: int | None = None,
    order: str | None = None,
) -> LundFile:
    """Parse a LUND generator file.

    Parameters
    ----------
    order
        ``None`` (default) auto-detects the four-momentum order; pass
        :data:`E_FIRST` or :data:`MOMENTA_FIRST` to force it.

    Raises
    ------
    ValueError
        If the file does not start with a header, or a track line has the wrong
        number of fields -- both of which usually mean the file was truncated
        mid-write, which matters because the original writes incrementally.
    """
    rows: dict[str, list[float]] = {k: [] for k in _TRACK_FIELDS}
    with Path(path).open() as fh:
        first = fh.readline()
        if not first.strip():
            raise ValueError(f"{path}: file is empty")
        header = dict(zip(HEADER_FIELDS, _floats(first, path, 0)))
        npart = int(header["npart"])
        line_no = 1
        for line in fh:
            line_no += 1
            if not line.strip():
                continue
            values = _floats(line, path, line_no)
            # The original repeats the 10-field header before every event, so a
            # line is a header when it has header width and starts with npart.
            if len(values) == len(HEADER_FIELDS) and int(values[0]) == npart:
                continue
            if len(values) != _N_TRACK_FIELDS:
                raise ValueError(
                    f"{path}:{line_no}: expected {_N_TRACK_FIELDS} fields in a "
                    f"track line, found {len(values)} -- the file looks truncated"
                )
            if not 1 <= int(values[0]) <= npart:
                raise ValueError(
                    f"{path}:{line_no}: track {int(values[0])} is outside the "
                    f"declared npart={npart}"
                )
            for key, value in zip(_TRACK_FIELDS, values):
                rows[key].append(value)

    tracks = {k: np.asarray(v, dtype=np.float64) for k, v in rows.items()}
    for key in _INT_KEYS:
        tracks[key] = tracks[key].astype(np.int64)

    if max_events is not None and len(tracks["track"]):
        event_of_row = np.cumsum(tracks["track"] == 1) - 1
        keep = event_of_row < max_events
        tracks = {k: v[keep] for k, v in tracks.items()}

    resolved = order or detect_order(tracks)
    _reorder_four_vectors(tracks, resolved)
    return LundFile(header=header, tracks=tracks, order=resolved)


def detect_order(tracks: dict[str, np.ndarray]) -> str:
    """Decide the four-momentum order from the mass shells.

    Read one way, ``E^2 - |p|^2`` is the particle mass squared -- zero for the
    electron, 0.0195 GeV^2 for a pion.  Read the other way it picks up ``-E^2``,
    which for a multi-GeV track is two orders of magnitude larger, so the
    comparison is decisive.  The electron alone is *not* enough: it is
    ultrarelativistic, so ``E ~ |p|`` and both readings give nearly zero.
    """
    f7, f8, f9, f10 = (tracks[k] for k in ("f7", "f8", "f9", "f10"))
    if not len(f7):
        return E_FIRST
    resid_e_first = np.abs(f7**2 - (f8**2 + f9**2 + f10**2)).mean()
    resid_momenta_first = np.abs(f10**2 - (f7**2 + f8**2 + f9**2)).mean()
    return E_FIRST if resid_e_first <= resid_momenta_first else MOMENTA_FIRST


def _reorder_four_vectors(tracks: dict[str, np.ndarray], order: str) -> None:
    """Rewrite ``f7..f10`` in place as canonical ``energy, px, py, pz``."""
    f7, f8, f9, f10 = (tracks[k] for k in ("f7", "f8", "f9", "f10"))
    if order == MOMENTA_FIRST:
        energy, px, py, pz = f10, f7, f8, f9
    else:
        energy, px, py, pz = f7, f8, f9, f10
    tracks["energy"], tracks["px"], tracks["py"], tracks["pz"] = energy, px, py, pz


def _floats(line: str, path: str | Path, line_no: int) -> list[float]:
    try:
        return [float(tok) for tok in line.split()]
    except ValueError as exc:
        raise ValueError(f"{path}:{line_no}: cannot parse {line.strip()!r}") from exc


def complete_events(path: str | Path) -> int:
    """How many whole events a possibly-truncated file contains.

    The original writes its LUND output incrementally and is routinely killed,
    so a partial trailing event is normal.  Both three-line (the original) and
    five-line (this port) layouts are recognised from the first header line.
    """
    lines = [ln for ln in Path(path).read_text().splitlines() if ln.strip()]
    if not lines:
        return 0
    npart = int(float(lines[0].split()[0]))
    for n_written in (npart, 2):
        block = n_written + 1
        whole, remainder = divmod(len(lines), block)
        if remainder == 0:
            return whole
    # Truncated tail: drop it and use the longest prefix that divides evenly.
    best = 0
    for n_written in (npart, 2):
        block = n_written + 1
        if len(lines) - (len(lines) % block) >= block:
            best = max(best, (len(lines) // block))
    return best