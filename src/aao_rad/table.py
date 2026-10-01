r"""Reader for the MAID07 multipole tables.

Each ``maid07-*.tbl`` file holds the 62 real/imaginary parts of the
``L_{0..5}^{\pm}``, ``E_{0..5}^{\pm}`` and ``M_{0..5}^{\pm}`` multipole
amplitudes on a regular ``(Q2, W)`` grid.  The layout of the original ASCII
files is, for every grid point::

    W=  1.08 Q2=   0.05000            <- 1 header line
       S_{L+}                         <- label
     <6 numbers>                      <- value line
     <6 numbers>
       S_{L-}
     ...
    (6 blocks x [1 label + 2 value lines] = 18 lines)

so 19 lines per grid point and 72 numbers, of which 10 are padding columns that
the Fortran reader threw away (``dumvar1..dumvar4``).

The parsed tables are cached as ``.npz`` files so that repeated runs (and
multiple processes) pay the parse cost only once.
"""

from __future__ import annotations

import hashlib
import logging
import os
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np

__all__ = [
    "MaidTable",
    "N_AMPS",
    "AMPS_PER_MULTIPOLE",
    "CHANNEL_FILES",
    "PADDING_COLUMNS",
    "TABLE_FORMAT_VERSION",
    "converted_path_for",
    "convert_all",
    "convert_table",
    "find_parms_dir",
    "load_table",
    "load_table_npz",
    "save_table_npz",
]

log = logging.getLogger(__name__)

#: Number of tabulated multipole components (real + imaginary parts).
N_AMPS = 62
#: Number of multipole partial waves (L = 0..5).
AMPS_PER_MULTIPOLE = 6

#: Lines consumed per grid point (1 header + 6 blocks x 3 lines).
_LINES_PER_POINT = 19
#: Numbers per grid point.
_NUMS_PER_POINT = 72
#: Zero-based offsets of the 10 padding columns inside each grid point.
#:
#: Each of the six blocks contributes 12 columns but the Fortran skipped a
#: leading run of them, because those partial waves are not tabulated:
#:
#:     block   columns   padding   waves actually tabulated
#:     S_{L+}      0-11        0      S0+..S5+   (6 waves, 12 values)
#:     S_{L-}     12-23        2      S1-..S5-   (5 waves, 10 values)
#:     E_{L+}     24-35        0      E0+..E5+   (6 waves, 12 values)
#:     E_{L-}     36-47        4      E2-..E5-   (4 waves,  8 values)
#:     M_{L+}     48-59        2      M1+..M5+   (5 waves, 10 values)
#:     M_{L-}     60-71        2      M1-..M5-   (5 waves, 10 values)
#:
#: 72 columns - 10 padding = 62 amplitudes, indexed SF1..SF62 as documented in
#: ``mpintp.inc``.
#:
#: NOTE: ``src/read_sf_file.f90`` skips **four** columns in the ``M_{L-}`` block
#: instead of two.  That makes it ask for 14 values where the file has 12, so a
#: Fortran list-directed read runs off the end of the block into the next grid
#: point's header and dies in its ``err=`` handler -- on the very first grid
#: point.  We follow the labeling in ``mpintp.inc`` instead; see
#: ``validation/dump_table.f90`` for the working reference reader.
PADDING_COLUMNS = (12, 13, 36, 37, 38, 39, 48, 49, 60, 61)
_PADDING = PADDING_COLUMNS

#: Fallback header pattern for tables that were re-wrapped and are no longer
#: column-aligned with the Fortran format.
_HEADER_RE = re.compile(
    r"^\s*W\s*=\s*(?P<w>[-+0-9.eEdD]+)\s*Q\s*2\s*=\s*(?P<q2>[-+0-9.eEdD]+)\s*$",
    re.IGNORECASE,
)

#: Default table file for each production channel (``epirea`` code).
CHANNEL_FILES = {
    1: "maid07-PPpi.tbl",  # p -> pi0 p
    2: "maid07-NPpi.tbl",  # n -> pi- p
    3: "maid07-PNpi.tbl",  # p -> pi+ n
    5: "maid07-NNpi.tbl",  # n -> pi0 n
}

#: Hard fallbacks for channels the Fortran driver never wired up.
_CHANNEL_FALLBACK = {2: "maid07-NPpi.tbl", 5: "maid07-PNpi.tbl"}


@dataclass(frozen=True)
class MaidTable:
    """A regular ``(Q2, W)`` grid of the 62 multipole components.

    Attributes
    ----------
    q2, w
        Grid coordinates in GeV^2 and GeV, shape ``(nq2,)`` and ``(nw,)``.
    amps
        Tabulated values, shape ``(62, nq2, nw)``.  Index ``i`` follows the
        MAID convention documented in ``mpintp.inc``::

            0..11   S0+ .. S5+
            12..21  S1- .. S5-
            22..33  E0+ .. E5+
            34..41  E2- .. E5-
            42..51  M1+ .. M5+
            52..61  M1- .. M5-
    uniform
        ``True`` when both axes are exactly equidistant, which enables the
        much faster closed-form grid-index computation during interpolation.
    source
        Path the table was loaded from (for provenance/reproducibility).
    """

    q2: np.ndarray
    w: np.ndarray
    amps: np.ndarray
    uniform: bool
    source: Path

    @property
    def q2_range(self) -> tuple[float, float]:
        return float(self.q2[0]), float(self.q2[-1])

    @property
    def w_range(self) -> tuple[float, float]:
        return float(self.w[0]), float(self.w[-1])

    def __repr__(self) -> str:  # pragma: no cover - debugging aid
        q2r, wr = self.q2_range, self.w_range
        return (
            f"MaidTable(nq2={self.q2.size}, nw={self.w.size}, "
            f"Q2={q2r[0]:g}..{q2r[1]:g}, W={wr[0]:g}..{wr[1]:g}, src={self.source.name})"
        )


def find_parms_dir(explicit: str | os.PathLike[str] | None = None) -> Path:
    """Locate the directory holding ``spp_tbl/maid07-*.tbl``.

    Resolution order:

    1. ``explicit``, if given.
    2. ``$AAO_RAD_PARMS``.
    3. ``$CLAS_PARMS`` (what the Fortran ``revinm`` routine used).
    4. ``<package>/parms`` and ``<repo root>/parms`` for in-tree use.

    Raises
    ------
    FileNotFoundError
        If no candidate directory contains ``spp_tbl``.
    """
    candidates: list[Path] = []
    if explicit is not None:
        candidates.append(Path(explicit))
    for var in ("AAO_RAD_PARMS", "CLAS_PARMS"):
        value = os.environ.get(var)
        if value:
            candidates.append(Path(value))
    # Walk up from this file looking for an in-tree parms/ directory.
    here = Path(__file__).resolve()
    candidates.extend(parent / "parms" for parent in here.parents[:4])
    candidates.append(Path.cwd() / "parms")

    seen: set[Path] = set()
    for cand in candidates:
        cand = cand.expanduser()
        if cand in seen:
            continue
        seen.add(cand)
        if (cand / "spp_tbl").is_dir():
            return cand
    raise FileNotFoundError(
        "Could not find the MAID tables. Unpack 'parms.tar.gz' and point the "
        "AAO_RAD_PARMS (or CLAS_PARMS) environment variable at the resulting "
        "'parms' directory, or pass parms_dir=... explicitly. Searched:\n  "
        + "\n  ".join(str(c) for c in seen)
    )


def _cache_dir() -> Path:
    root = os.environ.get("AAO_RAD_CACHE")
    base = Path(root).expanduser() if root else Path.home() / ".cache" / "aao_rad"
    base.mkdir(parents=True, exist_ok=True)
    return base


def _parse_maid_file(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Parse a MAID07 ``.tbl`` file into ``(q2, w, amps)``.

    The fast path assumes the well-formed 19-lines-per-point layout and
    validates it; anything unexpected falls back to a tolerant line-by-line
    parser so that tables with cosmetic differences still load.
    """
    text = path.read_text(errors="replace")
    lines = text.splitlines()

    headers = [i for i, line in enumerate(lines) if line.lstrip().startswith("W=")]
    if not headers:
        raise ValueError(f"{path}: no 'W=' grid headers found; not a MAID table?")

    fast_ok = len(headers) == len(lines) // _LINES_PER_POINT and all(
        headers[i] == i * _LINES_PER_POINT for i in range(len(headers))
    )

    value_text: str
    if fast_ok:
        keep: list[int] = []
        for i in range(len(lines)):
            if i % _LINES_PER_POINT in (0,):
                continue  # header
            if (i % _LINES_PER_POINT) % 3 == 1:
                continue  # block label
            keep.append(i)
        value_text = "\n".join(lines[i] for i in keep)
    else:
        log.debug("%s: non-standard layout, using tolerant parser", path.name)
        value_text = "\n".join(
            line for line in lines if line.strip() and not line.lstrip().startswith("W=")
        )

    values = _parse_floats(value_text, path)
    n_expected = len(headers) * _NUMS_PER_POINT
    if values.size != n_expected:
        raise ValueError(
            f"{path}: expected {_NUMS_PER_POINT} numbers per grid point "
            f"({n_expected} total) but parsed {values.size}"
        )

    values = values.reshape(len(headers), _NUMS_PER_POINT)
    amps = np.delete(values, _PADDING, axis=1).T.copy()  # (62, npoint)

    # Grid coordinates come from the header line itself ("W= 1.08 Q2= 0.05000").
    coords_w, coords_q2 = _parse_headers([lines[i] for i in headers], path)
    w = coords_w
    q2 = coords_q2
    # Each header contributes [W, Q2]; the Fortran read them as f4.2 and f7.5.
    # The file is ordered Q2-slowest, W-fastest (Fortran: do jvar1=Q2 outside,
    # do jvar2=W inside), so npoint = nq2 * nw and W varies within each block.
    nw = _count_unique(w)
    if nw <= 0 or len(headers) % nw:
        raise ValueError(
            f"{path}: {len(headers)} grid points are not divisible by the "
            f"{nw} distinct W values; the 19-lines-per-point layout is broken"
        )
    amps = amps.reshape(N_AMPS, len(headers) // nw, nw)

    # Rebuild the axes from the first Q2 row / first W column, then check the
    # whole block really is a mesh (this catches mis-aligned tables early).
    w_axis = w[:nw].astype(np.float64)
    q2_axis = q2.reshape(-1, nw)[:, 0].astype(np.float64)
    if not np.allclose(w.reshape(-1, nw), w_axis[None, :], atol=1e-6):
        raise ValueError(f"{path}: W is not constant within each Q2 block")
    if not np.allclose(q2.reshape(-1, nw), q2_axis[:, None], atol=1e-6):
        raise ValueError(f"{path}: Q2 is not constant within each W block")
    return q2_axis, w_axis, amps


def _parse_headers(lines: list[str], path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Extract ``(W, Q2)`` from the ``W=  1.08 Q2=   0.05000`` header lines.

    The number text is read first, because the header fields are *not* fixed
    width across tables: the Fortran's ``(A8, f4.2, A7, f7.5)`` only has room
    for two decimals of W, so a table with 0.001 W steps (``*-thresh.tbl``)
    widens the field and shifts every later column.  Slicing by column would
    then silently read garbage.  The Fortran fixed-column slice is kept as a
    fallback for a header that has lost its labels entirely.
    """
    w_out: list[float] = []
    q2_out: list[float] = []
    for line in lines:
        m = _HEADER_RE.match(line)
        if m is not None:
            w_out.append(float(m.group("w")))
            q2_out.append(float(m.group("q2")))
            continue
        try:
            w_out.append(float(line[8:12]))  # f4.2
            q2_out.append(float(line[19:26]))  # f7.5
        except ValueError:
            raise ValueError(
                f"{path}: unrecognised grid header {line!r}; expected "
                "'W= 1.08 Q2= 0.05000'"
            ) from None
    return np.asarray(w_out), np.asarray(q2_out)


def _count_unique(values: np.ndarray) -> int:
    """Number of distinct values, tolerating the f4.2/f7.5 rounding of repeats."""
    if values.size == 0:
        return 0
    return int(np.unique(np.round(values, 6)).size)


def _parse_floats(text: str, path: Path) -> np.ndarray:
    """Parse whitespace-separated floats, preferring the fast C parser."""
    try:
        return np.fromstring(text, dtype=np.float64, sep=" ")  # type: ignore[call-overload]
    except Exception:  # pragma: no cover - NumPy >= 2 deprecation fallback
        pass
    try:
        return np.array(text.split(), dtype=np.float64)
    except ValueError as exc:
        raise ValueError(f"{path}: could not parse numeric data ({exc})") from exc


def _is_uniform(axis: np.ndarray, tol: float = 1e-6) -> bool:
    if axis.size < 2:
        return False
    d = np.diff(axis)
    return bool(np.all(np.abs(d - d[0]) <= tol * max(1.0, abs(d[0]))))


def load_table(
    channel: int,
    *,
    parms_dir: str | os.PathLike[str] | None = None,
    table_dir: str | os.PathLike[str] | None = None,
    use_cache: bool = True,
) -> MaidTable:
    """Load the MAID07 table for a production ``channel``.

    Parameters
    ----------
    channel
        ``epirea`` code (1, 2, 3 or 5).  Channels 2 and 5 are not reachable
        from the original driver; they fall back to the closest available
        table with a warning.
    parms_dir
        Root parameter directory (contains ``spp_tbl/``).  Auto-detected when
        omitted.
    table_dir
        Explicit directory holding the ``.tbl`` file, bypassing ``spp_tbl``.
    use_cache
        Read/write the ``.npz`` cache of parsed tables.
    """
    if table_dir is not None:
        path = Path(table_dir).expanduser()
        if path.is_dir():
            path = path / CHANNEL_FILES.get(channel, CHANNEL_FILES[3])
    else:
        root = find_parms_dir(parms_dir)
        path = root / "spp_tbl" / CHANNEL_FILES.get(channel, CHANNEL_FILES[3])
    if not path.is_file():
        fallback = _CHANNEL_FALLBACK.get(channel)
        if fallback and (path.parent / fallback).is_file():
            log.warning(
                "no table for channel %s, falling back to %s", channel, fallback
            )
            path = path.parent / fallback
        else:
            raise FileNotFoundError(f"MAID table not found: {path}")

    cached = _cache_path(path) if use_cache else None
    if cached is not None and cached.is_file():
        try:
            return load_table_npz(cached)
        except Exception as exc:  # pragma: no cover - corrupt or stale cache
            log.debug("ignoring table cache %s (%s)", cached, exc)
            cached.unlink(missing_ok=True)

    # A deliberately pre-converted table is not a cache -- it is the intended
    # production input format -- so it is honoured even with use_cache=False.
    conv = converted_path_for(path)
    if conv.is_file():
        try:
            return load_table_npz(conv)
        except Exception as exc:  # pragma: no cover - corrupt or stale
            log.warning("ignoring unusable converted table %s (%s)", conv, exc)

    log.info("parsing MAID table %s", path)
    q2, w, amps = _parse_maid_file(path)
    uniform = _is_uniform(q2) and _is_uniform(w)
    table = MaidTable(q2=q2, w=w, amps=amps, uniform=uniform, source=path)

    if cached is not None:
        try:
            save_table_npz(table, cached)
        except OSError as exc:  # pragma: no cover - read-only cache dir
            log.debug("could not write table cache %s (%s)", cached, exc)
    return table


def _cache_path(path: Path) -> Path:
    stat = path.stat()
    key = hashlib.sha1(
        f"{path.resolve()}:{stat.st_size}:{int(stat.st_mtime)}:{N_AMPS}".encode()
    ).hexdigest()[:16]
    return _cache_dir() / f"{path.stem}-{key}.npz"


def iter_available_tables(parms_dir: str | os.PathLike[str] | None = None) -> list[Path]:
    """List the ``maid07-*.tbl`` files found in the parameter tree."""
    root = find_parms_dir(parms_dir)
    return sorted((root / "spp_tbl").glob("maid07-*.tbl"))


# ---------------------------------------------------------------------------
# Converted format
# ---------------------------------------------------------------------------
#: Bumped whenever the on-disk layout of the converted ``.npz`` changes.
TABLE_FORMAT_VERSION = 2


def save_table_npz(table: MaidTable, path: str | os.PathLike[str]) -> Path:
    """Write ``table`` in the converted format (see module docstring).

    The converted file is an uncompressed ``.npz`` holding five arrays:

    ``format_version``
        int scalar, :data:`TABLE_FORMAT_VERSION`.  Present so that a stale
        file written by an older release is rejected rather than misread.
    ``q2``, ``w``
        float64 grid axes, shape ``(nq2,)`` and ``(nw,)``, ascending.
    ``amps``
        float32 amplitudes, shape ``(62, nq2, nw)``, indexed as ``SF1..SF62``
        per ``mpintp.inc``.  The ASCII tables carry ~5 significant digits and
        the Fortran declares ``SF`` as ``REAL`` (single precision), so float32
        stores them exactly -- the same values the generator sees.
    ``source``
        the ``.tbl`` this was converted from, as a 0-d string array.

    Uncompressed so that the ~2.3 MB payload can be memory-mapped and paged in
    lazily instead of being inflated on every process start.
    """
    path = Path(path)
    tmp = path.with_name(f".{path.name}.tmp{os.getpid()}")
    path.parent.mkdir(parents=True, exist_ok=True)
    # Pass an open handle: np.savez(path, ...) would append a second ".npz".
    with tmp.open("wb") as fh:
        np.savez(
            fh,
            format_version=np.array(TABLE_FORMAT_VERSION, dtype=np.int32),
            q2=np.ascontiguousarray(table.q2, dtype=np.float64),
            w=np.ascontiguousarray(table.w, dtype=np.float64),
            amps=np.ascontiguousarray(table.amps, dtype=np.float32),
            source=np.array(table.source.name),
        )
    os.replace(tmp, path)
    return path


def load_table_npz(path: str | os.PathLike[str]) -> MaidTable:
    """Read the converted format written by :func:`save_table_npz`."""
    with np.load(Path(path), allow_pickle=False) as npz:
        version = int(npz["format_version"])
        if version != TABLE_FORMAT_VERSION:
            raise ValueError(
                f"{path}: converted table format version {version}, "
                f"expected {TABLE_FORMAT_VERSION}; re-run 'aao-rad-table convert'"
            )
        amps = npz["amps"]
        if amps.shape[0] != N_AMPS:
            raise ValueError(
                f"{path}: expected {N_AMPS} amplitudes, found {amps.shape[0]}"
            )
        source = Path(str(npz["source"])) if "source" in npz else Path(path)
        q2, w = npz["q2"], npz["w"]
        return MaidTable(
            q2=q2,
            w=w,
            amps=amps,
            uniform=_is_uniform(q2) and _is_uniform(w),
            source=source,
        )


def converted_path_for(tbl: Path, out_dir: str | os.PathLike[str] | None = None) -> Path:
    """Where the converted twin of ``tbl`` lives.

    Defaults to ``<spp_tbl>/../tables/<stem>.npz`` -- a sibling ``tables``
    directory next to ``spp_tbl`` -- so converted data never pollutes the
    original, read-only parameter tree.
    """
    if out_dir is not None:
        return Path(out_dir) / f"{tbl.stem}.npz"
    return tbl.parent.parent / "tables" / f"{tbl.stem}.npz"


def convert_table(
    tbl: str | os.PathLike[str],
    out_dir: str | os.PathLike[str] | None = None,
    *,
    force: bool = False,
) -> Path:
    """Convert one MAID ``.tbl`` to the fast format and return the output path."""
    tbl = Path(tbl)
    out = converted_path_for(tbl, out_dir)
    if out.is_file() and not force:
        log.debug("%s already converted", out)
        return out
    q2, w, amps = _parse_maid_file(tbl)
    table = MaidTable(
        q2=q2, w=w, amps=amps, uniform=_is_uniform(q2) and _is_uniform(w), source=tbl
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    save_table_npz(table, out)
    log.info(
        "converted %s -> %s (%d x %d x %d)", tbl.name, out, N_AMPS, q2.size, w.size
    )
    return out


def convert_all(
    parms_dir: str | os.PathLike[str] | None = None,
    out_dir: str | os.PathLike[str] | None = None,
    *,
    force: bool = False,
) -> list[Path]:
    """Convert every MAID table found in the parameter tree."""
    root = find_parms_dir(parms_dir)
    out: list[Path] = []
    for tbl in sorted((root / "spp_tbl").glob("maid07-*.tbl")):
        out.append(convert_table(tbl, out_dir, force=force))
    return out
