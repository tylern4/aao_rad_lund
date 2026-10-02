"""Tests for the MAID07 table reader and its converted ``.npz`` form.

The reader is the one piece of the port that can be checked *exactly* against
the original, so these tests are deliberately strict: the parsed amplitudes must
agree bit-for-bit with the values the Fortran reference reader produces.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from conftest import requires_tables

from aao_rad.table import (
    PADDING_COLUMNS,
    TABLE_FORMAT_VERSION,
    convert_table,
    converted_path_for,
    load_table,
    load_table_npz,
)

REPO = Path(__file__).resolve().parent.parent
REFERENCE_DUMP = REPO / "validation" / "fortran_table.npz"




#: The ``.tbl`` file behind the pi+ n channel.
PNPI_TBL = "maid07-PNpi.tbl"


def _convert(parms_dir: str):
    """Convert the pi+ table, forcing a rewrite."""
    tbl = Path(parms_dir) / "spp_tbl" / PNPI_TBL
    return convert_table(tbl, force=True)


@requires_tables
class TestParsedTable:
    def test_shape_and_dtype(self, pnpi_table):
        t = pnpi_table
        assert t.amps.shape == (62, len(t.q2), len(t.w))
        assert t.amps.dtype == np.float32
        assert np.isfinite(t.amps).all()

    def test_axes_are_increasing(self, pnpi_table):
        assert np.all(np.diff(pnpi_table.q2) > 0)
        assert np.all(np.diff(pnpi_table.w) > 0)

    def test_axis_extents_match_the_shipped_grid(self, pnpi_table):
        # mpintp.inc declares NVAR1=101 over Q^2 in [0, 5] and NVAR2=93 over
        # W in [1.08, 2.0].
        assert len(pnpi_table.q2) == 101
        assert len(pnpi_table.w) == 93
        assert pnpi_table.q2[0] == pytest.approx(0.0, abs=1e-6)
        assert pnpi_table.q2[-1] == pytest.approx(5.0, abs=1e-6)
        assert pnpi_table.w[0] == pytest.approx(1.08, abs=1e-6)
        assert pnpi_table.w[-1] == pytest.approx(2.0, abs=1e-6)

    def test_padding_columns_are_absent(self, pnpi_table):
        # 72 numbers per grid point on disk, 10 of which are the format's filler
        # columns.  If any leaked through, the amplitude count would be wrong.
        assert pnpi_table.amps.shape[0] == 72 - len(PADDING_COLUMNS)
        assert len(PADDING_COLUMNS) == 10
        assert len(set(PADDING_COLUMNS)) == 10

    def test_uniform_axes_are_detected(self, pnpi_table):
        assert pnpi_table.uniform

    def test_source_is_recorded(self, pnpi_table):
        assert pnpi_table.source.name == PNPI_TBL


@requires_tables
class TestConvertedFormat:
    def test_round_trip_is_exact(self, pnpi_table, parms_dir):
        # The ASCII carries ~5 significant digits and the Fortran declared SF as
        # REAL, so float32 is lossless; the converted file must be identical.
        _convert(parms_dir)
        again = load_table_npz(_convert(parms_dir))
        assert again.amps.dtype == np.float32
        assert np.array_equal(again.amps, pnpi_table.amps)
        assert np.array_equal(again.q2, pnpi_table.q2)
        assert np.array_equal(again.w, pnpi_table.w)

    def test_converted_path_is_a_sibling_of_the_ascii_tables(self, parms_dir):
        tbl = Path(parms_dir) / "spp_tbl" / PNPI_TBL
        p = converted_path_for(tbl)
        assert p.suffix == ".npz"
        assert p.parent.name == "tables"
        assert p.parent.parent == Path(parms_dir)

    def test_format_version_is_written(self, parms_dir):
        with np.load(_convert(parms_dir)) as data:
            assert int(data["format_version"]) == TABLE_FORMAT_VERSION

    def test_load_prefers_the_converted_file(
        self, parms_dir, pnpi_table, monkeypatch
    ):
        import aao_rad.table as tbl_mod

        path = _convert(parms_dir)
        assert path.is_file()

        # Make the ASCII parser explode: if the converted file is preferred,
        # the load still succeeds and the ASCII is never touched.
        def boom(*_a, **_k):
            raise AssertionError("ASCII table was re-parsed despite the .npz")

        monkeypatch.setattr(tbl_mod, "_parse_maid_file", boom)
        t = load_table(3, parms_dir=parms_dir, use_cache=False)
        assert np.array_equal(t.amps, pnpi_table.amps)
        # Provenance still names the .tbl the data came from.
        assert t.source.name == PNPI_TBL

    def test_explicit_out_dir_is_honoured(self, parms_dir, tmp_path):
        tbl = Path(parms_dir) / "spp_tbl" / PNPI_TBL
        out = convert_table(tbl, out_dir=tmp_path, force=True)
        assert out.parent == tmp_path
        assert np.array_equal(load_table_npz(out).amps, load_table_npz(_convert(parms_dir)).amps)


@requires_tables
@pytest.mark.skipif(
    not REFERENCE_DUMP.exists(),
    reason="Fortran reference dump not built; run validation/compare_table_to_fortran.py --save",
)
def test_matches_fortran_reference_bit_for_bit(pnpi_table):
    """The parsed amplitudes must equal the Fortran's, exactly.

    ``validation/dump_table.f90`` is a verbatim copy of the original
    ``read_sf_file.f90`` read sequence, corrected only in the ``M_{L-}`` block
    where the shipped file has 2 filler columns and the Fortran skips 4.
    """
    with np.load(REFERENCE_DUMP) as data:
        ref = data["amps"]
        points = data["points"]
    assert ref.shape == pnpi_table.amps.shape
    assert len(points) >= 10, "reference dump should cover many grid points"

    worst = 0.0
    worst_at = ""
    for iq2, iw in points:
        got = pnpi_table.amps[:, iq2, iw]
        want = ref[:, iq2, iw]
        denom = np.maximum(np.abs(want), 1e-30)
        rel = np.abs(got - want) / denom
        k = int(np.argmax(rel))
        if rel[k] > worst:
            worst = float(rel[k])
            worst_at = f"({iq2},{iw}) SF{k + 1}: fortran={want[k]:.8e} python={got[k]:.8e}"
    assert worst == 0.0, f"worst relative difference {worst:.3e} at {worst_at}"
