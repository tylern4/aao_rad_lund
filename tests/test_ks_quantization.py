"""KS must run at the precision the reference actually reports.

The Fortran writes its n-tuple with ``es16.8`` (src/aao_rad.f90:1082), so it
reports 8 significant digits and cannot express anything finer.  The port's npz is
float32.  On a continuous observable that difference is invisible, but on a
degenerate one it destroys the test.

``E_s`` at ``ebeam = 4.244`` GeV is a point mass: 97% of events sit at exactly
the beam energy.  4.244 is not representable in binary32, so the same physical
value is ``4.24399996`` in the Fortran's text and ``4.24399995803833`` in the
port's float32 -- 2e-9 apart.  Two atomic masses that close put the empirical
CDFs 97% apart, and KS reported **0.9739** for all 16 configurations at that beam
energy, which read as the two codes sampling disjoint physics.

Beam energies that *are* exact in binary32 (2.0, 6.0, 12.0) showed no such
artifact, and that is what identified the cause as float representation rather
than a port defect.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent.parent / "validation" / "perlmutter"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))

from compare_distributions import ks_test  # noqa: E402
from verify_statistics import FORT_SIG_DIGITS, ks_statistic, quantize_sig  # noqa: E402

N = 20_000
F32_4_244 = 4.24399995803833  # float32(4.244) widened to float64
F64_8DIG = 4.24399996  # what es16.8 prints


def point_mass(mode: float, n: int = N, seed: int = 0, frac: float = 0.97) -> np.ndarray:
    """97% at `mode`, with a thin tail below it -- the real shape at ebeam=4.244."""
    rng = np.random.default_rng(seed)
    n_mass = int(frac * n)
    tail = mode - np.abs(rng.normal(0.0, 0.02, n - n_mass))
    return np.concatenate([np.full(n_mass, mode), tail])


def test_the_real_case_collapses_once_quantized():
    a = point_mass(F64_8DIG, seed=1)
    b = point_mass(F32_4_244, seed=2)

    # The bug: a 2e-9 representation difference reported as a total disagreement.
    raw, _ = ks_test(a, b)
    assert raw > 0.9, f"expected the artifact to be visible, got {raw}"

    # The fix.
    assert ks_statistic(a, b) < 0.01, ks_statistic(a, b)


def test_quantization_is_at_the_reference_precision_not_the_finer_one():
    """Only quantizing to the precision the reference reports merges them.  At
    float64's full precision they stay apart, so this is a deliberate choice to
    compare at 8 digits, not a general rounding."""
    assert quantize_sig([F32_4_244])[0] == quantize_sig([F64_8DIG])[0]
    assert F32_4_244 != F64_8DIG
    assert quantize_sig([F32_4_244], digits=12)[0] != quantize_sig([F64_8DIG], digits=12)[0]


def test_a_real_difference_survives_quantization():
    """The fix must not blind the test to genuine disagreement."""
    a = point_mass(F64_8DIG, seed=1)
    b = point_mass(F64_8DIG - 0.05, seed=2)  # shifted by 5% of the mass
    assert ks_statistic(a, b) > 0.3, ks_statistic(a, b)


def test_exact_binary32_beam_energies_are_unaffected():
    """2.0, 6.0 and 12.0 are exact in float32, so they needed no rescue -- and
    must come through the quantized path unchanged."""
    for e in (2.0, 6.0, 12.0):
        assert quantize_sig([e])[0] == e


@pytest.mark.parametrize(
    "value, expected",
    [(0.0, 0.0), (-1.0, -1.0), (0.000123456789, 0.0001234568), (1e-30, 1e-30)],
)
def test_quantize_edge_cases(value, expected):
    assert quantize_sig([value])[0] == pytest.approx(expected, rel=1e-6)


def test_quantize_leaves_nan_and_inf_alone():
    out = quantize_sig([np.nan, np.inf, -np.inf, 1.5])
    assert np.isnan(out[0]) and np.isposinf(out[1]) and np.isneginf(out[2])
    assert out[3] == 1.5


def test_negative_values_round_correctly():
    """Log-of-absolute must not be used for the sign."""
    assert quantize_sig([-0.00587654321])[0] == pytest.approx(-0.0058765432, rel=1e-9)


def test_precision_is_documented_as_the_reference_precision():
    assert FORT_SIG_DIGITS == 8
