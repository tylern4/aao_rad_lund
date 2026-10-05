"""Percentages in the stats summary must not collapse to ``0.000%``.

Acceptance is ``mean(weight) / ceiling`` against a heavy-tailed integrand, so on
the validation grid it spanned 0.0000% to 0.117% -- five orders of magnitude.
The old fixed three-decimal format printed the bottom of that range as
``0.000%``, which reads as "no events were accepted" rather than "very few" and
hides exactly the variation one compares runs by.
"""

from __future__ import annotations

import pytest

from aao_rad.generate import _pct


@pytest.mark.parametrize(
    ("fraction", "expected"),
    [
        (0.0, "0.0000%"),
        (1.0, "100.0000%"),
        (0.5, "50.0000%"),
        (0.001, "0.1000%"),
        # The real cases from the validation grid's slow configurations.
        (0.00117, "0.1170%"),  # 0.117%, the fastest observed
        (2.0e-6, "2.000e-04%"),  # ~0.0002%
        (4.6e-7, "4.600e-05%"),  # ~0.000046%
        # Just below and just above the scientific-notation threshold (1e-3%).
        (9.9e-6, "9.900e-04%"),
        (1.0e-5, "0.0010%"),  # exactly 1e-3%, fixed form
    ],
)
def test_pct_keeps_small_values_visible(fraction, expected):
    assert _pct(fraction) == expected


def test_pct_small_values_are_not_all_zero():
    """The specific failure: distinct small acceptances must stay distinct."""
    rendered = [_pct(f) for f in (1e-7, 1e-6, 1e-5, 1e-4)]
    assert len(set(rendered)) == len(rendered), rendered
    assert not any(r.startswith("0.0000%") for r in rendered), rendered


def test_pct_zero_still_reads_as_exactly_zero():
    """A genuine zero is not a small number; it should not get an exponent."""
    assert _pct(0.0) == "0.0000%"
