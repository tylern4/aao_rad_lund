"""Shared fixtures and helpers for the test suite.

The tests need the real MAID07 parameter files, so they are located relative to
the repository root rather than to the installed package.  If ``parms/`` is
missing the table tests skip instead of failing, so the suite still runs
against a bare checkout.

``requires_tables`` is injected as a module-level global by the ``pytest``
plugin hook below, so test modules can simply ``from conftest import
requires_tables`` without needing ``tests`` to be a package.
"""

from __future__ import annotations

from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
PARMS = ROOT / "parms"
SPP_TBL = PARMS / "spp_tbl"


def _have_tables() -> bool:
    return SPP_TBL.is_dir() and any(SPP_TBL.glob("*.tbl"))


requires_tables = pytest.mark.skipif(
    not _have_tables(),
    reason="MAID07 parameter files are not unpacked (expected parms/spp_tbl/*.tbl)",
)


@pytest.fixture(scope="session")
def parms_dir() -> str:
    return str(PARMS)


@pytest.fixture(scope="session")
def pnpi_table(parms_dir: str):
    """The default table for the pi+ n channel."""
    from aao_rad.table import load_table

    return load_table("maid07-PNpi", parms_dir=parms_dir)


@pytest.fixture(scope="session")
def small_config():
    """A cheap configuration used by the end-to-end tests.

    The scattered-electron window is chosen so the sampled ``W`` stays inside
    the tabulated 1.08-2.0 GeV range, which is where the response is physical.
    """
    from aao_rad.config import GeneratorConfig

    return GeneratorConfig(
        channel=3,
        beam_energy=4.244,
        q2_min=0.2,
        q2_max=1.9,
        ep_min=1.6,
        ep_max=2.9,
        n_events=200,
        seed=20240501,
        batch_size=4096,
    )
