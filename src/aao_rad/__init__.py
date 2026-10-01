"""GPU-accelerated radiative pion electroproduction event generator.

A JAX port of the Fortran ``aao_rad`` generator.  The physics -- MAID07
multipole amplitudes, the AO response functions, the Mo & Tsai radiative kernel
and the beam/target kinematics -- follows the original algorithm, but every
kinematic point in a batch is evaluated together on the device instead of one
at a time on the host.

Quick start
-----------
>>> from aao_rad import GeneratorConfig, build_grid, generate_lund
>>> cfg = GeneratorConfig(channel=3, n_events=100_000, seed=7)  # doctest: +SKIP
>>> n = generate_lund(cfg, "events.lund")                        # doctest: +SKIP

Command line::

    aao-rad --experiment rgb --n-events 1000000 --output events.lund
"""

from __future__ import annotations

from .amplitudes import (
    cgln_amplitudes,
    helicity_amplitudes,
    legendre_polynomials,
    multipole_amplitudes,
)
from .config import EXPERIMENTS, GeneratorConfig
from .constants import ALPHA, M_E, M_N, M_PI0, M_PIP, PI, pi_masses
from .generate import (
    EVENT_COLUMNS,
    EVENT_DTYPE,
    EventGenerator,
    GenerationStats,
    build_grid,
    build_kinematics,
)
from .interpolate import InterpolationGrid
from .kinematics import hadronic_final_state
from .lund import LundWriter
from .motsa import motsa_sigma, non_radiative_sigma, spence
from .table import MaidTable, load_table
from .xsection import Response, response_functions

__version__ = "1.0.0"

__all__ = [
    "ALPHA",
    "EVENT_COLUMNS",
    "EVENT_DTYPE",
    "EXPERIMENTS",
    "EventGenerator",
    "GeneratorConfig",
    "GenerationStats",
    "InterpolationGrid",
    "LundWriter",
    "M_E",
    "M_N",
    "M_PI0",
    "M_PIP",
    "MaidTable",
    "PI",
    "Response",
    "__version__",
    "build_grid",
    "build_kinematics",
    "cgln_amplitudes",
    "generate",
    "generate_lund",
    "hadronic_final_state",
    "helicity_amplitudes",
    "legendre_polynomials",
    "load_table",
    "motsa_sigma",
    "multipole_amplitudes",
    "non_radiative_sigma",
    "pi_masses",
    "response_functions",
    "spence",
]


def generate(cfg: GeneratorConfig, parms_dir: str | None = None):
    """Generate events for ``cfg``; see :meth:`EventGenerator.generate`."""
    gen = EventGenerator(build_grid(cfg.channel, parms_dir=parms_dir, scheme=cfg.interp_scheme))
    return gen.generate(cfg)


def generate_lund(
    cfg: GeneratorConfig,
    path: str,
    parms_dir: str | None = None,
    *,
    chunk_events: int = 500_000,
    on_progress=None,
):
    """Stream events straight into a LUND file without holding them all in memory.

    Returns the :class:`GenerationStats` for the run.
    """
    import logging

    gen = EventGenerator(
        build_grid(cfg.channel, parms_dir=parms_dir, scheme=cfg.interp_scheme)
    )
    log = logging.getLogger(__name__)
    stats = None
    with LundWriter(path, cfg) as writer:
        for block, blk_stats in gen.stream(cfg, chunk_events=chunk_events):
            writer.write(block)
            if on_progress is not None:
                on_progress(writer.n_written, cfg.n_events)
            if blk_stats is not None:
                stats = blk_stats
    assert stats is not None
    log.info("wrote %d events to %s", stats.n_events, path)
    return stats
