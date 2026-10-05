"""stream() must hand the caller blocks it can actually hold on to.

The CLI accumulates every block a run yields, because the whole point of the
generator is that the n-tuple is written once at the end.  A run needs
thousands of iterations to reach 20k events at the original's sub-percent
acceptance, and the device-side buffer is (chunk_events, 32) float32 -- 64 MB at
the default chunk.  Yielding a *view* of that buffer therefore pinned 64 MB per
iteration, ~192 GB per run, which is what the 103 OOM kills on the 64-worker
Perlmutter arm and the ~100 GB of host memory per GPU process actually were.

What matters is not whether the yielded array is a view -- ``view()`` and
``reshape()`` always produce one -- but how much memory stays *reachable* from
it.  A view of a small owning array is harmless; a view of the 500k-row buffer
is the whole problem.  So these tests walk the ``base`` chain and measure what
would still be resident if the caller kept every block.
"""

from __future__ import annotations

import numpy as np
import pytest

from aao_rad.generate import EVENT_DTYPE, EventGenerator, build_grid

# Big enough that retaining one full buffer dwarfs what 200 events need, so the
# tests fail loudly if a whole buffer is handed out or kept alive.
CHUNK = 200_000

requires_tables = pytest.mark.skipif(
    not (__import__("pathlib").Path(__file__).resolve().parent.parent / "parms" / "spp_tbl").is_dir(),
    reason="MAID07 parameter files are not unpacked (expected parms/spp_tbl/*.tbl)",
)


def retained_bytes(arrays) -> int:
    """Bytes still resident if every array in ``arrays`` is kept alive.

    Follows each array to the root of its ``base`` chain -- the allocation that
    cannot be freed while the array lives -- and counts each root once.
    """
    roots: dict[int, np.ndarray] = {}
    for arr in arrays:
        root = arr
        # Stop at the first non-array base: device_get hands back a PyCapsule
        # underneath, which owns nothing we can measure.
        while isinstance(getattr(root, "base", None), np.ndarray):
            root = root.base
        if isinstance(root, np.ndarray):
            roots[id(root)] = root
    return sum(root.nbytes for root in roots.values())


@pytest.fixture(scope="module")
def generator(parms_dir):
    return EventGenerator(build_grid(3, parms_dir=parms_dir))


@requires_tables
def test_blocks_do_not_keep_the_whole_chunk_buffer_alive(generator, small_config):
    blocks = [
        block
        for block, _ in generator.stream(small_config, chunk_events=CHUNK, n_events=200)
    ]

    assert blocks, "stream() yielded nothing"
    whole_buffer = CHUNK * EVENT_DTYPE.itemsize
    kept = retained_bytes(blocks)
    assert kept < whole_buffer // 4, (
        f"keeping {len(blocks)} blocks retains {kept / 2**20:.0f} MB, but 200 "
        f"events need only {200 * EVENT_DTYPE.itemsize / 2**20:.3f} MB and one "
        f"full buffer is {whole_buffer / 2**20:.0f} MB"
    )


@requires_tables
def test_accumulating_every_block_costs_proportional_to_the_events(generator, small_config):
    """What the CLI actually does: keep every block, then write one npz."""
    blocks = [
        block
        for block, _ in generator.stream(small_config, chunk_events=CHUNK, n_events=200)
    ]

    assert sum(block.size for block in blocks) == 200
    # Scale the check by the chunk width: the old behaviour grew linearly with
    # it (~192 GB for a 20k-event run at the default 500k chunk), the fixed one
    # does not grow with it at all.
    per_chunk = retained_bytes(blocks) / CHUNK
    assert per_chunk * 500_000 < 64 * 2**20, (
        f"retaining blocks costs {per_chunk * 500_000 / 2**20:.0f} MB per 500k "
        "chunk, which extrapolates to the old ~192 GB per run"
    )


@requires_tables
def test_a_big_chunk_does_not_change_the_events(generator, small_config):
    """The chunk is a buffer size, not physics: the same seed gives the same
    events however wide the buffer is, so shrinking it is free if it ever needs
    to be."""

    def collect(chunk):
        return np.concatenate(
            [
                block
                for block, _ in generator.stream(
                    small_config, chunk_events=chunk, n_events=200
                )
            ]
        )

    assert np.array_equal(collect(CHUNK), collect(1_000))