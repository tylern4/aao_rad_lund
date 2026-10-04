"""Tests for the CPU slicing in the Perlmutter batch driver.

``taskset_prefix`` decides which hardware threads each worker process gets, and
getting it wrong is invisible until the scan runs: two processes sharing one core
just run slower, with no error anywhere.  The node layout is therefore worth
pinning down here, where the properties are checkable in a second.

A Perlmutter PM-CPU node is 128 cores with 2 hardware threads each, so its
affinity mask has 256 entries.  Which entries belong to one core depends on how
Slurm numbered the mask, and a probe on the node (``siblings of core 0 = [0,
128]``) says it is thread-major: entry *c* is thread *c mod 128* of core *c div
128*.  Nothing guarantees that, so the tests below run against both numberings --
the slicing must not depend on which one it gets.
"""

from __future__ import annotations

import argparse
import importlib.util
import sys
import threading
import time
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
DRIVER = ROOT / "validation" / "perlmutter" / "run_batch.py"


def _load_driver():
    spec = importlib.util.spec_from_file_location("aao_run_batch", DRIVER)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


run_batch = _load_driver()

NCORES = 128
MASK = list(range(2 * NCORES))
# How each mask entry maps to a physical core, per numbering scheme.
NUMBERINGS = {
    # Verified on the node: siblings of core 0 are [0, 128].
    "thread-major": {c: c % NCORES for c in MASK},
    # Siblings of core 0 would be [0, 1].
    "core-major": {c: c // 2 for c in MASK},
}


@pytest.fixture(params=sorted(NUMBERINGS), ids=sorted(NUMBERINGS))
def node(request, monkeypatch):
    """A simulated PM-CPU node: its sorted mask, core grouping, cpu -> core map."""
    cpu2core = NUMBERINGS[request.param]
    groups: dict[int, list[int]] = {}
    for cpu in MASK:
        groups.setdefault(cpu2core[cpu], []).append(cpu)
    # core_groups returns groups in order of first appearance in the sorted mask.
    ordered = [groups[k] for k in sorted(groups)]
    monkeypatch.setattr(run_batch, "allocated_cpus", lambda: list(MASK))
    monkeypatch.setattr(run_batch, "core_groups", lambda cpus: ordered)
    # taskset_prefix short-circuits when taskset is absent, which it is on macOS.
    monkeypatch.setattr(run_batch.shutil, "which", lambda name: "/usr/bin/taskset")
    return type("Node", (), {"cpu2core": cpu2core, "ngroups": len(ordered)})()


def slice_of(node, worker, cpus_per_job, jobs):
    prefix = run_batch.taskset_prefix(worker, cpus_per_job, jobs)
    assert prefix[0] == "taskset"
    return [int(c) for c in prefix[2].split(",")]


LAYOUTS = [(128, 1), (128, 2), (64, 4), (32, 8), (16, 16), (4, 32), (2, 128)]


@pytest.mark.parametrize(("jobs", "cpus_per_job"), LAYOUTS)
def test_each_worker_gets_the_cores_it_asked_for(node, jobs, cpus_per_job):
    """A worker pinned to N CPUs must span N distinct cores, not N/2.

    Cutting N consecutive entries of the mask gives N distinct cores under
    thread-major numbering but only N/2 under core-major, where consecutive
    entries are hyperthread siblings.  That is why the slicing cannot assume a
    numbering.  It also wasted half the node under thread-major at ``4 x 32``,
    where four workers' 32-entry blocks covered only 64 of the 128 cores.
    """
    want = min(cpus_per_job, node.ngroups)
    for worker in range(jobs):
        cores = {node.cpu2core[c] for c in slice_of(node, worker, cpus_per_job, jobs)}
        assert len(cores) == want, f"worker {worker} got {len(cores)} cores, wanted {want}"


@pytest.mark.parametrize(("jobs", "cpus_per_job"), LAYOUTS)
def test_every_worker_gets_the_cpu_count_it_asked_for(node, jobs, cpus_per_job):
    """The slice must hold exactly cpus_per_job CPUs, since that sets XLA's pool."""
    for worker in range(jobs):
        got = slice_of(node, worker, cpus_per_job, jobs)
        assert len(got) == min(cpus_per_job, 2 * node.ngroups)
        assert len(set(got)) == len(got), f"worker {worker} was given a CPU twice"


@pytest.mark.parametrize(("jobs", "cpus_per_job"), LAYOUTS)
def test_no_core_is_left_idle(node, jobs, cpus_per_job):
    """No core may sit unused while another is running two workers' threads."""
    used = set()
    for worker in range(jobs):
        used |= {node.cpu2core[c] for c in slice_of(node, worker, cpus_per_job, jobs)}
    assert len(used) == node.ngroups, f"{node.ngroups - len(used)} cores idle"


@pytest.mark.parametrize(("jobs", "cpus_per_job"), LAYOUTS)
def test_slices_partition_the_mask(node, jobs, cpus_per_job):
    """Where the workers ask for exactly the mask, no CPU is used twice.

    A duplicated CPU means two processes pinned to one hardware thread; a CPU
    outside the mask is a fatal error under a cgroup-enforced step.
    """
    seen = [c for worker in range(jobs) for c in slice_of(node, worker, cpus_per_job, jobs)]
    assert len(seen) == len(set(seen)), "a CPU was handed to more than one worker"
    if jobs * cpus_per_job == len(MASK):
        assert sorted(seen) == MASK


def test_unpinned_worker_gets_no_taskset(node):
    """``--cpus-per-job 0`` means the whole allocation, so no restriction at all."""
    assert run_batch.taskset_prefix(0, 0, 1) == []
    assert run_batch.taskset_prefix(3, -1, 4) == []


# The thread-count variables a pinned worker must be told the width of, and the
# value a worker must never keep when it has been given fewer cores than that.
THREAD_VARS = ["OMP_NUM_THREADS", "JAX_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"]


def _env_for(node, tmp_path, cpus_per_job, jobs=4, worker=0, code="py_cpu"):
    """The environment build_job hands to one child process."""
    grid = tmp_path / "grid"
    grid.mkdir(exist_ok=True)
    (grid / "cfg_000.txt").write_text("run card placeholder\n")
    run = {"run_id": "r0", "cfg_id": "cfg_000", "n_events": "10", "seed": "1"}
    _, env, _, _ = run_batch.build_job(
        run,
        code,
        tmp_path,
        tmp_path,
        sys.executable,
        cpus_per_job,
        jobs,
        worker,
        0,
    )
    return env


@pytest.mark.parametrize("cpus_per_job", [1, 2, 4])
def test_pinned_worker_is_told_its_slice_width(node, tmp_path, monkeypatch, cpus_per_job):
    """A worker pinned to N CPUs must be told N threads by every knob.

    The sbatch scripts export these for the whole node, and a pinned worker
    inherits them.  run_py_cpu.sbatch sets JAX_NUM_THREADS=128 for the node while
    handing each of its 128 workers 2 cores, so a worker that kept the inherited
    value would be told to use 128 threads on 2 cores.
    """
    monkeypatch.setenv("JAX_NUM_THREADS", "128")
    monkeypatch.setenv("OMP_NUM_THREADS", "128")
    env = _env_for(node, tmp_path, cpus_per_job)
    for var in THREAD_VARS:
        assert env[var] == str(cpus_per_job), f"{var}={env[var]}, wanted {cpus_per_job}"


def test_unpinned_worker_inherits_the_sbatch_thread_count(node, tmp_path, monkeypatch):
    """With no slice there is nothing to correct, so the sbatch value stands."""
    monkeypatch.setenv("JAX_NUM_THREADS", "128")
    monkeypatch.setenv("OMP_NUM_THREADS", "128")
    env = _env_for(node, tmp_path, 0)
    assert env["OMP_NUM_THREADS"] == "128"
    assert env["JAX_NUM_THREADS"] == "128"


def test_core_groups_partitions_the_real_mask():
    """Against this machine: every CPU lands in exactly one group, once."""
    cpus = run_batch.allocated_cpus()
    groups = run_batch.core_groups(cpus)
    flat = [c for group in groups for c in group]
    assert sorted(flat) == sorted(cpus)
    assert all(groups), "no group may be empty"


# --- seed stagger -------------------------------------------------------------
#
# The Fortran original seeds myran from unixtime at 1 s resolution and cannot be
# told a seed, so two runs of one configuration that start in the same second
# produce byte-identical output.  With every run in flight at once that is not an
# edge case: job 59269588 came back with 67 of 96 configurations having both
# seeds identical, which also throws away half of each configuration's
# statistics.  The stagger is the only thing standing between the layout and a
# silently degenerate seed axis, so it is pinned here.


def test_fortran_sibling_runs_are_pushed_ahead_of_each_other():
    """Two runs of one configuration: the second waits out the gap."""
    run_batch._stagger_last.clear()
    assert run_batch.stagger_seed("fortran", "cfg_000", gap=0.05) == pytest.approx(0, abs=0.05)
    delay = run_batch.stagger_seed("fortran", "cfg_000", gap=0.05)
    assert delay == pytest.approx(0.05, abs=0.03)


def test_unrelated_configurations_do_not_wait_on_each_other():
    """The gap is per configuration, so it does not serialise the whole batch."""
    run_batch._stagger_last.clear()
    assert run_batch.stagger_seed("fortran", "cfg_a", gap=5.0) == pytest.approx(0, abs=0.05)
    assert run_batch.stagger_seed("fortran", "cfg_b", gap=5.0) == pytest.approx(0, abs=0.05)


@pytest.mark.parametrize("code", ["py_cpu", "py_gpu"])
def test_the_port_is_never_delayed(code):
    """The port takes --seed from the manifest, so it has nothing to stagger."""
    run_batch._stagger_last.clear()
    assert run_batch.stagger_seed(code, "cfg_000", gap=5.0) == 0.0
    assert run_batch.stagger_seed(code, "cfg_000", gap=5.0) == 0.0


def test_a_forced_batch_of_siblings_still_gets_distinct_seconds():
    """All runs of one configuration launched at once, as the thread pool does.

    Each waits for its own slot, so the k-th is at least (k-1) gaps behind the
    first -- the property that matters is that no two land in the same second.
    """
    run_batch._stagger_last.clear()
    gap = 0.05
    delays = []
    lock = threading.Lock()

    def worker():
        d = run_batch.stagger_seed("fortran", "cfg_000", gap=gap)
        with lock:
            delays.append(d)

    threads = [threading.Thread(target=worker) for _ in range(4)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    starts = sorted(time.time() - d for d in delays)
    assert all(b - a >= gap * 0.9 for a, b in zip(starts, starts[1:])), starts
    assert {round(d, 3) for d in delays} == {round(i * gap, 3) for i in range(4)}


def test_gap_exceeds_the_clock_resolution_the_fortran_sees():
    """unixtime() has 1 s resolution, so anything under it is not a stagger."""
    assert run_batch.SEED_STAGGER_SECONDS >= 2.0


def test_run_all_actually_staggers_before_launching(tmp_path, monkeypatch):
    """The unit tests above would still pass if the dispatch path forgot to call
    it, and that is the failure that produced 67 identical seed pairs, so the
    call itself is pinned: every pending run is staggered, and always before its
    own execute()."""
    order: list[str] = []
    monkeypatch.setattr(run_batch, "stagger_seed", lambda code, cfg: order.append(f"stagger:{cfg}") or 0.0)
    monkeypatch.setattr(
        run_batch, "build_job", lambda run, *a, **k: ([f"run-{run['run_id']}"], {}, tmp_path, None)
    )
    monkeypatch.setattr(
        run_batch,
        "execute",
        lambda job: (order.append(f"exec:{job[0][0]}"), ("cfg_000_s0", 0, 0.0))[1],
    )
    monkeypatch.setattr(run_batch, "run_complete", lambda code, run_dir: False)

    args = argparse.Namespace(
        code="fortran", gpus=0, jobs=2, cpus_per_job=1, dry_run=False, python=sys.executable
    )
    runs = [
        {"run_id": "cfg_000_s0", "cfg_id": "cfg_000", "seed": "1"},
        {"run_id": "cfg_000_s1", "cfg_id": "cfg_000", "seed": "2"},
    ]
    assert run_batch.run_all(args, runs, tmp_path, tmp_path) == 0
    assert sorted(order) == [
        "exec:run-cfg_000_s0",
        "exec:run-cfg_000_s1",
        "stagger:cfg_000",
        "stagger:cfg_000",
    ], order
    first_exec = min(i for i, e in enumerate(order) if e.startswith("exec:"))
    assert first_exec > 0, order
    # Each run is staggered immediately before its own launch, so every exec is
    # directly preceded by a stagger (with 2 workers the two interleave).
    for i, entry in enumerate(order):
        if entry.startswith("exec:"):
            assert order[i - 1].startswith("stagger:"), order
