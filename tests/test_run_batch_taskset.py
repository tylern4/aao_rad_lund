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

import importlib.util
import sys
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
