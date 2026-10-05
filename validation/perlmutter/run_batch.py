#!/usr/bin/env python3
"""Run every manifest job for one code, in parallel, inside a Slurm allocation.

One driver serves all three implementations so batching, resumability and
bookkeeping stay identical between them:

    --code fortran   the original: one CWD per run (the MAID tables are opened
                     relative to CWD), card on stdin, AAO_TRIALS=0 so the
                     trial dump is off but the 32-variable event n-tuple is
                     still written
    --code py_cpu    the port on CPU (JAX_PLATFORMS=cpu), npz output
    --code py_gpu    the port on GPU; one worker per GPU, each pinned to its
                     device via CUDA_VISIBLE_DEVICES

The node layout is processes x threads, and each arm picks its own:

    fortran   one process per hardware thread -- the original is serial Fortran
              with no OpenMP directives, so one thread is all a run can use
    py_cpu    the same: one thread per run with every run in flight.  A run's
              cost grows faster than the number of cores it is given, so
              widening one only takes capacity away from a run that has none
    py_gpu    one process per GPU, each given an equal share of the host
              threads -- here the threads feed the device rather than saturate
              a core

Each run is an independent 20k-event batch, so runs are the unit of parallelism
and the only question a node layout has to answer is how many of them to have in
flight at once.

Runs are addressed by ``run_id`` from the manifest written by ``make_grid.py``
and land in ``<root>/<code>/<run_id>/``.  A run whose stdout carries the final
cross-section line is complete and skipped, so every batch script is
resumable.

The port is invoked with ``--n-events`` and ``--seed`` from the manifest --
the CLI flags override the card, and the card carries no seed, so both must be
explicit.  ``--ek-sampling fortran`` matches the original's RNG-limited photon
sampling, which is what the validated local comparison uses.

For the Fortran there is one extra hazard: it seeds ``myran`` from
``unixtime`` (1 s resolution), so two runs starting in the same second sample
identically.  The manifest orders the two runs of each configuration far
apart, and ``--check-collisions`` compares md5 checksums of every
same-configuration pair afterwards and re-runs any collided pair
sequentially so the seeds differ.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import shutil
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

SIGMA_MARKERS = {
    "fortran": "Integrated cross section",
    "py_cpu": "sigma (MC)",
    "py_gpu": "sigma (MC)",
}


def default_root() -> Path:
    env = os.environ.get("AAO_SCAN_ROOT")
    if env:
        return Path(env)
    scratch = os.environ.get("SCRATCH", "/tmp")
    return Path(scratch) / "aao_rad_scan"


def load_manifest(root: Path) -> list[dict]:
    with open(root / "grid" / "manifest.csv", newline="") as fh:
        return list(csv.DictReader(fh))


def run_complete(code: str, run_dir: Path) -> bool:
    """A run is complete when its final cross section was printed."""
    out = run_dir / "out.txt"
    if not out.is_file():
        return False
    try:
        text = out.read_text(errors="replace")
    except OSError:
        return False
    if SIGMA_MARKERS[code] not in text:
        return False
    if code == "fortran":
        return (run_dir / "aao_rad.ntuple").is_file()
    return (run_dir / "out.npz").is_file()


def allocated_cpus() -> list[int]:
    """The CPUs this process may actually run on.

    Slurm hands out an arbitrary slice of the node, not CPUs 0..n-1, so pinning
    to a guessed core can land outside the allocation -- which a cgroup-enforced
    step treats as a fatal error.  Read the real mask instead.
    """
    try:
        cpus = sorted(os.sched_getaffinity(0))
    except AttributeError:  # not Linux
        cpus = list(range(os.cpu_count() or 1))
    return cpus or [0]


def core_groups(cpus: list[int]) -> list[list[int]]:
    """Group the allocation's CPUs by the physical core each belongs to.

    A PM-CPU node's affinity mask carries one entry per hardware thread, so two
    entries that sit next to each other in the mask are usually the two
    hyperthreads of a single core.  Knowing the grouping is what lets a worker be
    handed N *cores* instead of N neighbouring mask entries.  Falls back to one
    CPU per group where the kernel does not expose the topology.
    """
    groups: list[list[int]] = []
    index: dict[int, int] = {}
    for cpu in cpus:
        try:
            raw = Path(
                f"/sys/devices/system/cpu/cpu{cpu}/topology/thread_siblings_list"
            ).read_text()
        except OSError:
            return [[c] for c in cpus]
        siblings: set[int] = set()
        for part in raw.strip().split(","):
            part = part.strip()
            if not part:
                continue
            if "-" in part:
                lo, hi = part.split("-")
                siblings.update(range(int(lo), int(hi) + 1))
            else:
                siblings.add(int(part))
        key = min(siblings) if siblings else cpu
        if key not in index:
            index[key] = len(groups)
            groups.append([cpu])
        else:
            groups[index[key]].append(cpu)
    return groups


def taskset_prefix(worker: int, cpus_per_job: int, jobs: int) -> list[str]:
    """Give one worker process its own slice of the allocation.

    XLA sizes its CPU thread pool from ``sched_getaffinity``, so restricting a
    process to N CPUs is what makes it use N threads.  The node layout is
    therefore a choice of processes x threads: ``--jobs 1`` with no restriction
    gives one process the whole node's threads, and ``--jobs 4
    --cpus-per-job 32`` gives four processes 32 threads each.  Slurm hands out
    an arbitrary slice of the node, so the slices are cut from the real mask
    rather than from cores numbered 0..n-1.

    The CPUs are handed out round-robin over *cores* -- the first hardware thread
    of every core, then the second of every core -- and each worker takes the
    next ``cpus_per_job`` of them.  This used to cut contiguous runs of the mask
    instead, which happens to be right only when consecutive entries belong to
    different cores.  A probe on the node says they do: the mask is thread-major,
    with the siblings of core 0 being entries 0 and 128.  So on this node type
    contiguous cutting did give each worker the cores it asked for, and the real
    cost of the change is that it no longer depends on that holding.  Number the
    mask core-major instead, where consecutive entries *are* siblings, and every
    layout halves: ``--cpus-per-job 1`` puts two workers on each of 64 cores, and
    ``--cpus-per-job 2`` gives a worker both threads of one core when it was
    promised two cores.

    Round-robin keeps each worker spread over distinct cores under either
    numbering, and never leaves a core idle while another is doubled up.  Where
    there are more CPUs than cores the cores are shared deliberately -- that is
    what using all 256 hardware threads means.

    The sibling pairing comes from core_groups, so nothing here assumes how Slurm
    numbered the mask.
    """
    if cpus_per_job <= 0 or shutil.which("taskset") is None:
        return []
    groups = core_groups(allocated_cpus())
    if not groups:
        return []
    depth = max(len(g) for g in groups)
    order = [g[i] for i in range(depth) for g in groups if len(g) > i]
    n = min(cpus_per_job, len(order))
    start = (worker * n) % len(order)
    picked = [order[(start + i) % len(order)] for i in range(n)]
    return ["taskset", "-c", ",".join(str(c) for c in picked)]


# The Fortran original seeds myran from unixtime at 1 s resolution and has no way
# to be told a seed, so two runs of one configuration that start in the same
# second produce byte-identical output.  That is not hypothetical: with every run
# in flight at once, 67 of 96 configurations came back with both seeds identical
# (job 59269588).  A second seed that is really the first seed also throws away
# half of every configuration's statistics, so it has to be prevented rather than
# detected -- check_collisions can only report it, and re-running the losers
# sequentially costs ~40 minutes each.
#
# The port is immune (it takes --seed from the manifest), so this only ever
# delays a Fortran run, and only the second and later runs of any one
# configuration.  The gap is comfortably more than the 1 s the clock resolves.
SEED_STAGGER_SECONDS = 2.0
_stagger_lock = threading.Lock()
_stagger_last: dict[str, float] = {}

# Events between progress lines in a port run's out.txt.  See build_job.
PROGRESS_EVERY = 2000

# Run id the compilation warm-up writes to.  It is not a manifest run id, so it
# can never be mistaken for one; see warmup_cache.
WARMUP_RUN_ID = "_warmup"


def stagger_seed(code: str, cfg_id: str, gap: float = SEED_STAGGER_SECONDS) -> float:
    """Block until this run of ``cfg_id`` may start without sharing a unixtime
    second with a sibling run of the same configuration.  Returns the delay
    actually taken, so a caller (or a test) can see what happened."""
    if code != "fortran":
        return 0.0
    with _stagger_lock:
        previous = _stagger_last.get(cfg_id)
        claimed = time.time() if previous is None else max(time.time(), previous + gap)
        _stagger_last[cfg_id] = claimed
    delay = claimed - time.time()
    if delay > 0:
        time.sleep(delay)
    return delay


def build_job(
    run: dict,
    code: str,
    root: Path,
    repo: Path,
    python: str,
    cpus_per_job: int,
    jobs: int,
    worker: int,
    gpus: int,
) -> tuple[list[str], dict, Path, Path | None]:
    """Assemble (command, env, run_dir, stdin_path) for one run."""
    run_dir = root / code / run["run_id"]
    run_dir.mkdir(parents=True, exist_ok=True)

    card = run_dir / "run_card.txt"
    if not card.is_file():
        shutil.copy(root / "grid" / f"{run['cfg_id']}.txt", card)

    env = dict(os.environ)

    if code == "fortran":
        binary = repo / "bin" / "aao_rad_lund"
        if not binary.is_file():
            sys.exit(f"fortran binary missing: {binary} (build the repo first)")
        link = run_dir / "spp_tbl"
        if not link.exists():
            link.symlink_to(repo / "parms" / "spp_tbl")
        env["AAO_TRIALS"] = "0"
        # The original is serial Fortran with no OpenMP directives, so a core is
        # all it can use; one process per core is the only way to fill the node.
        cmd = taskset_prefix(worker, cpus_per_job or 1, jobs) + [str(binary)]
        return cmd, env, run_dir, card

    env["PYTHONPATH"] = os.pathsep.join(
        [str(repo / "src"), env["PYTHONPATH"]] if env.get("PYTHONPATH") else [str(repo / "src")]
    )
    # One XLA compilation cache shared by every port run in this scan root.  XLA
    # keeps no compiled code between processes, so without this each of the 192
    # runs compiles its own copy of the same two modules.  It lives under the scan
    # root rather than a run directory because the warm-up that populates it and
    # the runs that then reuse it do not have to be on the same node.
    env.setdefault("AAO_JAX_CACHE_DIR", str(root / "jax_cache"))
    # A worker given a slice of the allocation is told the same width in thread
    # counts that its taskset mask gives it in cores.  XLA sizes its pool from the
    # affinity mask, so taskset is what really bounds the port's own loops; the
    # thread-count variables matter for the OpenMP-backed pieces (BLAS under
    # numpy/scipy, used to read and interpolate the tables), which would otherwise
    # serialise, and for anything reading them.  All of them are set, not just
    # OMP_NUM_THREADS: run_py_cpu.sbatch exports JAX_NUM_THREADS for the whole node
    # (128), and a worker pinned to 2 CPUs that inherits it is being told to use
    # more threads than it has cores.  An unpinned worker inherits the sbatch
    # script's values untouched.
    if cpus_per_job:
        env["OMP_NUM_THREADS"] = str(cpus_per_job)
        env["JAX_NUM_THREADS"] = str(cpus_per_job)
        env["OPENBLAS_NUM_THREADS"] = str(cpus_per_job)
        env["MKL_NUM_THREADS"] = str(cpus_per_job)
    else:
        env.setdefault("OMP_NUM_THREADS", "1")
    if gpus:
        vis = os.environ.get("CUDA_VISIBLE_DEVICES", "")
        devices = vis.split(",") if vis else [str(i) for i in range(gpus)]
        env["CUDA_VISIBLE_DEVICES"] = devices[worker % len(devices)]
    else:
        env["JAX_PLATFORMS"] = "cpu"

    cmd = taskset_prefix(worker, cpus_per_job, jobs) + [
        python,
        "-m",
        "aao_rad.cli",
        "--legacy-input",
        str(card),
        "--parms",
        str(repo / "parms"),
        "--n-events",
        run["n_events"],
        "--seed",
        run["seed"],
        "--ek-sampling",
        "fortran",
        # Without this a run that has stopped making progress is indistinguishable
        # from one that is merely slow: out.txt stays empty from the pilot all the
        # way through the sampling loop.  One line per 2000 events is enough to
        # tell the two apart while a node is saturated.
        "--progress",
        str(PROGRESS_EVERY),
        "--format",
        "npz",
        "-o",
        str(run_dir / "out.npz"),
    ]
    return cmd, env, run_dir, None


def execute(job: tuple[list[str], dict, Path, Path | None]) -> tuple[str, int, float]:
    cmd, env, run_dir, stdin_path = job
    t0 = time.time()
    stdin_fh = open(stdin_path, "rb") if stdin_path else subprocess.DEVNULL
    try:
        # cwd is the run's own directory, and that is not cosmetic: the Fortran
        # opens the MAID tables relative to the CWD and writes its n-tuple to the
        # compile-time constant path "aao_rad.ntuple", so without a private CWD
        # 128 parallel runs would read the tables fine and then overwrite each
        # other's n-tuple.
        with open(run_dir / "out.txt", "wb") as out:
            proc = subprocess.run(
                cmd, env=env, stdin=stdin_fh, stdout=out, stderr=subprocess.STDOUT, cwd=run_dir
            )
    finally:
        if stdin_path:
            stdin_fh.close()
    return run_dir.name, proc.returncode, time.time() - t0


def md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def warmup_cache(args, runs: list[dict], root: Path, repo: Path, gpus: int) -> int:
    """Compile the jitted modules once, serially, before the fan-out.

    Every port run needs the same compiled modules.  Launched together, each one
    compiles its own, and on a saturated node that compile is what runs away:
    186 copies of ``jit_sample`` took 11 minutes apiece before contending with
    each other.  Compiling once up front costs a few seconds and leaves the
    fan-out with nothing to compile -- the results land in the shared
    ``AAO_JAX_CACHE_DIR`` and every later process loads them.

    The warm-up asks for very few events.  The compiled program does not depend
    on the event count, but the output does, so the warm-up runs under its own
    id and never writes into a real run's directory -- a short run left in a real
    directory would satisfy run_complete() and silently stand in for the
    full-length one.  The Fortran binary is serial C-like code with nothing to
    compile, so this is skipped for it.
    """
    if not args.warmup_events or args.code == "fortran" or not runs:
        return 0

    warm = dict(runs[0])
    warm["run_id"] = WARMUP_RUN_ID
    warm["n_events"] = str(args.warmup_events)

    print(f"warming the XLA cache: one {warm['n_events']}-event run, one process")
    run_dir = root / args.code / WARMUP_RUN_ID
    if run_dir.is_dir():
        shutil.rmtree(run_dir, ignore_errors=True)
    cmd, env, run_dir, _ = build_job(
        warm, args.code, root, repo, args.python, args.cpus_per_job, 1, 0, gpus
    )
    t0 = time.time()
    proc = subprocess.run(
        cmd, env=env, stdin=subprocess.DEVNULL, cwd=run_dir,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
    )
    dt = time.time() - t0
    if proc.returncode != 0:
        # Not fatal: the fan-out will just compile the modules itself.
        print(f"  warm-up exited {proc.returncode} after {dt:.0f}s, continuing anyway")
        print("  " + (proc.stdout or "").strip().replace("\n", "\n  ")[-800:])
    else:
        print(f"  warm-up ok in {dt:.0f}s, cache at {env['AAO_JAX_CACHE_DIR']}")
    return 0


def run_all(args, runs: list[dict], root: Path, repo: Path) -> int:
    gpus = args.gpus
    workers = args.jobs or (gpus if gpus else os.cpu_count() or 1)
    pending = [r for r in runs if not run_complete(args.code, root / args.code / r["run_id"])]
    if args.limit:
        # Smoke-test aid: run only the first N pending runs, so a short debug job
        # exercises the real launch path instead of the whole 192-run grid.
        pending = pending[: args.limit]
    print(
        f"{args.code}: {len(runs) - len(pending)}/{len(runs)} complete already, "
        f"{len(pending)} to run, {workers} workers"
    )

    if args.dry_run:
        for run in pending[:3]:
            job = build_job(
                run, args.code, root, repo, args.python, args.cpus_per_job, workers, 0, gpus
            )
            print("would run:", " ".join(job[0]))
        return 0

    warmup_cache(args, pending, root, repo, gpus)

    failures: list[str] = []
    t0 = time.time()

    # build_job reads a run's slice of the node out of `worker`: which core
    # taskset pins it to, and which GPU it is handed.  Passing a constant there
    # gave every concurrent run the same slice -- `taskset -c <one CPU>` and
    # `CUDA_VISIBLE_DEVICES=<device 0>` -- so 128 processes time-sliced a single
    # core and three of the four GPUs sat idle.  That is what job 59339555 was:
    # 30 minutes of wall shared 128 ways is 14.1 CPU-seconds per run, it measured
    # 13.6, the node was 99% idle, and every thread sat in futex_wait_queue.
    #
    # Slots come off a free list rather than a counter, so "no two runs in flight
    # share a slice" holds by construction.  A counter would nearly do it -- the
    # mapping round-robins -- but a slow run holding index 5 would collide with
    # index 133 once the other workers had churned that far, and the whole point
    # is that this failure is invisible until it is very expensive.
    slots = list(range(workers))
    slots_lock = threading.Lock()

    def launch(run: dict) -> tuple[str, int, float]:
        # Before build_job, so the wait is not charged to the run's own timing.
        stagger_seed(args.code, run["cfg_id"])
        with slots_lock:
            worker = slots.pop()
        try:
            return execute(
                build_job(
                    run,
                    args.code,
                    root,
                    repo,
                    args.python,
                    args.cpus_per_job,
                    workers,
                    worker,
                    gpus,
                )
            )
        finally:
            with slots_lock:
                slots.append(worker)

    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(launch, run): run["run_id"] for run in pending}
        done = 0
        for fut in as_completed(futures):
            name, rc, dt = fut.result()
            done += 1
            if rc != 0:
                failures.append(name)
            if done % 20 == 0 or rc != 0:
                print(
                    f"  [{done}/{len(pending)}] {name}: {'ok' if rc == 0 else f'EXIT {rc}'} ({dt:.0f}s)"
                )

    print(
        f"{args.code}: {len(pending) - len(failures)} finished in {time.time() - t0:.0f}s, "
        f"{len(failures)} failures"
    )
    for name in failures:
        print(f"  FAILED: {name}")
    return 1 if failures else 0


def check_collisions(runs: list[dict], root: Path, repo: Path) -> int:
    """md5-compare the runs of each configuration; re-run collided pairs one at
    a time so their unixtime seeds differ."""
    pairs: dict[str, list[dict]] = {}
    for run in runs:
        pairs.setdefault(run["cfg_id"], []).append(run)

    collided = []
    for cfg_id, group in sorted(pairs.items()):
        checksums = []
        for run in group:
            ntp = root / "fortran" / run["run_id"] / "aao_rad.ntuple"
            if ntp.is_file():
                checksums.append((run["run_id"], md5(ntp)))
        if len(checksums) >= 2 and len({h for _, h in checksums}) < len(checksums):
            collided.append(cfg_id)
            print(f"seed collision in {cfg_id}: " + ", ".join(f"{n}={h[:8]}" for n, h in checksums))

    if not collided:
        print("no seed collisions")
        return 0

    print(f"re-running {len(collided)} collided configurations sequentially")
    for cfg_id in collided:
        for run in pairs[cfg_id]:
            run_dir = root / "fortran" / run["run_id"]
            (run_dir / "aao_rad.ntuple").unlink(missing_ok=True)
            (run_dir / "out.txt").unlink(missing_ok=True)
        for run in pairs[cfg_id]:
            job = build_job(run, "fortran", root, repo, sys.executable, 0, 1, 0, 0)
            _, rc, dt = execute(job)
            print(f"  {run['run_id']}: exit {rc} ({dt:.0f}s)")
            time.sleep(2.0)
    return 0


def main() -> int:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--code", required=True, choices=tuple(SIGMA_MARKERS))
    p.add_argument("--root", type=Path, default=None)
    p.add_argument(
        "--repo",
        type=Path,
        default=None,
        help="repo root holding src/, bin/ and parms/ (default: two levels up)",
    )
    p.add_argument(
        "--python", default=sys.executable, help="python for the port runs (a venv with jax)"
    )
    p.add_argument(
        "--jobs",
        type=int,
        default=0,
        help="parallel workers (default: cpu count, or --gpus for py_gpu)",
    )
    p.add_argument("--gpus", type=int, default=0, help="for --code py_gpu: one worker per GPU")
    p.add_argument(
        "--cpus-per-job",
        type=int,
        default=0,
        help="CPUs (and so XLA threads) per worker process, cut from the "
        "allocation's real affinity mask; 0 leaves the process the whole "
        "allocation, so --jobs 1 means one process using every thread. The "
        "serial Fortran ignores this and gets one core per process.",
    )
    p.add_argument(
        "--check-collisions",
        action="store_true",
        help="fortran: md5-compare same-configuration pairs, re-run collisions",
    )
    p.add_argument(
        "--warmup-events",
        type=int,
        default=2000,
        help="for the port arms: compile the jitted modules once, serially, with a "
        "throwaway run of this many events before the fan-out, so 192 processes "
        "do not each compile their own copy. 0 skips it (default: %(default)s)",
    )
    p.add_argument(
        "--limit",
        type=int,
        default=0,
        help="run at most this many of the pending runs, 0 for all of them; a "
        "smoke-test aid so a short debug job exercises the real launch path",
    )
    p.add_argument("--dry-run", action="store_true")
    args = p.parse_args()

    root = args.root or default_root()
    repo = args.repo or Path(__file__).resolve().parent.parent.parent
    runs = load_manifest(root)

    if args.check_collisions:
        return check_collisions(runs, root, repo)
    return run_all(args, runs, root, repo)


if __name__ == "__main__":
    raise SystemExit(main())
