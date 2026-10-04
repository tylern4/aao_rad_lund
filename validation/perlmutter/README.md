# Perlmutter scan: Fortran vs JAX-CPU vs JAX-GPU

Batch scripts for NERSC Perlmutter that sample the same 96-configuration grid
with all three implementations, then verify the statistics at scale:

* **Fortran** — the original `aao_rad` (gfortran, this repo's Makefile)
* **JAX on CPU** — the port, `JAX_PLATFORMS=cpu`
* **JAX on GPU** — the port, one worker per A100

Everything runs on `$SCRATCH` (`/tmp` is not shared between Perlmutter's
login and compute nodes and is routinely scrubbed).

## The grid

`make_grid.py` writes 96 run cards: beam energies **2, 4.244, 6, 8, 10 and
12 GeV** crossed with

| axis | values |
|---|---|
| channel | pi+ n (`epirea=3`), pi0 p (`epirea=1`) |
| beam polarisation | on, off |
| explicit photon cut (`delta`) | 0.005, 0.05 GeV |
| target thickness | 5.0, 2.5 cm |

with 20 000 events per run and two runs per configuration (192 run jobs per
implementation).  Q^2 stays in [0.2, 1.9] GeV^2 — below the 90-degree elastic
clamp down to 2 GeV — and the E' window is derived per beam energy so the
whole (Q^2, E') rectangle keeps W in [1.20, 1.98] GeV, inside the MAID07
table and below both kinematic clamps `aao_rad.f90` applies.  Neither code
therefore silently samples a clipped window, which would make the cross-code
comparison meaningless.  The 4.244 GeV point anchors the grid to the
validated local configuration.

## Usage

On a Perlmutter login node, from a clone of this repo:

```bash
validation/perlmutter/submit_all.sh
```

That submits everything to project `m3792` (override with `AAO_ACCOUNT`). It
creates the python venv with `jax[cuda12]` (one-time, on the login node
— compute nodes have no internet), builds the Fortran binary, writes the
grid, and submits four jobs:

| job | nodes | what |
|---|---|---|
| `run_fortran.sbatch` | 1 PM-CPU, exclusive | 192 runs, one process per hardware thread, all in flight, then an md5 seed-collision check |
| `run_py_cpu.sbatch` | 1 PM-CPU, exclusive | 192 runs from `nproc` processes × 1 thread, all in flight |
| `run_py_gpu.sbatch` | 1 PM-GPU, exclusive | 192 runs from 4 processes, one per A100, 32 host threads each |
| `run_verify.sbatch` | 1 PM-CPU | verification, after the three above (`--dependency=afterok`) |

None of the three asks Slurm for a node size. `--cpus-per-task=128` looks like
a whole node and is not: a probe with no directives at all (`debug` 59300339)
came back with **256 CPUs and 487 GB** on a PM-CPU node, so the 128 an explicit
request bought was 64 cores' worth. The old `128 processes × 2 threads` did not
fit in that either — 256 slots on 64 cores, so every core ran two processes and
four threads. Asking for nothing gets the whole node, `--exclusive` keeps a
co-tenant off it, and the scripts read `nproc` at run time rather than assuming
how many cores that is. On a PM-GPU node the same probe (59300338) returned
128 CPUs, 229 GB and `CUDA_VISIBLE_DEVICES=0,1,2,3` — which is also the only
request that satisfies the `gpu_shared_ss11` QoS's 32-cores-per-GPU minimum for
4 GPUs, and the 112 GB of the previous 2-GPU attempt is what its 2 processes
exhausted (job 59269609, `OUT_OF_MEMORY`, 95 of 192 runs killed).

Individual arms can be (re)submitted on their own, and every batch is resumable
— completed runs are recognised by their final cross-section line and skipped.

Through the `nersc` CLI (`nersc submit <remote path>`), which takes only a
remote path and no other `sbatch` flags, each script carries its own
`#SBATCH --account` and `#SBATCH --output` directives:

```bash
nersc submit $SCRATCH/aao_rad_lund/validation/perlmutter/run_py_gpu.sbatch
nersc jobs --user tylern --command squeue
nersc job <jobid>
nersc cancel <jobid>
```

Slurm's log then lands as `slurm_<jobid>.out` next to the script, whereas
`submit_all.sh` passes `-o` and puts every arm's output in
`$SCRATCH/aao_rad_scan/logs/`.

### Processes and threads

Each run is an independent 20k-event batch, and `generate.stream` sizes its
vectorised chunk as `max(min(chunk_events, n_events), batch_size)` — so a run
costs one fixed 65,536-trial pass whatever it keeps.  Runs are therefore the
unit of parallelism, and the only question a node layout has to answer is how
many to have in flight.  `--jobs` sets the processes and `--cpus-per-job` cuts
each one a slice of the allocation (XLA sizes its thread pool from
`sched_getaffinity`, so the slice is what sets the thread count):

| arm | layout | reason |
|---|---|---|
| fortran | `nproc` processes × 1 thread | the original is serial Fortran with no OpenMP, so a thread is all a run can use |
| py_cpu | `nproc` processes × 1 thread | a second thread buys 1.2× of work for 2× the cores, and there are more runs than cores |
| py_gpu | 4 processes × `nproc`/4 threads | one per A100, splitting the node's affinity entries evenly; here the threads feed a device |

`calibrate.sbatch` measures the CPU choice rather than assuming it: it times one
real run at 1, 2, 4 … 128 threads and projects the wall time of every candidate
layout.  Run it before the scan and set `run_py_cpu.sbatch`'s `--jobs` /
`--cpus-per-job` to the winner:

```bash
nersc submit $SCRATCH/aao_rad_lund/validation/perlmutter/calibrate.sbatch
```

It times 300-event runs and scales them to the grid's 20 000, because the
1-thread point is the whole point of the measurement and a full-size run there
takes hours — long past the 30 minutes the `debug` queue allows.  Each label is
appended to `$SCRATCH/aao_rad_scan/calib_results.txt` as it finishes and the
projection reads that file, so a job killed by its wall clock still leaves usable
numbers.  A label that exceeds `AAO_CALIB_LABEL_TIMEOUT` (600 s) is recorded as
`TIMEOUT` and dropped rather than guessed at.

The first table below was measured on `nid004960` with `--cpus-per-task=128`,
i.e. half a node, and it is kept here only to be corrected — **its thread counts
are not cores.** The mask carries one entry per hardware thread, and the script
used to pick a slice by striding it, which on a hyperthreaded node reaches a
core's second thread before it reaches a second core: checked on the node,
`sorted(sched_getaffinity(0))[::128]` is entries 0 and 128, which are both
threads of core 0. So "2 threads" was one core's two hyperthreads and "4 threads"
was two cores. `--cpus-per-job N` never does that — it round-robins over cores,
so N threads is N different cores — which means the measurement and the layout
chosen from it were different quantities, and the `128 processes × 2 threads`
recommendation that came out of it was 256 slots on 128 cores, twice
oversubscribed. `one_process` now takes its slice from `run_batch.core_groups`,
so the two are the same thing by construction.

| threads asked for | 1 | 2 | 4 | 8 | 16 | 32 | 64 | 128 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 300-event run (s) | 35.0 | 28.8 | 20.2 | 13.9 | 10.6 | 8.9 | 8.4 | 8.7 |
| cores those threads actually spanned | 1 | 1 | 2 | 4 | 8 | 16 | 32 | 64 |
| speedup vs 1 thread | 1.0× | 1.2× | 1.7× | 2.5× | 3.3× | 3.9× | 4.2× | 4.0× |

The last row is the one that matters for a layout: read against the row above
it, a run's cost in cores was 35.0, 28.8, 40.4, 55.6, 84.8, 142, 269, 557
core-seconds — the second *hyperthread* is nearly free (35.0 → 28.8 for the same
core) while every *extra core* costs more than it returns. Widen a run only by
using both threads of a core it already has, never by taking a core from another
run.

Both CPU arms therefore give every run one thread and keep every run in flight:
the grid's 192 runs are independent, so there are always more of them than the
node has cores, and a wider run can only take capacity away from a run that has
none. `run_fortran.sbatch` cannot use a second thread at all (serial Fortran,
no OpenMP). `run_py_cpu.sbatch` defaults to `AAO_CPU_JOBS=$(nproc)
AAO_CPU_CPUS=1`; a changed node or grid wants those re-derived rather than
carried over. Re-run `calibrate.sbatch` for the corrected table.

## Layout on `$SCRATCH`

```
$SCRATCH/aao_rad_scan/
  grid/           cfg_NNN.txt run cards + manifest.csv (one row per run job)
  fortran/cfg_NNN_sK/   aao_rad.ntuple (32+15 columns), out.txt, aao_rad.lund
  py_cpu/cfg_NNN_sK/    out.npz (lossless event n-tuple), out.txt
  py_gpu/cfg_NNN_sK/    out.npz, out.txt
  jax_cache/      persistent JAX compilation cache
  logs/           Slurm job output
  verify/         verification CSVs + summary.txt
```

Each run gets its own directory, and the run's CWD *is* that directory — not
cosmetic: the Fortran opens the MAID tables relative to the CWD
(`maid_lee.f90`) and writes its n-tuple to the compile-time constant path
`aao_rad.ntuple`, so without a private CWD parallel runs would read the tables
fine and then overwrite each other's n-tuple.  `AAO_TRIALS=0` keeps the 0.5 GB
trial dump off while the 32-variable event n-tuple is still written.

## Verification

`verify_statistics.py` (also runnable standalone: `python verify_statistics.py
--root $SCRATCH/aao_rad_scan`) checks two levels:

**Like-for-like** — the two runs of one configuration within one
implementation.  The run-to-run cross-section scatter, pooled across all 96
configurations into one fractional sigma per implementation (the heavy photon
tail makes the naive `1/sqrt(n_trials)` error optimistic, so the scatter is
measured, not assumed), plus KS distances per observable against the 95%
noise floor `1.36*sqrt(2/N)`.

**Cross-code** — Fortran vs JAX-CPU, Fortran vs JAX-GPU and JAX-CPU vs
JAX-GPU per configuration: the cross-section ratio with an error bar from
each implementation's pooled scatter (deviation in sigmas), and KS per
observable against `1.36*sqrt(1/N_A + 1/N_B)`.

CPU and GPU share seeds per configuration, so any CPU-vs-GPU difference is
pure backend floating-point, not sampling; the Fortran seeds `myran` from
`unixtime` (1 s resolution) and its two runs are simply independent samples.
Because the two runs of one configuration are ordered far apart in the work
list they cannot start in the same second, and `--check-collisions` md5-
compares every same-configuration pair afterwards as a guard.

Both sigma estimators are the same one the sigma tables quote — mean trial
weight over all trials times phase space (the Fortran's `sig_sum`, the
port's `sigma (MC)`) — in micro-barns on both sides.

Outputs in `verify/`:

| file | contents |
|---|---|
| `sigmas.csv` | one row per run: cross section, trial count |
| `like_for_like.csv` | per implementation x configuration: scatter z, worst KS |
| `cross_code.csv` | per pair x configuration: sigma ratio, z, worst KS |
| `ks_full.csv` | every configuration x observable x level KS with its floor |
| `per_observable.csv` | per pair x observable: median/max KS, fraction within floor |
| `summary.txt` | the overall agreement verdict |

## Notes

* The port runs with `--ek-sampling fortran`, matching the original's
  RNG-limited photon sampling — the same setting the local validation uses.
* `--n-events` and `--seed` are passed from the manifest: the CLI flags
  override the card, and the card carries no seed.
* One driver (`run_batch.py`) serves all three arms so batching,
  resumability and bookkeeping stay identical between them.
* Wall-time estimates: Fortran ~2-3 h (128 single-core processes, two waves),
  JAX-CPU and JAX-GPU depend on the thread scaling `calibrate.sbatch` measures,
  verify ~20 min; the `#SBATCH --time` headers have headroom.
