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

All three run in the `regular` QoS — see below for why not `shared`.
| `run_verify.sbatch` | 1 PM-CPU | verification, after the three above (`--dependency=afterok`) |

None of the three asks Slurm for a node size. `--cpus-per-task=128` looks like
a whole node and is not: a probe with no directives at all (`debug` 59300339)
came back with **256 CPUs and 487 GB** on a PM-CPU node, so the 128 an explicit
request bought was 64 cores' worth. The old `128 processes × 2 threads` did not
fit in that either — 256 slots on 64 cores, so every core ran two processes and
four threads. Asking for nothing gets the whole node, `--exclusive` keeps a
co-tenant off it, and the scripts read `nproc` at run time rather than assuming
how many cores that is. On a PM-GPU node the same probe (59300338) returned
128 CPUs, 229 GB and `CUDA_VISIBLE_DEVICES=0,1,2,3`; the 112 GB of the previous
2-GPU attempt is what its 2 processes exhausted (job 59269609, `OUT_OF_MEMORY`,
95 of 192 runs killed).

All four jobs use the **`regular`** QoS, not `shared`. Both allow 2 days, but a
GPU job in `shared` is capped at `gres/gpu=2,node=1`, so it can never hold the
node's four A100s, and it additionally demands exactly 32 cores per GPU — which
is why `--exclusive`, meaning a whole node and so 128 cores, is refused outright
under `shared` (`Requested node configuration is not available`). `regular` has
no per-job TRES cap:

```
regular   2-00:00:00   (no cap)
shared    2-00:00:00   node=1
gpu_shared 2-00:00:00  gres/gpu=2,node=1
debug     00:30:00     node=8
gpu_debug 00:30:00     node=8
```

`debug` also grants a whole 4-GPU node, but 30 minutes is far too short for a
192-run arm, which is why it is kept for `calibrate.sbatch`.

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

The corrected measurement, whole node (256 hardware threads, 128 cores), each
thread count on a **distinct core** — steady state of two 300-event runs, so the
compile is not in the number:

| threads (= cores) | 1 | 2 | 4 | 8 | 16 | 32 | 64 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 300-event run (s) | 35.1 | 19.9 | 12.9 | 9.7 | 8.2 | 7.2 | 7.1 |
| speedup vs 1 thread | 1.0× | 1.8× | 2.7× | 3.6× | 4.3× | 4.9× | 4.9× |
| core-seconds per run | 35.1 | 39.7 | 51.4 | 77.3 | 131 | 232 | 456 |

The last row decides the layout. Speedup saturates near 5×, but the *cost* does
not: a second core costs more than the 1.8× it returns, so widening a run is
strictly worse than giving that core to another run. (The second *hyperthread* of
a core is a different matter and nearly free — 35.0 s → 28.8 s for the same core
in the flawed table below — but `taskset_prefix` never hands a run the sibling of
a core another run already has, so it is not on offer.)

Both CPU arms therefore give every run one thread and keep every run in flight:
the grid's 192 runs are independent, so there are always more of them than the
node's 128 cores, and a wider run can only take capacity away from a run that has
none. `run_fortran.sbatch` cannot use a second thread at all (serial Fortran, no
OpenMP). The projection now models wave quantisation *and* slot contention, and
agrees:

```
 procs  threads  waves  slots  projected
     1 x 64        192     64    1521.1m
     2 x 64         96    128     760.5m
     4 x 32         48    128     386.1m
     8 x 16         24    128     219.2m
    16 x 8          12    128     128.8m
    32 x 4           6    128      85.7m
    64 x 2           3    128      66.2m
   128 x 1           2    128      78.0m
   192 x 1           1    192      58.5m
```

`run_py_cpu.sbatch` defaults to `AAO_CPU_JOBS=$(nproc) AAO_CPU_CPUS=1`; workers
beyond the run count are idle, so the count does not need pinning to 192. A
changed node or grid wants both re-derived rather than carried over.

### The superseded table

The first table was measured on `nid004960` with `--cpus-per-task=128`, i.e.
half a node, and **its thread counts are not cores.** The mask carries one entry
per hardware thread, and the script used to pick a slice by striding it, which on
a hyperthreaded node reaches a core's second thread before it reaches a second
core: checked on the node, `sorted(sched_getaffinity(0))[::128]` is entries 0
and 128, which are both threads of core 0. So "2 threads" was one core's two
hyperthreads and "4 threads" was two cores. `--cpus-per-task` never does that —
`taskset_prefix` round-robins over cores, so N threads is N different cores —
which means the measurement and the layout chosen from it were different
quantities, and the `128 processes × 2 threads` recommendation that came out of
it was 256 slots on 64 cores. `one_process` now takes its slice from
`run_batch.core_groups`, so the two are the same thing by construction.

| threads asked for | 1 | 2 | 4 | 8 | 16 | 32 | 64 | 128 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 300-event run (s) | 35.0 | 28.8 | 20.2 | 13.9 | 10.6 | 8.9 | 8.4 | 8.7 |
| cores those threads actually spanned | 1 | 1 | 2 | 4 | 8 | 16 | 32 | 64 |

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

### The reference fails on a quarter of this grid

`verify_statistics.py` drops any run that could not produce a usable answer,
because averaging those in silently is worse than dropping them: they reported a
Fortran-vs-GPU sigma ratio spanning **-10.4 to +5.4** with a |z| of 24.8 and
78% of configurations inside 2 sigma, which reads as a catastrophic port failure
and is entirely an artifact of the reference.

The cause is in the original, not the port.  `integer*4 ntries` and
`integer*4 ngeom` (`src/aao_rad.f90:154,158`) are 32-bit.  A configuration that
needs more than ~2.1e9 trials wraps the counter, and since `ntries` is a *divisor*
in the cross section the printed value goes **negative** rather than merely
becoming imprecise — e.g. `cfg_088_s1` reports `-0.150309518` where the port
reports `0.00687858`.  Such runs also tend to stop short of their event quota.

On this grid, **46 of 192 Fortran runs are unusable**: 29 overflowed counters, 22
non-positive cross sections, 18 short n-tuple records.  (A run can carry more
than one.)  That leaves 146 usable runs and **61 of 96 configurations with both
seeds valid**, so cross-code agreement is measured on 61 configurations, not 96.
The remaining 35 have no valid cross-code reference and are covered only by the
CPU-vs-GPU comparison.

`missm-2` is deliberately *not* treated as a defect.  It looks like one — a
numerical guard printing to stdout — but it is per-event, inside the generation
loop (`src/aao_rad.f90:1666-1668`): whenever `csthcm**2` rounds just past 1 it
clamps `snthcm` to 1e-7 and continues.  **167 of the 192** runs print it,
including runs with a full 20,000-record n-tuple and a healthy positive cross
section.  Counting it as a defect discarded 167 good runs and left 1 usable
configuration of 96.

The defect checks live in `verify_statistics.defects()` and are Fortran-only by
construction — the port has no 32-bit counters and writes its own n-tuple.
Record counting uses newlines rather than dividing the byte size by 800, so it
does not depend on `es16.8` continuing to print exactly 16 characters.

With those runs excluded, the Fortran's own run-to-run scatter is **0.137%**,
against 0.077% for JAX-CPU and 0.070% for JAX-GPU — the same order, which is the
point: the port's residual noise is the reference's, not extra noise of its own.

### Two measurement artifacts that had to be fixed first

Both were defects in the *verification*, not in the port, and both produced
alarming-looking numbers until they were identified.

**`missm-2` is not a defect.** See above — counting it discarded 167 good runs
and left 1 usable configuration of 96.

**KS must run at the precision the reference reports.** KS reported **0.9739** on
`E_s` for all 16 configurations at `ebeam = 4.244` GeV, which reads as the two
codes sampling disjoint physics.  It was float representation: `E_s` at that
beam energy is a point mass (97% of events at exactly 4.244), the Fortran writes
`es16.8` so it reports `4.24399996`, and the port's npz is float32 so it stores
`4.24399995803833`.  Two atomic masses 2e-9 apart put the empirical CDFs 97%
apart.  What identified it: the affected configurations were *exactly* the 16 at
4.244 GeV, while 2.0, 6.0 and 12.0 — which are exact in binary32 — showed no such
artifact.  Quantizing to 8 significant digits drops the worst KS from 0.9739 to
0.041.

### The remaining disagreement is real, and it is energy-dependent

With both artifacts removed, the cross-code cross section still does not agree
to within Monte Carlo noise.  Comparing the two-seed means, Fortran minus
JAX-GPU, on the 61 configurations with a valid reference:

| `ebeam` (GeV) | n | median difference |
|---|---|---|
| 2.0 | 8 | +0.487% |
| 4.244 | 13 | +0.114% |
| 6.0 | 12 | +0.065% |
| 8.0 | 8 | +0.197% |
| 10.0 | 10 | +0.620% |
| 12.0 | 10 | +0.734% |

Overall median +0.28%, mean +0.39%, sd 0.59%, `|d|` p90 1.19%, max 2.34% — against
a same-code run-to-run median of 0.07%.  So the offset at 10 and 12 GeV is roughly
ten times the noise: **this is a systematic, not sampling.**

It is *not* polarisation (median +0.279% polarised vs +0.282% unpolarised).  It
is smallest at 6 GeV and grows toward both ends of the energy range, and every one
of the eight largest disagreements sits at 10 or 12 GeV.  The `ep` acceptance
window scales with the beam energy (0.28–0.68 GeV at 2 GeV up to 10.28–10.68 GeV at
12 GeV), so the extremes probe the hadronic amplitude furthest outside the region
where it is best constrained — the MAID interpolation is the first place to look.
**This is unresolved.**  The port is not wrong by a bug-sized margin; it is off by
a few tenths of a percent in a way that tracks where in kinematics the
configuration sits.

The same offset shows up in the distributions: 151 of 976 fortran-vs-GPU
observable checks sit outside the KS noise floor (15.5%, against ~5% expected at
a 1.36-sigma threshold), with the worst at `mm^2` (KS 0.041).  By contrast
JAX-CPU vs JAX-GPU — the two port backends, sharing seeds, so pure floating point
— is 0 of 592 over floor with a cross-section ratio of 1.00049 and a range of
0.99933 to 1.00343.

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
