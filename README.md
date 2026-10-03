# aao_rad

Radiative pion electroproduction event generator for CLAS12, ported from the
Fortran `aao_rad` to JAX.

The physics is unchanged: MAID07 multipole amplitudes, the AO response
functions, the Mo & Tsai radiative kernel, and the same importance-sampling
scheme in `cos(theta_k)` and photon energy. What changed is how the integral is
evaluated. The original was a scalar accept/reject loop that drew one point,
interpolated the response tables, and moved on. This version evaluates a whole
batch of trial points at once on the device, turns the kinematic `go to` chain
into boolean masks, and replaces the hand-tuned `sigr_max` prompt with a
ceiling estimated from a pilot batch.

The original Fortran is still in `src/` and still builds; it is only needed for
validation. See [Building the Fortran reference](#building-the-fortran-reference).

---

## Contents

- [Install](#install)
- [Parameter files](#parameter-files)
- [Running it](#running-it)
- [Performance and tuning](#performance-and-tuning)
- [Output](#output)
- [Validation](#validation)
  - [The shared-key trap](#the-shared-key-trap)
  - [A second, output-only bug: `asym_p`](#a-second-output-only-bug-asym_p)
- [Differences from the Fortran original](#differences-from-the-fortran-original)
- [Bugs found in the Fortran original](#bugs-found-in-the-fortran-original)
- [Development](#development)
- [Building the Fortran reference](#building-the-fortran-reference)

---

## Install

```sh
pip install .            # or: uv pip install .
```

Requires Python >= 3.10, JAX >= 0.4.28 and NumPy. Extras cover the optional
pieces: `.[yaml]` for YAML run cards (JSON and TOML need no parser), `.[dev]`
for the tests, the linter and the validation plots.

JAX picks the backend itself. On a CUDA machine with `jax[cuda12]` installed it
runs on the GPU with no further configuration; otherwise it runs on the CPU.
Check what you got with:

```sh
python -c "import jax; print(jax.devices())"
```

## Parameter files

The MAID07 response tables ship as fixed-width Fortran `.tbl` text in
`parms/spp_tbl/`. They are parsed once and cached as an uncompressed `.npz` in
`parms/tables/` — same directory, so there is nothing to install separately:

```sh
tar -xvf parms.tar.gz          # if you do not already have parms/spp_tbl
aao-rad-table --list           # what is available
aao-rad-table -c 3             # grid, shape and first values
```

The `.npz` holds `q2`, `w`, `amps (62, nq2, nw)`, `source` and a
`format_version`. Amplitudes are float32 because the `.tbl` files carry ~5
significant digits and the Fortran declared `SF` as single precision, so float32
stores them exactly. To force a re-parse after replacing a `.tbl`:

```sh
aao-rad-table --rebuild-cache -c 3
```

The tables are found automatically. To override, set `AAO_RAD_PARMS` (or
`CLAS_PARMS`, which is what the Fortran used), or pass `--parms DIR`:

```sh
export AAO_RAD_PARMS=/path/to/parms
```

Only MAID07 (`theory=7`) is implemented. It is the only model the original
driver could actually reach — see `src/maid_lee.f90`.

## Running it

```sh
aao-rad --experiment rgb --n-events 1000000 --output events.lund
```

`--experiment` takes `rgb` or `default`, the two presets the old driver script
used. Everything is overridable:

```sh
aao-rad -n 5e6 --channel 1 --beam-energy 4.8 \
        --q2min 0.9 --q2max 2.5 --emin 2.66 --emax 3.29 \
        --mm-cut 0.2 --seed 7 --output pi0.lund
```

`aao-rad --dump-config` prints the fully resolved configuration, which is the
easiest way to see what a preset or a config file actually resolved to.

From Python:

```python
from aao_rad import GeneratorConfig, generate_lund

cfg = GeneratorConfig(
    channel=3,
    beam_energy=4.244,
    q2_min=0.2, q2_max=1.9,
    ep_min=1.6, ep_max=2.9,
    n_events=1_000_000,
    seed=7,
)
stats = generate_lund(cfg, "events.lund")
print(stats.summary())
```

`generate_lund` streams, so peak memory does not grow with `n_events`. For the
raw array:

```python
from aao_rad import EventGenerator, build_grid, GeneratorConfig

grid = build_grid(3)                       # cached after the first call
events, stats = EventGenerator(grid).generate(cfg)
events["q2"], events["w"], events["mm2"]   # 32-column structured array
```

Old run cards still work:

```sh
aao-rad --legacy-input test.inp --output test.lund
```

The original `stdin` prompt sequence is parsed by
`GeneratorConfig.from_legacy_input`, minus the two `sigr_max` prompts — the
acceptance ceiling is now estimated on the device.

### Configuration files

`--config FILE` takes JSON, TOML, or YAML. Any subset of `GeneratorConfig`'s
fields; unknown keys are an error rather than a silent no-op. A `--config` file
is overridden by command-line flags, and both override the preset.

## Performance and tuning

The cost is dominated by the response-function evaluation: one table
interpolation plus small sums over six partial waves, per trial point. That is
perfectly parallel, which is why the batched version exists.

```sh
aao-rad -n 1e7 --batch-size 1048576 --progress 1000000 -o big.lund
```

- `--batch-size` is the number of trial points per vectorised step. It must be
  at least 1024. Too small and the device is underfilled; too large and you pay
  the compile cost and hold a big buffer. Somewhere in `2^18`–`2^20` is the
  usual compromise; each distinct value compiles once and is then reused for
  every step. (Heuristic — no GPU was available while this was written, so only
  the CPU numbers below are measured.)
- `--chunk-events` bounds how many events are held in memory at once (default
  500 000). A billion-event run uses the same memory as a million-event one.
- `--progress N` logs every `N` events.

The end-of-run summary is worth reading:

```
sampling accept. : 0.035%
event yield      : 0.035%
weight max / mean: 0.0288 / 8.6e-06
above ceiling    : 3 trials, 0.02% of the cross section
sigma (MC)       : 0.0218 micro-barn
sigma (accepted) : 0.0216 micro-barn
```

`sampling accept.` is the importance-sampler efficiency — the fraction of trial
points that became an event. `event yield` is the same number over the trimmed
run, so the two should agree exactly; they do not if the last batch was cut
short, which is the only difference between them.

The two cross sections are independent estimators and should also agree:
`sigma_mc` sums the trial weight over *every* trial, while `sigma_accepted`
sums `accepted x ceiling`. They share no code, so a disagreement means the
acceptance logic is wrong.

**The acceptance ceiling.** Accepting `u * ceiling < weight` is
`min(1, weight/ceiling)`, which samples events as `weight` only while the
ceiling exceeds every weight in the run. The integrand has a heavy tail near the
beam line and no finite maximum, so the ceiling is estimated from a pilot batch,
scaled by `--weight-max-margin` (default 1.5), and raised again with a warning
if anything still exceeds it. Acceptance is `mean(weight) / ceiling`, so the
margin is a straight throughput dial: 2.0 accepts about a third fewer trials
than 1.5. `above ceiling` reports the residual as a fraction of the cross
section; it should be well below 1%.

Neither cross-section estimator depends on the ceiling, only the per-event
distributions do.

What is measured, on CPU: 20000 events took 570M trials in 328 s before the
ceiling was managed and 57M trials in 33 s after — a 10x throughput gain from
the same physics, with the two independent cross-section estimators agreeing to
0.8% and a residual ceiling bias of 0.02%.

## Output

**LUND** (`--format lund`, default) writes the CLAS12 generator format: a header
declaring 4 particles and three track lines per event, in the order
`px py pz E`, with the vertex in cm.

**NPZ** (`--format npz`) writes a lossless `.npz` holding the same 32 variables
as a structured array, plus the run metadata. Use this when you want to
re-analyse without re-parsing a text file.

The 32 columns, in order, are the original n-tuple:

| | | | |
|---|---|---|---|
| `es` | `ep` | `theta` | `w` |
| `w_real` | `ppx` `ppy` `ppz` | `eprot` | `ppix` `ppiy` `ppiz` |
| `epi` | `csthcm` | `phicm` | `mm2` |
| `eg` | `cstk` | `phik` | `qx` `qz` |
| `q0` | `csthe` | `egx` `egy` `egz` | `vx` `vy` `vz` |
| `q2` | `e_hel` | `asym_p` | |

`theta`, `phicm` and `phik` are in degrees; `csthcm`, `cstk` and `csthe` are
cosines. `es` and `ep` are *post*-radiation-loss, while `w`, `qx`, `qz`, `q0`
and `csthe` use the *pre*-exit electron, exactly as `aao_rad.f90` records them.

Note the original wrote the four-momentum in CLAS12's `E px py pz` order but
filled `px py pz E`; see
[Bugs found in the Fortran original](#bugs-found-in-the-fortran-original). This
port writes `px py pz E`, which is what the LUND standard and CLAS12's own
reader expect.

## Validation

`validation/` holds scripts that compare this port against the Fortran, plus the
plots. They are not part of the installed package.

| script | what it checks |
|---|---|
| `compare_table_to_fortran.py` | the parsed tables against the Fortran reader — **bit-for-bit** |
| `compare_sigma_stages.py` | every intermediate stage of `sigma()` against a Fortran driver |
| `compare_dsigma_to_fortran.py` | the amplitude and response stages of `dsigma()` |
| `compare_sigma_points.py` | end-to-end event weights, region by region |
| `compare_trials.py` | the *trial* measure: acceptance and mean weight, separately |
| `compare_distributions.py` | 16 observables, KS distances and histograms -> CSV + plots |
| `compare_mm2.py` | survivor `mm^2`/`W_real` and the `A(ek)` acceptance profile |
| `compare_sigr.py` | the raw cross section, recovered from the trial dump |
| `weight_breakdown.py` | the cross section split by branch and importance region |
| `compare_references.py` | two Fortran n-tuples from one run card — proves a change to the reference was surgical |

The Fortran side needs the instrumented `aao_rad.f90` (the version in this
branch), which dumps a 32-variable n-tuple plus the per-trial measure to units
13 and 14. Set `AAO_TRIALS=0` to skip the trial dump (it is ~0.5 GB per 8000
events and costs about as much wall time as the physics). See
[Building the Fortran reference](#building-the-fortran-reference).

Reproducing the headline number:

```sh
# 1. build and run the reference
cd build && cmake --build . --target aao_rad && cd ..
mkdir -p /tmp/frun_mm2 && cd /tmp/frun_mm2
ln -sf "$OLDPWD/parms" spp_tbl
../../build/aao_rad < "$OLDPWD/validation/maid07_pi+.inp"

# 2. compare the trial measures.  --tdump is the stride of the Fortran's trial
#    dump (tdump = 97 at aao_rad.f90:394, applied by the mod() guard at :967),
#    so P(reach) is n_dumped * tdump / ntries.
cd -
python validation/compare_trials.py --trials /tmp/frun_mm2/aao_rad.trials \
       --fortran-ntries 136336902 --fortran-sigma 1.54079217e-02 \
       --tdump 97 --seed 707 --n 16777216
```

which reports

```
mean(weight | reached the weight stage)
  fortran         : 5.25834763e-06
  port            : 5.21834409e-06   ratio 0.992392

P(reach the weight stage)
  fortran         : 0.507793   [713,981 / 136,336,902]
  port            : 0.508174
  ratio           : 1.000751

sigma = P(reached) * mean(weight | reached) * PS
  fortran         : 0.015407922 micro-barn
  port            : 0.015302186 micro-barn
  ratio           : 0.993138
  from P(reached) : +0.0751%
  from mean(w)    : -0.7608%
```

Splitting the cross section that way is the point of the exercise. The accepted
n-tuple cannot distinguish a sampling-density difference from a cross-section
difference, because accepted events are importance-sampled *by the weight*; the
trial dump has no importance sampling in it and can. The two factors also have
very different uncertainties, and conflating them is how a real bug hides:

| factor | port, 6 seeds x 16.7M trials | uncertainty |
|---|---|---|
| `P(reach)` | +0.075% | 0.02% — a binomial over 16.7M trials |
| `mean(weight)` | −0.76% | **0.44%** — see below |
| `sigma` | **−0.20% ± 0.18%** | dominated by `mean(weight)` |

`P(reach)` is measured to 0.02%. `mean(weight)` is not: the raw cross section
`sigma_r` spans ten decades, so the variance of its sample mean converges far
more slowly than `1/N` and the statistic is worth about 0.4% no matter how many
trials are thrown at it. Across six independent seeds the port gives
`sigma = 0.0152668 ± 0.0000067` micro-barn, i.e. the port-to-port scatter alone
is 0.44%. Quoting a single seed's `sigma` ratio to three decimal places is
therefore meaningless; the scan is in the history of this branch and the error
bar is the number to compare against.

Against the three reference runs available, the port sits at

| reference | events | `sigma` | port ratio |
|---|---|---|---|
| `frun` (pre-fix) | 8000 | 1.5297985e-2 | **0.998 (−0.20% ± 0.18%)** |
| `frun_fixed` (post-fix) | 8000 | 1.5298942e-2 | **0.998 (−0.21% ± 0.18%)** |
| `frun_mm2` | 2000 | 1.5407922e-2 | 0.991 (−0.92% ± 0.18%) |

The first two rows are the *same* run card and the same 8000-event target,
from builds that differ only by the [`asym_p` fix](#a-second-output-only-bug-asym_p),
and they differ by **0.006%**. `myran` seeds from `unixtime`, so these are
independent samples, not the same events twice — which makes the agreement
meaningful rather than tautological. It is the control the fix needs: `asym_p`
enters the n-tuple only, so a correct fix must leave the cross section untouched,
and it does. It also calibrates how much of any gap between two 8000-event
references is run-to-run noise rather than physics.

The widest gap between two Fortran references is 0.72% (`frun` against
`frun_mm2`), which is larger than the port's own scatter and comparable to the
discrepancy being measured. The
honest statement is therefore that **the port and the Fortran agree to within
the Fortran's own run-to-run scatter**; pinning the residual below that needs a
reference with enough events that its `sigma` is itself precise, which is the
main thing still missing from this validation.

### Where the distributions differ

Every KS distance has to be read against the two-sample 5% critical value,
`1.36 * sqrt(1/N_py + 1/N_f77)`. It is dominated by the *reference* sample, not
the port's, so a useful reference run needs at least as many events as the one
it is compared against — otherwise the floor is the thing being measured.

The trial-level comparison is the sharper instrument, because it conditions on
reaching the weight stage and so compares ~8.5M port trials against ~714k
Fortran rows — a noise floor of **0.0020**:

| variable | KS | mean ratio | mean dz |
|---|---|---|---|
| `es` | 0.00012 | 1.00000 | +0.49 |
| `q2` | 0.00131 | 0.99958 | −0.48 |
| `ep` | 0.00106 | 1.00014 | +0.69 |
| `ek` | 0.00054 | 0.99966 | −0.33 |
| `cstk` | 0.00056 | 0.99988 | −0.22 |
| `phik` | 0.00077 | — (mean is 0 on both sides) | −0.19 |

All six are below the floor. `phik` is a signed angle, so its mean is
statistically zero in both samples and a *ratio* of the two is the quotient of
two noise values; the script prints a z-score instead and suppresses the ratio.

The importance-region mix agrees to better than 0.5% in every region
(1.0034, 1.0015, 0.9974, 0.9979, 1.0006), which matters because the regions are
what `mcfac` compensates for: if the port sampled them in different proportions
the weights would be wrong by construction.

At the event level, 200k port events against the 8000-event reference gives a
noise floor of `1.36 * sqrt(1/200000 + 1/8000) = 0.0154`. This table was measured
against `frun`, the **pre-fix** reference, and is kept as-is because
`validation/compare_references.py` shows the fix left all 32 other event columns
unchanged between `frun` and `frun_fixed` — so every row except `asym_p` still
holds against the current build. Only the `asym_p` row is stale, and it is
marked:

| observable | KS D | | observable | KS D |
|---|---|---|---|---|
| `E_gamma` | 0.0083 | | `E'` | 0.0143 |
| `cos(theta_k)` | 0.0077 | | `E_pion` | 0.0151 |
| `phi_k` | 0.0081 | | `q0` | 0.0150 |
| `Q^2` | 0.0091 | | `cos(theta*)` | 0.0103 |
| `theta_e` | 0.0085 | | `phi*` | 0.0133 |
| `cos(theta_e)` | 0.0085 | | `W` | 0.0125 |
| `W_real` | 0.0129 | | `mm^2` | 0.0359 |
| | | | `asym_p` | ~~0.0254~~ *(pre-fix)* |
| | | | `E_s` | 0.9459 |

Thirteen of the sixteen sit at or below the floor. The three that do not are
`mm^2` at 2.3x, `asym_p` at 1.6x — the latter now explained and closed, see below
— and `E_s`, which is degenerate. None of the three can account for the
cross-section residual above: `asym_p` and `E_s` never enter the weight, and
`mm^2` is cut identically on both sides.

Per bin (40 bins each), which is the stricter test — a KS distance can hide a
compensating pair of errors that per-bin pulls cannot:

| observable | bins > 3σ | max abs z | ratio p05 / p50 / p95 |
|---|---|---|---|
| `phi*` | 0/40 | 2.33 | 0.900 / 0.996 / 1.156 |
| `cos(theta*)` | 0/40 | 2.02 | 0.889 / 1.009 / 1.176 |
| `Q^2` | 0/40 | 1.96 | 0.858 / 0.996 / 1.122 |
| `theta_e` | 0/40 | 1.86 | 0.849 / 0.999 / 1.072 |
| `cos(theta_k)` | 0/40 | 2.50 | 0.897 / 0.998 / 1.146 |
| `phi_k` | 0/40 | 2.69 | 0.893 / 0.993 / 1.121 |
| `W` | 0/40 | 2.26 | 0.827 / 1.002 / 1.484 |
| `E_gamma` | 1/40 | 3.08 | 0.296 / 1.025 / 2.672 |
| `mm^2` | 1/40 | 3.79 | 0.900 / 1.049 / 1.267 |
| `W_real` | 0/40 | 2.73 | 0.807 / 0.981 / 1.226 |
| `E_pion` | 0/40 | 2.76 | 0.848 / 1.000 / 1.281 |
| `asym_p` | ~~**5/40**~~ *(pre-fix)* | ~~**5.30**~~ *(pre-fix)* | ~~0.153 / **0.709** / 1.059~~ *(pre-fix)* |

`phi*` is the one to look at: it was 12 of 12 bins over 3σ with max abs z = 59
and a median ratio of 0.761 before the key fix, and is now indistinguishable.
That is the visible symptom of the bug described below.

`asym_p` was the last row still over 3σ, and it was the *reference* that was
wrong: `ntp(32)` held an unrelated trial's asymmetry for the whole radiative
branch. Those struck-through figures are the pre-fix symptom, kept for the
record. Note also that an aggregate column mixes the soft and radiative
branches, which is a weaker instrument than splitting them — so the numbers
that actually close this are the branch-split ones in
[the `asym_p` write-up](#a-second-output-only-bug-asym_p).

`ref-empty` bins in `validation/distributions_bins.csv` are a binning artifact
for the two concentrated observables, not a disagreement: `E_gamma` has a mean
of 0.014 GeV in a range plotted over 0.6, and `E_s` puts 97% of its events in a
`1e-4` GeV window, so uniform bins leave the dense region nearly empty in an
8000-event reference. The `outside` column is worth checking before trusting
any single observable, since a mis-set plot range silently drops events. That
is not hypothetical: `phi_k` was initially plotted over `[-30, 390]` when it
spans `[-180, 180]`, which discarded 31.5% of events while the histogram still
looked plausible.

`E_s` at D = 0.95 is a degenerate comparison, not a 95% discrepancy. The beam
energy loss is sampled as `xs**(1/targs)` with `targs` in radiation lengths
(0.0077 for a 5 cm target), so `1/targs ~ 130` and `eloss` spans about eighty
orders of magnitude:
97% of events land inside a `1e-4` GeV window at the top of the range, and the
whole variable spans less than 0.4% of the beam energy. The means agree to
`3e-4`. The port's closed form was checked directly against a NumPy
transcription of the Fortran's own accept/reject loop and the two agree. A KS
statistic on a distribution with that dynamic range measures the shape of a
numerically negligible tail.

The one substantive lesson from this comparison is in the next section.

### What is still open

- **The reference limits the cross-section comparison, not the port.** The
  widest gap between two Fortran references is 0.72%, and it brackets the
  residual being measured. A single long reference run (or several) settles it;
  nothing on the port side needs changing.
- **`mm^2`** is the last observable still above the floor: 2.3x on KS, and
  1 bin at 3.8σ. It is the same quantity on both sides — the cut is applied to
  the pre-exit electron energy
  but the recorded column comes from the final state, which uses the post-exit
  one, and the Fortran does exactly the same (it updates `ep` at
  `aao_rad.f90:1023` and only then calls `missm` again at `:1034`) — and the two
  means agree to 0.1%. The port's recorded maximum reaches 1.111 where the
  reference stops at 1.080 because the `eloss`-driven overshoot past the cut
  ceiling affects about one event in 10^4-10^5, and the port sample is five
  times larger.
- **Event-by-event comparison is impossible.** `myran` seeds from `unixtime`,
  so two Fortran runs disagree with each other. Everything above is
  distributional as a result.
- **No GPU was available** while this was written, so the batching claims are
  unmeasured here; only the CPU numbers are.

### The shared-key trap

The single largest discrepancy found during the port was one line of code:

```python
csthcm    = 2.0 * u(k_cm) - 1.0
phicm_deg = 360.0 * u(k_cm)      # same key as csthcm
```

`cos(theta*)` and `phi*` are two independent `myran` calls in the Fortran
(`aao_rad.f90:794-795`). Drawing both from one JAX key makes
`phi* = 180*(cos(theta*) + 1)` **exactly** — verified
`max |phi* - 180(cos(theta*)+1)| = 0.0` — which confines the pion decay
direction to a curve on the sphere instead of covering it.

The instructive part is how well it hid. Both marginals stay perfectly uniform
under the constraint, so neither histogram moved: `phi*`'s KS distance was 0.26
before the fix and 0.0133 after, but that is the *aggregate*, and the
cross-section residual stayed at +0.6% throughout because the two errors
cancelled — the survivor fraction came out +2.2% and the mean weight −1.6%.
The per-bin view is what separates it: 12 of 12 bins over 3σ at max abs z = 59
before, 0 of 40 at max abs z = 2.3 after.

What did expose it was plotting the photon-energy acceptance profile
`A(ek)` (`validation/compare_mm2.py`). Sharing the key correlates the decay
direction with `cos(theta*)`, and `cos(theta*)` is one of the terms in the
`ek_max` window, so it moved the *acceptance* as a function of `ek` by up to
+31% mid-range and −84% at the top while leaving `sigma` nearly intact. After
the fix `A(ek)` agrees within ~1% across the well-populated range.

The general lesson, and the reason `tests/test_kinematics.py` has a
`TestPhotonDecayAngles` class asserting *joint* independence rather than
marginal uniformity: with a batched RNG the failure mode is not "this variable
is wrong" but "these variables are correlated", and marginals cannot see it.
The regression test checks the correlation and the mean of `phi*` within each
`cos(theta*)` decile; both marginals still pass with the bug in place, which is
exactly the point.

### A second, output-only bug: `asym_p`

After the shared-key fix, `asym_p` was the last observable still disagreeing:
5 of 40 per-bin ratios over 3σ, max abs z = 5.30, and a p50 ratio of 0.709. This
one was **in the Fortran, not the port**, and it is worth writing down because
the way it hid is instructive.

`asym_p` is produced by `dsigma`, which the original calls from two places:

- the soft branch, `if (ek .lt. delta)` → `call dsigma(th0, qsq, epw, ...)` at
  `aao_rad.f90:807`, where `epw` is the driver's `sqrt(w_sq)`;
- inside `real function sigma(...)` at `aao_rad.f90:1366`, with
  `epw = sqrt(mf2)` and `mf2 = uu - 2*ek*(u0 - pu*csthk)`.

Those are two *different* hadronic masses — the second is lower by the radiated
photon energy — and `asym_p` depends on `W`, so the branches genuinely need
separate evaluations. The port does both and selects between them on `ek`.

The original's problem was scoping. `asym_p` is declared twice: once in the
driver and once, separately, inside `sigma()`. No `COMMON` block connects them,
so the driver's copy was written **only** by the soft branch. For every
`ek >= delta` event — the radiative branch, 24% of the accepted sample —
`ntp(32) = asym_p` recorded whatever the last *soft trial* had left behind, not
the event's own value. Trials outnumber accepted events by ~3e4, so the stale
values came from an unrecorded trial with unrelated kinematics.

Isolating it took evaluating the port's response at the *reference rows' own
kinematics*, split by branch. Each row is compared against the port's `asym`
for the same event, at that branch's `W`, so the only difference left is which
`W` (and, before the fix, which event) produced the recorded number. The floor
is `1.36 * sqrt(2/n)` for an n-vs-n comparison:

| reference | branch | rows | KS | x floor | std ref | std port | median abs diff |
|---|---|---|---|---|---|---|---|
| before fix | soft | 6039 | 0.0015 | 0.06x | 0.0462 | 0.0462 | 7.9e-6 |
| before fix | radiative | 1962 | 0.1172 | **2.70x** | **0.0628** | 0.0433 | 4.7e-2 |
| after fix | soft | 6038 | 0.0015 | 0.06x | 0.0453 | 0.0453 | 7.9e-6 |
| after fix | radiative | 1963 | 0.0234 | **0.54x** | 0.0445 | 0.0440 | 5.7e-3 |

The soft rows reproduce the port pointwise to 7.9e-6 — float32 noise — which
confirms the port's soft branch and pins the two `W` conventions. The radiative
rows matched neither, and their spread was **wider than the port's**, which is
the signature of a value drawn from a different point in kinematic space rather
than of a mis-convention: a wrong `W` would have been wrong by a systematic
amount, not drawn from a wider pool. The first hypothesis was wrong and worth
recording — stale-in-a-loop would imply radiative rows duplicating earlier rows,
and they do not (2 of 1962 coincide with any accepted soft row), because the
value comes from the last soft *trial*, not the last soft *event*.

After the fix the radiative rows sit at 0.54x the floor, i.e. indistinguishable
from the port, and `validation/compare_references.py` confirms that all 32
event columns of the reference are unchanged between the two runs.

The fix is three lines in the instrumented `aao_rad.f90` — promote `asym_p` to
`COMMON /radasy/` shared by the driver and `sigma()`, and initialise it to zero
on entry to `sigma()` so the three early returns before `dsigma` (`mf2` below
threshold, `ffac <= 0`, `gfac <= 0`) cannot leave it stale either. No physics
changes: `dsigma` already computed the value, it was being discarded.

`tests/test_kinematics.py::TestBeamAsymmetryBranch` pins the port side: it
checks that each branch's `asym_p` equals `response_functions` evaluated at that
branch's own `W`, that the *other* branch's `W` gives a materially different
answer (so the test cannot pass with the branch selection dropped), and that
`asym_p` is not a repeated value. Reverting the port's branch selection makes
the radiative half fail.

## Differences from the Fortran original

Deliberate, and each one is a config option you can turn off:

- **Precise physical constants.** `alpha = 1/137.035999084` rather than
  `1/137`, `pi` rather than the Fortran's single-precision `3.14159`, `m_pi`
  from PDG rather than `0.1395`. Worth ~0.1% in `sig_r`.
- **Saturation past the table edge.** `w_max="clamp"` (default) reproduces the
  original's saturating lookup; passing a number rejects those trials instead.
  Rejecting is more physical but drops the cross section by ~30% on a window
  with a third of its weight above `W = 2`.
- **Truncated photon-energy sampling.** `ek_sampling="truncated"` (default)
  draws `ek` from the exponential truncated to `[0, ek_max]` instead of
  redrawing, which is the same distribution and roughly 2x more efficient.
  `"fortran"` reproduces the original's `myran` transform, including its 1e-9
  resolution cap.
- **Branch-free rejection.** The original's `go to` chain becomes masks, so a
  rejected trial occupies a slot rather than being redrawn in place. This
  changes throughput, not the sampled distribution.
- **Estimated acceptance ceiling.** `sigr_max` is gone. The ceiling is
  estimated from a pilot batch, scaled by `--weight-max-margin`, and raised
  again if a later batch exceeds it. See
  [Performance and tuning](#performance-and-tuning).
- **Table cache.** `.tbl` parsed once into `.npz`.
- **`spence` in float64.** The Fortran's is a 100-step Riemann sum in float32,
  which is not reproducible in float64 and not worth reproducing.
- **Four LUND track lines written by default**, where the original declared four
  particles and wrote two. `--two-tracks` restores the original behaviour.

## Bugs found in the Fortran original

Found while porting; each is annotated in the source at the corresponding place.
Line numbers are given for the instrumented `aao_rad.f90` in this branch, which
the validation dumps shift around — search for the quoted statement instead.

- `read_sf_file.f90` skips the wrong number of columns for the `M_{L-}` row, so
  the shipped table reader mis-parses. `validation/dump_table.f90` is the
  corrected reader, and the port's parser is validated bit-for-bit against it.
- **`asym_p` is never computed for the radiative branch.** It is a local of
  `sigma()`, so the `dsigma` call inside that function wrote a value that died
  with the call and the driver's `ntp(32)` kept whatever the last *soft trial*
  had left in it. Since the radiative branch is the one `sigma()` serves, every
  event with `ek >= delta` — 24% of the sample — recorded an unrelated trial's
  beam asymmetry. The physics was always computed; only the output column was
  wrong. Fixed in this branch by promoting `asym_p` to `COMMON /radasy/` and
  initialising it so the three early returns before `dsigma` cannot leave it
  stale either. See [the asym_p write-up](#a-second-output-only-bug-asym_p).
- `interp.f90` has its `STOP` statements commented out, so an out-of-range
  lookup silently returns garbage instead of failing.
- `multipole_amps.f90` clamps `W > 2` to `2` *and* writes the clamped value back
  into a `COMMON` block, so a later caller sees a wrong `w`. The port
  reproduces the clamp on the lookup only, evaluating `nu_cm`, `qv_mag_cm`,
  `ppi_mag_cm` and `ekin` at the unclamped `W` as `maid_lee.f90` does.
- `sigr_max = sigr_max * fmcall` (`:546`) zeroes a user-supplied ceiling when
  `fmcall = 0`.
- `sig_tot = sig_tot - sigr_max` (`:1028`) subtracts the *ceiling* on the
  exit-loss rejection path instead of `sigr`, biasing `sig_sum`.
- `delphi` (`:714-718`) is computed geometrically, bounded to `[pi/9, 2pi]`, and
  then unconditionally overwritten with `pi/9` at `:720`, so the whole
  computation is dead code.
- The region-5 band test at `:750` uses `.or.` where the surrounding logic implies
  `.and.`, so a trial near either band is rejected.
- The LUND writer emits `px py pz E` into fields the format defines as
  `E px py pz`, which is why the original's output will not load.
- `missm()` computes `wreal` from the post-entrance energy `es` but `mm2` from
  the nominal `ebeam`, so the two are not consistent with each other. The port
  reproduces this — it is load-bearing for the missing-mass cut's `ek`
  dependence, not a typo to be tidied away.

## Development

```sh
pip install -e ".[dev]"        # pytest, ruff, and the validation plotting deps
pytest                         # 118 tests
ruff check src tests validation
```

`uv.lock` is there if you prefer `uv sync && uv run pytest`; it is not the only
supported route.

The suite needs `parms/spp_tbl/*.tbl`; without it the table-dependent tests skip
rather than fail, so a bare checkout still runs.

Layout:

| file | contents |
|---|---|
| `config.py` | `GeneratorConfig`, presets, `ep_window_for_w` |
| `table.py` | `.tbl` parsing, `.npz` conversion and cache |
| `interpolate.py` | batched linear and spline interpolation |
| `amplitudes.py` | multipole, helicity and CGLN amplitude sums |
| `xsection.py` | the AO response functions |
| `motsa.py` | the Mo & Tsai radiative kernel and the soft-photon folding |
| `kinematics.py` | the hadronic two-body decay and the missing mass |
| `generate.py` | sampling, the integrand, and the batching/streaming driver |
| `lund.py` | LUND and NPZ writers |

One note for reading `generate.py`: the integration region factors
(`mcfac`, `mpfac`) look wrong until you check them against the region-sampling
block in `aao_rad.f90` (`:711`–`:758`). They are reproduced exactly, and
`compare_sigma_points.py` checks the reconstruction to six digits per region.

## Building the Fortran reference

Only needed for validation. The instrumented sources are in `src/`; the Fortran
in this branch differs from upstream only by clearly-marked `write(13,...)` /
`write(14,...)` validation dumps.

```sh
mkdir -p build && cd build && cmake .. && cmake --build . --target aao_rad
```

The Fortran locates `spp_tbl/` via `CLAS_PARMS`, so run it from a directory that
has one:

```sh
mkdir -p /tmp/frun && cd /tmp/frun
ln -s "$REPO/parms" spp_tbl
"$REPO/build/aao_rad" < "$REPO/validation/maid07_pi+.inp"
```

This writes `aao_rad.ntuple` (47 columns: the 32 n-tuple fields plus 15
validation mirrors, one row per accepted event), `aao_rad.trials` (12 columns,
one row per 97th trial that reached the weight stage), and `aao_rad.lund`.

`AAO_TRIALS=0` skips the trial dump. It is ~0.5 GB per 8000 events and the
formatting costs about as much wall time as the physics, so long runs that only
need the event n-tuple should set it. Nothing in it touches the random
sequence.

The trial dump's column 10 is the Fortran's `sigr` *after* the region,
multipole and Jacobian factors, i.e. the full trial weight — not the bare cross
section. `validation/compare_sigr.py` recovers the raw cross section from it as
`weight / (mcfac * mpfac * jacob)`, which is what separates a physics difference
from a geometry one.

Two traps in the n-tuple columns, both of which cost time:

- **`soft_ek` is only assigned in the radiative branch** (`aao_rad.f90:934`, past
  the `go to 28` at `:914` that the soft branch takes). Soft events therefore record
  whatever `ek` the last radiative *trial* held. Do not use it to identify the
  branch — and note the corollary, which is the surprise: **75% of accepted
  events are soft**, not the ~10% the `E_gamma` mean of 0.014 GeV suggests,
  because the radiative tail is heavy enough to carry the mean. The reliable
  branch indicator is `ntp(17) = eg < delta`.
- The trial dump is written *after* the missing-mass test but *before* the
  accept, and the n-tuple records the post-exit `missm` call at `:1034`, not the
  pre-exit one at `:941`. The `mm2` in the trial file and the `mm2` in the
  n-tuple are therefore different quantities by construction.

The `validation/dump_*.f90` drivers lift the code under test verbatim out of
`aao_rad.f90` at build time, so they cannot drift from it:

```sh
validation/build_dump_sigma.sh    # sigma(), 32 output columns
validation/build_dump_soft.sh     # the soft-photon branch
```