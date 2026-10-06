#!/usr/bin/env python3
"""Generate the Perlmutter scan grid: one run card per configuration + manifest.

The grid spans beam energies from 2 to 12 GeV crossed with channel (pi+ n and
pi0 p), beam polarisation, the explicit-photon-energy cut and the target
thickness -- a full cross product, so every setting is scanned against every
other one.

The ``E'`` window is derived per beam energy so the whole ``(Q^2, E')``
rectangle stays inside the MAID07 table (``W`` kept in [1.20, 1.98] GeV,
comfortably within the table's [1.08, 2.0]) *and* below both kinematic clamps
the original applies (``aao_rad.f90`` clamps ``q2_max`` at the 90-degree
elastic limit and ``ep_max`` at the W-threshold limit; the port clamps both
identically, ``generate.py``).  Keeping the windows inside the clamps means
neither code silently samples a different rectangle, which would turn the
cross-code comparison into a useless diagnostic.

Outputs (under ``--root``):
    grid/cfg_NNN.txt   the run card for configuration NNN (the Fortran card
                       format: 22 whitespace-separated values, read from
                       stdin by ``aao_rad`` and by ``--legacy-input``)
    grid/manifest.csv  one row per *run job* (configuration x seed); every
                       batch script walks this file.

The same manifest drives the Fortran, JAX-on-CPU and JAX-on-GPU batches, so
the three implementations sample identical run cards.  The seeds only matter
to the port (the Fortran seeds ``myran`` from ``unixtime``); CPU and GPU use
the *same* seeds per configuration, which is deliberate -- any CPU-vs-GPU
difference is then pure backend floating-point, not sampling.
"""

from __future__ import annotations

import argparse
import csv
import math
from itertools import product
from pathlib import Path

# The Fortran hard-codes mp = .938 GeV (see the port's constants module
# docstring), so the window derivation uses that value: it is what
# aao_rad.f90's clamps actually evaluate.
MP = 0.938
M_PI = {3: 0.13957039, 1: 0.1349768}  # pi+, pi0
MEL_WG_MARGIN = 0.0005  # wg = mp + m_pi + margin, as in aao_rad.f90

# Fixed kinematics, matching the validated reference card
# (/tmp/frun_fixed/run_card.txt): Q^2 in [0.2, 1.9] is below the 90-degree
# elastic clamp for every beam energy down to 2 GeV (2.55 GeV^2 there).
Q2_MIN, Q2_MAX = 0.2, 1.9

# W window kept inside the MAID07 table's [1.08, 2.0], with margin so the
# rounded E' limits stay inside.
W_LO, W_HI = 1.20, 1.98

ENERGIES = (2.0, 4.244, 6.0, 8.0, 10.0, 12.0)
CHANNELS = (3, 1)  # pi+ n (the validated channel) first, then pi0 p
POLARIZED = (1, 0)
DELTAS = (0.005, 0.05)
TARGETS = (5.0, 2.5)

# The rest of the card: identical to the validated reference card.
THEORY = 7
NPART = 4
REGIONS = (0.20, 0.12, 0.20, 0.20)
MM_CUT = 0.2
R_TARG = 0.486
VERTEX = (0.3, 0.03, 0.0)
FMCALL = 1.0

N_EVENTS = 20_000
SEEDS_PER_CFG = 2
SEED_BASE = 8_675_309

CHANNEL_NAME = {3: "pi+ n", 1: "pi0 p"}


def ep_window(ebeam: float, channel: int) -> tuple[float, float]:
    """E' window that keeps W in [W_LO, W_HI] across the whole (Q^2, E') rectangle.

    W^2 = mp^2 + 2*mp*nu - Q^2 with nu = ebeam - E', so the maximum W sits at
    (Q^2_min, E'_min) and the minimum W at (Q^2_max, E'_max).  Rounded
    *inwards* to 0.01 GeV so the rounded window is strictly inside the
    unrounded one.
    """
    nu_max = (W_HI**2 - MP**2 + Q2_MIN) / (2.0 * MP)
    nu_min = (W_LO**2 - MP**2 + Q2_MAX) / (2.0 * MP)
    ep_min = math.ceil((ebeam - nu_max) * 100.0) / 100.0
    ep_max = math.floor((ebeam - nu_min) * 100.0) / 100.0
    return ep_min, ep_max


def elastic_q2_limit(ebeam: float) -> float:
    """The 90-degree elastic Q^2 clamp applied by aao_rad.f90:322-323."""
    s = 0.5
    return 4.0 * ebeam**2 * s / (1.0 + 2.0 * ebeam * s / MP)


def threshold_ep_limit(ebeam: float, channel: int) -> float:
    """The W-threshold E' clamp applied by aao_rad.f90:336-338."""
    wg = MP + M_PI[channel] + MEL_WG_MARGIN
    return ebeam - (wg**2 + Q2_MIN - MP**2) / (2.0 * MP)


def w_of(ebeam: float, q2: float, ep: float) -> float:
    return math.sqrt(MP**2 + 2.0 * MP * (ebeam - ep) - q2)


def card_text(cfg: dict) -> str:
    """The 22-value legacy card, in the Fortran's prompt order.

    ``n_events`` is read from ``cfg`` rather than the module constant so the card
    and the manifest can never disagree: a configuration whose event quota was
    reduced to keep the reference's 32-bit counters inside range must have *both*
    the run card and the manifest updated, or the Fortran would generate the
    original quota while the verifier expected the reduced one and reported every
    run as short.
    """
    lines = [
        str(THEORY),
        str(cfg["polarized"]),
        " ".join(f"{r:g}" for r in REGIONS),
        str(NPART),
        str(cfg["channel"]),
        f"{MM_CUT:g}",
        f"{cfg['t_target']:g}",
        f"{R_TARG:g}",
        *[f"{v:g}" for v in VERTEX],
        f"{cfg['ebeam']:g}",
        f"{Q2_MIN:g} {Q2_MAX:g}",
        f"{cfg['ep_min']:g} {cfg['ep_max']:g}",
        f"{cfg['delta']:g}",
        str(cfg.get("n_events", N_EVENTS)),
        f"{FMCALL:g}",
    ]
    return "\n".join(lines) + "\n"


def check_config(cfg: dict) -> None:
    """Raise if the card would be clamped by either code or leave the table."""
    e = cfg["ebeam"]
    tag = cfg["cfg_id"]
    if not 0.0 < cfg["ep_min"] < cfg["ep_max"] < e:
        raise ValueError(f"{tag}: invalid E' window {cfg['ep_min']}..{cfg['ep_max']}")
    if cfg["q2_max"] >= elastic_q2_limit(e):
        raise ValueError(
            f"{tag}: Q^2_max {cfg['q2_max']} hits the elastic clamp "
            f"{elastic_q2_limit(e):.3f} at {e:g} GeV"
        )
    if cfg["ep_max"] >= threshold_ep_limit(e, cfg["channel"]):
        raise ValueError(
            f"{tag}: E'_max {cfg['ep_max']} hits the W-threshold clamp "
            f"{threshold_ep_limit(e, cfg['channel']):.3f} at {e:g} GeV"
        )
    corners = {
        (Q2_MIN, cfg["ep_min"]): "max W",
        (Q2_MAX, cfg["ep_max"]): "min W",
    }
    for (q2, ep), label in corners.items():
        w = w_of(e, q2, ep)
        if not W_LO <= w <= W_HI:
            raise ValueError(f"{tag}: {label} = {w:.3f} leaves [{W_LO}, {W_HI}]")


def load_overrides(path: Path | None) -> dict[str, int]:
    """``cfg_id,n_events`` rows restricting a configuration's event quota.

    Only a *reduction* is accepted.  The quota exists to keep ``integer*4
    ntries`` (src/aao_rad.f90:154) inside 2**31 for configurations whose trials
    per event run to 3e5; raising it would reintroduce the overflow this is for,
    and a quota above the default is never useful.
    """
    if path is None:
        return {}
    out: dict[str, int] = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            cfg_id, raw = row["cfg_id"].strip(), row["n_events"].strip()
            n = int(raw)
            if n <= 0:
                raise ValueError(f"{path}: {cfg_id} has non-positive n_events {n}")
            if n > N_EVENTS:
                raise ValueError(
                    f"{path}: {cfg_id} asks for {n} events, more than the "
                    f"standard quota of {N_EVENTS}"
                )
            out[cfg_id] = n
    return out


def build_configs(overrides: dict[str, int] | None = None) -> list[dict]:
    overrides = overrides or {}
    cfgs = []
    idx = 0
    for ebeam, channel, pol, delta, targ in product(ENERGIES, CHANNELS, POLARIZED, DELTAS, TARGETS):
        ep_min, ep_max = ep_window(ebeam, channel)
        cfg_id = f"cfg_{idx:03d}"
        cfg = {
            "cfg_id": cfg_id,
            "ebeam": ebeam,
            "channel": channel,
            "channel_name": CHANNEL_NAME[channel],
            "polarized": pol,
            "delta": delta,
            "t_target": targ,
            "q2_min": Q2_MIN,
            "q2_max": Q2_MAX,
            "ep_min": ep_min,
            "ep_max": ep_max,
            "n_events": overrides.get(cfg_id, N_EVENTS),
        }
        check_config(cfg)
        cfgs.append(cfg)
        idx += 1
    return cfgs


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--root",
        default=None,
        help="scan root (default: $AAO_SCAN_ROOT or $SCRATCH/aao_rad_scan)",
    )
    p.add_argument(
        "--n-events-override",
        default=None,
        type=Path,
        help=(
            "CSV of cfg_id,n_events restricting those configurations' event "
            "quota; applied to the run cards and the manifest together so the "
            "Fortran and the verifier always agree on the quota.  Defaults to "
            "<root>/grid/n_events_override.csv when that exists, because every "
            "sbatch script re-runs this with only --root: an override passed as "
            "a flag alone would be silently dropped the next time any job "
            "started, and the quota it set would survive only in the file it "
            "rewrote."
        ),
    )
    args = p.parse_args()

    import os

    root = (
        args.root
        or os.environ.get("AAO_SCAN_ROOT")
        or os.path.join(os.environ.get("SCRATCH", "/tmp"), "aao_rad_scan")
    )
    grid_dir = f"{root}/grid"
    os.makedirs(grid_dir, exist_ok=True)

    override_path = args.n_events_override
    if override_path is None:
        default_path = Path(grid_dir) / "n_events_override.csv"
        override_path = default_path if default_path.is_file() else None
    overrides = load_overrides(override_path)
    if overrides:
        print(
            f"event quotas overridden for {len(overrides)} configuration(s) "
            f"from {override_path}"
        )

    cfgs = build_configs(overrides)

    for cfg in cfgs:
        path = os.path.join(grid_dir, f"{cfg['cfg_id']}.txt")
        with open(path, "w") as fh:
            fh.write(card_text(cfg))

    manifest = os.path.join(grid_dir, "manifest.csv")
    with open(manifest, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(
            [
                "run_id",
                "cfg_id",
                "seed",
                "ebeam",
                "channel",
                "channel_name",
                "polarized",
                "delta",
                "t_target",
                "q2_min",
                "q2_max",
                "ep_min",
                "ep_max",
                "mm_cut",
                "n_events",
                "theory",
                "regions",
                "r_targ",
                "vx",
                "vy",
                "vz",
            ]
        )
        n = 0
        for cfg in cfgs:
            for si in range(SEEDS_PER_CFG):
                w.writerow(
                    [
                        f"{cfg['cfg_id']}_s{si}",
                        cfg["cfg_id"],
                        SEED_BASE + 1000 * int(cfg["cfg_id"].split("_")[1]) + si,
                        cfg["ebeam"],
                        cfg["channel"],
                        cfg["channel_name"],
                        cfg["polarized"],
                        cfg["delta"],
                        cfg["t_target"],
                        cfg["q2_min"],
                        cfg["q2_max"],
                        cfg["ep_min"],
                        cfg["ep_max"],
                        MM_CUT,
                        cfg["n_events"],
                        THEORY,
                        " ".join(f"{r:g}" for r in REGIONS),
                        R_TARG,
                        *VERTEX,
                    ]
                )
                n += 1

    print(f"grid: {len(cfgs)} configurations x {SEEDS_PER_CFG} seeds = {n} run jobs")
    print(f"cards + manifest under {grid_dir}")
    for e in ENERGIES:
        ep_min, ep_max = ep_window(e, CHANNELS[0])
        print(
            f"  E = {e:5.3f} GeV:  Q^2 {Q2_MIN}-{Q2_MAX} GeV^2,  "
            f"E' {ep_min:.2f}-{ep_max:.2f} GeV,  "
            f"W {w_of(e, Q2_MAX, ep_max):.3f}-{w_of(e, Q2_MIN, ep_min):.3f} GeV"
        )
    print(
        f"settings axes: channel {CHANNEL_NAME[3]}/{CHANNEL_NAME[1]}, "
        f"polarized {POLARIZED}, delta {DELTAS}, target {TARGETS} cm"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
