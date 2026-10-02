"""Command line interface.

Two entry points are installed:

``aao-rad``
    Generate events.
``aao-rad-table``
    Inspect and convert the MAID tables.

The flag names of the original ``aao_rad`` driver script (``--trig``,
``--experiment``, ``--q2min``, ``--q2max``, ``--channel``, ``--unpolarized``,
``--seed``) are all still accepted, so existing shell scripts keep working.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from . import __version__
from .config import EXPERIMENTS, GeneratorConfig

log = logging.getLogger("aao_rad")

_EPILOG = """\
examples:
  aao-rad --experiment rgb --n-events 1000000 --output rgb.lund
  aao-rad --experiment default --n-events 5e6 --format npz --output rgb.npz
  aao-rad --config run.yaml --output run.lund
  aao-rad --legacy-input test.inp --output test.lund
  echo 7 | aao-rad --n-events 1000            # interactive, like the Fortran

The weight maximum used by the rejection sampler is estimated on the device,
so the original's 'sigr_max' prompt is gone and can no longer be got wrong.
"""


def _positive_int(text: str) -> int:
    """Accept ``1000000``, ``1e6`` and ``5M``."""
    t = text.strip().rstrip("kKmM").replace("_", "")
    mult = 1
    if t[-1:] in "kK":
        mult = 1_000
    elif t[-1:] in "mM":
        mult = 1_000_000
    return int(float(t) * mult)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="aao-rad",
        description="GPU-accelerated radiative pion electroproduction event generator",
        epilog=_EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--version", action="version", version=f"aao-rad {__version__}")

    out = p.add_argument_group("output")
    out.add_argument("-o", "--output", default="aao_rad.lund", help="output file (default: %(default)s)")
    out.add_argument(
        "-f", "--format", choices=("lund", "npz"), default="lund",
        help="lund = CLAS12 generator file, npz = lossless n-tuple (default: %(default)s)",
    )
    out.add_argument(
        "--chunk-events", type=_positive_int, default=500_000,
        help="events held in memory at once when streaming (default: %(default)s)",
    )
    out.add_argument(
        "--two-tracks", action="store_true",
        help="write only the electron and hadron lines, as the original did",
    )
    out.add_argument("-q", "--quiet", action="store_true", help="only report errors")

    run = p.add_argument_group("run control")
    run.add_argument("-n", "--n-events", "--trig", dest="n_events", type=_positive_int,
                     default=1_000_000, help="number of events to generate (default: %(default)s)")
    run.add_argument("--seed", type=int, default=None, help="PRNG seed (default: random)")
    run.add_argument("--batch-size", type=int, default=65_536,
                     help="trial points per vectorised step (default: %(default)s)")
    run.add_argument("--progress", type=_positive_int, default=0, metavar="N",
                     help="log progress every N events (0 = only at the end)")

    sel = p.add_argument_group("event selection")
    sel.add_argument("--experiment", choices=sorted(EXPERIMENTS), default=None,
                     help="preset for a known running condition")
    sel.add_argument("--channel", "-c", type=int, choices=(1, 3), default=None,
                     help="1 = pi0 p, 3 = pi+ n")
    sel.add_argument("--unpolarized", "--unpolarised", dest="unpolarized",
                     action="store_true", help="unpolarized beam")
    sel.add_argument("--q2min", type=float, help="minimum Q^2 [GeV^2]")
    sel.add_argument("--q2max", type=float, help="maximum Q^2 [GeV^2]")
    sel.add_argument("--emin", type=float, help="minimum scattered electron energy [GeV]")
    sel.add_argument("--emax", type=float, help="maximum scattered electron energy [GeV]")
    sel.add_argument("--beam-energy", type=float, help="incident electron energy [GeV]")
    sel.add_argument("--mm-cut", type=float, help="half-width of the (mm^2 - mm_exp^2) cut [GeV^2]")
    sel.add_argument("--delta", type=float, help="minimum photon energy treated explicitly [GeV]")
    sel.add_argument("--ntracks", type=int, choices=(2, 4), help="particles per event")

    tgt = p.add_argument_group("target")
    tgt.add_argument("--target-length", type=float, help="target length [cm]")
    tgt.add_argument("--target-radius", type=float, help="target radius [cm]")
    tgt.add_argument("--beam-x", type=float, help="beam x position [cm]")
    tgt.add_argument("--beam-y", type=float, help="beam y position [cm]")
    tgt.add_argument("--beam-z", type=float, help="beam z position [cm]")

    adv = p.add_argument_group("advanced")
    adv.add_argument("--parms", default=None,
                     help="directory holding spp_tbl/ (default: $AAO_RAD_PARMS or $CLAS_PARMS)")
    adv.add_argument("--interp", choices=("linear", "spline"), default=None,
                     help="table interpolation scheme (default: linear, as the original)")
    adv.add_argument("--ek-sampling", choices=("truncated", "fortran"), default=None,
                     help="photon energy sampling; 'fortran' reproduces the original RNG-limited range")
    adv.add_argument("--cos-step", type=float, help="csrng: width of the narrow angular bands")
    adv.add_argument("--weight-max-margin", type=float,
                     help="safety factor on the estimated acceptance ceiling "
                          "(default: %(default)s)")
    adv.add_argument("--k-exp", type=float, help="slope of the photon energy variable")
    adv.add_argument("--region", type=float, nargs=4, metavar=("R1", "R2", "R3", "R4"),
                     help="sizes of the four importance sampling regions")

    src = p.add_argument_group("input sources")
    src.add_argument("--config", type=Path, help="JSON/TOML/YAML configuration file")
    src.add_argument("--legacy-input", type=Path, metavar="FILE",
                     help="read the original interactive prompt sequence from FILE ('-' = stdin)")
    src.add_argument("--dump-config", action="store_true",
                     help="print the fully resolved configuration and exit")
    return p


def config_from_args(args: argparse.Namespace) -> GeneratorConfig:
    """Merge (in increasing priority) defaults, preset, config file, CLI flags."""
    cfg = GeneratorConfig()

    if args.experiment:
        cfg = cfg.merge(**EXPERIMENTS[args.experiment])

    if args.config:
        cfg = GeneratorConfig.from_file(args.config)

    if args.legacy_input:
        text = (
            sys.stdin.read()
            if str(args.legacy_input) == "-"
            else Path(args.legacy_input).read_text()
        )
        legacy = GeneratorConfig.from_legacy_input(text)
        # Keep anything the .inp file did not set (seed, batch size, ...).
        keep = {
            "seed": cfg.seed,
            "batch_size": cfg.batch_size,
            "interp_scheme": cfg.interp_scheme,
            "verbose_every": cfg.verbose_every,
        }
        cfg = legacy.merge(**keep)

    cli = {
        "n_events": args.n_events,
        "seed": args.seed,
        "batch_size": args.batch_size,
        "verbose_every": args.progress,
        "channel": args.channel,
        "q2_min": args.q2min,
        "q2_max": args.q2max,
        "ep_min": args.emin,
        "ep_max": args.emax,
        "beam_energy": args.beam_energy,
        "missing_mass_cut": args.mm_cut,
        "min_photon_energy": args.delta,
        "n_tracks": args.ntracks,
        "target_length_cm": args.target_length,
        "target_radius_cm": args.target_radius,
        "beam_x_cm": args.beam_x,
        "beam_y_cm": args.beam_y,
        "beam_z_cm": args.beam_z,
        "interp_scheme": args.interp,
        "ek_sampling": args.ek_sampling,
        "cos_step": args.cos_step,
        "k_exp": args.k_exp,
        "weight_max_margin": args.weight_max_margin,
        "regions": tuple(args.region) if args.region else None,
        "write_tracks": False if args.two_tracks else None,
    }
    if args.unpolarized:
        cli["polarized_beam"] = False
    return cfg.merge(**cli)


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.WARNING if args.quiet else logging.INFO,
        format="%(message)s",
        stream=sys.stderr,
    )

    try:
        cfg = config_from_args(args)
    except (ValueError, FileNotFoundError, ImportError) as exc:
        parser.error(str(exc))
        return 2  # unreachable, keeps type checkers happy

    if args.dump_config:
        import json

        print(json.dumps(cfg.to_dict(), indent=2, sort_keys=True))
        return 0

    # Imported lazily so that --help and --dump-config do not pay for jax.
    import jax

    from .generate import EventGenerator, build_grid
    from .lund import LundWriter

    log.info(
        "device: %s | channel %d | %.3f GeV beam | Q2 %.2f-%.2f | E' %.2f-%.2f",
        jax.devices()[0].platform.upper() + ": " + str(jax.devices()[0].device_kind),
        cfg.channel, cfg.beam_energy, cfg.q2_min, cfg.q2_max, cfg.ep_min, cfg.ep_max,
    )

    try:
        grid = build_grid(
            cfg.channel, parms_dir=args.parms, scheme=cfg.interp_scheme
        )
    except FileNotFoundError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    log.info("table: %s", _describe_table(grid, cfg, args))
    gen = EventGenerator(grid)

    stats_summary = ""
    try:
        if args.format == "lund":
            with LundWriter(args.output, cfg) as writer:
                for block, stats in gen.stream(cfg, chunk_events=args.chunk_events):
                    writer.write(block)
                    if stats is not None:
                        stats_summary = stats.summary()
        else:
            from .lund import write_npz

            blocks = []
            stats = None
            for block, blk_stats in gen.stream(cfg, chunk_events=args.chunk_events):
                blocks.append(block)
                if blk_stats is not None:
                    stats = blk_stats
            write_npz(args.output, iter(blocks))
            if stats is not None:
                stats_summary = stats.summary()
    except RuntimeError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    log.info("wrote %s events to %s", f"{cfg.n_events:,}", args.output)
    if not args.quiet and stats_summary:
        print(stats_summary, file=sys.stderr)
    return 0


def _describe_table(grid, cfg: GeneratorConfig, args) -> str:
    from .table import load_table

    try:
        t = load_table(cfg.channel, parms_dir=args.parms)
    except FileNotFoundError:  # pragma: no cover - already reported
        return "<unavailable>"
    q2 = t.q2_range
    w = t.w_range
    return (
        f"{t.source.name} [{t.amps.shape[1]}x{t.amps.shape[2]}] "
        f"Q2={q2[0]:g}..{q2[1]:g} W={w[0]:g}..{w[1]:g} {cfg.interp_scheme}"
    )


# ---------------------------------------------------------------------------
# aao-rad-table
# ---------------------------------------------------------------------------
def table_main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="aao-rad-table",
        description="Inspect and convert the MAID07 multipole tables",
    )
    parser.add_argument("--parms", default=None, help="directory holding spp_tbl/")
    parser.add_argument("--channel", "-c", type=int, default=3, help="production channel")
    parser.add_argument("--list", action="store_true", help="list available tables and exit")
    parser.add_argument("--convert", type=Path, metavar="OUT.npz",
                        help="write the parsed table to a compressed .npz and exit")
    parser.add_argument("--values", type=int, default=0, metavar="N",
                        help="print the first N tabulated values")
    parser.add_argument("--rebuild-cache", action="store_true",
                        help="ignore and overwrite the .npz parse cache")
    args = parser.parse_args(argv)

    from .table import find_parms_dir, iter_available_tables, load_table

    try:
        if args.list:
            for path in iter_available_tables(args.parms):
                print(path)
            return 0
        print(f"parms: {find_parms_dir(args.parms)}")
        table = load_table(args.channel, parms_dir=args.parms, use_cache=not args.rebuild_cache)
    except FileNotFoundError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    print(table)
    q2, w = table.q2_range, table.w_range
    print(f"  Q^2 grid : {table.q2.size} points, {q2[0]:g} .. {q2[1]:g} GeV^2")
    print(f"  W   grid : {table.w.size} points, {w[0]:g} .. {w[1]:g} GeV")
    print(f"  uniform  : {table.uniform}")

    if args.values:
        print("  first values (S0+ re/im, S1+ re/im, ...):")
        for k in range(min(args.values, table.amps.shape[0])):
            print(f"    sf{k + 1:<3d} {table.amps[k, 0, 0]: .6e}")

    if args.convert:
        np_save = __import__("numpy").savez_compressed
        np_save(
            args.convert,
            q2=table.q2, w=table.w, amps=table.amps, uniform=table.uniform,
        )
        print(f"wrote {args.convert}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
