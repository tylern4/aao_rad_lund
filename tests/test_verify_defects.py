"""The verifier must refuse to average in a reference that could not count.

The Fortran declares ``integer*4 ntries`` / ``integer*4 ngeom``, so a
configuration needing more than 2**31 trials wraps the counter its cross
section divides by and prints a *negative* sigma -- while still looking like a
finished run.  On this grid that happened to 46 of 192 runs, and including them
reported a fortran-vs-py_gpu sigma ratio spanning -10.4 to +5.4 with a |z| of
24.8: indistinguishable from a catastrophic port failure, and entirely an
artifact of the reference.

These tests pin the gate using the actual log signatures, including the one
that matters most -- an overflowed run that still prints a plausible sigma.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "validation" / "perlmutter"))

from verify_statistics import (  # noqa: E402
    FORT_RECORD_BYTES,
    count_ntuple_records,
    defect_summary,
    defects,
    parse_sigma,
)

N_EVENTS = 20_000


def ntuple_bytes(records: int) -> bytes:
    """A fixed-width Fortran record: 47 items of '1x,es16.8' plus a newline."""
    one = (" " + " ".join("0.12345678E+00" for _ in range(47))).ljust(FORT_RECORD_BYTES - 1)
    assert len(one) + 1 == FORT_RECORD_BYTES
    return (one + "\n").encode() * records


def write_run(tmp_path: Path, text: str, records: int | None = N_EVENTS) -> Path:
    run_dir = tmp_path / "run"
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / "out.txt").write_text(text)
    if records is not None:
        (run_dir / "aao_rad.ntuple").write_bytes(ntuple_bytes(records))
    return run_dir


CLEAN = """\
  ntries, ngeom =   152553895   111206141
  ntries, nevent, mcall_max:    211748729       20000           2
  Integrated cross section (MC, numerical) =   6.58308296E-03   6.58196583E-03  mu-barns
"""

# integer*4 wrapped past 2**31, so the counter sigma divides by is negative.
OVERFLOWED = """\
  ntries, ngeom =   2117483647   111206141
  ntries, ngeom =  -382487988 -1006273173
  ntries, nevent, mcall_max:   -382487988        8800           1
  Integrated cross section (MC, numerical) =  -1.50309518E-01  -1.50309518E-01  mu-barns
  missm-2: snthcm =   0.00000000
"""


def test_a_clean_fortran_run_has_no_defects(tmp_path):
    run_dir = write_run(tmp_path, CLEAN)
    assert defects("fortran", run_dir, CLEAN, N_EVENTS) == []


def test_a_clean_run_parses(tmp_path):
    run_dir = write_run(tmp_path, CLEAN)
    info = parse_sigma("fortran", run_dir, N_EVENTS)
    # sig_sum is the second capture group, per parse_sigma.
    assert info["sigma"] == pytest.approx(6.58196583e-03)
    assert info["trials"] == 211748729
    assert info["defects"] == []


def test_counter_overflow_is_caught(tmp_path):
    """The wrapped ntries prints with a leading minus the unsigned regex misses."""
    run_dir = write_run(tmp_path, OVERFLOWED, records=8800)
    found = defects("fortran", run_dir, OVERFLOWED, N_EVENTS)
    assert "counter_overflow" in found
    assert "non_positive_sigma" in found


def test_overflow_is_caught_even_when_the_ntuple_is_full(tmp_path):
    """The dangerous case: a full 20,000-record n-tuple and a sigma that still
    looks like a number.  Only the counter reveals the wrap."""
    text = CLEAN.replace("211748729", "-382487988").replace(
        "6.58308296E-03", "6.58308296E-03"
    )
    run_dir = write_run(tmp_path, text)
    found = defects("fortran", run_dir, text, N_EVENTS)
    assert "counter_overflow" in found
    # A positive sigma means the sign check cannot be what caught it.
    assert "non_positive_sigma" not in found


def test_short_ntuple_is_caught(tmp_path):
    run_dir = write_run(tmp_path, CLEAN, records=9577)
    found = defects("fortran", run_dir, CLEAN, N_EVENTS)
    assert any(d.startswith("short_ntuple") for d in found), found
    assert any("9577" in d for d in found), found


def test_missm_is_not_a_defect(tmp_path):
    """``missm-2`` must NOT exclude a run.

    It is a per-event clamp inside the generation loop (src/aao_rad.f90:1666-1668):
    when ``csthcm**2`` rounds just past 1 it prints the warning, sets snthcm to
    1e-7 and continues with the next event.  167 of the 192 real runs hit it,
    including runs with a full 20,000-record n-tuple and a healthy positive cross
    section.  Treating it as a defect discarded 167 good runs and left exactly one
    usable configuration out of 96.
    """
    text = CLEAN + "  missm-2: snthcm =   0.00000000\n" * 500
    run_dir = write_run(tmp_path, text)
    assert defects("fortran", run_dir, text, N_EVENTS) == []
    info = parse_sigma("fortran", run_dir, N_EVENTS)
    assert info["sigma"] > 0.0
    assert info["defects"] == []


def test_missm_does_not_mask_a_real_defect(tmp_path):
    """It must not shield a run that genuinely overflowed."""
    text = OVERFLOWED + "  missm-2: snthcm =   0.00000000\n" * 500
    run_dir = write_run(tmp_path, text, records=8800)
    found = defects("fortran", run_dir, text, N_EVENTS)
    assert "counter_overflow" in found
    assert "missm-2" not in found


def test_a_non_positive_sigma_is_flagged_on_its_own_merits(tmp_path):
    """No counter overflow, but a cross section of zero is still unusable -- the
    sign check must not be gated on some other defect being present."""
    text = CLEAN.replace("6.58196583E-03", "0.00000000E+00")
    run_dir = write_run(tmp_path, text)
    found = defects("fortran", run_dir, text, N_EVENTS)
    assert found == ["non_positive_sigma"], found


def test_the_port_is_never_flagged(tmp_path):
    """The port has no 32-bit counters and writes the n-tuple itself, so none of
    these checks may fire on it -- otherwise the gate would eat the good arm."""
    run_dir = tmp_path / "port"
    run_dir.mkdir()
    (run_dir / "out.txt").write_text(
        "sigma (MC)       : 0.00648899 micro-barn\n"
        "throughput       : 70 events/s (7,092,083 trials/s)\n"
    )
    assert defects("py_gpu", run_dir, "", N_EVENTS) == []


def test_defect_summary_counts_and_excludes_clean(tmp_path):
    clean = write_run(tmp_path / "a", CLEAN)
    bad = write_run(tmp_path / "b", OVERFLOWED, records=8800)
    results = {
        "fortran": {
            "clean_run": {"defects": defects("fortran", clean, CLEAN, N_EVENTS)},
            "bad_run": {"defects": defects("fortran", bad, OVERFLOWED, N_EVENTS)},
        }
    }
    summary = defect_summary(results)
    assert summary["fortran"]["counter_overflow"] == 1
    assert summary["fortran"]["clean"] == 1


def test_defect_summary_buckets_the_short_ntuple_detail(tmp_path):
    """Each short run must not get its own report line: the counts are per defect
    kind, with the exact record counts kept in the per-run detail."""
    runs = {}
    for i, n in enumerate((19032, 18084, 9582)):
        d = write_run(tmp_path / f"s{i}", CLEAN, records=n)
        runs[f"run{i}"] = {"defects": defects("fortran", d, CLEAN, N_EVENTS)}
    summary = defect_summary({"fortran": runs})
    assert summary["fortran"] == {"short_ntuple": 3}
    # The detail is still recoverable per run.
    assert "9582" in runs["run2"]["defects"][0]


def test_defect_counts_can_exceed_the_excluded_run_count(tmp_path):
    """A run may carry several defects, so the per-defect total can exceed the
    number of runs dropped.  The report prints both numbers, so make sure the
    distinction is real rather than a coincidence of the fixture."""
    d = write_run(tmp_path, OVERFLOWED, records=8800)
    results = {"fortran": {"run0": {"defects": defects("fortran", d, OVERFLOWED, N_EVENTS)}}}
    summary = defect_summary(results)
    n_runs = sum(1 for r in results["fortran"].values() if r["defects"])
    assert sum(v for k, v in summary["fortran"].items() if k != "clean") > n_runs


def test_parse_sigma_reports_defects_through_to_the_caller(tmp_path):
    """collect() reads defects off parse_sigma, so they have to travel together."""
    run_dir = write_run(tmp_path, OVERFLOWED, records=8800)
    info = parse_sigma("fortran", run_dir, N_EVENTS)
    assert info["defects"], "defects must survive parse_sigma"
    assert info["sigma"] is not None  # the whole problem: a sigma is still printed


def test_record_count_does_not_depend_on_the_byte_width(tmp_path):
    """The count must survive a record that is not 800 bytes wide.

    Dividing the file size by 800 is only right while es16.8 prints to exactly 16
    characters.  If it ever stops, every run reads as short, the reference arm
    drops out entirely, and the report calls that missing data rather than a
    defect in its own check.
    """
    ntp = tmp_path / "odd.ntuple"
    ntp.write_bytes(b"0.5 1.5 2.5\n" * 4242)  # narrow records, no relation to 800
    assert ntp.stat().st_size % FORT_RECORD_BYTES != 0  # not divisible by the old width
    assert count_ntuple_records(ntp) == 4242


def test_a_full_record_in_the_real_format_counts_as_one(tmp_path):
    ntp = tmp_path / "real.ntuple"
    ntp.write_bytes(ntuple_bytes(20_000))
    assert ntp.stat().st_size == 20_000 * FORT_RECORD_BYTES
    assert count_ntuple_records(ntp) == 20_000


def test_a_truncated_last_record_still_counts(tmp_path):
    """A run killed mid-write leaves a partial final line; it is still a record,
    and it is the short_ntuple check's job to notice the run fell short."""
    ntp = tmp_path / "cut.ntuple"
    ntp.write_bytes(ntuple_bytes(9) + b" 0.1")  # 9 whole records, then a partial one
    assert count_ntuple_records(ntp) == 10


def test_a_hundred_record_run_reads_as_short_not_as_missing(tmp_path):
    """The failure mode this guards: 100 records must be 100 short, not 'no file'."""
    run_dir = write_run(tmp_path, CLEAN, records=100)
    found = defects("fortran", run_dir, CLEAN, N_EVENTS)
    assert found == ["short_ntuple(100<20000)"], found
