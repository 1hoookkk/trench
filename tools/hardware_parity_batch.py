"""Hardware parity batch harness — runs hardware_parity_check.py against
every registered (capture, cartridge, morph, Q) triple and reports pass /
fail as a table. Exit code 0 only if all triples pass their tolerances.

This is the gate-pyramid rung that catches engine math drift against
real hardware. Unlike every other existing gate (hedz_cascade,
talking_hedz_parity, render_diff_harness) which verify internal
consistency, this one fails when the engine diverges from captured
E-mu hardware.

To add a new triple: append to the TRIPLES table. Each entry points at
a 44.1 kHz mono WAV capture of a P2K filter ringing out after brief
excitation (impulse-response style) with known (morph, Q).

Usage:

    python tools/hardware_parity_batch.py          # run all, print table
    python tools/hardware_parity_batch.py --plot   # also write PNGs
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path
from typing import NamedTuple

REPO = Path(__file__).resolve().parent.parent
PARITY_CHECK = REPO / "tools" / "hardware_parity_check.py"


class Triple(NamedTuple):
    label: str
    capture: Path
    cartridge: Path
    morph: float
    q: float


# Only captures where we have ground-truth (morph, Q) labels in the
# filename. Benchmark directory one-off captures lack position info and
# are deferred until positions are recorded alongside the WAVs.
CAPTURE_ROOT = Path(r"C:\Users\hooki\trenchwork_clean")
P2K_ROOT = REPO / "trench-juce" / "cartridges" / "p2k"

TRIPLES: list[Triple] = [
    Triple("P2k_013 M0 Q0",
           CAPTURE_ROOT / "P2k_013_m0_q0.wav",
           P2K_ROOT / "P2k_013.json", 0.0, 0.0),
    Triple("P2k_013 M50 Q0",
           CAPTURE_ROOT / "P2k_013_m50_q0.wav",
           P2K_ROOT / "P2k_013.json", 0.5, 0.0),
    Triple("P2k_013 M100 Q0",
           CAPTURE_ROOT / "P2k_013_m100_q0.wav",
           P2K_ROOT / "P2k_013.json", 1.0, 0.0),
    Triple("P2k_013 M50 Q50",
           CAPTURE_ROOT / "P2k_013_m50_q50.wav",
           P2K_ROOT / "P2k_013.json", 0.5, 0.5),
]


# Tolerances chosen to pass the current engine (which nulls against
# hardware) and fail any material drift. mean and max are in dB; corr
# is Pearson over the 60 Hz – 16 kHz compare band.
TOL_MEAN_DB = 2.0
TOL_MAX_DB  = 50.0
TOL_CORR    = 0.95


def run_one(tr: Triple, plot_dir: Path | None) -> tuple[bool, dict]:
    argv = [
        sys.executable, str(PARITY_CHECK),
        "--capture",   str(tr.capture),
        "--cartridge", str(tr.cartridge),
        "--morph", str(tr.morph),
        "--q",     str(tr.q),
        "--input", "impulse",
        "--tol-mean-db", str(TOL_MEAN_DB),
        "--tol-max-db",  str(TOL_MAX_DB),
        "--tol-corr",    str(TOL_CORR),
    ]
    if plot_dir is not None:
        plot_dir.mkdir(parents=True, exist_ok=True)
        png = plot_dir / (tr.label.replace(" ", "_") + ".png")
        argv += ["--plot", str(png)]

    proc = subprocess.run(argv, capture_output=True, text=True)
    # The checker prints JSON as its stdout payload; exit code is the verdict.
    try:
        payload = json.loads(proc.stdout)
    except json.JSONDecodeError:
        payload = {"stdout": proc.stdout.strip(), "stderr": proc.stderr.strip()}
    ok = proc.returncode == 0
    return ok, payload


def format_row(tr: Triple, payload: dict, ok: bool) -> str:
    mean = payload.get("mean_err_db", "?")
    mx   = payload.get("max_err_db",  "?")
    corr = payload.get("correlation", "?")
    tag  = "PASS" if ok else "FAIL"
    return f"  [{tag}] {tr.label:<22} mean={mean!s:>6}  max={mx!s:>6}  corr={corr!s:>7}"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--plot", action="store_true",
                    help="write per-triple PNG overlay plots to parity_plots/")
    ap.add_argument("--plot-dir", type=Path, default=REPO / "parity_plots")
    args = ap.parse_args(argv)

    plot_dir = args.plot_dir if args.plot else None

    missing = [tr for tr in TRIPLES
               if not tr.capture.exists() or not tr.cartridge.exists()]
    if missing:
        print("!! Missing inputs:")
        for tr in missing:
            if not tr.capture.exists():
                print(f"   capture   not found: {tr.capture}")
            if not tr.cartridge.exists():
                print(f"   cartridge not found: {tr.cartridge}")
        return 2

    print(f"hardware parity batch — {len(TRIPLES)} triples")
    print(f"  tolerances: mean≤{TOL_MEAN_DB} dB  max≤{TOL_MAX_DB} dB  corr≥{TOL_CORR}")
    print()
    failed = 0
    for tr in TRIPLES:
        ok, payload = run_one(tr, plot_dir)
        print(format_row(tr, payload, ok))
        if not ok:
            failed += 1

    print()
    print(f"  {len(TRIPLES) - failed} / {len(TRIPLES)} passed")
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
