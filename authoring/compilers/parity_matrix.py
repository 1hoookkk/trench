"""4×4 capture × corner-position parity matrix.

For each hardware capture (M0_Q0, M100_Q0, M0_Q100, M100_Q100) test against
each predicted corner position of the cartridge (M0_Q0, M100_Q0, M0_Q100,
M100_Q100). Prints a grid of (correlation, mean_dB_error).

Use to diagnose whether a diverging corner is:
  - mislabeled on the hardware side (off-diagonal beats diagonal)
  - axis-flipped runtime side (every capture prefers the opposite morph/q)
  - actual coefficient error (no position matches cleanly)

Usage:
    python tools/parity_matrix.py <cartridge.json> <capture_stem> <capture_dir>

    <capture_stem> is the prefix before _M*_Q*.wav (e.g. "P2k_026")
"""
from __future__ import annotations

import subprocess
import sys
import json
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
PARITY_CHECK = REPO / "tools" / "hardware_parity_check.py"

CORNERS = [
    ("M0_Q0",     0.0, 0.0),
    ("M100_Q0",   1.0, 0.0),
    ("M0_Q100",   0.0, 1.0),
    ("M100_Q100", 1.0, 1.0),
]


def run_one(capture: Path, cartridge: Path, morph: float, q: float) -> dict | None:
    cmd = [
        sys.executable, str(PARITY_CHECK),
        "--capture", str(capture),
        "--cartridge", str(cartridge),
        "--morph", str(morph),
        "--q", str(q),
        "--input", "pink",
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    try:
        start = proc.stdout.index("{")
        return json.loads(proc.stdout[start:])
    except (ValueError, json.JSONDecodeError):
        return None


def main() -> int:
    if len(sys.argv) != 4:
        print("usage: parity_matrix.py <cartridge.json> <capture_stem> <capture_dir>")
        return 2

    cartridge = Path(sys.argv[1])
    capture_stem = sys.argv[2]
    capture_dir = Path(sys.argv[3])

    print(f"Cartridge: {cartridge}")
    print(f"Captures:  {capture_dir}/{capture_stem}_M*_Q*.wav\n")

    # Header
    print(f"{'capture ↓ / cart →':<22}", end="")
    for label, _, _ in CORNERS:
        print(f"{label:>18}", end="")
    print()
    print("-" * (22 + 18 * 4))

    for cap_label, _, _ in CORNERS:
        cap = capture_dir / f"{capture_stem}_{cap_label}.wav"
        if not cap.is_file():
            print(f"{cap_label:<22}  (missing capture)")
            continue
        print(f"{cap_label:<22}", end="")
        best_corr = -1.0
        best_label = ""
        for cart_label, morph, q in CORNERS:
            r = run_one(cap, cartridge, morph, q)
            if r is None:
                cell = "      err       "
            else:
                corr = r["correlation"]
                mean = r["mean_err_db"]
                cell = f"{corr:.3f}/{mean:5.1f}dB"
                if corr > best_corr:
                    best_corr = corr
                    best_label = cart_label
            print(f"{cell:>18}", end="")
        print(f"   → best: {best_label}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
