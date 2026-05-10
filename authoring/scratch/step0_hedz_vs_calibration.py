"""Step 0 diagnostic: render baked `hedz` vs `Talking_Hedz` calibration.

Produces a 4-corner overlay plot to measure the gap between what the
runtime ships today (trench-core/src/hedz_rom.rs, compiled from the
`hedz` Morph Designer grid via tools/bake_hedz_const.py) and the
reverse-engineered calibration truth (docs/calibration/Talking_Hedz.json).

Loss: log-spaced response magnitude in dB, 512 bins, 20 Hz -> 20 kHz.
"""
from __future__ import annotations

import json
import math
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

REPO = Path(__file__).resolve().parents[2]
HEDZ_ROM = REPO / "trench-core/src/hedz_rom.rs"
CAL_PATH = REPO / "docs/calibration/Talking_Hedz.json"
OUT_PATH = REPO / "forge/scratch/step0_hedz_vs_calibration.png"

HEDZ_SR = 44100.0
IMPULSE_LEN = 8192
FREQ_MIN = 20.0
FREQ_MAX = 20000.0
N_BINS = 512
CORNER_LABELS = ("M0_Q0", "M100_Q0", "M0_Q100", "M100_Q100")


def parse_hedz_rom(path: Path) -> list[list[tuple[float, ...]]]:
    txt = path.read_text()
    m = re.search(
        r"pub const HEDZ_CORNERS: \[CornerData; NUM_CORNERS\] = \[(.*)\];",
        txt, flags=re.DOTALL,
    )
    assert m, "HEDZ_CORNERS block not found"
    block = m.group(1)
    stage_re = re.compile(r"\[([-0-9eE\.\+ ,]+)\]\s*,\s*//\s*stage")
    flat = [
        tuple(float(x) for x in mm.group(1).split(","))
        for mm in stage_re.finditer(block)
    ]
    assert len(flat) == 24, f"expected 24 stages, got {len(flat)}"
    return [flat[i * 6:(i + 1) * 6] for i in range(4)]


def parse_boosts(path: Path) -> list[float]:
    txt = path.read_text()
    m = re.search(r"HEDZ_BOOSTS: \[f64; NUM_CORNERS\] = \[(.*?)\];", txt)
    assert m, "HEDZ_BOOSTS not found"
    return [float(x.strip()) for x in m.group(1).split(",") if x.strip()]


def run_kernel_impulse(stages, length: int, boost: float) -> np.ndarray:
    """Mirrors thumbnail.py run_impulse. c0..c4 kernel form, DF2T."""
    state = [(0.0, 0.0) for _ in stages]
    out = np.zeros(length, dtype=np.float64)
    for n in range(length):
        sig = 1.0 if n == 0 else 0.0
        for i, (c0, c1, c2, c3, c4) in enumerate(stages):
            w1, w2 = state[i]
            y = c0 * sig + w1
            state[i] = (c1 * sig - c3 * y + w2, c2 * sig - c4 * y)
            sig = y
        out[n] = sig * boost
    return out


def run_biquad_cascade(biquads, length: int) -> np.ndarray:
    """Standard direct-form II transposed. b0..b2 / a1..a2 per stage."""
    state = [(0.0, 0.0) for _ in biquads]
    out = np.zeros(length, dtype=np.float64)
    for n in range(length):
        sig = 1.0 if n == 0 else 0.0
        for i, (b0, b1, b2, a1, a2) in enumerate(biquads):
            w1, w2 = state[i]
            y = b0 * sig + w1
            state[i] = (b1 * sig - a1 * y + w2, b2 * sig - a2 * y)
            sig = y
        out[n] = sig
    return out


def mag_db_on_grid(impulse: np.ndarray, sr: float,
                   grid: np.ndarray) -> np.ndarray:
    spec = np.fft.rfft(impulse)
    mag = np.maximum(np.abs(spec), 1e-12)
    db = 20.0 * np.log10(mag)
    freqs = np.fft.rfftfreq(impulse.size, d=1.0 / sr)
    return np.interp(grid, freqs, db)


def build_cal_biquads(cal: dict, label: str) -> list[tuple[float, ...]]:
    corner = cal["corners"][label]
    sr = cal["sample_rate"]
    biquads = []
    for s in corner["stages"]:
        omega = 2.0 * math.pi * s["pole_freq_hz"] / sr
        r = s["radius"]
        a1 = -2.0 * r * math.cos(omega)
        a2 = r * r
        b0 = 1.0
        b1 = s["zeros"]["b1"]
        b2 = s["zeros"]["b2"]
        biquads.append((b0, b1, b2, a1, a2))
    return biquads


def main() -> int:
    hedz_corners = parse_hedz_rom(HEDZ_ROM)
    hedz_boosts = parse_boosts(HEDZ_ROM)
    cal = json.loads(CAL_PATH.read_text())
    cal_sr = cal["sample_rate"]
    cal_boost = cal["boost"]

    grid = np.logspace(math.log10(FREQ_MIN), math.log10(FREQ_MAX), N_BINS)

    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True, sharey=True)
    axes = axes.flatten()

    gaps = {}
    for i, label in enumerate(CORNER_LABELS):
        hedz_imp = run_kernel_impulse(hedz_corners[i], IMPULSE_LEN, hedz_boosts[i])
        hedz_db = mag_db_on_grid(hedz_imp, HEDZ_SR, grid)

        cal_biquads = build_cal_biquads(cal, label)
        cal_imp = run_biquad_cascade(cal_biquads, IMPULSE_LEN) * cal_boost
        cal_db = mag_db_on_grid(cal_imp, cal_sr, grid)

        ax = axes[i]
        ax.semilogx(grid, hedz_db, label=f"hedz @ {int(HEDZ_SR)} Hz",
                    color="#00aaff", linewidth=1.6)
        ax.semilogx(grid, cal_db, label=f"calibration @ {cal_sr} Hz",
                    color="#ff6644", linestyle="--", linewidth=1.4)
        ax.set_title(label)
        ax.set_xlim(FREQ_MIN, FREQ_MAX)
        ax.set_ylim(-60, 80)
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(loc="lower right", fontsize=8)
        ax.set_xlabel("Hz")
        ax.set_ylabel("dB")

        rms = float(np.sqrt(np.mean((hedz_db - cal_db) ** 2)))
        gaps[label] = rms
        ax.text(0.02, 0.97, f"RMS gap: {rms:.1f} dB",
                transform=ax.transAxes, fontsize=9, va="top")

    fig.suptitle("step 0: hedz (baked) vs Talking_Hedz calibration", y=1.00)
    plt.tight_layout()
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(OUT_PATH, dpi=100, bbox_inches="tight")
    print(f"wrote {OUT_PATH}")
    for label, rms in gaps.items():
        print(f"  {label}: RMS dB gap = {rms:.2f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
