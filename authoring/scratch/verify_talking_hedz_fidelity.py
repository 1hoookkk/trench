"""Real fidelity check for talking_hedz.compiled.json.

Renders:
  - calibration cascade at its authored sample rate (39062.5 Hz)
  - the pill's cascade at host sample rate (44100 Hz)
Both use identical monic stage kernels (c0=1) with a per-corner boost that
pins the cascade peak to calibration's reported cascade_peak_db. The only
structural difference between them is the sample rate the poles/zeros land
at — so the delta between the two curves is the SR-warping error of the
pill, i.e. the genuine fidelity loss of rebuilding at 44.1 kHz.

Reports RMS dB and peak-error dB per corner across the full 20 Hz - 20 kHz
log-spaced grid.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

REPO = Path(__file__).resolve().parents[2]
CAL = REPO / "docs/calibration/Talking_Hedz.json"
OUT_PLOT = REPO / "forge/scratch/talking_hedz_fidelity.png"
OUT_TXT = REPO / "forge/scratch/talking_hedz_fidelity.txt"

HOST_SR = 44100
IMPULSE_LEN = 32768
FMIN, FMAX, NBINS = 20.0, 20000.0, 1024
CORNERS = ("M0_Q0", "M100_Q0", "M0_Q100", "M100_Q100")


def stage_monic(s, sr):
    omega = 2.0 * math.pi * s["pole_freq_hz"] / sr
    a1 = -2.0 * s["radius"] * math.cos(omega)
    a2 = s["radius"] * s["radius"]
    return (1.0, s["zeros"]["b1"], s["zeros"]["b2"], a1, a2)


def render(stages, n, boost):
    state = [(0.0, 0.0) for _ in stages]
    out = np.zeros(n, dtype=np.float64)
    for i in range(n):
        x = 1.0 if i == 0 else 0.0
        for k, (c0, c1, c2, c3, c4) in enumerate(stages):
            w1, w2 = state[k]
            y = c0 * x + w1
            state[k] = (c1 * x - c3 * y + w2, c2 * x - c4 * y)
            x = y
        out[i] = x * boost
    return out


def mag_db(imp, sr, grid):
    spec = np.fft.rfft(imp)
    f = np.fft.rfftfreq(imp.size, d=1.0 / sr)
    return np.interp(grid, f, 20.0 * np.log10(np.maximum(np.abs(spec), 1e-12)))


def boost_to_peak(stages, sr, target_db, grid):
    imp = render(stages, IMPULSE_LEN, 1.0)
    raw_peak = float(mag_db(imp, sr, grid).max())
    return 10.0 ** ((target_db - raw_peak) / 20.0)


def main():
    cal = json.loads(CAL.read_text())
    cal_sr = float(cal["sample_rate"])
    grid = np.logspace(math.log10(FMIN), math.log10(FMAX), NBINS)

    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True, sharey=True)
    lines = [
        f"Talking Hedz fidelity — calibration @ {cal_sr} Hz vs pill @ {HOST_SR} Hz",
        f"Both cascades peak-normalized to calibration cascade_peak_db per corner.",
        f"{'corner':<12} {'rms_dB':>10} {'max_err_dB':>12} {'peak_cal':>10} {'peak_pill':>10}",
    ]
    for i, label in enumerate(CORNERS):
        stages_cal = [stage_monic(s, cal_sr) for s in cal["corners"][label]["stages"]]
        stages_pill = [stage_monic(s, HOST_SR) for s in cal["corners"][label]["stages"]]
        target_db = float(cal["corners"][label]["cascade_peak_db"])

        b_cal = boost_to_peak(stages_cal, cal_sr, target_db, grid)
        b_pill = boost_to_peak(stages_pill, HOST_SR, target_db, grid)

        db_cal = mag_db(render(stages_cal, IMPULSE_LEN, b_cal), cal_sr, grid)
        db_pill = mag_db(render(stages_pill, IMPULSE_LEN, b_pill), HOST_SR, grid)

        delta = db_pill - db_cal
        rms = float(np.sqrt(np.mean(delta ** 2)))
        max_err = float(np.max(np.abs(delta)))
        peak_cal = float(db_cal.max())
        peak_pill = float(db_pill.max())

        lines.append(
            f"{label:<12} {rms:>10.2f} {max_err:>12.2f} "
            f"{peak_cal:>10.2f} {peak_pill:>10.2f}"
        )

        ax = axes.flat[i]
        ax.semilogx(grid, db_cal, color="#808080", lw=1.4, ls="--",
                    label=f"calibration @ {int(cal_sr)} Hz")
        ax.semilogx(grid, db_pill, color="#ff6644", lw=1.4,
                    label=f"pill @ {HOST_SR} Hz")
        ax.set_title(f"{label}  (RMS {rms:.2f} dB, max |err| {max_err:.2f} dB)")
        ax.set_xlim(FMIN, FMAX)
        ax.set_ylim(-60, 50)
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(loc="lower right", fontsize=8)
        ax.set_xlabel("Hz")
        ax.set_ylabel("dB")

    text = "\n".join(lines)
    print(text)
    OUT_TXT.write_text(text)

    fig.suptitle("Talking Hedz fidelity — native vs host SR")
    plt.tight_layout()
    plt.savefig(OUT_PLOT, dpi=100, bbox_inches="tight")


if __name__ == "__main__":
    main()
