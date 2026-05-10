"""Offline prototype: Talking Hedz with internal SR resampling.

Signal path:
  host SR impulse (44.1k)
    --resample_poly-->  native SR (39062.5)
    --cascade@native--> filtered
    --resample_poly-->  host SR
    --FFT-->            magnitude response at host SR

Reference:
  Direct cascade render at native SR on a native impulse → its magnitude
  response at native SR, plotted on the same Hz axis.

Target: RMS dB delta between the two curves under 1 dB (limited by the
resampler's stopband/passband). If this holds, the engine-side path is
worth building. If it doesn't, we know offline before writing Rust.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
from fractions import Fraction

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal

REPO = Path(__file__).resolve().parents[2]
CAL = REPO / "docs/calibration/Talking_Hedz.json"
OUT_PLOT = REPO / "forge/scratch/talking_hedz_resample.png"
OUT_TXT = REPO / "forge/scratch/talking_hedz_resample.txt"

HOST_SR = 44100
NATIVE_SR = 39062.5
IMPULSE_LEN_HOST = 32768
FMIN, FMAX, NBINS = 20.0, 19000.0, 1024  # cap below native Nyquist 19531
CORNERS = ("M0_Q0", "M100_Q0", "M0_Q100", "M100_Q100")


def stage_monic(s, sr):
    omega = 2.0 * math.pi * s["pole_freq_hz"] / sr
    return (
        1.0, s["zeros"]["b1"], s["zeros"]["b2"],
        -2.0 * s["radius"] * math.cos(omega),
        s["radius"] * s["radius"],
    )


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
    imp = render(stages, 32768, 1.0)
    raw = float(mag_db(imp, sr, grid).max())
    return 10.0 ** ((target_db - raw) / 20.0)


def main():
    cal = json.loads(CAL.read_text())
    grid = np.logspace(math.log10(FMIN), math.log10(FMAX), NBINS)

    # FFT-based resampler handles the 44100 / 39062.5 ratio cleanly offline.
    native_len = int(round(IMPULSE_LEN_HOST * NATIVE_SR / HOST_SR))

    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True, sharey=True)
    lines = [
        f"Talking Hedz with internal resampling (FFT) — "
        f"host {HOST_SR} Hz, native {NATIVE_SR} Hz",
        f"{'corner':<12} {'rms_dB':>10} {'max_err_dB':>12} "
        f"{'native_peak':>12} {'resamp_peak':>12}",
    ]

    for i, label in enumerate(CORNERS):
        target_db = float(cal["corners"][label]["cascade_peak_db"])
        stages_native = [stage_monic(s, NATIVE_SR)
                          for s in cal["corners"][label]["stages"]]

        # Reference: cascade at native SR on native impulse
        b_native = boost_to_peak(stages_native, NATIVE_SR, target_db, grid)
        imp_native_ref = render(stages_native, native_len, b_native)
        db_native = mag_db(imp_native_ref, NATIVE_SR, grid)

        # Resampling path: host impulse -> native -> cascade -> host
        imp_host = np.zeros(IMPULSE_LEN_HOST, dtype=np.float64)
        imp_host[0] = 1.0
        imp_at_native = signal.resample(imp_host, native_len)
        filtered = render(stages_native, native_len, b_native)
        imp_back_to_host = signal.resample(filtered, IMPULSE_LEN_HOST)
        db_resamp = mag_db(imp_back_to_host, HOST_SR, grid)

        delta = db_resamp - db_native
        rms = float(np.sqrt(np.mean(delta ** 2)))
        max_err = float(np.max(np.abs(delta)))

        lines.append(
            f"{label:<12} {rms:>10.3f} {max_err:>12.3f} "
            f"{float(db_native.max()):>12.2f} {float(db_resamp.max()):>12.2f}"
        )

        ax = axes.flat[i]
        ax.semilogx(grid, db_native, color="#808080", lw=1.4, ls="--",
                    label=f"native @ {NATIVE_SR} Hz (reference)")
        ax.semilogx(grid, db_resamp, color="#00aaff", lw=1.4,
                    label=f"resampled, host {HOST_SR} Hz")
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

    fig.suptitle("Talking Hedz with internal resampling (offline prototype)")
    plt.tight_layout()
    plt.savefig(OUT_PLOT, dpi=100, bbox_inches="tight")


if __name__ == "__main__":
    main()
