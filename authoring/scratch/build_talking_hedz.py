"""Build Talking Hedz compiled-v1 pill, 1:1 from calibration.

Source: docs/calibration/Talking_Hedz.json (ROM-RE pole/zero truth plus
per-corner cascade_peak_db summaries).

Strategy:
  - Each stage: monic numerator + pole denominator at host SR.
      c0 = 1.0
      c1 = zeros.b1
      c2 = zeros.b2
      c3 = a1 = -2 r cos(2 pi f_pole / host_sr)
      c4 = a2 = r^2
  - Per-corner boost empirically pins the cascade peak magnitude to the
    calibration's reported cascade_peak_db. Stage zero locations (shape)
    are already preserved by the monic form — boost only rescales.

This sidesteps the b0-per-stage ambiguity in the calibration's packed
format and produces a cartridge whose cascade response curve matches the
calibration ROM truth in both shape and absolute peak level per corner.
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
OUT_PILL = REPO / "forge/scratch/talking_hedz.compiled.json"
OUT_PLOT = REPO / "forge/scratch/talking_hedz_cascade.png"
OUT_VERIFY = REPO / "forge/scratch/talking_hedz_verify.txt"

HOST_SR = 44100
IMPULSE_LEN = 32768
FMIN, FMAX, NBINS = 20.0, 20000.0, 1024
CORNERS = ("M0_Q0", "M100_Q0", "M0_Q100", "M100_Q100")


def stage_kernel_monic(s: dict, sr: int) -> tuple[float, ...]:
    omega = 2.0 * math.pi * s["pole_freq_hz"] / sr
    a1 = -2.0 * s["radius"] * math.cos(omega)
    a2 = s["radius"] * s["radius"]
    return (1.0, s["zeros"]["b1"], s["zeros"]["b2"], a1, a2)


def render_impulse(stages, length: int, boost: float) -> np.ndarray:
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


def mag_db(imp: np.ndarray, sr: int, grid: np.ndarray) -> np.ndarray:
    spec = np.fft.rfft(imp)
    mag = np.maximum(np.abs(spec), 1e-12)
    db = 20.0 * np.log10(mag)
    f = np.fft.rfftfreq(imp.size, d=1.0 / sr)
    return np.interp(grid, f, db)


def main() -> int:
    cal = json.loads(CAL.read_text())
    grid = np.logspace(math.log10(FMIN), math.log10(FMAX), NBINS)

    keyframes = []
    reports = []
    stages_by_corner = {}
    dbs_by_corner = {}

    for label in CORNERS:
        stages = [stage_kernel_monic(s, HOST_SR) for s in cal["corners"][label]["stages"]]
        stages_by_corner[label] = stages

        # Measure peak at monic (boost=1)
        imp_monic = render_impulse(stages, IMPULSE_LEN, 1.0)
        peak_monic_db = float(mag_db(imp_monic, HOST_SR, grid).max())

        # Target peak = calibration cascade_peak_db
        target_db = float(cal["corners"][label]["cascade_peak_db"])
        corner_boost = 10.0 ** ((target_db - peak_monic_db) / 20.0)

        morph = 1.0 if label.startswith("M100") else 0.0
        q = 1.0 if label.endswith("Q100") else 0.0
        keyframes.append({
            "label": label,
            "morph": morph,
            "q": q,
            "boost": corner_boost,
            "stages": [
                {"c0": s[0], "c1": s[1], "c2": s[2], "c3": s[3], "c4": s[4]}
                for s in stages
            ],
        })

        # Verify final (with boost)
        imp_final = render_impulse(stages, IMPULSE_LEN, corner_boost)
        db_final = mag_db(imp_final, HOST_SR, grid)
        dbs_by_corner[label] = db_final
        measured_final = float(db_final.max())
        reports.append((label, target_db, peak_monic_db, corner_boost,
                         measured_final, measured_final - target_db))

    pill = {
        "format": "compiled-v1",
        "name": "talking_hedz",
        "sampleRate": HOST_SR,
        "provenance": {
            "source": "docs/calibration/Talking_Hedz.json",
            "authored_sr": cal["sample_rate"],
            "rebuilt_at_sr": HOST_SR,
            "strategy": "monic stage numerators; per-corner boost pins "
                        "cascade peak to calibration cascade_peak_db",
            "tool": "forge/scratch/build_talking_hedz.py",
        },
        "keyframes": keyframes,
    }
    OUT_PILL.parent.mkdir(parents=True, exist_ok=True)
    OUT_PILL.write_text(json.dumps(pill, indent=2))

    lines = [
        f"Talking Hedz cascade peak match (host sr {HOST_SR} Hz):",
        f"{'corner':<12} {'target_dB':>10} {'monic_dB':>10} "
        f"{'boost':>10} {'final_dB':>10} {'delta_dB':>10}",
    ]
    for label, target, monic, b, final, delta in reports:
        lines.append(
            f"{label:<12} {target:>10.2f} {monic:>10.2f} "
            f"{b:>10.4f} {final:>10.2f} {delta:>+10.3f}"
        )
    verify_text = "\n".join(lines)
    OUT_VERIFY.write_text(verify_text)
    print(verify_text)

    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True, sharey=True)
    for i, label in enumerate(CORNERS):
        ax = axes.flat[i]
        ax.semilogx(grid, dbs_by_corner[label], color="#ff6644",
                    linewidth=1.4, label=f"talking_hedz pill @ {HOST_SR} Hz")
        ax.axhline(cal["corners"][label]["cascade_peak_db"], ls="--", lw=1.0,
                   color="#808080",
                   label=f"cal cascade_peak_db = "
                         f"{cal['corners'][label]['cascade_peak_db']:.2f} dB")
        ax.set_title(label)
        ax.set_xlim(FMIN, FMAX)
        ax.set_ylim(-60, 40)
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(loc="lower right", fontsize=8)
        ax.set_xlabel("Hz")
        ax.set_ylabel("dB")

    fig.suptitle(
        f"Talking Hedz compiled-v1 pill — rebuilt from calibration at "
        f"{HOST_SR} Hz"
    )
    plt.tight_layout()
    plt.savefig(OUT_PLOT, dpi=100, bbox_inches="tight")
    print(f"wrote {OUT_PILL.name}, {OUT_PLOT.name}, {OUT_VERIFY.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
