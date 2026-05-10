"""Build hedz and P2K Talking Hedz as workbench pills + overlay comparison.

Outputs (forge/scratch/):
  - hedz.pill.json               — compiled-v1 pill, transcribed from hedz_rom.rs
  - talking_hedz_p2k.pill.json   — compiled-v1 pill, built from the calibration
                                    pole/zero truth rebuilt at 44.1 kHz
  - hedz_vs_p2k.png              — 4-corner overlay of response magnitude
  - hedz_{corner}.wav / p2k_{corner}.wav — impulse responses, per corner
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
from scipy.io import wavfile

REPO = Path(__file__).resolve().parents[2]
HEDZ_ROM = REPO / "trench-core/src/hedz_rom.rs"
CAL = REPO / "docs/calibration/Talking_Hedz.json"
OUT = REPO / "forge/scratch"

SR = 44100
IMPULSE_LEN = 16384
N_BINS = 512
FMIN, FMAX = 20.0, 20000.0
CORNERS = ("M0_Q0", "M100_Q0", "M0_Q100", "M100_Q100")


def parse_hedz_rom() -> tuple[dict, dict]:
    txt = HEDZ_ROM.read_text()
    m = re.search(
        r"HEDZ_CORNERS: \[CornerData; NUM_CORNERS\] = \[(.*)\];",
        txt, flags=re.DOTALL,
    )
    stage_re = re.compile(r"\[([-0-9eE\.\+ ,]+)\]\s*,\s*//\s*stage")
    flat = [
        tuple(float(x) for x in mm.group(1).split(","))
        for mm in stage_re.finditer(m.group(1))
    ]
    assert len(flat) == 24
    mb = re.search(r"HEDZ_BOOSTS: \[f64; NUM_CORNERS\] = \[(.*?)\];", txt)
    boosts_flat = [float(x.strip()) for x in mb.group(1).split(",") if x.strip()]
    corners = {label: flat[i * 6:(i + 1) * 6] for i, label in enumerate(CORNERS)}
    return corners, dict(zip(CORNERS, boosts_flat))


def pole_zero_to_kernel(pole_freq_hz, radius, zero_b1, zero_b2, sr):
    omega = 2.0 * math.pi * pole_freq_hz / sr
    a1 = -2.0 * radius * math.cos(omega)
    a2 = radius * radius
    return (1.0, zero_b1, zero_b2, a1, a2)


def build_p2k_corners(cal: dict, sr: float) -> tuple[dict, dict]:
    corners, boosts = {}, {}
    for label in CORNERS:
        stages = []
        for s in cal["corners"][label]["stages"]:
            stages.append(pole_zero_to_kernel(
                s["pole_freq_hz"], s["radius"],
                s["zeros"]["b1"], s["zeros"]["b2"], sr,
            ))
        corners[label] = stages
        boosts[label] = float(cal["boost"])
    return corners, boosts


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


def mag_db_on_grid(imp: np.ndarray, sr: float, grid: np.ndarray) -> np.ndarray:
    spec = np.fft.rfft(imp)
    mag = np.maximum(np.abs(spec), 1e-12)
    db = 20.0 * np.log10(mag)
    f = np.fft.rfftfreq(imp.size, d=1.0 / sr)
    return np.interp(grid, f, db)


def compiled_v1(name: str, corners: dict, boosts: dict, sr: int,
                source_note: str) -> dict:
    kfs = []
    for label in CORNERS:
        morph = 1.0 if label.startswith("M100") else 0.0
        q = 1.0 if label.endswith("Q100") else 0.0
        stages = [
            {"c0": s[0], "c1": s[1], "c2": s[2], "c3": s[3], "c4": s[4]}
            for s in corners[label]
        ]
        kfs.append({
            "label": label, "morph": morph, "q": q,
            "boost": boosts[label], "stages": stages,
        })
    return {
        "format": "compiled-v1",
        "name": name,
        "sampleRate": sr,
        "provenance": {
            "source": "forge/scratch/build_hedz_and_p2k.py",
            "note": source_note,
        },
        "keyframes": kfs,
    }


def save_wav(imp: np.ndarray, sr: int, path: Path) -> None:
    peak = float(np.abs(imp).max())
    norm = imp / max(peak, 1e-12) * 0.5
    wavfile.write(str(path), sr, (norm * 32767).astype(np.int16))


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    hedz_corners, hedz_boosts = parse_hedz_rom()
    cal = json.loads(CAL.read_text())
    p2k_corners, p2k_boosts = build_p2k_corners(cal, SR)

    (OUT / "hedz.pill.json").write_text(json.dumps(compiled_v1(
        "hedz", hedz_corners, hedz_boosts, SR,
        "Transcribed from trench-core/src/hedz_rom.rs (MD save compiled "
        "via pyruntime/heritage_coeffs Type 1/2/3 recipes).",
    ), indent=2))

    (OUT / "talking_hedz_p2k.pill.json").write_text(json.dumps(compiled_v1(
        "talking_hedz_p2k", p2k_corners, p2k_boosts, SR,
        "Built from docs/calibration/Talking_Hedz.json pole/zero truth; "
        f"biquads reconstructed at {SR} Hz (calibration sr was "
        f"{cal['sample_rate']}); b0=1, a1=-2rcos(w), a2=r^2, "
        "b1/b2=zeros.b1/b2.",
    ), indent=2))

    grid = np.logspace(math.log10(FMIN), math.log10(FMAX), N_BINS)
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True, sharey=True)
    for i, label in enumerate(CORNERS):
        ax = axes.flat[i]
        h_imp = render_impulse(hedz_corners[label], IMPULSE_LEN, hedz_boosts[label])
        p_imp = render_impulse(p2k_corners[label], IMPULSE_LEN, p2k_boosts[label])
        ax.semilogx(grid, mag_db_on_grid(h_imp, SR, grid),
                    label="hedz (MD save)", color="#00aaff", linewidth=1.6)
        ax.semilogx(grid, mag_db_on_grid(p_imp, SR, grid),
                    label="P2K Talking Hedz (ROM RE)", color="#ff6644", linewidth=1.6)
        ax.set_title(label)
        ax.set_xlim(FMIN, FMAX)
        ax.set_ylim(-60, 40)
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(loc="lower right", fontsize=8)
        ax.set_xlabel("Hz")
        ax.set_ylabel("dB")

        save_wav(h_imp, SR, OUT / f"hedz_{label}.wav")
        save_wav(p_imp, SR, OUT / f"p2k_{label}.wav")

    fig.suptitle("hedz (MD save) vs P2K Talking Hedz (ROM RE) — 44.1 kHz")
    plt.tight_layout()
    plt.savefig(OUT / "hedz_vs_p2k.png", dpi=100, bbox_inches="tight")
    print(f"wrote {OUT / 'hedz.pill.json'}")
    print(f"wrote {OUT / 'talking_hedz_p2k.pill.json'}")
    print(f"wrote {OUT / 'hedz_vs_p2k.png'}")
    print("impulse WAVs per corner written alongside")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
