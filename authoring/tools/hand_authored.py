"""forge/hand_authored.py — hand-authored cartridge editor.

Procedural extrapolation isn't getting us there. Edit the STAGES_M0 and
STAGES_M100 lists directly: each stage is one biquad with explicit pole
frequency, pole radius, zero frequency, zero radius, and c0 attenuation.

No formulas. No calibration. Whatever you write IS what gets baked.

Workflow:
  1. Edit STAGES_M0 / STAGES_M100 below.
  2. Run:  python forge/hand_authored.py
  3. Plot opens. Audition WAVs at cartridges/.../audition/_hand_authored_*.wav.
  4. Iterate.

Each stage: {pole_hz, pole_r, zero_hz, zero_r, c0}.
- pole_hz: where the resonance peak lives (Hz, < 19531 = Nyquist at 39062.5 SR)
- pole_r:  pole radius. 0.98-0.99 is the user-specified range. Higher = sharper.
- zero_hz: where the notch carves (Hz). Position relative to pole_hz controls
           peak-notch character: zero below pole = bandpass-ish, zero above = bandstop-ish.
- zero_r:  zero radius. 1.0 = on unit circle (infinite notch depth). < 1.0 = interior
           (finite notch). Heritage uses 0.5-0.85 typically.
- c0:      per-stage attenuation. < 1.0. Compounds across cascade — hot c0 = louder.

The runtime LERPs every coefficient from M0 → M100 over the morph axis. M0_Q* and
M100_Q* are duplicated (Q is a runtime knob, not authored).
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


REPO = Path(__file__).resolve().parent.parent
PNG_PATH = REPO / "cartridges" / "factory" / "generated" / "qlaw" / "_hand_authored.png"
JSON_PATH = REPO / "cartridges" / "factory" / "generated" / "qlaw" / "_hand_authored.json"
AUTHORING_SR = 39062.5

# ─── EDIT BELOW ───────────────────────────────────────────────────────────────

CARTRIDGE_NAME = "_hand_authored"

# M0 frame — surgical violent end. Six bandpass features distributed across the band.
# Pole frequencies log-spaced ~factor 3 to keep peaks visually separated.
STAGES_M0 = [
    # pole_hz, pole_r, zero_hz, zero_r, c0
    # Each stage: zero on unit circle ~2 semis below pole = clean LF rolloff
    # per stage AND deep notch below each peak. c0 tapered so LF contributes
    # less to the cumulative shelf. No poles below ~250 Hz: heritage Anchor
    # at 200 Hz adds LF shelf that the rest of the design fights uphill.
    {"pole_hz":   250, "pole_r": 0.990, "zero_hz":   220, "zero_r": 0.92, "c0": 0.50},
    {"pole_hz":   600, "pole_r": 0.990, "zero_hz":   530, "zero_r": 0.92, "c0": 0.60},
    {"pole_hz":  1400, "pole_r": 0.990, "zero_hz":  1240, "zero_r": 0.92, "c0": 0.70},
    {"pole_hz":  3200, "pole_r": 0.990, "zero_hz":  2850, "zero_r": 0.92, "c0": 0.80},
    {"pole_hz":  7200, "pole_r": 0.990, "zero_hz":  6400, "zero_r": 0.92, "c0": 0.90},
    {"pole_hz": 14500, "pole_r": 0.985, "zero_hz": 13000, "zero_r": 0.92, "c0": 0.95},
]

# M100 frame — wider, broader, less surgical. Same six features but lower Q,
# wider zero offsets, interior zeros (less surgical notches), pole frequencies
# shifted to a different distribution. Each stage c0 independent.
STAGES_M100 = [
    {"pole_hz":   400, "pole_r": 0.980, "zero_hz":   320, "zero_r": 0.85, "c0": 0.55},
    {"pole_hz":   900, "pole_r": 0.980, "zero_hz":   720, "zero_r": 0.85, "c0": 0.65},
    {"pole_hz":  2000, "pole_r": 0.980, "zero_hz":  1600, "zero_r": 0.85, "c0": 0.75},
    {"pole_hz":  4500, "pole_r": 0.980, "zero_hz":  3600, "zero_r": 0.85, "c0": 0.85},
    {"pole_hz":  9000, "pole_r": 0.980, "zero_hz":  7200, "zero_r": 0.85, "c0": 0.90},
    {"pole_hz": 16000, "pole_r": 0.980, "zero_hz": 13000, "zero_r": 0.85, "c0": 0.95},
]

# ─── END EDIT ─────────────────────────────────────────────────────────────────


def stage_coefs(s: dict, sr: float = AUTHORING_SR) -> dict[str, float]:
    """Hand-authored stage dict → compiled-v1 (c0..c4). Numerator scales with c0
    so zero locations are preserved under attenuation: num = c0·(1 + b1z⁻¹ + b2z⁻²)."""
    pole_hz = float(s["pole_hz"])
    pole_r  = float(s["pole_r"])
    zero_hz = float(s["zero_hz"])
    zero_r  = float(s["zero_r"])
    c0      = float(s["c0"])
    if zero_hz >= sr * 0.5:
        zero_hz = sr * 0.49
    theta_p = 2.0 * math.pi * pole_hz / sr
    theta_z = 2.0 * math.pi * zero_hz / sr
    b1 = -2.0 * zero_r * math.cos(theta_z)
    b2 = zero_r * zero_r
    return {
        "c0": c0,
        "c1": c0 * b1,
        "c2": c0 * b2,
        "c3": -2.0 * pole_r * math.cos(theta_p),
        "c4": pole_r * pole_r,
    }


def pad_passthrough(stages: list[dict]) -> list[dict]:
    out = list(stages)
    while len(out) < 12:
        out.append({"c0": 1.0, "c1": 0.0, "c2": 0.0, "c3": 0.0, "c4": 0.0})
    return out


def cascade_db(stages: list[dict], freqs: np.ndarray, sr: float = AUTHORING_SR) -> np.ndarray:
    w = 2.0 * np.pi * freqs / sr
    e1 = np.exp(-1j * w); e2 = np.exp(-2j * w)
    mag = np.ones_like(freqs)
    for s in stages:
        c0, c1, c2, c3, c4 = s["c0"], s["c1"], s["c2"], s["c3"], s["c4"]
        if c3 == 0.0 and c4 == 0.0:
            continue
        num = c0 + c1*e1 + c2*e2
        den = 1.0 + c3*e1 + c4*e2
        mag *= np.abs(num/den)
    return 20*np.log10(np.maximum(mag, 1e-12))


def stage_db(s: dict, freqs: np.ndarray, sr: float = AUTHORING_SR) -> np.ndarray:
    """Single-stage magnitude in dB."""
    w = 2.0 * np.pi * freqs / sr
    e1 = np.exp(-1j * w); e2 = np.exp(-2j * w)
    c0, c1, c2, c3, c4 = s["c0"], s["c1"], s["c2"], s["c3"], s["c4"]
    if c3 == 0.0 and c4 == 0.0:
        return np.zeros_like(freqs)
    num = c0 + c1*e1 + c2*e2
    den = 1.0 + c3*e1 + c4*e2
    return 20*np.log10(np.maximum(np.abs(num/den), 1e-12))


def write_cartridge(stages_m0_compiled: list[dict], stages_m100_compiled: list[dict]) -> None:
    cart = {
        "format": "compiled-v1",
        "name": CARTRIDGE_NAME,
        "sampleRate": AUTHORING_SR,
        "recipe": "hand_authored",
        "provenance": {
            "mode": "hand_authored",
            "stages_m0_source": STAGES_M0,
            "stages_m100_source": STAGES_M100,
        },
        "keyframes": [
            {"label": "M0_Q0",    "morph": 0.0, "q": 0.0, "boost": 1.0, "stages": stages_m0_compiled},
            {"label": "M0_Q100",  "morph": 0.0, "q": 1.0, "boost": 1.0, "stages": stages_m0_compiled},
            {"label": "M100_Q0", "morph": 1.0, "q": 0.0, "boost": 1.0, "stages": stages_m100_compiled},
            {"label": "M100_Q100","morph": 1.0,"q": 1.0, "boost": 1.0, "stages": stages_m100_compiled},
        ],
    }
    JSON_PATH.parent.mkdir(parents=True, exist_ok=True)
    JSON_PATH.write_text(json.dumps(cart, indent=2), encoding="utf-8")


def main() -> int:
    stages_m0   = pad_passthrough([stage_coefs(s) for s in STAGES_M0])
    stages_m100 = pad_passthrough([stage_coefs(s) for s in STAGES_M100])

    freqs = np.geomspace(20.0, 18000.0, 4096)
    db_m0   = cascade_db(stages_m0,   freqs)
    db_m100 = cascade_db(stages_m100, freqs)
    peak_m0   = float(np.max(db_m0))
    peak_m100 = float(np.max(db_m100))

    print(f"M0   stages: {len(STAGES_M0)}, cascade peak {peak_m0:+.2f} dB")
    print(f"M100 stages: {len(STAGES_M100)}, cascade peak {peak_m100:+.2f} dB")

    fig, axes = plt.subplots(2, 2, figsize=(16.0, 9.0), dpi=140, sharex=True)
    fig.patch.set_facecolor("#16191a")
    fig.suptitle(f"HAND_AUTHORED  ::  {CARTRIDGE_NAME}\n"
                 "left: cascade   right: per-stage breakdown",
                 color="#ffba00", fontweight="bold", fontsize=13)

    cascade_panels = [
        (axes[0,0], f"M0 cascade ({peak_m0:+.1f} dB)",   db_m0,   STAGES_M0,   "#ff5555"),
        (axes[1,0], f"M100 cascade ({peak_m100:+.1f} dB)", db_m100, STAGES_M100, "#ff7777"),
    ]
    breakdown_panels = [
        (axes[0,1], "M0 per-stage",   stages_m0[:len(STAGES_M0)],   STAGES_M0,   "#ff5555"),
        (axes[1,1], "M100 per-stage", stages_m100[:len(STAGES_M100)], STAGES_M100, "#ff7777"),
    ]

    all_db = np.concatenate([db_m0, db_m100])
    y_lo = float(np.percentile(all_db, 0.5)) - 5.0
    y_hi = float(np.max(all_db)) + 5.0

    for ax, title, db, stages_meta, color in cascade_panels:
        ax.set_facecolor("#0d0f10")
        ax.semilogx(freqs, db, color=color, lw=2.2)
        ax.set_title(title, color=color, fontsize=10, fontweight="bold")
        ax.set_ylim(y_lo, y_hi)
        ax.grid(True, which="both", alpha=0.18, color="white")
        ax.tick_params(colors="white", labelsize=8)
        for spine in ax.spines.values(): spine.set_color("#3a3e42")
        ax.set_ylabel("dB", color="white", fontsize=8)
        ax.axhline(0.0, color="#444", lw=0.5, alpha=0.5)
        for s in stages_meta:
            ax.axvline(s["pole_hz"], color=color, lw=0.4, alpha=0.4, linestyle=":")
            ax.axvline(s["zero_hz"], color="#22ddff", lw=0.4, alpha=0.3, linestyle="--")

    # Per-stage breakdown — colored by stage index, individual transfer functions
    stage_colors = ["#22ddff", "#88ee99", "#ffeb3b", "#ff9933", "#ff4488", "#cc66ff"]
    for ax, title, compiled_stages, meta_stages, _ in breakdown_panels:
        ax.set_facecolor("#0d0f10")
        for i, (cs, m) in enumerate(zip(compiled_stages, meta_stages)):
            db = stage_db(cs, freqs)
            ax.semilogx(freqs, db, color=stage_colors[i % len(stage_colors)], lw=1.4,
                        label=f"S{i}: pole {m['pole_hz']:.0f} Hz, c0={m['c0']:.2f}")
        ax.set_title(title, color="white", fontsize=10, fontweight="bold")
        ax.set_ylim(y_lo, y_hi)
        ax.grid(True, which="both", alpha=0.18, color="white")
        ax.tick_params(colors="white", labelsize=8)
        for spine in ax.spines.values(): spine.set_color("#3a3e42")
        ax.set_ylabel("dB", color="white", fontsize=8)
        ax.axhline(0.0, color="#444", lw=0.5, alpha=0.5)
        ax.legend(facecolor="#22262a", edgecolor="#444", labelcolor="white",
                  fontsize=7, loc="lower left")

    for ax in axes[-1, :]:
        ax.set_xlabel("Hz  (cascade: red dotted = pole, cyan dashed = zero)", color="white", fontsize=8)

    fig.tight_layout()
    PNG_PATH.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(PNG_PATH, facecolor=fig.get_facecolor())
    plt.close(fig)

    write_cartridge(stages_m0, stages_m100)
    print(f"plot      → {PNG_PATH}")
    print(f"cartridge → {JSON_PATH}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
