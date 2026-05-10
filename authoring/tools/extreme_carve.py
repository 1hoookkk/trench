"""forge/extreme_carve.py — prototype hostile-mode authoring.

Smarter than EMU silicon, calibrated by EMU silicon.

Design intent:
  - 6 active stages with poles distributed across the band (acoustic intent,
    not heritage role labels): 100 Hz, 380, 1200, 2800, 5600, 11000.
  - Pole radius = 0.998 across all stages (above heritage's 0.95–0.997 ceiling).
  - Zeros on the unit circle (zero_r = 1.0) at M0 — full cancellation, infinite
    notch depth, surgical carving silicon couldn't render in fixed-point DF2T.
  - Each zero offset +5 semitones above its pole, creating peak-notch pairs
    distributed across the band. Six peak-notch pairs = a multi-band scalpel.
  - c0 bisected per stage so cascade peak at M0 = +24 dB (calibrated by
    heritage discipline: hot but not a volume bomb).
  - M100: zeros relax to zero_r = 0.85, c0 unchanged. M0 → M100 = the carve
    softening from surgical to broad.

Output:
  cartridges/factory/generated/qlaw/_extreme_carve_proof.png  (4-panel plot)
  cartridges/factory/generated/qlaw/_extreme_carve_proof.json (cartridge)
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
CALIB_DIR = REPO / "reference" / "calibration"
OUT_DIR = REPO / "cartridges" / "factory" / "generated" / "qlaw"
PNG_PATH = OUT_DIR / "_extreme_carve_proof.png"
JSON_PATH = OUT_DIR / "_extreme_carve_proof.json"

AUTHORING_SR = 39062.5

# Pole frequencies — acoustic intent, NOT heritage roles.
# Real morph: pole frequencies SHIFT between M0 and M100.
# M0 = low growl cluster (80 Hz → 6 kHz); M100 = high shimmer (400 Hz → 15 kHz).
# Same six features but the whole spectrum slides up. Like Talking_Hedz where
# anchor drops 21% and HF stage rises 9% — pole movement IS the morph.
POLE_FREQS_HZ_M0   = [80.0,  200.0, 500.0,  1200.0, 2800.0, 6000.0]
POLE_FREQS_HZ_M100 = [400.0, 800.0, 1600.0, 3000.0, 5500.0, 11000.0]

POLE_RADIUS_M0 = 0.99        # tight resonance end of range
POLE_RADIUS_M100 = 0.98      # softer end of range — both stay 0.98-0.99
ZERO_OFFSET_M0_SEMITONES = -3
ZERO_OFFSET_M100_SEMITONES = -7
ZERO_RADIUS_M0 = 0.92        # interior: notches stay deep but pole peak survives
ZERO_RADIUS_M100 = 0.85      # softer notches at M100
TARGET_CASCADE_PEAK_DB_M0 = 24.0   # hostile but not a volume bomb


def stage_coefs(pole_hz: float, pole_r: float, zero_hz: float, zero_r: float,
                c0: float, sr: float = AUTHORING_SR) -> dict[str, float]:
    """One biquad stage's c0..c4 from (pole, zero, c0).

    The renderer evaluates num = c0 + c1·z⁻¹ + c2·z⁻². For c0 to act as a true
    scalar multiplier on the numerator polynomial (so zero locations are
    preserved regardless of attenuation), c1 and c2 must scale with c0:
        num = c0·(1 + b1·z⁻¹ + b2·z⁻²)
            = c0 + c0·b1·z⁻¹ + c0·b2·z⁻²
    so c1_stored = c0·b1, c2_stored = c0·b2 with b1 = -2·zero_r·cos(θ_z),
    b2 = zero_r². Without this scaling, attenuating c0 shifts the zeros out
    of their intended positions and breaks per-stage character.
    """
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


def cascade_peak_db(stages: list[dict], freqs: np.ndarray, sr: float) -> float:
    """Max magnitude of the multiplied cascade response, in dB."""
    w = 2.0 * np.pi * freqs / sr
    e1 = np.exp(-1j * w)
    e2 = np.exp(-2j * w)
    mag = np.ones_like(freqs)
    for s in stages:
        c0, c1, c2, c3, c4 = s["c0"], s["c1"], s["c2"], s["c3"], s["c4"]
        if c3 == 0.0 and c4 == 0.0:
            continue
        num = c0 + c1 * e1 + c2 * e2
        den = 1.0 + c3 * e1 + c4 * e2
        mag *= np.abs(num / den)
    return float(np.max(20.0 * np.log10(np.maximum(mag, 1e-12))))


def cascade_db(stages: list[dict], freqs: np.ndarray, sr: float) -> np.ndarray:
    w = 2.0 * np.pi * freqs / sr
    e1 = np.exp(-1j * w)
    e2 = np.exp(-2j * w)
    mag = np.ones_like(freqs)
    for s in stages:
        c0, c1, c2, c3, c4 = s["c0"], s["c1"], s["c2"], s["c3"], s["c4"]
        if c3 == 0.0 and c4 == 0.0:
            continue
        num = c0 + c1 * e1 + c2 * e2
        den = 1.0 + c3 * e1 + c4 * e2
        mag *= np.abs(num / den)
    return 20.0 * np.log10(np.maximum(mag, 1e-12))


def build_stages(pole_freqs_hz: list[float], pole_radius: float,
                 zero_radius: float, zero_offset_semitones: float,
                 c0_uniform: float, sr: float = AUTHORING_SR) -> list[dict]:
    """6 stages with shared (pole_r, zero_r, offset, c0) at given pole_freqs."""
    stages = []
    for pole_hz in pole_freqs_hz:
        zero_hz = pole_hz * (2.0 ** (zero_offset_semitones / 12.0))
        if zero_hz >= sr * 0.5:
            zero_hz = sr * 0.49
        if zero_hz <= 0.0:
            zero_hz = 1.0
        s = stage_coefs(
            pole_hz=pole_hz, pole_r=pole_radius,
            zero_hz=zero_hz, zero_r=zero_radius,
            c0=c0_uniform, sr=sr,
        )
        stages.append(s)
    while len(stages) < 12:
        stages.append({"c0": 1.0, "c1": 0.0, "c2": 0.0, "c3": 0.0, "c4": 0.0})
    return stages


def calibrate_c0_for_target_peak(pole_freqs_hz: list[float], pole_radius: float,
                                 zero_radius: float, zero_offset_semitones: float,
                                 target_peak_db: float,
                                 sr: float = AUTHORING_SR) -> float:
    """Bisection on c0 against the M0 frame to land cascade peak at target."""
    freqs = np.geomspace(20.0, 18000.0, 8192)
    lo, hi = 1e-6, 1.0
    for _ in range(80):
        mid = math.sqrt(lo * hi)
        stages = build_stages(pole_freqs_hz, pole_radius, zero_radius,
                              zero_offset_semitones, mid, sr)
        peak = cascade_peak_db(stages, freqs, sr)
        if peak > target_peak_db:
            hi = mid
        else:
            lo = mid
        if abs(peak - target_peak_db) < 0.1:
            return mid
    return math.sqrt(lo * hi)


def heritage_render(name: str, freqs: np.ndarray) -> tuple[np.ndarray, float, int]:
    """Load heritage M0_Q0, render via same composite math."""
    fp = CALIB_DIR / f"{name}.json"
    d = json.loads(fp.read_text(encoding="utf-8"))
    sr = float(d.get("sample_rate", AUTHORING_SR))
    boost = float(d.get("boost", 1.0))
    stages_in = d["corners"]["M0_Q0"]["stages"]
    stages = []
    for s in stages_in:
        radius = float(s.get("radius", 0.0))
        pole_hz = float(s.get("pole_freq_hz", 0.0))
        if radius == 0.0 or pole_hz == 0.0:
            stages.append({"c0": 1.0, "c1": 0.0, "c2": 0.0, "c3": 0.0, "c4": 0.0})
            continue
        theta = 2.0 * math.pi * pole_hz / sr
        z = s["zeros"]
        stages.append({
            "c0": 1.0 + float(s.get("val1", 0.0)),
            "c1": float(z.get("b1", 0.0)),
            "c2": float(z.get("b2", 0.0)),
            "c3": -2.0 * radius * math.cos(theta),
            "c4": radius * radius,
        })
    n_active = len(stages)
    while len(stages) < 12:
        stages.append({"c0": 1.0, "c1": 0.0, "c2": 0.0, "c3": 0.0, "c4": 0.0})
    db = cascade_db(stages, freqs, sr) + 20.0 * math.log10(boost)
    return db, sr, n_active


def write_cartridge(stages_m0: list[dict], stages_m100: list[dict],
                    c0_m0: float, c0_m100: float) -> None:
    cart = {
        "format": "compiled-v1",
        "name": "_extreme_carve_proof",
        "sampleRate": AUTHORING_SR,
        "recipe": "extreme_carve",
        "provenance": {
            "mode": "extreme_carve_prototype",
            "pole_freqs_hz_m0": POLE_FREQS_HZ_M0,
            "pole_freqs_hz_m100": POLE_FREQS_HZ_M100,
            "pole_radius_m0": POLE_RADIUS_M0,
            "pole_radius_m100": POLE_RADIUS_M100,
            "zero_offset_m0_semitones": ZERO_OFFSET_M0_SEMITONES,
            "zero_offset_m100_semitones": ZERO_OFFSET_M100_SEMITONES,
            "zero_radius_m0": ZERO_RADIUS_M0,
            "zero_radius_m100": ZERO_RADIUS_M100,
            "c0_m0": c0_m0,
            "c0_m100": c0_m100,
            "target_cascade_peak_db": TARGET_CASCADE_PEAK_DB_M0,
        },
        "keyframes": [
            {"label": "M0_Q0",   "morph": 0.0, "q": 0.0, "boost": 1.0, "stages": stages_m0},
            {"label": "M0_Q100", "morph": 0.0, "q": 1.0, "boost": 1.0, "stages": stages_m0},
            {"label": "M100_Q0", "morph": 1.0, "q": 0.0, "boost": 1.0, "stages": stages_m100},
            {"label": "M100_Q100","morph": 1.0,"q": 1.0, "boost": 1.0, "stages": stages_m100},
        ],
    }
    JSON_PATH.parent.mkdir(parents=True, exist_ok=True)
    JSON_PATH.write_text(json.dumps(cart, indent=2), encoding="utf-8")


def main() -> int:
    freqs = np.geomspace(20.0, 18000.0, 4096)

    # Calibrate M0 and M100 INDEPENDENTLY. Runtime LERPs c0 between corners,
    # so each corner can carry its own attenuation. Sharing c0 across morph
    # endpoints chokes whichever corner doesn't need the heavy attenuation.
    c0_m0 = calibrate_c0_for_target_peak(
        POLE_FREQS_HZ_M0, POLE_RADIUS_M0, ZERO_RADIUS_M0,
        ZERO_OFFSET_M0_SEMITONES, TARGET_CASCADE_PEAK_DB_M0,
    )
    c0_m100 = calibrate_c0_for_target_peak(
        POLE_FREQS_HZ_M100, POLE_RADIUS_M100, ZERO_RADIUS_M100,
        ZERO_OFFSET_M100_SEMITONES, TARGET_CASCADE_PEAK_DB_M0,
    )
    stages_m0 = build_stages(
        POLE_FREQS_HZ_M0, POLE_RADIUS_M0, ZERO_RADIUS_M0,
        ZERO_OFFSET_M0_SEMITONES, c0_m0,
    )
    stages_m100 = build_stages(
        POLE_FREQS_HZ_M100, POLE_RADIUS_M100, ZERO_RADIUS_M100,
        ZERO_OFFSET_M100_SEMITONES, c0_m100,
    )

    db_m0   = cascade_db(stages_m0,   freqs, AUTHORING_SR)
    db_m100 = cascade_db(stages_m100, freqs, AUTHORING_SR)
    peak_m0   = float(np.max(db_m0))
    peak_m100 = float(np.max(db_m100))

    print(f"calibrated c0 M0   = {c0_m0:.6f}  (val1 = {c0_m0 - 1:+.4f})")
    print(f"calibrated c0 M100 = {c0_m100:.6f}  (val1 = {c0_m100 - 1:+.4f})")
    print(f"cascade peak M0    = {peak_m0:+.2f} dB  (target {TARGET_CASCADE_PEAK_DB_M0:+.1f})")
    print(f"cascade peak M100  = {peak_m100:+.2f} dB  (target {TARGET_CASCADE_PEAK_DB_M0:+.1f})")

    fig, axes = plt.subplots(2, 1, figsize=(13.0, 8.5), dpi=140, sharex=True)
    fig.patch.set_facecolor("#16191a")
    fig.suptitle(
        "EXTREME_CARVE  ::  morph trajectory M0 → M100\n"
        f"pole_r {POLE_RADIUS_M0:.3f}→{POLE_RADIUS_M100:.3f},  "
        f"zero_r {ZERO_RADIUS_M0:.2f}→{ZERO_RADIUS_M100:.2f},  "
        f"offset {ZERO_OFFSET_M0_SEMITONES}→{ZERO_OFFSET_M100_SEMITONES} semis,  "
        f"c0  M0={c0_m0:.3f}  M100={c0_m100:.3f}",
        color="#ffba00", fontweight="bold", fontsize=13,
    )

    m0_sub  = f"M0   ::  pole_r={POLE_RADIUS_M0:.3f}, zero_r={ZERO_RADIUS_M0:.2f}, offset {ZERO_OFFSET_M0_SEMITONES} semis, c0={c0_m0:.3f},  cascade peak {peak_m0:+.1f} dB"
    m100_sub = f"M100 ::  pole_r={POLE_RADIUS_M100:.3f}, zero_r={ZERO_RADIUS_M100:.2f}, offset {ZERO_OFFSET_M100_SEMITONES} semis, c0={c0_m100:.3f},  cascade peak {peak_m100:+.1f} dB"
    panels = [
        (axes[0], m0_sub,  db_m0,   "#ff5555"),
        (axes[1], m100_sub, db_m100, "#ff7777"),
    ]
    # Share y-range across both panels.
    all_db = np.concatenate([db_m0, db_m100])
    y_lo = float(np.percentile(all_db, 0.5)) - 5.0
    y_hi = float(np.max(all_db)) + 5.0

    for ax, title, db, color in panels:
        ax.set_facecolor("#0d0f10")
        ax.semilogx(freqs, db, color=color, lw=2.2)
        ax.set_title(title, color=color, fontsize=11, fontweight="bold")
        ax.set_ylim(y_lo, y_hi)
        ax.grid(True, which="both", alpha=0.18, color="white")
        ax.tick_params(colors="white", labelsize=9)
        for spine in ax.spines.values():
            spine.set_color("#3a3e42")
        ax.set_ylabel("dB", color="white", fontsize=9)
        ax.axhline(0.0, color="#444", lw=0.5, alpha=0.5)
        # Mark the pole frequencies so peaks are obvious
        # Mark this panel's pole frequencies (different per panel now)
        pf_list = POLE_FREQS_HZ_M0 if "M0" in title and "M100" not in title else POLE_FREQS_HZ_M100
        for pf in pf_list:
            ax.axvline(pf, color=color, lw=0.4, alpha=0.35, linestyle=":")

    axes[1].set_xlabel("Hz", color="white", fontsize=9)

    fig.tight_layout()
    PNG_PATH.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(PNG_PATH, facecolor=fig.get_facecolor())
    plt.close(fig)

    write_cartridge(stages_m0, stages_m100, c0_m0, c0_m100)

    print(f"plot     → {PNG_PATH}")
    print(f"cartridge → {JSON_PATH}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
