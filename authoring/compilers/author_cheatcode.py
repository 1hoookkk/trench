"""Cheatcode body generator for the playable TRENCH plugin.

Emits `compiled-v1` keyframe JSON in plain-DF2T kernel form — the shape the
shipping runtime actually decodes (`trench-core/src/cartridge.rs`,
`cascade.rs`). Numerator = [c0, c1, c2]; Denominator = [1, c3, c4].
Passthrough stage = (1, 0, 0, 0, 0), matching `PASSTHROUGH_COEFFS`.

Previous version of this file imported `pyruntime` from a sibling repo and
produced minifloat-encoded coefficients (passthrough = (2, 1, 2, 1, 1)) that
the shipping loader rejects. This rewrite is self-contained, uses runtime-
native math, and writes JSON the plugin will load without complaint.

Preserves the five candidates Gemini drafted — Scissors / Sub_Tear / Chase /
Vocal_Rip / Collapse — by translating the intent (pole/notch frequency sets,
Q stress via radius push) directly into stable DF2T coefficients.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

# Shipping runtime uses host-SR kernel coefficients (per
# memory:compiled_v1_sample_rate — hardware-parity proven at 44100 Hz).
SR = 44100.0

REPO = Path(__file__).resolve().parents[1]
OUT_DIR = REPO / "looperator" / "cheatcode_presets"

# Corner labels in the order trench-core expects for bilinear interpolation
# (see cartridge.rs::from_keyframe_json). Q-axis first, then morph.
CORNER_LABELS = ("M0_Q0", "M100_Q0", "M0_Q100", "M100_Q100")

# Gemini's design: 4 meaningful stages + 2 passthrough per corner = 6 total.
# Runtime caps active stages at NUM_STAGES=6, so 6-stage form (not 12) is
# exactly what the cascade consumes — no padding gymnastics.
NUM_ACTIVE = 4
NUM_STAGES = 6

PASSTHROUGH = [1.0, 0.0, 0.0, 0.0, 0.0]


def clamp(v: float, lo: float, hi: float) -> float:
    return max(lo, min(hi, v))


def stage_coeffs(
    pole_hz: float,
    pole_r: float,
    *,
    zero_hz: float | None = None,
    zero_r: float = 1.0,
    gain: float = 1.0,
) -> list[float]:
    """Plain-DF2T biquad → [c0, c1, c2, c3, c4].

    - Pole pair at (pole_hz, pole_r) → denominator [1, a1, a2].
    - Optional zero pair at (zero_hz, zero_r) → numerator [b0, b1, b2].
    - Omitting zero_hz gives an all-pole resonator (b1=b2=0).
    """
    pole_hz = clamp(pole_hz, 20.0, SR * 0.48)
    pole_r = clamp(pole_r, 0.0, 0.999)  # never on/outside unit circle

    theta = 2.0 * math.pi * pole_hz / SR
    a1 = -2.0 * pole_r * math.cos(theta)
    a2 = pole_r * pole_r

    if zero_hz is None:
        b0, b1, b2 = gain, 0.0, 0.0
    else:
        zero_hz = clamp(zero_hz, 20.0, SR * 0.48)
        phi = 2.0 * math.pi * zero_hz / SR
        b0 = gain
        b1 = -2.0 * zero_r * math.cos(phi) * gain
        b2 = (zero_r * zero_r) * gain

    return [b0, b1, b2, a1, a2]


def build_corner(
    p1_hz: float,
    p2_hz: float,
    n1_hz: float,
    n2_hz: float,
    q_boost: float,
) -> list[list[float]]:
    """One 6-stage corner: two all-pole resonators, two pole+zero notches.

    q_boost pushes pole radii closer to the unit circle for the Q100 corners.
    Stays below 0.995 so we never blow up.

    Per-stage gains come from Gemini's design (via pyruntime val1): throat/
    bite poles peak at b0 ≈ 0.68; notch stages use b0 ≈ 0.085 so the unit-
    circle null dominates over numerator scale.
    """
    GAIN_POLE = 0.68   # = 1 + (-1 + 0.8 * 0.85)
    GAIN_ZERO = 0.085  # = 1 + (-1 + 0.1 * 0.85)
    OFFSET_SEMITONES = 2.5

    stages: list[list[float]] = [
        # 1. Throat pole (massive low-band resonance)
        stage_coeffs(p1_hz, min(0.995, 0.95 + q_boost), gain=GAIN_POLE),
        # 2. Bite pole
        stage_coeffs(p2_hz, min(0.995, 0.92 + q_boost), gain=GAIN_POLE),
        # 3. Unit-circle null co-located with a pole at n1 → surgical notch
        stage_coeffs(
            n1_hz, min(0.98, 0.85 + q_boost),
            zero_hz=n1_hz, zero_r=1.0, gain=GAIN_ZERO,
        ),
        # 4. Offset notch: zero at n2 * 2^(2.5/12), pole at n2 → phase tear
        stage_coeffs(
            n2_hz, min(0.95, 0.80 + q_boost),
            zero_hz=n2_hz * (2.0 ** (OFFSET_SEMITONES / 12.0)),
            zero_r=1.0, gain=GAIN_ZERO,
        ),
        # 5-6: passthrough (unused active-stage budget)
        PASSTHROUGH,
        PASSTHROUGH,
    ]
    assert len(stages) == NUM_STAGES
    return stages


def build_cheat(name: str, m0_freqs: list[float], m1_freqs: list[float]) -> dict:
    """Four keyframes (M0_Q0, M100_Q0, M0_Q100, M100_Q100) from two frame recipes.

    Q-axis stress = +0.03 radius push per Gemini's original q_boost.
    """
    frames = {
        "M0_Q0":    (m0_freqs, 0.0),
        "M100_Q0":  (m1_freqs, 0.0),
        "M0_Q100":  (m0_freqs, 0.03),
        "M100_Q100":(m1_freqs, 0.03),
    }
    keyframes = []
    for label in CORNER_LABELS:
        freqs, q_boost = frames[label]
        keyframes.append({
            "label": label,
            "boost": 4.0,
            "stages": [
                {"c0": s[0], "c1": s[1], "c2": s[2], "c3": s[3], "c4": s[4]}
                for s in build_corner(*freqs, q_boost=q_boost)
            ],
        })
    return {
        "format": "compiled-v1",
        "name": name,
        "provenance": "author_cheatcode.py (plain-DF2T, runtime-native)",
        "sampleRate": SR,
        "stages": 12,
        "keyframes": keyframes,
    }


def audit_cartridge(cart: dict) -> list[str]:
    """Return a list of problems; empty list = clean."""
    issues: list[str] = []
    for kf in cart["keyframes"]:
        for si, s in enumerate(kf["stages"]):
            c0, c1, c2, c3, c4 = s["c0"], s["c1"], s["c2"], s["c3"], s["c4"]
            # Stability: |pole radius|² < 1  ⇔  |c4| < 1
            if abs(c4) >= 1.0:
                issues.append(
                    f"{kf['label']} stage {si+1}: c4={c4:.4f} ≥ 1 (unstable pole)"
                )
            # Finite
            for i, v in enumerate((c0, c1, c2, c3, c4)):
                if not math.isfinite(v):
                    issues.append(f"{kf['label']} stage {si+1}: c{i} not finite")
    return issues


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    candidates = [
        # name                    m0_freqs                m1_freqs
        ("Cheat_01_Scissors",     [400, 1200, 800, 2400], [1200, 3200, 200, 600]),
        ("Cheat_02_Sub_Tear",     [80, 2000, 600, 4000],  [120, 2400, 80, 2000]),
        ("Cheat_03_Chase",        [300, 900, 450, 1350],  [1000, 3000, 1500, 4500]),
        ("Cheat_04_Vocal_Rip",    [600, 2400, 200, 8000], [600, 2400, 2000, 400]),
        ("Cheat_05_Collapse",     [200, 6000, 50, 12000], [1000, 1500, 800, 1800]),
    ]

    print(f"Authoring cheatcode presets → {OUT_DIR}")
    print(f"  Encoding: plain-DF2T, SR={SR} Hz, passthrough=(1,0,0,0,0)\n")

    total_issues = 0
    for name, m0, m1 in candidates:
        cart = build_cheat(name, m0, m1)
        issues = audit_cartridge(cart)
        out = OUT_DIR / f"{name}.json"
        out.write_text(json.dumps(cart, indent=2), encoding="utf-8")
        status = "✓ clean" if not issues else f"✗ {len(issues)} issue(s)"
        print(f"  {name}: {status}")
        for issue in issues:
            print(f"    - {issue}")
        total_issues += len(issues)

    print(f"\n{len(candidates)} presets written; {total_issues} total issues.")
    return 0 if total_issues == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
