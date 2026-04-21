"""Bake `Small Talk Ah-Ee` vocal body → compiled-v1 cartridge JSON.

Tier-1 authoring of a 5-active-stage vocal body that morphs between an
"Ah" formant configuration (M0, throat-open) and an "Ee" configuration
(M100, mouth-closed). Slot 6 is forced passthrough to satisfy the
"highest stage hard-converted" rule under the frozen 6-active topology.

## Mechanics

- Stages 1-3: Type 3 vocal sections with morphing freq bytes carrying
  the F1/F2/F3 formant trajectory. MorphLP table entry override on the
  numerator words (w0, w1, w4) gives per-corner zero variation
  (deep at M0, near-passthrough at M100) — the mathematical source of
  "talking" character.
- Stages 4-5: Type 3 anchor stages. Freq byte > 0xDB; per-stage shift
  set to -32 so the split-code frequency-compression hack clamps the
  pole frequency to exactly 220 (semitone 57). High-formant anchors
  stay locked across the morph regardless of UI position.
- Stage 6: Type 0 passthrough.

## Q-axis

Degenerate. The MorphLP zero values are corner-specific in morph but
collapsed in Q (Q0 and Q100 share the same numerator at each morph
endpoint). The pole skeleton from the byte-derived denominator is also
morph-only; Q is reflected in the runtime gain-offset alone.

## Output

- `cartridges/engine/character/small_talk_ah_ee.json` — engine-side.
- `trench-juce/plugin/assets/cartridges/small_talk_ah_ee.json` — plugin
  asset mirror (the JUCE side ships from this directory).

## Provenance

Doctrine paste 2026-04-21 + MorphLP table sanitized via
`tools/morphlp_zeros.py`. The Type 3 firmware compile recipe is
identical to the proven Talking Hedz path (`tools/bake_hedz_const.py`).

Stdlib only (matches the discipline of bake_hedz_const.py).
"""
from __future__ import annotations

import argparse
import json
import struct
from pathlib import Path

from morphlp_zeros import morphlp_corner_words

REPO = Path(__file__).resolve().parents[1]
ENGINE_OUT = REPO / "cartridges" / "engine" / "character" / "small_talk_ah_ee.json"
PLUGIN_OUT = REPO / "trench-juce" / "plugin" / "assets" / "cartridges" / "small_talk_ah_ee.json"

NUM_PHYSICAL_STAGES = 12
NUM_ACTIVE_SLOTS = 6
SR_DECLARED = 39062.5
SR_COMPILE = 44100  # see memory: compiled-v1 sampleRate is stale metadata
BOOST = 8.0  # Type 3 output is quiet without numerator boost; Talking Hedz uses 4.0
             # but our 5-stage with anchored highs runs ~6 dB lower at M100. 8.0
             # lands M100 peak near -10 dBFS for impulse, leaves audition headroom.

# ---------------------------------------------------------------------------
# Firmware compile (mirror of tools/bake_hedz_const.py — keep in sync)
# ---------------------------------------------------------------------------

SR_FAMILY = {44100: 0, 48000: 1, 96000: 2, 192000: 3}
FW_BASE = [18, 18, 4, 1]
FW_SCALE = [220, 220, 200, 177]
COMBINE_K = 4.0


def f32(x: float) -> float:
    return struct.unpack("<f", struct.pack("<f", x))[0]


def mf_decode(packed: int) -> float:
    packed &= 0xFFFF
    if packed == 0xFFFF:
        return 1.0
    if packed == 0x0000:
        return 0.0
    u = packed + 1
    exp = ((u >> 12) & 0xF) - 15
    mant = u & 0xFFF
    if exp < -14:
        return mant / 134217728.0
    scaled = (mant | 0x1000) / 8192.0
    return scaled * (2.0 ** exp)


def fw_words_to_kernel(w0: int, w1: int, w2: int, w3: int, w4: int) -> tuple[float, float, float, float, float]:
    d0 = mf_decode(w0)
    d1 = mf_decode(w1)
    d2 = mf_decode(w2)
    d3 = mf_decode(w3)
    d4 = mf_decode(w4)
    return (
        d0 * COMBINE_K + d1,
        d1,
        d2 * COMBINE_K + d3,
        d3,
        d4,
    )


def fw_freq_value(packed_freq: int, sr: int = SR_COMPILE) -> int:
    idx = SR_FAMILY.get(sr, 0)
    return (FW_SCALE[idx] * packed_freq) // 128 + FW_BASE[idx]


def fw_radius(freq_value: int) -> int:
    return (freq_value * 124) // 256 + 118


def fw_gain_offset(packed_gain: int, global_shift: int = 0) -> int:
    raw = (packed_gain - 64) // 2 + global_shift
    return max(-32, min(31, raw))


def type3_freq_compression(freq_value: int, shift: int) -> int:
    if freq_value > 0xDB and shift < 0:
        freq_value = (((freq_value - 220) * (shift + 32)) >> 5) + 220
    return freq_value


def type3_vocal_kernel(
    freq_packed: int,
    gain_packed: int,
    shift: int,
    sr: int,
    morphlp_idx: int | None,
    corner: str,
) -> tuple[float, float, float, float, float]:
    """Type 3 vocal kernel with optional MorphLP numerator override.

    When `morphlp_idx` is None the numerator words match the stock Type 3
    pattern (passthrough-shaped numerator from `bake_hedz_const.py`). When
    set, the (w0, w1, w4) numerator words are replaced with the MorphLP
    table values for the given corner — driving the per-corner zero
    dissolution that produces the "talking" character.
    """
    idx = SR_FAMILY.get(sr, 0)
    fv = fw_freq_value(freq_packed, sr)
    go = fw_gain_offset(gain_packed, shift)
    fv_compressed = type3_freq_compression(fv, shift)
    rad = fw_radius(fv_compressed)
    if morphlp_idx is not None:
        v1, v2, v3 = morphlp_corner_words(morphlp_idx, corner)
        w0 = v1
        w1 = v2
        w4 = v3
    else:
        base = FW_BASE[idx]
        w0 = base << 8
        w1 = ((base * 124) // 256 + 150) << 8
        if idx < 2:
            w4 = (fv_compressed - 18) * (-12) + (-8192)
            w4 &= 0xFFFF
        else:
            w4 = 0xE000
    w2 = fv_compressed << 8
    w3 = max(0, min(255, rad - go)) << 8
    return fw_words_to_kernel(w0, w1, w2, w3, w4)


PASSTHROUGH = (1.0, 0.0, 0.0, 0.0, 0.0)


# ---------------------------------------------------------------------------
# Body spec — Small Talk Ah-Ee
# ---------------------------------------------------------------------------

# Each section: (type_id, low_freq_byte, low_gain_byte, high_freq_byte,
# high_gain_byte, per_stage_shift, morphlp_idx_or_None).
#
# Formant trajectory (adult-male phonetics):
#   Ah (/ɑ/, M0):  F1≈730   F2≈1090  F3≈2440
#   Ee (/i/, M100): F1≈270   F2≈2290  F3≈3010
# Byte values are first-pass; expect audition-driven retune.
#
# Stages 4-5 are anchors: freq bytes > 0xDB and shift = -32 → split-code
# clamp pins the pole frequency at exactly 220. They give the dual-formant
# "anchor points" that make the motion parse as biomechanical vocal cords.
#
# MorphLP table index choices:
#   - moving formants (1-3): idx=4 — strong M0 divergence, M100 nearer to
#     dissolution, mid-deep zeros. Gives clear vocal motion without going
#     to extreme cancellation.
#   - anchors (4-5): idx=8 — saturated zone where M0/M100 zeros are
#     converging. Anchors sit still both spectrally and zero-wise.
SECTIONS = [
    # (type, lo_f, lo_g, hi_f, hi_g, shift, morphlp_idx)
    # MorphLP override DISABLED: Type 3 (FUN_1802c3600) and CPhantomMorphLP
    # (FUN_1802c59b0) are different firmware compile paths with their own
    # numerator/denominator polynomial coupling. Forcing MorphLP zeros onto
    # Type 3 poles broke gain balance by 90+ dB across all 16 indices.
    # `tools/morphlp_zeros.py` stays in place for future MorphLP-native bodies.
    (3,  64,  80,  32,  80,   0, None),  # Stage 1 — F1
    (3,  80,  90, 160,  90,   0, None),  # Stage 2 — F2 (talking crossover)
    (3, 140,  80, 170,  80,   0, None),  # Stage 3 — F3
    (3, 224,  96, 224,  96, -32, None),  # Stage 4 — anchor F4 (split-code clamp)
    (3, 240,  96, 240,  96, -32, None),  # Stage 5 — anchor F5 (split-code clamp)
    (0,   0,   0,   0,   0,   0, None),  # Stage 6 — forced passthrough
]


def compile_corner(morph: float, corner_label: str, sr: int) -> list[tuple[float, ...]]:
    """Compile the 6 active stages for one corner. Q-axis is degenerate so
    `corner_label` selects M0 vs M100 zeros from the MorphLP table; the
    Q distinction (Q0 vs Q100) does not affect the numerator."""
    if morph == 0.0:
        zero_corner = "M0"
    elif morph == 1.0:
        zero_corner = "M100"
    else:
        raise ValueError("compile_corner expects morph endpoint (0.0 or 1.0)")
    stages: list[tuple[float, ...]] = []
    for type_id, lo_f, lo_g, hi_f, hi_g, shift, mlp_idx in SECTIONS:
        if type_id == 0:
            stages.append(PASSTHROUGH)
            continue
        freq = round(lo_f + morph * (hi_f - lo_f))
        gain = round(lo_g + morph * (hi_g - lo_g))
        coeffs = type3_vocal_kernel(freq, gain, shift, sr, mlp_idx, zero_corner)
        stages.append(tuple(f32(c) for c in coeffs))
    while len(stages) < NUM_PHYSICAL_STAGES:
        stages.append(PASSTHROUGH)
    return stages


def build_cartridge() -> dict:
    m0_stages = compile_corner(0.0, "M0", SR_COMPILE)
    m100_stages = compile_corner(1.0, "M100", SR_COMPILE)
    keyframes = [
        {"label": "M0_Q0",     "morph": 0.0, "q": 0.0, "boost": BOOST,
         "stages": [stage_dict(s) for s in m0_stages]},
        {"label": "M100_Q0",   "morph": 1.0, "q": 0.0, "boost": BOOST,
         "stages": [stage_dict(s) for s in m100_stages]},
        {"label": "M0_Q100",   "morph": 0.0, "q": 1.0, "boost": BOOST,
         "stages": [stage_dict(s) for s in m0_stages]},
        {"label": "M100_Q100", "morph": 1.0, "q": 1.0, "boost": BOOST,
         "stages": [stage_dict(s) for s in m100_stages]},
    ]
    return {
        "format": "compiled-v1",
        "name": "Small Talk Ah-Ee",
        "provenance": "tier1-vocal-morphlp-2026-04-21",
        "sampleRate": SR_DECLARED,
        "stages": NUM_PHYSICAL_STAGES,
        "keyframes": keyframes,
    }


def stage_dict(coeffs: tuple[float, ...]) -> dict:
    c0, c1, c2, c3, c4 = coeffs
    return {"c0": c0, "c1": c1, "c2": c2, "c3": c3, "c4": c4}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Bake Small Talk Ah-Ee body.")
    parser.add_argument("--engine-out", type=Path, default=ENGINE_OUT)
    parser.add_argument("--plugin-out", type=Path, default=PLUGIN_OUT)
    args = parser.parse_args(argv)

    cart = build_cartridge()
    payload = json.dumps(cart, indent=2)

    args.engine_out.parent.mkdir(parents=True, exist_ok=True)
    args.engine_out.write_text(payload + "\n", encoding="utf-8")
    print(f"wrote {args.engine_out}")

    if args.plugin_out is not None:
        args.plugin_out.parent.mkdir(parents=True, exist_ok=True)
        args.plugin_out.write_text(payload + "\n", encoding="utf-8")
        print(f"wrote {args.plugin_out}")

    print(f"sections: {len(SECTIONS)} authored, "
          f"{sum(1 for s in SECTIONS if s[0] != 0)} active vocal, "
          f"{sum(1 for s in SECTIONS if s[0] == 0)} passthrough")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
