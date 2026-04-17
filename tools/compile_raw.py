#!/usr/bin/env python3
"""Compile an internal raw stage surface into a compiled-v1 cartridge.

Input format: raw-stage-v1 JSON with either:
  - `corners: {M0_Q0: {...}, ...}`
  - `keyframes: [{label, stages, ...}, ...]`

Each stage is authored in musical units:
  - `allpole` / `resonator`: pole frequency + radius + stage gain. A textbook
    2-pole resonator with no numerator zeros: `b = (g, 0, 0)`. This is the
    firmware-native pole-only shape (Morpheus cubes per E-mu RE).
  - `zero_forced`: pole frequency + radius + stage gain, zero on unit circle
    at the pole frequency
  - `zero_forced_offset`: pole frequency + radius + stage gain, unit-circle
    zero offset by semitones
  - `explicit_zero`: pole frequency + radius + stage gain + explicit zero
    frequency/radius
  - `passthrough`

Output: compiled-v1 cartridge JSON that loads through
`trench_core::Cartridge::from_json`.

Stdlib only. No RBJ derivation. This is a deterministic format bridge from an
internal authoring surface to kernel-form coefficients.
"""

from __future__ import annotations

import argparse
import json
import math
import struct
import sys
from pathlib import Path

AUTHORING_SR = 39062.5
NYQUIST_LIMIT_HZ = AUTHORING_SR * 0.5
CORNER_ORDER = ("M0_Q0", "M0_Q100", "M100_Q0", "M100_Q100")
CORNER_MORPH_Q = {
    "M0_Q0": (0.0, 0.0),
    "M0_Q100": (0.0, 1.0),
    "M100_Q0": (1.0, 0.0),
    "M100_Q100": (1.0, 1.0),
}
PASSTHROUGH = {"c0": 1.0, "c1": 0.0, "c2": 0.0, "c3": 0.0, "c4": 0.0}
ACTIVE_STAGE_LIMIT = 6
OUTPUT_STAGE_COUNT = 12


class CompileError(Exception):
    pass


def f32(value: float) -> float:
    return struct.unpack("<f", struct.pack("<f", float(value)))[0]


def require_finite(value: object, field: str) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError) as exc:
        raise CompileError(f"{field} must be numeric") from exc
    if not math.isfinite(out):
        raise CompileError(f"{field} must be finite")
    return out


def require_range(value: float, field: str, low: float, high: float, high_open: bool = False) -> float:
    if value < low:
        raise CompileError(f"{field} must be >= {low}")
    if high_open:
        if value >= high:
            raise CompileError(f"{field} must be < {high}")
    elif value > high:
        raise CompileError(f"{field} must be <= {high}")
    return value


def pole_a1(freq_hz: float, radius: float) -> float:
    theta = 2.0 * math.pi * freq_hz / AUTHORING_SR
    return -2.0 * radius * math.cos(theta)


def validate_pole_stage(stage: dict, stage_path: str) -> tuple[float, float, float]:
    freq_hz = require_range(
        require_finite(stage.get("pole_freq_hz"), f"{stage_path}.pole_freq_hz"),
        f"{stage_path}.pole_freq_hz",
        0.0,
        NYQUIST_LIMIT_HZ,
        high_open=True,
    )
    radius = require_range(
        require_finite(stage.get("radius"), f"{stage_path}.radius"),
        f"{stage_path}.radius",
        0.0,
        1.0,
        high_open=True,
    )
    stage_gain = require_finite(stage.get("stage_gain", 1.0), f"{stage_path}.stage_gain")
    return freq_hz, radius, stage_gain


def encode_stage(stage: dict, stage_path: str) -> dict:
    kind = str(stage.get("kind", "allpole"))
    if kind == "passthrough":
        return PASSTHROUGH.copy()

    if kind in ("allpole", "resonator", "zero_forced", "zero_forced_offset", "explicit_zero"):
        pole_freq_hz, radius, stage_gain = validate_pole_stage(stage, stage_path)
    else:
        raise CompileError(
            f"{stage_path}.kind must be one of passthrough, allpole, resonator, "
            f"zero_forced, zero_forced_offset, explicit_zero"
        )

    a1 = pole_a1(pole_freq_hz, radius)
    val1 = stage_gain - 1.0

    # Firmware-faithful pole-only encoding. Morpheus cubes in the E-mu
    # firmware (per the RE decode) set val1 = val2 = val3 = 0, which
    # produces b = (1, 0, 0): a textbook 2-pole resonator with no zeros
    # in the numerator. With stage_gain g applied via val1 = g-1, this
    # becomes b = (g, 0, 0). Prior to 2026-04-17 this kind wrote
    # b1_target = a1 and b2_target = r² (numerator mirroring
    # denominator), a form not documented anywhere in the firmware RE
    # notes and one that produced massive DC gain instead of a
    # resonance peak at the authored pole frequency.
    if kind in ("allpole", "resonator"):
        b1_target = 0.0
        b2_target = 0.0
    elif kind == "zero_forced":
        b1_target = -2.0 * math.cos(2.0 * math.pi * pole_freq_hz / AUTHORING_SR)
        b2_target = 1.0
    elif kind == "zero_forced_offset":
        offset = require_finite(stage.get("offset_semitones"), f"{stage_path}.offset_semitones")
        zero_freq_hz = pole_freq_hz * (2.0 ** (offset / 12.0))
        zero_freq_hz = max(20.0, min(AUTHORING_SR * 0.48, zero_freq_hz))
        b1_target = -2.0 * math.cos(2.0 * math.pi * zero_freq_hz / AUTHORING_SR)
        b2_target = 1.0
    else:
        zero_freq_hz = require_range(
            require_finite(stage.get("zero_freq_hz"), f"{stage_path}.zero_freq_hz"),
            f"{stage_path}.zero_freq_hz",
            0.0,
            NYQUIST_LIMIT_HZ,
            high_open=True,
        )
        zero_radius = require_range(
            require_finite(stage.get("zero_radius"), f"{stage_path}.zero_radius"),
            f"{stage_path}.zero_radius",
            0.0,
            1.0,
        )
        b1_target = -2.0 * zero_radius * math.cos(2.0 * math.pi * zero_freq_hz / AUTHORING_SR)
        b2_target = zero_radius * zero_radius

    val2 = b1_target - a1
    val3 = radius * radius - b2_target

    a1 = f32(a1)
    radius = f32(radius)
    val1 = f32(val1)
    val2 = f32(val2)
    val3 = f32(val3)

    c0 = 1.0 + val1
    c1 = a1 + val2
    c2 = radius * radius - val3
    c3 = a1
    c4 = radius * radius
    return {"c0": c0, "c1": c1, "c2": c2, "c3": c3, "c4": c4}


def normalize_corner_map(raw_doc: dict) -> dict[str, dict]:
    if "corners" in raw_doc:
        corners = raw_doc["corners"]
        if not isinstance(corners, dict):
            raise CompileError("'corners' must be an object")
        out = {}
        for label in CORNER_ORDER:
            if label not in corners:
                raise CompileError(f"raw doc missing corner '{label}'")
            corner = corners[label]
            if isinstance(corner, list):
                out[label] = {"stages": corner}
            elif isinstance(corner, dict):
                out[label] = corner
            else:
                raise CompileError(
                    f"corner '{label}' must be an object or a stage list"
                )
        return out

    if "keyframes" in raw_doc:
        keyframes = raw_doc["keyframes"]
        if not isinstance(keyframes, list):
            raise CompileError("'keyframes' must be a list")
        out = {}
        for idx, keyframe in enumerate(keyframes):
            if not isinstance(keyframe, dict):
                raise CompileError(f"keyframe[{idx}] must be an object")
            label = keyframe.get("label")
            if label not in CORNER_ORDER:
                raise CompileError(
                    f"keyframe[{idx}].label must be one of {', '.join(CORNER_ORDER)}"
                )
            out[label] = keyframe
        missing = [label for label in CORNER_ORDER if label not in out]
        if missing:
            raise CompileError(f"raw doc missing keyframes: {missing}")
        return out

    raise CompileError("raw doc must have 'corners' or 'keyframes'")


def compile_raw(raw_doc: dict) -> dict:
    raw_format = raw_doc.get("format")
    if raw_format != "raw-stage-v1":
        raise CompileError(
            f"unsupported raw format: {raw_format!r}, expected 'raw-stage-v1'"
        )

    sample_rate_fields = []
    for field in ("authoring_sample_rate_hz", "sampleRate"):
        if field in raw_doc:
            sample_rate_fields.append((field, require_finite(raw_doc[field], field)))

    if len(sample_rate_fields) == 2:
        (_, first_rate), (_, second_rate) = sample_rate_fields
        if abs(first_rate - second_rate) > 1e-9:
            raise CompileError(
                "authoring_sample_rate_hz and sampleRate must match when both are present"
            )

    sample_rate = sample_rate_fields[0][1] if sample_rate_fields else AUTHORING_SR
    if abs(sample_rate - AUTHORING_SR) > 1e-9:
        raise CompileError(
            f"authoring sample rate must be {AUTHORING_SR} for the raw authoring surface"
        )

    global_boost = require_finite(raw_doc.get("boost", 1.0), "boost")
    corner_map = normalize_corner_map(raw_doc)

    keyframes = []
    for label in CORNER_ORDER:
        corner = corner_map[label]
        stages = corner.get("stages", [])
        if not isinstance(stages, list):
            raise CompileError(f"{label}.stages must be a list")
        if len(stages) > ACTIVE_STAGE_LIMIT:
            raise CompileError(
                f"{label}.stages has {len(stages)} entries, expected <= {ACTIVE_STAGE_LIMIT}"
            )

        compiled_stages = [
            encode_stage(stage, f"{label}.stages[{idx}]")
            for idx, stage in enumerate(stages)
        ]
        while len(compiled_stages) < ACTIVE_STAGE_LIMIT:
            compiled_stages.append(PASSTHROUGH.copy())
        while len(compiled_stages) < OUTPUT_STAGE_COUNT:
            compiled_stages.append(PASSTHROUGH.copy())

        boost = require_finite(corner.get("boost", global_boost), f"{label}.boost")
        morph, q = CORNER_MORPH_Q[label]
        keyframes.append(
            {
                "label": label,
                "morph": morph,
                "q": q,
                "boost": boost,
                "stages": compiled_stages,
            }
        )

    return {
        "format": "compiled-v1",
        "name": raw_doc.get("name", "unnamed_body"),
        "provenance": "compile_raw",
        "sampleRate": AUTHORING_SR,
        "keyframes": keyframes,
    }


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(
        description="Compile an internal raw-stage-v1 surface into compiled-v1 JSON."
    )
    parser.add_argument(
        "raw",
        type=Path,
        help="Path to raw-stage-v1 JSON file, or '-' for stdin.",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=None,
        help="Output cartridge path (stdout if omitted).",
    )
    args = parser.parse_args(argv)

    try:
        if str(args.raw) == "-":
            raw_doc = json.load(sys.stdin)
        else:
            with args.raw.open() as f:
                raw_doc = json.load(f)
    except json.JSONDecodeError as exc:
        print(f"raw parse error: {exc}", file=sys.stderr)
        return 1
    except OSError as exc:
        print(f"raw read error: {exc}", file=sys.stderr)
        return 1

    try:
        cart = compile_raw(raw_doc)
    except CompileError as exc:
        print(f"compile error: {exc}", file=sys.stderr)
        return 2

    out_json = json.dumps(cart, indent=2)
    if args.out is None:
        sys.stdout.write(out_json + "\n")
    else:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(out_json + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
