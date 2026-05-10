#!/usr/bin/env python3
"""Generate original 6-stage 2-frame bodies from bark-space sonic targets."""

from __future__ import annotations

import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

_TOOLS = Path(__file__).resolve().parent
if str(_TOOLS) not in sys.path:
    sys.path.insert(0, str(_TOOLS))

from compile_raw import compile_raw  # noqa: E402


ROOT = Path(__file__).resolve().parent.parent
CAL_INDEX = ROOT / "docs" / "calibration" / "index.json"
SONIC_TABLES = ROOT / "docs" / "sonic_tables" / "tables.json"


@dataclass(frozen=True)
class Experiment:
    name: str
    start_vowel: str
    end_vowel: str
    anchor_landmark: str
    air_landmark: str
    boost: float
    notes: str


EXPERIMENTS: tuple[Experiment, ...] = (
    Experiment(
        "bark_ee_to_ah",
        "ee",
        "ah",
        "chest_resonance",
        "air_band",
        1.8,
        "vowel arc from tight front vowel to open dark mouth",
    ),
    Experiment(
        "bark_oo_to_ee_belch",
        "oo",
        "ee",
        "room_mode_medium",
        "presence_peak",
        1.9,
        "dark rounded opening into aggressive upper-mid bark step",
    ),
    Experiment(
        "bark_er_to_aw_horn",
        "er",
        "aw",
        "chest_resonance",
        "ear_canal_resonance",
        1.7,
        "mid-brass throat with wider bark migration than the vowel pair",
    ),
)


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _bark_from_hz(freq_hz: float) -> float:
    return 26.81 / (1.0 + 1960.0 / max(freq_hz, 1e-6)) - 0.53


def _hz_from_bark(bark: float) -> float:
    return 1960.0 / (26.81 / (bark + 0.53) - 1.0)


def _nearest_semitone(freq_hz: float, semitone_table: list[float]) -> float:
    return min(semitone_table, key=lambda value: abs(value - freq_hz))


def _snap_bark(freq_hz: float, semitone_table: list[float]) -> float:
    return _nearest_semitone(max(40.0, min(freq_hz, 18000.0)), semitone_table)


def _blend_log(freq_a_hz: float, freq_b_hz: float, t: float) -> float:
    a = math.log(max(freq_a_hz, 1e-6))
    b = math.log(max(freq_b_hz, 1e-6))
    return math.exp(a * (1.0 - t) + b * t)


def _interior_zero(pole_freq_hz: float, zero_bark_offset: float, semitone_table: list[float]) -> float:
    zero_bark = _bark_from_hz(pole_freq_hz) + zero_bark_offset
    return _snap_bark(_hz_from_bark(zero_bark), semitone_table)


def _stage(kind: str, pole_freq_hz: float, radius: float, stage_gain: float, *, zero_freq_hz: float | None = None, zero_radius: float = 1.0, offset_semitones: float | None = None) -> dict[str, Any]:
    out: dict[str, Any] = {
        "kind": kind,
        "pole_freq_hz": pole_freq_hz,
        "radius": radius,
        "stage_gain": stage_gain,
    }
    if zero_freq_hz is not None:
        out["zero_freq_hz"] = zero_freq_hz
        out["zero_radius"] = zero_radius
    if offset_semitones is not None:
        out["offset_semitones"] = offset_semitones
    return out


def _dark_frame(
    low_body: float,
    low_mid: float,
    emergence: float,
    bite: float,
    presence: float,
    cap: float,
    semitone_table: list[float],
) -> list[dict[str, Any]]:
    return [
        _stage("allpole", low_body, 0.972, 0.94),
        _stage(
            "explicit_zero",
            low_mid,
            0.954,
            0.78,
            zero_freq_hz=_snap_bark(low_mid * 1.55, semitone_table),
            zero_radius=0.84,
        ),
        _stage("allpole", emergence, 0.966, 0.88),
        _stage(
            "explicit_zero",
            bite,
            0.956,
            0.76,
            zero_freq_hz=_snap_bark(bite * 1.12, semitone_table),
            zero_radius=0.80,
        ),
        _stage(
            "explicit_zero",
            presence,
            0.938,
            0.62,
            zero_freq_hz=_snap_bark(presence * 1.08, semitone_table),
            zero_radius=0.74,
        ),
        _stage("zero_forced_offset", cap, 0.904, 0.48, offset_semitones=2.0),
    ]


def _bright_frame(
    low_body: float,
    pivot: float,
    emergence: float,
    bite: float,
    presence: float,
    cap: float,
    semitone_table: list[float],
) -> list[dict[str, Any]]:
    return [
        _stage("allpole", low_body, 0.960, 0.82),
        _stage(
            "explicit_zero",
            pivot,
            0.948,
            0.74,
            zero_freq_hz=_snap_bark(pivot * 0.74, semitone_table),
            zero_radius=0.82,
        ),
        _stage("allpole", emergence, 0.972, 0.96),
        _stage(
            "explicit_zero",
            bite,
            0.966,
            0.86,
            zero_freq_hz=_snap_bark(bite * 0.90, semitone_table),
            zero_radius=0.80,
        ),
        _stage(
            "explicit_zero",
            presence,
            0.950,
            0.74,
            zero_freq_hz=_snap_bark(presence * 0.96, semitone_table),
            zero_radius=0.74,
        ),
        _stage("zero_forced_offset", cap, 0.910, 0.52, offset_semitones=3.0),
    ]


def _build_body(exp: Experiment, semitone_table: list[float], tables: dict[str, Any]) -> dict[str, Any]:
    vowels = tables["vowels"]
    landmarks = {entry["name"]: entry for entry in tables["landmarks"]["entries"]}

    v0 = vowels[exp.start_vowel]
    v1 = vowels[exp.end_vowel]
    anchor = landmarks[exp.anchor_landmark]
    air = landmarks[exp.air_landmark]

    f1_0 = _snap_bark(float(v0["f1"]), semitone_table)
    f2_0 = _snap_bark(float(v0["f2"]), semitone_table)
    f3_0 = _snap_bark(float(v0["f3"]), semitone_table)
    f4_0 = _snap_bark(float(v0["f4"]), semitone_table)
    f1_1 = _snap_bark(float(v1["f1"]), semitone_table)
    f2_1 = _snap_bark(float(v1["f2"]), semitone_table)
    f3_1 = _snap_bark(float(v1["f3"]), semitone_table)
    f4_1 = _snap_bark(float(v1["f4"]), semitone_table)

    anchor_hz = float(anchor["freq_hz"])
    air_hz = float(air["freq_hz"])
    low_body_0 = _snap_bark(max(anchor_hz * 2.4, f1_0 * 0.92), semitone_table)
    low_body_1 = _snap_bark(max(anchor_hz * 5.4, f1_1 * 1.22), semitone_table)
    low_mid_0 = _snap_bark(_blend_log(low_body_0, f2_0, 0.34), semitone_table)
    low_mid_1 = _snap_bark(_blend_log(low_body_1, f2_1, 0.40), semitone_table)
    emergence_0 = _snap_bark(_blend_log(f1_0, f2_0, 0.26), semitone_table)
    emergence_1 = _snap_bark(_blend_log(f1_1, f2_1, 0.56), semitone_table)
    bite_0 = _snap_bark(_blend_log(f2_0, f3_0, 0.46), semitone_table)
    bite_1 = _snap_bark(_blend_log(f2_1, f3_1, 0.58), semitone_table)
    presence_0 = _snap_bark(min(_blend_log(f3_0, max(f4_0, air_hz), 0.18), air_hz * 0.78), semitone_table)
    presence_1 = _snap_bark(min(_blend_log(f3_1, max(f4_1, air_hz), 0.28), air_hz * 0.92), semitone_table)
    cap_0 = _snap_bark(min(max(presence_0 * 1.9, 5400.0), air_hz * 0.62), semitone_table)
    cap_1 = _snap_bark(min(max(presence_1 * 2.0, 6800.0), air_hz * 0.82), semitone_table)

    frames = {
        "M0": _dark_frame(low_body_0, low_mid_0, emergence_0, bite_0, presence_0, cap_0, semitone_table),
        "M100": _bright_frame(low_body_1, low_mid_1, emergence_1, bite_1, presence_1, cap_1, semitone_table),
    }

    return {
        "format": "raw-stage-v1",
        "name": exp.name,
        "authoring_sample_rate_hz": 39062.5,
        "boost": max(1.12, exp.boost - 0.55),
        "notes": exp.notes,
        "frames": frames,
    }


def main() -> int:
    cal = _load_json(CAL_INDEX)
    tables = _load_json(SONIC_TABLES)
    semitone_table = list(cal["re_tables"]["semitone_table"]["values"])

    raw_dir = ROOT / "cartridges" / "engine" / "_source" / "bark_experiments"
    out_dir = ROOT / "cartridges" / "engine" / "generated" / "bark_experiments"
    raw_dir.mkdir(parents=True, exist_ok=True)
    out_dir.mkdir(parents=True, exist_ok=True)

    written: list[Path] = []
    for exp in EXPERIMENTS:
        raw_doc = _build_body(exp, semitone_table, tables)
        raw_path = raw_dir / f"{exp.name}.raw.json"
        raw_path.write_text(json.dumps(raw_doc, indent=2) + "\n", encoding="utf-8")
        compiled = compile_raw(raw_doc)
        compiled["provenance"] = "bark-semitone-sonic-experiment"
        compiled_path = out_dir / f"{exp.name}.json"
        compiled_path.write_text(json.dumps(compiled, indent=2) + "\n", encoding="utf-8")
        written.append(compiled_path)

    manifest = {"generated_files": len(written), "entries": [{"file": p.name} for p in written]}
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f"generated {len(written)} files in {out_dir}")
    for path in written:
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
