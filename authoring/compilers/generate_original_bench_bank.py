#!/usr/bin/env python3
"""Generate a small original 6-stage bench bank via compile_raw.py.

This is a reset path for body work: original authored stage pairs only,
no E-mu seed cloning, no MorphLP helper path, no authored Q corners.
"""

from __future__ import annotations

import json
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

_TOOLS = Path(__file__).resolve().parent
if str(_TOOLS) not in sys.path:
    sys.path.insert(0, str(_TOOLS))

from compile_raw import compile_raw  # noqa: E402


@dataclass(frozen=True)
class StageSpec:
    role: str
    kind: str
    pole_m0: float
    pole_m1: float
    zero_m0: float
    zero_m1: float
    radius_m0: float
    radius_m1: float
    gain_m0: float
    gain_m1: float
    zero_radius: float = 1.0


@dataclass(frozen=True)
class BodySpec:
    name: str
    boost: float
    notes: str
    stages: tuple[StageSpec, ...]


BODIES: tuple[BodySpec, ...] = (
    BodySpec(
        name="bench_peak_belch_a",
        boost=2.80,
        notes="anchor-first skeleton: low shelf opens, mids build, boundary trims at the top",
        stages=(
            StageSpec("anchor", "explicit_zero", 198.0, 170.0, 7280.0, 7060.0, 0.940, 0.940, 0.42, 0.40, 0.72),
            StageSpec("mass", "explicit_zero", 910.0, 980.0, 4620.0, 4720.0, 0.980, 0.992, 0.58, 0.55, 0.71),
            StageSpec("carve", "explicit_zero", 1540.0, 1670.0, 4860.0, 5030.0, 0.978, 0.993, 0.59, 0.56, 0.72),
            StageSpec("pressure", "explicit_zero", 2360.0, 2820.0, 5260.0, 5460.0, 0.994, 0.998, 0.60, 0.57, 0.73),
            StageSpec("detail", "explicit_zero", 4580.0, 4880.0, 8380.0, 7240.0, 0.949, 0.983, 0.58, 0.53, 0.56),
            StageSpec("boundary", "explicit_zero", 9280.0, 8560.0, 9800.0, 9200.0, 0.975, 0.964, 0.53, 0.50, 0.96),
        ),
    ),
    BodySpec(
        name="bench_vowel_arc_b",
        boost=2.30,
        notes="anchor-first vowel arc: low shelf opens, mid formants arc, boundary trims top",
        stages=(
            StageSpec("anchor", "explicit_zero", 212.0, 188.0, 7560.0, 7320.0, 0.940, 0.940, 0.38, 0.40, 0.68),
            StageSpec("mass", "explicit_zero", 760.0, 940.0, 4560.0, 4700.0, 0.981, 0.994, 0.56, 0.53, 0.70),
            StageSpec("carve", "explicit_zero", 1360.0, 1720.0, 4880.0, 5060.0, 0.982, 0.994, 0.57, 0.54, 0.72),
            StageSpec("pressure", "explicit_zero", 2480.0, 3020.0, 5320.0, 5520.0, 0.990, 0.997, 0.60, 0.57, 0.73),
            StageSpec("detail", "explicit_zero", 4280.0, 4520.0, 8140.0, 6880.0, 0.950, 0.976, 0.57, 0.53, 0.54),
            StageSpec("boundary", "explicit_zero", 8640.0, 7820.0, 9100.0, 8300.0, 0.969, 0.958, 0.50, 0.50, 0.96),
        ),
    ),
    BodySpec(
        name="bench_razor_fold_c",
        boost=2.55,
        notes="anchor-first razor fold: low shelf, sharp-Q mid razor, boundary trims top",
        stages=(
            StageSpec("anchor", "explicit_zero", 184.0, 162.0, 7920.0, 7640.0, 0.940, 0.940, 0.40, 0.38, 0.70),
            StageSpec("mass", "explicit_zero", 1080.0, 1220.0, 4720.0, 5020.0, 0.988, 0.995, 0.58, 0.54, 0.70),
            StageSpec("carve", "explicit_zero", 2080.0, 2440.0, 5120.0, 5360.0, 0.990, 0.997, 0.60, 0.56, 0.72),
            StageSpec("pressure", "explicit_zero", 2920.0, 3440.0, 5500.0, 5740.0, 0.992, 0.998, 0.60, 0.56, 0.73),
            StageSpec("detail", "explicit_zero", 5680.0, 6320.0, 8520.0, 7220.0, 0.962, 0.984, 0.58, 0.52, 0.58),
            StageSpec("boundary", "explicit_zero", 10120.0, 8920.0, 10800.0, 9400.0, 0.982, 0.972, 0.50, 0.46, 0.96),
        ),
    ),
)


def _stage_dict(spec: StageSpec, *, hi_frame: bool) -> dict[str, Any]:
    pole = spec.pole_m1 if hi_frame else spec.pole_m0
    radius = spec.radius_m1 if hi_frame else spec.radius_m0
    gain = spec.gain_m1 if hi_frame else spec.gain_m0

    if spec.kind in ("allpole", "lowpass2"):
        return {
            "role": spec.role,
            "kind": spec.kind,
            "pole_freq_hz": pole,
            "radius": radius,
            "stage_gain": gain,
        }
    if spec.kind == "zero_forced_offset":
        offset = 9.0 if hi_frame else 5.0
        return {
            "role": spec.role,
            "kind": "zero_forced_offset",
            "pole_freq_hz": pole,
            "radius": radius,
            "stage_gain": gain,
            "offset_semitones": offset,
        }
    return {
        "role": spec.role,
        "kind": "explicit_zero",
        "pole_freq_hz": pole,
        "radius": radius,
        "stage_gain": gain,
        "zero_freq_hz": spec.zero_m1 if hi_frame else spec.zero_m0,
        "zero_radius": spec.zero_radius,
    }


def _raw_body(spec: BodySpec) -> dict[str, Any]:
    frame_grid = {
        "M0": False,
        "M100": True,
    }
    frames: dict[str, list[dict[str, Any]]] = {}
    for label, hi_frame in frame_grid.items():
        frames[label] = [_stage_dict(stage, hi_frame=hi_frame) for stage in spec.stages]
    return {
        "format": "raw-stage-v1",
        "name": spec.name,
        "authoring_sample_rate_hz": 39062.5,
        "boost": spec.boost,
        "notes": spec.notes,
        "frames": frames,
    }


def main() -> int:
    root = Path(__file__).resolve().parent.parent
    raw_dir = root / "cartridges" / "engine" / "_source" / "bench_raw_bank"
    compiled_dir = root / "cartridges" / "engine" / "generated" / "bench_raw_bank"
    raw_dir.mkdir(parents=True, exist_ok=True)
    compiled_dir.mkdir(parents=True, exist_ok=True)

    written: list[Path] = []
    for spec in BODIES:
        raw_doc = _raw_body(spec)
        raw_path = raw_dir / f"{spec.name}.raw.json"
        compiled_path = compiled_dir / f"{spec.name}.json"
        raw_path.write_text(json.dumps(raw_doc, indent=2) + "\n", encoding="utf-8")
        compiled = compile_raw(raw_doc)
        compiled["provenance"] = "original-bench-bank"
        compiled_path.write_text(json.dumps(compiled, indent=2) + "\n", encoding="utf-8")
        written.append(compiled_path)

    manifest = {
        "generated_files": len(written),
        "entries": [{"file": p.name} for p in written],
    }
    (compiled_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f"generated {len(written)} files in {compiled_dir}")
    for path in written:
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
