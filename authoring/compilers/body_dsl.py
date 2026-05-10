"""Minimum-surface authoring DSL for TRENCH 6-section bodies."""
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


# Role-implied defaults, reverse-engineered from bench_peak_belch_a /
# bench_vowel_arc_b / bench_razor_fold_c (empirical, not invented).
ROLE_DEFAULTS: dict[str, dict[str, float]] = {
    "anchor":   {"zero_ratio": 35.0, "zero_radius": 0.70},
    "mass":     {"zero_ratio":  5.0, "zero_radius": 0.70},
    "carve":    {"zero_ratio":  3.0, "zero_radius": 0.72},
    "pressure": {"zero_ratio":  2.2, "zero_radius": 0.73},
    "detail":   {"zero_ratio":  1.7, "zero_radius": 0.56},
    "boundary": {"zero_ratio":  1.05, "zero_radius": 0.96},
}

NYQUIST_LIMIT_HZ = 39062.5 * 0.5 - 10.0


def _pair(value: Any) -> tuple[float, float]:
    if isinstance(value, (tuple, list)) and len(value) == 2:
        return float(value[0]), float(value[1])
    v = float(value)
    return v, v


@dataclass(frozen=True)
class Slot:
    role: str
    pole: Any
    Q: Any
    gain: Any
    zero: Any = None
    zero_r: float | None = None
    kind: str = "explicit_zero"

    def resolved(self) -> dict[str, Any]:
        pole_lo, pole_hi = _pair(self.pole)
        q_lo, q_hi = _pair(self.Q)
        gain_lo, gain_hi = _pair(self.gain)
        defaults = ROLE_DEFAULTS.get(self.role, {"zero_ratio": 3.0, "zero_radius": 0.7})
        if self.zero is None:
            ratio = defaults["zero_ratio"]
            zero_lo = min(pole_lo * ratio, NYQUIST_LIMIT_HZ)
            zero_hi = min(pole_hi * ratio, NYQUIST_LIMIT_HZ)
        else:
            zero_lo, zero_hi = _pair(self.zero)
        zero_r = float(self.zero_r) if self.zero_r is not None else defaults["zero_radius"]
        return {
            "role": self.role,
            "kind": self.kind,
            "pole_lo": pole_lo, "pole_hi": pole_hi,
            "q_lo": q_lo, "q_hi": q_hi,
            "gain_lo": gain_lo, "gain_hi": gain_hi,
            "zero_lo": zero_lo, "zero_hi": zero_hi,
            "zero_radius": zero_r,
        }


@dataclass(frozen=True)
class Body:
    name: str
    boost: float
    slots: tuple
    notes: str = ""


def _stage_dict(slot: dict, *, hi: bool) -> dict[str, Any]:
    return {
        "role": slot["role"],
        "kind": slot["kind"],
        "pole_freq_hz": slot["pole_hi"] if hi else slot["pole_lo"],
        "radius": slot["q_hi"] if hi else slot["q_lo"],
        "stage_gain": slot["gain_hi"] if hi else slot["gain_lo"],
        "zero_freq_hz": slot["zero_hi"] if hi else slot["zero_lo"],
        "zero_radius": slot["zero_radius"],
    }


def body_to_raw(body: Body) -> dict[str, Any]:
    resolved = [slot.resolved() for slot in body.slots]
    return {
        "format": "raw-stage-v1",
        "name": body.name,
        "authoring_sample_rate_hz": 39062.5,
        "boost": body.boost,
        "notes": body.notes,
        "frames": {
            "M0":   [_stage_dict(s, hi=False) for s in resolved],
            "M100": [_stage_dict(s, hi=True)  for s in resolved],
        },
    }


def build(body: Body, out_root: Path) -> tuple[Path, Path]:
    """Write raw + compiled JSON. Returns (raw_path, compiled_path)."""
    raw = body_to_raw(body)
    compiled = compile_raw(raw)
    compiled["provenance"] = "body-dsl"
    raw_dir = out_root / "raw"
    compiled_dir = out_root / "compiled"
    raw_dir.mkdir(parents=True, exist_ok=True)
    compiled_dir.mkdir(parents=True, exist_ok=True)
    raw_path = raw_dir / f"{body.name}.raw.json"
    compiled_path = compiled_dir / f"{body.name}.json"
    raw_path.write_text(json.dumps(raw, indent=2) + "\n", encoding="utf-8")
    compiled_path.write_text(json.dumps(compiled, indent=2) + "\n", encoding="utf-8")
    return raw_path, compiled_path
