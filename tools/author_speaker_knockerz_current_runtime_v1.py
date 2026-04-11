from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from pyruntime.analysis import body_profile
from pyruntime.body import Body
from pyruntime.corner import CornerArray, CornerName, CornerState
from pyruntime.freq_response import cascade_response_db, freq_points
from pyruntime.stage_math import resonator_with_zero


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "cartridges" / "p2k" / "P2k_006.json"
OUT_DIR = ROOT / "cartridges" / "candidates"
OUT_KEYFRAME = OUT_DIR / "Speaker_Knockerz_current_runtime_v1.keyframe.json"
OUT_COMPILED = OUT_DIR / "Speaker_Knockerz_current_runtime_v1.compiled.json"
PEAK_FREQS = freq_points(16384)


def rewrite_zero(corner: CornerState, stage_index: int, zero_hz: float, zero_radius: float) -> CornerState:
    stages = list(corner.stages)
    src = stages[stage_index]
    stages[stage_index] = resonator_with_zero(
        src.freq_hz(),
        src.r,
        src.val1,
        zero_hz,
        zero_radius,
    )
    return CornerState(stages=stages, boost=corner.boost)


def build_candidate(base: Body) -> Body:
    corners: dict[CornerName, CornerState] = {
        name: base.corners.corner(name) for name in CornerName
    }

    # One sparse numerator move only:
    # keep the full P2k_006 pole scaffold, but lift the C-corner low pole's zero
    # high enough to stop M100_Q0 from turning into a 1.35 kHz-led open state.
    corners[CornerName.C] = rewrite_zero(
        corners[CornerName.C],
        stage_index=1,
        zero_hz=5000.0,
        zero_radius=0.40,
    )

    return Body(
        name="Speaker Knockerz current runtime v1",
        corners=CornerArray(
            corners[CornerName.A],
            corners[CornerName.B],
            corners[CornerName.C],
            corners[CornerName.D],
        ),
        boost=base.boost,
        preset_schema=base.preset_schema,
    )


def peak_trough(body: Body, morph: float, q: float) -> dict[str, float]:
    db = cascade_response_db(body.corners.interpolate(morph, q), PEAK_FREQS)
    peak_i = int(np.argmax(db))
    trough_i = int(np.argmin(db))
    return {
        "peak_hz": float(PEAK_FREQS[peak_i]),
        "peak_db": float(db[peak_i]),
        "trough_hz": float(PEAK_FREQS[trough_i]),
        "trough_db": float(db[trough_i]),
    }


def summarize(body: Body) -> dict:
    profile = body_profile(body)
    return {
        "morph_distance": profile["morph_distance"],
        "crossings": profile["crossings"],
        "spectral_tilt_db": profile["spectral_tilt_db"],
        "midpoint_audit": profile["midpoint_audit"],
        "states": {
            "M0_Q0": peak_trough(body, 0.0, 0.0),
            "M50_Q50": peak_trough(body, 0.5, 0.5),
            "M100_Q0": peak_trough(body, 1.0, 0.0),
            "M100_Q100": peak_trough(body, 1.0, 1.0),
        },
    }


def main() -> None:
    base = Body.from_json(str(SOURCE))
    candidate = build_candidate(base)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    OUT_KEYFRAME.write_text(candidate.to_json(), encoding="utf-8")
    OUT_COMPILED.write_text(
        candidate.to_compiled_json(
            provenance="codex-speaker-knockerz-current-runtime-v1"
        ),
        encoding="utf-8",
    )

    report = {
        "source": str(SOURCE),
        "candidate_keyframe": str(OUT_KEYFRAME),
        "candidate_compiled": str(OUT_COMPILED),
        "edit": {
            "corner": "M100_Q0",
            "stage_index": 1,
            "zero_hz": 5000.0,
            "zero_radius": 0.40,
        },
        "raw_p2k_006": summarize(base),
        "speaker_knockerz_candidate": summarize(candidate),
    }
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
