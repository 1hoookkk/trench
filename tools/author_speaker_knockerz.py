"""Author Speaker Knockerz for the current 6-authored-stage runtime.

Scaffold: P2k_006
Method: keep the full pole skeleton, rewrite numerator placement only.

This candidate makes one sparse move:
- M100_Q0 stage 1: keep the 158.8 Hz pole, move its zero to 5.0 kHz / r=0.40

Reason:
- raw P2k_006 keeps M100_Q100 bass-centered, but M100_Q0 is led by a 1.2 kHz peak
- the prior 12-stage rewrite drifted into HF dominance, which is a regression
- this edit restores bass-centric low-Q opening without touching the high-Q open state
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pyruntime.analysis import body_profile, dense_midpoint_audit, shipping_gate
from pyruntime.body import Body
from pyruntime.constants import SR
from pyruntime.corner import CornerArray, CornerName, CornerState
from pyruntime.freq_response import cascade_response_db, freq_points
from pyruntime.stage_math import resonator_with_zero


SOURCE = ROOT / "cartridges" / "p2k" / "P2k_006.json"
OUT_DIR = ROOT / "cartridges" / "candidates"
OUT_KEYFRAME = OUT_DIR / "Speaker_Knockerz_current_runtime.keyframe.json"
OUT_COMPILED = OUT_DIR / "Speaker_Knockerz_current_runtime.compiled.json"
PEAK_FREQS = freq_points(65536, SR)

EDIT = {
    "corner": "M100_Q0",
    "stage_index": 1,
    "zero_hz": 5000.0,
    "zero_radius": 0.40,
}


def rewrite_zero(
    corner: CornerState,
    stage_index: int,
    zero_hz: float,
    zero_radius: float,
) -> CornerState:
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

    corners[CornerName.C] = rewrite_zero(
        corners[CornerName.C],
        stage_index=EDIT["stage_index"],
        zero_hz=EDIT["zero_hz"],
        zero_radius=EDIT["zero_radius"],
    )

    return Body(
        name="Speaker Knockerz",
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
    db = cascade_response_db(body.corners.interpolate(morph, q), PEAK_FREQS, SR)
    peak_i = int(np.argmax(db))
    trough_i = int(np.argmin(db))
    return {
        "peak_hz": float(PEAK_FREQS[peak_i]),
        "peak_db": float(db[peak_i]),
        "trough_hz": float(PEAK_FREQS[trough_i]),
        "trough_db": float(db[trough_i]),
        "sub_le_60_db": float(np.max(db[PEAK_FREQS <= 60.0])),
    }


def summarize(body: Body) -> dict:
    profile = body_profile(body)
    gate = shipping_gate(body, "Speaker Knockerz")
    dense = dense_midpoint_audit(body, peak_limit_db=40.0, n_morph=21, n_q=3)
    return {
        "morph_distance": profile["morph_distance"],
        "crossings": profile["crossings"],
        "spectral_tilt_db": profile["spectral_tilt_db"],
        "midpoint_audit": profile["midpoint_audit"],
        "dense_midpoint_audit": dense,
        "speaker_knockerz_gate": gate,
        "states": {
            "M0_Q0": peak_trough(body, 0.0, 0.0),
            "M0.5_Q0.5": peak_trough(body, 0.5, 0.5),
            "M1_Q0": peak_trough(body, 1.0, 0.0),
            "M1_Q1": peak_trough(body, 1.0, 1.0),
        },
    }


def main() -> None:
    base = Body.from_json(str(SOURCE))
    candidate = build_candidate(base)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    OUT_KEYFRAME.write_text(candidate.to_json(), encoding="utf-8")
    OUT_COMPILED.write_text(
        candidate.to_compiled_json(
            provenance="codex-speaker-knockerz-current-runtime"
        ),
        encoding="utf-8",
    )

    report = {
        "source": str(SOURCE),
        "candidate_keyframe": str(OUT_KEYFRAME),
        "candidate_compiled": str(OUT_COMPILED),
        "candidate_name": candidate.name,
        "edit": EDIT,
        "raw_p2k_006": summarize(base),
        "speaker_knockerz_candidate": summarize(candidate),
    }
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
