"""Deterministic shipping backbone translation.

This script translates fixed P2K source bodies into the four shipping body names
without inventing new topology or solving new coefficients.

Translation only:
- preserve source corner/stage content
- rename body
- write keyframe + compiled outputs
- report gate/profile measurements for visibility
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pyruntime.analysis import body_profile, shipping_gate
from pyruntime.body import Body


TRANSLATION_MAP = [
    ("Speaker Knockerz", ROOT / "cartridges" / "p2k" / "P2k_006.json"),
    ("Aluminum Siding", ROOT / "cartridges" / "p2k" / "P2k_026.json"),
    ("Small Talk Ah-Ee", ROOT / "cartridges" / "p2k" / "P2k_013.json"),
    ("Cul-De-Sac", ROOT / "cartridges" / "p2k" / "P2k_010.json"),
]

OUT_DIR = ROOT / "cartridges" / "translation"


def slug(name: str) -> str:
    return (
        name.lower()
        .replace("-", " ")
        .replace(" ", "_")
        .replace("__", "_")
    )


def translate_one(public_name: str, source_path: Path) -> dict:
    source = Body.from_json(str(source_path))
    translated = Body(
        name=public_name,
        corners=source.corners,
        boost=source.boost,
        preset_schema=source.preset_schema,
    )

    stem = slug(public_name)
    keyframe_path = OUT_DIR / f"{stem}.keyframe.json"
    compiled_path = OUT_DIR / f"{stem}.compiled.json"

    keyframe_path.write_text(translated.to_json(), encoding="utf-8")
    compiled_path.write_text(
        translated.to_compiled_json(provenance="translation-p2k-backbone"),
        encoding="utf-8",
    )

    profile = body_profile(translated)
    gate = shipping_gate(translated, public_name)

    return {
        "name": public_name,
        "source": str(source_path),
        "keyframe": str(keyframe_path),
        "compiled": str(compiled_path),
        "midpoint_passed": profile["midpoint_audit"]["passed"],
        "morph_distance_mean_db": profile["morph_distance"]["mean_rms_db"],
        "shipping_gate_passed": gate["passed"],
        "shipping_gate_failures": gate["failures"],
    }


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    report = [translate_one(name, src) for name, src in TRANSLATION_MAP]
    report_path = OUT_DIR / "translation_report.json"
    report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(report_path)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
