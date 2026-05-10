#!/usr/bin/env python3
"""Annotate generated packs with calibration-guided similarity scores.

This does NOT synthesize from calibration bodies. It scores already-generated
cartridges against high-level calibration descriptors from docs/calibration.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np

_TOOLS = Path(__file__).resolve().parent
if str(_TOOLS) not in sys.path:
    sys.path.insert(0, str(_TOOLS))

from cube_authoring.thumbnail import magnitude_spectrum_db, run_impulse  # noqa: E402


REPO = Path(__file__).resolve().parent.parent
CAL_INDEX = REPO / "docs" / "calibration" / "index.json"
SR = 44100
IR_LEN = 4096
FREQ_MIN = 40.0
FREQ_MAX = 20000.0


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _curve(kf: dict[str, Any]) -> tuple[np.ndarray, np.ndarray]:
    stages = [
        (
            float(s["c0"]),
            float(s["c1"]),
            float(s["c2"]),
            float(s["c3"]),
            float(s["c4"]),
        )
        for s in kf["stages"]
    ]
    ir = run_impulse(stages, IR_LEN, float(kf.get("boost", 1.0)))
    freqs, db = magnitude_spectrum_db(ir, SR)
    mask = (freqs >= FREQ_MIN) & (freqs <= FREQ_MAX)
    return freqs[mask], db[mask]


def _dominant_peak_hz(freqs: np.ndarray, db: np.ndarray) -> float:
    if len(db) == 0:
        return 0.0
    idx = int(np.argmax(db))
    return float(freqs[idx])


def _cart_metrics(cart: dict[str, Any]) -> dict[str, float]:
    keyframes = {kf["label"]: kf for kf in cart["keyframes"]}
    f00, db00 = _curve(keyframes["M0_Q0"])
    f10, db10 = _curve(keyframes["M100_Q0"])
    peak00 = float(np.max(db00)) if len(db00) else 0.0
    min00 = float(np.min(db00)) if len(db00) else 0.0
    return {
        "m0q0_peak_db": peak00,
        "m0q0_cancellation_db": peak00 - min00,
        "morph_freq_delta": abs(_dominant_peak_hz(f10, db10) - _dominant_peak_hz(f00, db00)),
    }


def _distance(candidate: dict[str, float], ref: dict[str, Any]) -> float:
    peak = abs(candidate["m0q0_peak_db"] - float(ref.get("m0q0_peak_db", 0.0))) / 30.0
    cancel = abs(candidate["m0q0_cancellation_db"] - float(ref.get("m0q0_cancellation_db", 0.0))) / 60.0
    morph = abs(candidate["morph_freq_delta"] - float(ref.get("morph_freq_delta", 0.0))) / 8000.0
    return peak * 0.35 + cancel * 0.35 + morph * 0.30


def _selection_score(distance: float) -> float:
    return round(max(0.0, 1.0 - distance), 4)


def _annotate_entry(entry: dict[str, Any], source_dir: Path, refs: list[dict[str, Any]]) -> dict[str, Any]:
    cart = _load_json(source_dir / str(entry["file"]))
    metrics = _cart_metrics(cart)
    ranked = sorted(
        ((_distance(metrics, ref), ref) for ref in refs),
        key=lambda item: item[0],
    )
    best_distance, best_ref = ranked[0]
    out = dict(entry)
    out["calibration_match"] = str(best_ref["name"])
    out["calibration_strategy"] = str(best_ref["strategy"])
    out["calibration_distance"] = round(float(best_distance), 4)
    out["selection_score"] = _selection_score(float(best_distance))
    out["candidate_m0q0_peak_db"] = round(metrics["m0q0_peak_db"], 2)
    out["candidate_m0q0_cancellation_db"] = round(metrics["m0q0_cancellation_db"], 2)
    out["candidate_morph_freq_delta"] = round(metrics["morph_freq_delta"], 1)
    return out


def _index_markdown(manifest: dict[str, Any]) -> str:
    lines = ["# Calibration-Annotated Pack Index", ""]
    lines.append(f"- entries: {len(manifest['entries'])}")
    lines.append("")
    lines.append("| file | score | match | strategy | peak | cancellation | morph_delta |")
    lines.append("|---|---|---|---|---|---|---|")
    for e in manifest["entries"]:
        lines.append(
            f"| `{e['file']}` | `{e['selection_score']}` | `{e['calibration_match']}` | "
            f"`{e['calibration_strategy']}` | `{e['candidate_m0q0_peak_db']}` | "
            f"`{e['candidate_m0q0_cancellation_db']}` | `{e['candidate_morph_freq_delta']}` |"
        )
    lines.append("")
    return "\n".join(lines)


def main() -> int:
    ap = argparse.ArgumentParser(description="Annotate a generated pack manifest with calibration-guided scores")
    ap.add_argument("--manifest", type=Path, required=True, help="source manifest.json")
    ap.add_argument("--source-dir", type=Path, required=True, help="directory containing cartridge JSONs")
    ap.add_argument("--out-manifest", type=Path, required=True, help="destination annotated manifest")
    args = ap.parse_args()

    manifest = _load_json(args.manifest)
    refs = list(_load_json(CAL_INDEX)["skins"])
    entries = list(manifest.get("entries", []))
    if not entries:
        raise SystemExit("manifest has no entries")

    annotated = [_annotate_entry(entry, args.source_dir, refs) for entry in entries]
    annotated.sort(key=lambda e: (-float(e["selection_score"]), str(e["file"])))

    out_doc = dict(manifest)
    out_doc["entries"] = annotated
    out_doc["scoring"] = "calibration_index_v1"
    args.out_manifest.parent.mkdir(parents=True, exist_ok=True)
    args.out_manifest.write_text(json.dumps(out_doc, indent=2) + "\n", encoding="utf-8")
    args.out_manifest.with_suffix(".md").write_text(_index_markdown(out_doc), encoding="utf-8")
    print(f"wrote {args.out_manifest}")
    print(f"wrote {args.out_manifest.with_suffix('.md')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
