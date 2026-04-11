"""Label overnight-role-loop candidates with corpus-grounded phoneme modules.

This is compute-only: no LLM calls.

Inputs:
  - vault/_overnight/<run_id>/manifest.json (from tools/run_overnight_role_loop.py)
  - vault/_phonemes/p2k_corner_atlas.jsonl + p2k_phoneme_inventory.json

Outputs (written into the run folder):
  - phoneme_labels.json
  - phoneme_labels.md

Each candidate gets, per corner:
  - nearest phoneme cluster id (PH_###)
  - mouth_band (chest/throat/mouth/teeth/air)
  - stress (calm/strained/violent)
  - optional closest vowel hint (from docs/sonic_tables/tables.json)
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from scipy.signal import find_peaks

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pyruntime.body import Body
from pyruntime.constants import SR, TWO_PI
from pyruntime.corner import CornerName
from pyruntime.freq_response import cascade_response_db, freq_points


CORNER_LABELS: list[tuple[CornerName, str]] = [
    (CornerName.A, "M0_Q0"),
    (CornerName.B, "M0_Q100"),
    (CornerName.C, "M100_Q0"),
    (CornerName.D, "M100_Q100"),
]


@dataclass(frozen=True)
class VowelHint:
    key: str
    dist_st: float
    f1_hz: float
    f2_hz: float
    confidence: str
    label: str


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _safe_log10_hz(hz: float) -> float:
    return float(math.log10(max(1.0, float(hz))))


def _mouth_band(centroid_hz: float) -> str:
    c = float(centroid_hz)
    if c < 400.0:
        return "chest"
    if c < 1000.0:
        return "throat"
    if c < 2500.0:
        return "mouth"
    if c < 7000.0:
        return "teeth"
    return "air"


def _stress_tag(dynamic_range_db: float, peak_db: float) -> str:
    # Same heuristic as tools/build_phoneme_atlas.py
    score = float(dynamic_range_db) + max(0.0, float(peak_db)) * 0.35
    if score < 110.0:
        return "calm"
    if score < 170.0:
        return "strained"
    return "violent"


def _st_distance(a_hz: float, b_hz: float) -> float:
    a = max(1.0, float(a_hz))
    b = max(1.0, float(b_hz))
    return abs(12.0 * math.log(a / b, 2.0))


def _pole_zero_from_kernel(c0: float, c1: float, c2: float, c3: float, c4: float) -> dict[str, float]:
    # Mirrors tools/profile_p2k_skin.py
    r = math.sqrt(abs(float(c4)))
    if r > 0.001:
        cos_t = max(-1.0, min(1.0, -float(c3) / (2.0 * r)))
        pole_hz = math.acos(cos_t) * SR / TWO_PI
    else:
        pole_hz = 0.0

    if abs(float(c0)) > 0.001:
        zr = math.sqrt(abs(float(c2) / float(c0)))
        if zr > 0.001:
            zcos = max(-1.0, min(1.0, -(float(c1) / float(c0)) / (2.0 * zr)))
            zero_hz = math.acos(zcos) * SR / TWO_PI
        else:
            zero_hz = 0.0
    else:
        zero_hz = 0.0
        zr = 0.0

    return {
        "pole_hz": float(pole_hz),
        "pole_r": float(r),
        "zero_hz": float(zero_hz),
        "zero_r": float(zr),
    }


def _response_stats(enc, *, freqs: np.ndarray) -> dict[str, float]:
    db = cascade_response_db(enc, freqs, SR)
    peak_db = float(np.max(db))
    notch_db = float(np.min(db))
    dyn_range = peak_db - notch_db

    linear = 10 ** (db / 20.0)
    total = float(np.sum(linear))
    centroid_hz = float(np.sum(freqs * linear) / total) if total > 0.0 else 0.0

    peaks, _ = find_peaks(db, prominence=3.0)
    notches, _ = find_peaks(-db, prominence=3.0)

    return {
        "peak_db": float(round(peak_db, 1)),
        "notch_db": float(round(notch_db, 1)),
        "dynamic_range_db": float(round(dyn_range, 1)),
        "centroid_hz": float(round(centroid_hz, 1)),
        "num_peaks": float(len(peaks)),
        "num_notches": float(len(notches)),
    }


def _vowel_hint_from_stages(stages: list[dict[str, float]], vowels: dict[str, Any]) -> VowelHint | None:
    # Same heuristic as tools/build_phoneme_atlas.py
    candidates = [s for s in stages if 250.0 <= float(s.get("pole_hz", 0.0)) <= 3500.0]
    if len(candidates) < 2:
        return None

    candidates.sort(key=lambda s: (-float(s.get("pole_r", 0.0)), float(s.get("pole_hz", 0.0))))
    a = float(candidates[0]["pole_hz"])
    b = float(candidates[1]["pole_hz"])
    f1_hz = min(a, b)
    f2_hz = max(a, b)

    best_key = ""
    best_dist = float("inf")
    for key, v in vowels.items():
        d1 = _st_distance(f1_hz, float(v["f1"]))
        d2 = _st_distance(f2_hz, float(v["f2"]))
        dist = float(math.sqrt(d1 * d1 + d2 * d2))
        if dist < best_dist:
            best_dist = dist
            best_key = str(key)

    if not best_key:
        return None

    if best_dist <= 4.0:
        conf = "high"
    elif best_dist <= 8.0:
        conf = "medium"
    elif best_dist <= 14.0:
        conf = "low"
    else:
        conf = "none"

    label = str(vowels.get(best_key, {}).get("label", best_key))
    return VowelHint(
        key=best_key,
        dist_st=float(round(best_dist, 2)),
        f1_hz=float(round(f1_hz, 1)),
        f2_hz=float(round(f2_hz, 1)),
        confidence=conf,
        label=label,
    )


def _candidate_corner_vector(*, stages: list[dict[str, float]], response: dict[str, float]) -> np.ndarray:
    feats: list[float] = []
    for s in stages[:6]:
        feats.append(_safe_log10_hz(float(s.get("pole_hz", 0.0))))
        feats.append(float(s.get("pole_r", 0.0)))
        feats.append(_safe_log10_hz(float(s.get("zero_hz", 0.0))))
        feats.append(float(s.get("zero_r", 0.0)))

    feats.extend(
        [
            _safe_log10_hz(float(response.get("centroid_hz", 0.0))),
            float(response.get("dynamic_range_db", 0.0)),
            float(response.get("peak_db", 0.0)),
            float(response.get("notch_db", 0.0)),
            float(response.get("num_peaks", 0.0)),
            float(response.get("num_notches", 0.0)),
        ]
    )
    return np.array(feats, dtype=np.float64)


def _load_phoneme_centroids(atlas_path: Path) -> dict[str, np.ndarray]:
    if not atlas_path.is_file():
        raise FileNotFoundError(f"Missing phoneme atlas: {atlas_path}")

    by: dict[str, list[np.ndarray]] = {}
    for line in atlas_path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        r = json.loads(line)
        key = str(r.get("phoneme", ""))
        vec = r.get("vector")
        if not key or not isinstance(vec, list):
            continue
        by.setdefault(key, []).append(np.array(vec, dtype=np.float64))

    centroids: dict[str, np.ndarray] = {}
    for key, vs in by.items():
        centroids[key] = np.mean(np.vstack(vs), axis=0)
    return centroids


def _nearest_phoneme(vec: np.ndarray, centroids: dict[str, np.ndarray]) -> tuple[str, float]:
    best_key = ""
    best_dist = float("inf")
    for key, mu in centroids.items():
        d = float(np.linalg.norm(vec - mu))
        if d < best_dist:
            best_dist = d
            best_key = key
    return best_key, best_dist


def label_run(*, run_dir: Path, phonemes_dir: Path) -> Path:
    manifest_path = run_dir / "manifest.json"
    if not manifest_path.is_file():
        raise FileNotFoundError(f"Missing manifest: {manifest_path}")

    inventory_path = phonemes_dir / "p2k_phoneme_inventory.json"
    atlas_path = phonemes_dir / "p2k_corner_atlas.jsonl"
    inv = _load_json(inventory_path)

    schema = inv.get("vector_schema", {})
    mean = np.array(schema.get("zscore_mean", []), dtype=np.float64)
    std = np.array(schema.get("zscore_std", []), dtype=np.float64)
    if mean.size == 0 or std.size == 0:
        raise ValueError(f"Phoneme inventory missing zscore stats: {inventory_path}")

    centroids = _load_phoneme_centroids(atlas_path)
    if not centroids:
        raise ValueError(f"No phoneme centroids loaded from: {atlas_path}")

    tables = _load_json(Path("docs/sonic_tables/tables.json"))
    vowels = tables.get("vowels", {})

    freqs = freq_points(sr=SR)

    manifest = _load_json(manifest_path)
    out: dict[str, Any] = {"run": str(run_dir), "generated_from": str(manifest_path), "candidates": []}

    for body_key, rows in manifest.items():
        if not isinstance(rows, list):
            continue
        for row in rows:
            path = Path(str(row.get("path", ""))).expanduser()
            if not path.is_file():
                continue
            body = Body.from_json(str(path))

            corners: dict[str, Any] = {}
            for cn, label in CORNER_LABELS:
                enc = body.corners.corner(cn).encode()
                stages: list[dict[str, float]] = []
                for s in enc[:6]:
                    pz = _pole_zero_from_kernel(s.c0, s.c1, s.c2, s.c3, s.c4)
                    stages.append({k: float(round(v, 6)) for k, v in pz.items()})

                resp = _response_stats(enc, freqs=freqs)
                vec_raw = _candidate_corner_vector(stages=stages, response=resp)
                vec = (vec_raw - mean) / np.where(std < 1e-12, 1.0, std)
                phoneme_key, phoneme_dist = _nearest_phoneme(vec, centroids)

                vh = _vowel_hint_from_stages(stages, vowels=vowels) if isinstance(vowels, dict) else None

                corners[label] = {
                    "phoneme": phoneme_key,
                    "phoneme_dist": float(round(phoneme_dist, 4)),
                    "mouth_band": _mouth_band(resp["centroid_hz"]),
                    "stress": _stress_tag(resp["dynamic_range_db"], resp["peak_db"]),
                    "response": resp,
                }
                if vh is not None:
                    corners[label]["closest_vowel"] = {
                        "key": vh.key,
                        "label": vh.label,
                        "dist_st": vh.dist_st,
                        "f1_hz": vh.f1_hz,
                        "f2_hz": vh.f2_hz,
                        "confidence": vh.confidence,
                    }

            out["candidates"].append(
                {
                    "body_key": body_key,
                    "name": str(row.get("name", body.name)),
                    "archetype": row.get("archetype"),
                    "seed": row.get("seed"),
                    "path": str(path),
                    "metrics": row.get("metrics", {}),
                    "corners": corners,
                }
            )

    json_out = run_dir / "phoneme_labels.json"
    json_out.write_text(json.dumps(out, indent=2), encoding="utf-8")

    md_lines: list[str] = []
    md_lines.append(f"# Overnight Candidate Phoneme Labels · {run_dir.name}")
    md_lines.append("")
    md_lines.append(f"- candidates: `{len(out['candidates'])}`")
    md_lines.append(f"- phonemes: `{phonemes_dir}`")
    md_lines.append("")
    for item in out["candidates"]:
        md_lines.append(f"## {item['name']}")
        md_lines.append(f"- path: `{item['path']}`")
        md_lines.append(f"- body_key: `{item['body_key']}` · archetype: `{item.get('archetype')}` · seed: `{item.get('seed')}`")
        md_lines.append("")
        md_lines.append("| Corner | Phoneme | Band | Stress | Vowel hint |")
        md_lines.append("|---|---|---|---|---|")
        corners = item.get("corners", {})
        for _, label in CORNER_LABELS:
            c = corners.get(label, {})
            ph = c.get("phoneme", "PH_???")
            band = c.get("mouth_band", "?")
            stress = c.get("stress", "?")
            vh = c.get("closest_vowel")
            if isinstance(vh, dict) and vh.get("confidence") in ("high", "medium"):
                vtxt = f"{vh.get('label')} ({vh.get('dist_st')} st)"
            else:
                vtxt = "—"
            md_lines.append(f"| {label} | `{ph}` | `{band}` | `{stress}` | {vtxt} |")
        md_lines.append("")

    md_out = run_dir / "phoneme_labels.md"
    md_out.write_text("\n".join(md_lines).strip() + "\n", encoding="utf-8")
    return md_out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True, help="Run directory under vault/_overnight/<run_id>")
    parser.add_argument("--phonemes", type=Path, default=Path("vault/_phonemes"), help="Phoneme atlas directory")
    args = parser.parse_args()

    md_out = label_run(run_dir=args.run.resolve(), phonemes_dir=args.phonemes.resolve())
    print(md_out)


if __name__ == "__main__":
    main()
