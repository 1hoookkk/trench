"""Build a corpus-grounded phoneme atlas from P2K profiles.

Outputs:
  - vault/_phonemes/p2k_corner_atlas.jsonl
  - vault/_phonemes/p2k_phoneme_inventory.json

Each corner (M0_Q0, M0_Q100, M100_Q0, M100_Q100) is tagged with:
  - a learned cluster id (blanket "phoneme" module)
  - stress tag (calm/strained/violent)
  - mouth-band tag (chest/throat/mouth/teeth/air)
  - optional closest-vowel hint from docs/sonic_tables/tables.json

This is intentionally non-judgmental: clustering is the "module" assignment;
vowel hints are just hints. Novel territory is marked explicitly.
"""

from __future__ import annotations

import argparse
import glob
import json
import math
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage


CORNER_ORDER = ("M0_Q0", "M0_Q100", "M100_Q0", "M100_Q100")


@dataclass(frozen=True)
class VowelHint:
    key: str
    dist_st: float
    f1_hz: float
    f2_hz: float
    confidence: str


def _load_json(path: Path) -> dict[str, Any]:
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


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
    # Thresholds are corpus-derived (P2k_* corners) after SR fixes.
    # Median range ~97 dB, P90 ~161 dB, max ~353 dB.
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


def _vowel_hint_from_stages(stages: list[dict[str, Any]], vowels: dict[str, Any]) -> VowelHint | None:
    """Pick two likely formant carriers and match to vowel table.

    Heuristic:
      - consider poles in [250, 3500] Hz (speech zone)
      - pick the two highest-radius poles as candidates
      - map lower->F1, higher->F2
    """
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
            best_key = key

    if best_key == "":
        return None

    # Confidence is purely geometric (no claim of perceptual match).
    if best_dist <= 4.0:
        conf = "high"
    elif best_dist <= 8.0:
        conf = "medium"
    elif best_dist <= 14.0:
        conf = "low"
    else:
        conf = "none"

    return VowelHint(
        key=best_key,
        dist_st=round(float(best_dist), 2),
        f1_hz=round(float(f1_hz), 1),
        f2_hz=round(float(f2_hz), 1),
        confidence=conf,
    )


def _corner_vector(corner: dict[str, Any]) -> np.ndarray:
    """Numeric vector for clustering (stage-order stable)."""
    stages = corner["stages"]
    response = corner["response"]
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


def _zscore(x: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    mean = np.mean(x, axis=0)
    std = np.std(x, axis=0)
    std = np.where(std < 1e-12, 1.0, std)
    return (x - mean) / std, mean, std


def _cluster_ids(xz: np.ndarray, max_clusters: int) -> np.ndarray:
    # Ward linkage tends to produce compact, interpretable clusters.
    z = linkage(xz, method="ward")
    labels = fcluster(z, t=int(max_clusters), criterion="maxclust")
    return labels.astype(int)


def _cluster_key(cluster_index: int) -> str:
    return f"PH_{cluster_index:03d}"


def _filter_type_from_name(name: str) -> int | None:
    # "P2k_013" -> 13
    try:
        if "_" not in name:
            return None
        return int(name.split("_", 1)[1])
    except Exception:
        return None


def _load_filter_names(names_path: Path) -> dict[int, str]:
    if not names_path.is_file():
        return {}
    d = _load_json(names_path)
    out: dict[int, str] = {}
    # Accept:
    # - {"names": {"13": "Phaser 2", ...}, ...} (repo default)
    # - {"13": "Phaser 2", ...}
    # - [{"filterType": 13, "name": "..."}]
    if isinstance(d, dict):
        name_map = d.get("names", d)
        if not isinstance(name_map, dict):
            return {}
        for k, v in name_map.items():
            try:
                out[int(k)] = str(v)
            except Exception:
                continue
        return out
    if isinstance(d, list):
        for row in d:
            try:
                out[int(row["filterType"])] = str(row["name"])
            except Exception:
                continue
    return out


def build_atlas(*, profiles_dir: Path, out_dir: Path, clusters: int) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)

    tables = _load_json(Path("docs/sonic_tables/tables.json"))
    vowels = tables.get("vowels", {})

    filter_names = _load_filter_names(Path("datasets/p2k_filter_names.json"))

    profile_paths = sorted(Path(p).resolve() for p in glob.glob(str(profiles_dir / "P2k_*_profile.json")))
    if not profile_paths:
        raise SystemExit(f"No profiles found in: {profiles_dir}")

    records: list[dict[str, Any]] = []
    vectors: list[np.ndarray] = []

    for path in profile_paths:
        prof = _load_json(path)
        body_id = str(prof.get("name", path.stem.replace("_profile", "")))
        filter_type = _filter_type_from_name(body_id)
        filter_name = filter_names.get(filter_type, "") if filter_type is not None else ""

        for corner_key in CORNER_ORDER:
            corner = prof["corners"][corner_key]
            response = corner["response"]

            centroid_hz = float(response["centroid_hz"])
            dyn_range = float(response["dynamic_range_db"])
            peak_db = float(response["peak_db"])
            notch_db = float(response["notch_db"])

            band = _mouth_band(centroid_hz)
            stress = _stress_tag(dyn_range, peak_db)
            vh = _vowel_hint_from_stages(corner["stages"], vowels=vowels)

            record: dict[str, Any] = {
                "id": body_id,
                "filter_type": filter_type,
                "filter_name": filter_name,
                "corner": corner_key,
                "mouth_band": band,
                "stress": stress,
                "response": response,
                "stages": corner["stages"],
            }
            if vh is not None:
                record["closest_vowel"] = {
                    "key": vh.key,
                    "label": str(vowels.get(vh.key, {}).get("label", vh.key)),
                    "dist_st": vh.dist_st,
                    "f1_hz": vh.f1_hz,
                    "f2_hz": vh.f2_hz,
                    "confidence": vh.confidence,
                }

            v = _corner_vector(corner)
            record["_vector_raw"] = v.tolist()
            records.append(record)
            vectors.append(v)

    x = np.vstack(vectors)
    xz, mean, std = _zscore(x)
    labels = _cluster_ids(xz, max_clusters=clusters)

    # Attach cluster ids + normalized vectors (for vector DB ingestion)
    for record, lab, vz in zip(records, labels.tolist(), xz.tolist()):
        record["phoneme"] = _cluster_key(int(lab))
        record["vector"] = [round(float(v), 6) for v in vz]
        record.pop("_vector_raw", None)

    # Inventory summary
    inventory: dict[str, Any] = {
        "source": "P2K profiles (vault/_profiles/P2k_*_profile.json)",
        "clusters": int(clusters),
        "vector_schema": {
            "per_stage": ["log10(pole_hz)", "pole_r", "log10(zero_hz)", "zero_r"],
            "global": ["log10(centroid_hz)", "dynamic_range_db", "peak_db", "notch_db", "num_peaks", "num_notches"],
            "stage_count": 6,
            "dims": int(xz.shape[1]),
            "zscore_mean": [round(float(v), 9) for v in mean.tolist()],
            "zscore_std": [round(float(v), 9) for v in std.tolist()],
        },
        "phonemes": {},
    }

    by_cluster: dict[str, list[dict[str, Any]]] = {}
    for r in records:
        by_cluster.setdefault(str(r["phoneme"]), []).append(r)

    for key in sorted(by_cluster.keys()):
        rs = by_cluster[key]
        # Aggregate hints: which vowel keys appear with high/med confidence?
        vowel_votes: dict[str, int] = {}
        for r in rs:
            cv = r.get("closest_vowel")
            if not isinstance(cv, dict):
                continue
            conf = str(cv.get("confidence", "none"))
            if conf not in ("high", "medium"):
                continue
            vk = str(cv.get("key", ""))
            if vk:
                vowel_votes[vk] = vowel_votes.get(vk, 0) + 1

        # Representative examples: choose 3 closest to cluster mean vector.
        idxs = [i for i, r in enumerate(records) if r["phoneme"] == key]
        sub = xz[np.array(idxs)]
        mu = np.mean(sub, axis=0)
        dists = np.linalg.norm(sub - mu, axis=1)
        best = np.argsort(dists)[:3].tolist()
        reps = []
        for j in best:
            r = rs[j]
            reps.append({"id": r["id"], "corner": r["corner"], "filter_name": r.get("filter_name", "")})

        # Mark novel if no high/medium vowel votes at all.
        is_novel = len(vowel_votes) == 0

        inventory["phonemes"][key] = {
            "count": int(len(rs)),
            "novel": bool(is_novel),
            "vowel_votes": vowel_votes,
            "representatives": reps,
        }

    atlas_path = out_dir / "p2k_corner_atlas.jsonl"
    with open(atlas_path, "w", encoding="utf-8") as f:
        for r in records:
            f.write(json.dumps(r, ensure_ascii=False) + "\n")

    inv_path = out_dir / "p2k_phoneme_inventory.json"
    with open(inv_path, "w", encoding="utf-8") as f:
        json.dump(inventory, f, indent=2, ensure_ascii=False)

    print(f"Wrote {atlas_path}")
    print(f"Wrote {inv_path}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--profiles", default="vault/_profiles", help="Directory containing P2k_*_profile.json files")
    parser.add_argument("--out", default="vault/_phonemes", help="Output directory")
    parser.add_argument("--clusters", type=int, default=32, help="Maximum cluster count (blanket phoneme modules)")
    args = parser.parse_args()

    build_atlas(
        profiles_dir=Path(args.profiles),
        out_dir=Path(args.out),
        clusters=int(args.clusters),
    )


if __name__ == "__main__":
    main()
