"""Heritage splice — cross-breed top P2K bodies into novel morph trajectories.

Takes the highest-scoring heritage bodies and splices corners between them.
Each splice creates a new 4-corner body with heritage-quality corners but
a novel morph gesture. MCP analysis ranks the results.

Usage:
    python tools/heritage_splice.py
    python tools/heritage_splice.py --top 8 --out vault/_splice_candidates
"""
import argparse
import json
import os
import sys
from itertools import combinations
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from pyruntime.body import Body
from pyruntime.splice import SpliceMode, splice_corners
from pyruntime.freq_response import cascade_response_db, freq_points
from pyruntime.constants import SR

import numpy as np

CARTRIDGE_DIR = Path(__file__).parent.parent / "cartridges" / "p2k"
FREQS = freq_points()
_VOCAL_MASK = (FREQS >= 200) & (FREQS <= 5000)

# Heritage bodies ranked by MCP composite score
HERITAGE_RANKED = [
    "P2k_012", "P2k_022", "P2k_025", "P2k_007",
    "P2k_024", "P2k_026", "P2k_020", "P2k_021",
    "P2k_030", "P2k_010",
]

MODES = [
    SpliceMode.REST_TO_MORPHED,
    SpliceMode.CROSS,
]


def _talkingness(db: np.ndarray) -> float:
    vocal_db = db[_VOCAL_MASK]
    if len(vocal_db) < 3:
        return 0.0
    mean_db = float(np.mean(db))
    peaks = 0
    for i in range(1, len(vocal_db) - 1):
        if vocal_db[i] > vocal_db[i - 1] and vocal_db[i] > vocal_db[i + 1]:
            if vocal_db[i] > mean_db + 3.0:
                peaks += 1
    return min(1.0, peaks / 4.0)


def _ridge_prominence(db: np.ndarray) -> float:
    pk = float(np.max(db))
    mn = float(np.mean(db))
    return max(0.0, min(1.0, (pk - mn) / 40.0))


def _coherence(talk_vals: list[float]) -> float:
    if not talk_vals or max(talk_vals) == 0:
        return 0.0
    return max(0.0, 1.0 - float(np.std(talk_vals)) / max(talk_vals))


def _centroid(db: np.ndarray) -> float:
    lin = 10.0 ** (db / 20.0)
    denom = float(np.sum(lin))
    if denom <= 0:
        return 1000.0
    return float(np.sum(FREQS * lin) / denom)


def score_body(body: Body) -> dict:
    """Quick 5-point morph sweep at Q=0.5. Returns composite + components."""
    morph_points = [0.0, 0.25, 0.5, 0.75, 1.0]
    talks = []
    ridges = []
    centroids = []
    peak_dbs = []

    for m in morph_points:
        enc = body.corners.interpolate(m, 0.5)
        db = cascade_response_db(enc, FREQS, SR)
        talks.append(_talkingness(db))
        ridges.append(_ridge_prominence(db))
        centroids.append(_centroid(db))
        peak_dbs.append(float(np.max(db)))

    talk = float(np.mean(talks))
    ridge = float(np.mean(ridges))
    coher = _coherence(talks)
    peak = max(peak_dbs)

    # Gesture commitment: monotonicity of centroid trajectory
    diffs = [centroids[i + 1] - centroids[i] for i in range(len(centroids) - 1)]
    signs = [1 if d > 0 else -1 if d < 0 else 0 for d in diffs]
    if all(s >= 0 for s in signs) or all(s <= 0 for s in signs):
        monotonicity = 1.0
    else:
        changes = sum(1 for i in range(1, len(signs)) if signs[i] != signs[i - 1])
        monotonicity = max(0.0, 1.0 - changes / 3.0)

    centroid_range = max(centroids) - min(centroids)
    gesture = min(1.0, centroid_range / 2000.0) * monotonicity

    composite = 0.35 * talk + 0.25 * ridge + 0.25 * coher + 0.15 * gesture

    return {
        "composite": composite,
        "talk": talk,
        "ridge": ridge,
        "coher": coher,
        "gesture": gesture,
        "peak_db": peak,
        "centroid_range": centroid_range,
        "monotonicity": monotonicity,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--top", type=int, default=6, help="Use top N heritage bodies")
    parser.add_argument("--out", type=str, default="vault/_splice_candidates",
                        help="Output directory")
    parser.add_argument("--threshold", type=float, default=0.50,
                        help="Minimum composite score to keep")
    args = parser.parse_args()

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Load heritage bodies
    sources = HERITAGE_RANKED[:args.top]
    bodies = {}
    for name in sources:
        path = CARTRIDGE_DIR / f"{name}.json"
        if path.exists():
            bodies[name] = Body.from_json(str(path))
            print(f"  loaded {name}")

    # Generate all splices
    candidates = []
    pairs = list(combinations(bodies.keys(), 2))
    total = len(pairs) * len(MODES) * 2  # x2 for A→B and B→A

    print(f"\n  generating {total} splice candidates from {len(bodies)} bodies...")

    for name_a, name_b in pairs:
        ba = bodies[name_a]
        bb = bodies[name_b]
        for mode in MODES:
            for src_a, src_b, label_a, label_b in [
                (ba, bb, name_a, name_b),
                (bb, ba, name_b, name_a),
            ]:
                splice_name = f"{label_a}_x_{label_b}_{mode.value}"
                corners = splice_corners(src_a.corners, src_b.corners, mode)
                spliced = Body(name=splice_name, corners=corners, boost=4.0)

                metrics = score_body(spliced)

                if metrics["peak_db"] > 50.0:
                    continue
                if metrics["composite"] < args.threshold:
                    continue

                candidates.append((splice_name, spliced, metrics))

    # Sort by composite
    candidates.sort(key=lambda x: -x[2]["composite"])

    # Report
    print(f"\n  {'name':<45} {'comp':>5} {'talk':>5} {'ridge':>5} {'coher':>5} {'gest':>5} {'peak':>5}")
    print(f"  {'─' * 45} {'─' * 5} {'─' * 5} {'─' * 5} {'─' * 5} {'─' * 5} {'─' * 5}")

    saved = 0
    for name, body, m in candidates[:30]:
        print(f"  {name:<45} {m['composite']:5.2f} {m['talk']:5.2f} {m['ridge']:5.2f} "
              f"{m['coher']:5.2f} {m['gesture']:5.2f} {m['peak_db']:5.1f}")

        # Save top candidates
        if saved < 20:
            out_path = out_dir / f"{name}.json"
            with open(out_path, "w") as f:
                f.write(body.to_compiled_json(provenance="heritage-splice"))
            saved += 1

    print(f"\n  {len(candidates)} candidates passed threshold ({args.threshold})")
    print(f"  {saved} saved to {out_dir}/")


if __name__ == "__main__":
    main()
