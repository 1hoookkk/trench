"""Mass candidate generator for all 4 shipping bodies.

Usage: python tools/mass_generate.py <body_name> [--count 50]

Bodies: speaker_knockerz, ear_bender, nite_mode, ah_ay_ee

Each body has a search space defined by:
- Per-stage type (1/2/3) ranges
- Per-stage freq ranges (0-127 packed)
- Per-stage gain ranges (0-127 packed)
- Frame A (low morph) and Frame B (high morph) independently

The heritage compiler (type1/2/3_compile) produces coefficients with
real zeros — the only path that produces stable bodies.
"""
import argparse
import json
import math
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np

from pyruntime.designer_compile import legacy_sections_to_four_corner, compile_four_corner_to_body
from pyruntime.body import Body
from pyruntime.analysis import midpoint_audit, morph_trajectory_distance
from pyruntime.validator import validate as validate_corners
from pyruntime.freq_response import cascade_response_db, freq_points
from pyruntime.constants import SR

OUT_BASE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "vault")


# =========================================================================
# SEARCH SPACES
#
# Each body defines per-stage constraints for Frame A and Frame B.
# Format: list of 6 stages, each stage is:
#   {
#     "types": [1,2,3],           # allowed types
#     "a_freq": (lo, hi),         # Frame A freq range (packed 0-127)
#     "a_gain": (lo, hi),         # Frame A gain range
#     "b_freq": (lo, hi),         # Frame B freq range
#     "b_gain": (lo, hi),         # Frame B gain range
#   }
# =========================================================================

SEARCH_SPACES = {
    "speaker_knockerz": {
        "desc": "Bass sculpt. Sub cluster sweeps to mids. Presence anchors.",
        "stages": [
            # S0-S3: sub cluster at rest, sweeps to mids at open
            {"types": [1, 2], "a_freq": (2, 12),  "a_gain": (64, 80), "b_freq": (40, 70), "b_gain": (56, 68)},
            {"types": [1, 2], "a_freq": (5, 16),  "a_gain": (62, 78), "b_freq": (45, 72), "b_gain": (58, 70)},
            {"types": [1, 2], "a_freq": (8, 22),  "a_gain": (60, 74), "b_freq": (50, 76), "b_gain": (58, 68)},
            {"types": [1, 2], "a_freq": (12, 28), "a_gain": (58, 70), "b_freq": (55, 82), "b_gain": (56, 66)},
            # S4-S5: presence/bite anchors
            {"types": [1, 2, 3], "a_freq": (82, 100), "a_gain": (50, 64), "b_freq": (88, 108), "b_gain": (60, 72)},
            {"types": [1, 2, 3], "a_freq": (92, 110), "a_gain": (48, 62), "b_freq": (96, 115), "b_gain": (58, 70)},
        ],
    },
    "ear_bender": {
        "desc": "HF wall fold. 4 HF stages fold down to 800Hz. Complex morph.",
        "stages": [
            # S0: low anchor stays put
            {"types": [1, 2], "a_freq": (25, 45), "a_gain": (60, 72), "b_freq": (20, 40), "b_gain": (58, 68)},
            # S1-S4: HF wall at rest, folds to mids at open
            {"types": [1, 2, 3], "a_freq": (95, 115), "a_gain": (58, 70), "b_freq": (55, 80), "b_gain": (60, 72)},
            {"types": [1, 2, 3], "a_freq": (100, 120), "a_gain": (56, 68), "b_freq": (50, 75), "b_gain": (58, 70)},
            {"types": [1, 2, 3], "a_freq": (105, 122), "a_gain": (56, 68), "b_freq": (45, 72), "b_gain": (60, 72)},
            {"types": [1, 2, 3], "a_freq": (108, 125), "a_gain": (54, 66), "b_freq": (40, 68), "b_gain": (58, 70)},
            # S5: ceiling suppressor
            {"types": [1, 3], "a_freq": (115, 127), "a_gain": (48, 60), "b_freq": (60, 90), "b_gain": (50, 64)},
        ],
    },
    "nite_mode": {
        "desc": "Collapsing ceiling. Dark shelf muffler. Not a lowpass.",
        "stages": [
            # S0-S1: low shelf anchors (stay dark)
            {"types": [1, 2], "a_freq": (10, 30), "a_gain": (62, 72), "b_freq": (8, 25), "b_gain": (56, 66)},
            {"types": [1, 2], "a_freq": (15, 40), "a_gain": (60, 70), "b_freq": (12, 35), "b_gain": (54, 64)},
            # S2-S3: mid presence — alive at rest, collapses at open
            {"types": [1, 2], "a_freq": (55, 80), "a_gain": (62, 74), "b_freq": (30, 55), "b_gain": (50, 62)},
            {"types": [1, 2, 3], "a_freq": (65, 90), "a_gain": (60, 72), "b_freq": (35, 60), "b_gain": (48, 60)},
            # S4-S5: HF ceiling — present at rest, crushed at open
            {"types": [1, 2, 3], "a_freq": (90, 115), "a_gain": (58, 70), "b_freq": (70, 95),  "b_gain": (40, 56)},
            {"types": [1, 3],    "a_freq": (100, 125), "a_gain": (56, 68), "b_freq": (80, 105), "b_gain": (38, 54)},
        ],
    },
    "ah_ay_ee": {
        "desc": "Vowel transition aa->ay->ee. Formant morph.",
        "stages": [
            # F1: 700Hz (aa) -> 400Hz (ee) — F1 drops
            {"types": [1, 2], "a_freq": (58, 68), "a_gain": (64, 76), "b_freq": (38, 50), "b_gain": (62, 74)},
            # F2: 1100Hz (aa) -> 2300Hz (ee) — F2 rises dramatically
            {"types": [1, 2], "a_freq": (72, 82), "a_gain": (62, 74), "b_freq": (90, 102), "b_gain": (62, 74)},
            # F3: 2500Hz (aa) -> 3000Hz (ee) — F3 slight rise
            {"types": [1, 2], "a_freq": (88, 98), "a_gain": (58, 70), "b_freq": (92, 104), "b_gain": (58, 70)},
            # F4: nasal anti-formant anchor
            {"types": [2, 3], "a_freq": (80, 95), "a_gain": (50, 62), "b_freq": (82, 98), "b_gain": (50, 62)},
            # F5: HF coloring
            {"types": [1, 2], "a_freq": (100, 115), "a_gain": (52, 64), "b_freq": (102, 118), "b_gain": (52, 64)},
            # Sub anchor
            {"types": [1, 2], "a_freq": (5, 18), "a_gain": (60, 70), "b_freq": (5, 18), "b_gain": (58, 68)},
        ],
    },
}


def generate_candidate(body_name: str, idx: int, rng: random.Random) -> dict:
    """Generate one candidate from the search space."""
    space = SEARCH_SPACES[body_name]
    sections = []

    for stage_spec in space["stages"]:
        typ = rng.choice(stage_spec["types"])
        a_freq = rng.randint(*stage_spec["a_freq"])
        a_gain = rng.randint(*stage_spec["a_gain"])
        b_freq = rng.randint(*stage_spec["b_freq"])
        b_gain = rng.randint(*stage_spec["b_gain"])
        sections.append((typ, a_freq, a_gain, b_freq, b_gain))

    safe_name = f"{body_name}_{idx:03d}"
    template = legacy_sections_to_four_corner(safe_name, sections)
    body = compile_four_corner_to_body(template, boost=4.0)

    # Gates
    audit = midpoint_audit(body)
    morph = morph_trajectory_distance(body)
    issues = validate_corners(body.corners)
    errors = [i for i in issues if i.severity == "error"]

    # Q diff
    freqs = freq_points()
    enc_q0 = body.corners.interpolate(0.5, 0.0)
    enc_q1 = body.corners.interpolate(0.5, 1.0)
    db_q0 = cascade_response_db(enc_q0, freqs, SR)
    db_q1 = cascade_response_db(enc_q1, freqs, SR)
    q_diff = float(np.sqrt(np.mean((db_q0 - db_q1) ** 2)))

    passed = (
        audit["passed"]
        and morph["mean_rms_db"] >= 3.0
        and len(errors) == 0
    )

    return {
        "name": safe_name,
        "passed": passed,
        "midpoint_peak_db": audit["worst_peak_db"],
        "morph_rms_db": round(morph["mean_rms_db"], 1),
        "q_diff_rms_db": round(q_diff, 1),
        "errors": len(errors),
        "sections": sections,
        "body": body,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("body", choices=list(SEARCH_SPACES.keys()))
    parser.add_argument("--count", type=int, default=50)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    out_dir = os.path.join(OUT_BASE, f"{args.body}_candidates")
    os.makedirs(out_dir, exist_ok=True)

    rng = random.Random(args.seed)
    results = []
    passing = []

    for i in range(args.count):
        r = generate_candidate(args.body, i, rng)
        results.append(r)
        status = "PASS" if r["passed"] else "FAIL"
        print(f"[{status}] {r['name']:30s}  mid={r['midpoint_peak_db']:+.1f}dB  morph={r['morph_rms_db']:.1f}dB  Q={r['q_diff_rms_db']:.1f}dB")

        if r["passed"]:
            path = os.path.join(out_dir, f"{r['name']}.json")
            with open(path, "w") as f:
                f.write(r["body"].to_compiled_json(provenance=f"mass-gen-{args.body}"))
            passing.append({**{k: v for k, v in r.items() if k != "body"}, "path": path})

    # Sort passing by morph distance (most motion first)
    passing.sort(key=lambda x: -x["morph_rms_db"])

    print(f"\n{'='*70}")
    print(f"{len(passing)}/{args.count} passed all gates for {args.body}")
    if passing:
        print(f"\nTop 10 by morph motion:")
        for r in passing[:10]:
            secs = r["sections"]
            types_str = "".join(str(s[0]) for s in secs)
            print(f"  {r['name']:30s}  morph={r['morph_rms_db']:.1f}dB  types={types_str}  mid={r['midpoint_peak_db']:+.1f}dB")

    # Save summary
    summary_path = os.path.join(out_dir, "_summary.json")
    with open(summary_path, "w") as f:
        json.dump({
            "body": args.body,
            "total": args.count,
            "passed": len(passing),
            "top_10": passing[:10],
        }, f, indent=2)
    print(f"\nSummary: {summary_path}")


if __name__ == "__main__":
    main()
