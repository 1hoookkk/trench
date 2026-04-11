"""Generate 10 Speaker Knockerz candidates from spec constraints.

Spec: Bass sculpt. Cluster Sweep strategy.
- 4 poles tight at 80-200 Hz
- 2 poles at 2-4 kHz
- Interior zeros carving 400-800 Hz mud
- Morph pulls low cluster into mids
- Q sharpens

Each candidate varies placement within the spec bounds.
All are compiled, gated, and reported.
"""
import json
import math
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np

from pyruntime.forge_generator import (
    Blueprint, Stage, blueprint_to_body,
    solve_c4_surface, apply_gain_budget, cascade_peak_db,
    solve_zero_network,
)
from pyruntime.analysis import midpoint_audit, morph_trajectory_distance
from pyruntime.validator import validate as validate_corners
from pyruntime.freq_response import cascade_response_db, freq_points
from pyruntime.constants import SR
from pyruntime.corner import CornerName

OUT_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "vault", "sk_candidates")
os.makedirs(OUT_DIR, exist_ok=True)


def make_candidate(idx: int, seed: int) -> dict:
    """Build one Speaker Knockerz candidate."""
    rng = random.Random(seed)

    # Spec constraints:
    # S0-S3: 4 poles tight at 80-200 Hz (sub cluster)
    # S4-S5: 2 poles at 2-4 kHz (presence/bite anchors)
    # Morph: sub cluster sweeps UP into mids (300-1200 Hz range)
    # Q: sharpens all poles (radius increases)

    # Sub cluster: spread across 60-250 Hz with moderate radii
    # BassBox 303 reference: S0=79Hz r=0.962, not 0.997
    sub_base = rng.uniform(55, 90)
    sub_spread = rng.uniform(30, 80)
    s0_hz = sub_base
    s1_hz = sub_base + sub_spread * rng.uniform(0.8, 1.2)
    s2_hz = sub_base + sub_spread * rng.uniform(1.5, 2.5)
    s3_hz = sub_base + sub_spread * rng.uniform(2.5, 4.0)

    # Presence anchors: 2-4 kHz, moderate radius
    s4_hz = rng.uniform(2000, 3500)
    s5_hz = rng.uniform(3500, 5000)

    # Radii: moderate — BassBox range, not P2K ceiling
    sub_r = rng.uniform(0.960, 0.985)
    pres_r = rng.uniform(0.950, 0.980)

    anchor = [
        Stage(s0_hz, sub_r),
        Stage(s1_hz, sub_r - rng.uniform(0, 0.008)),
        Stage(s2_hz, sub_r - rng.uniform(0, 0.010)),
        Stage(s3_hz, sub_r - rng.uniform(0.005, 0.020)),
        Stage(s4_hz, pres_r),
        Stage(s5_hz, pres_r - rng.uniform(0, 0.015)),
    ]

    # Morph: sub cluster fans out to 200-1500 Hz (wide spread to avoid pileup)
    morph_targets = sorted([
        rng.uniform(200, 500),
        rng.uniform(350, 800),
        rng.uniform(500, 1000),
        rng.uniform(700, 1500),
    ])
    morph = [
        (morph_targets[0] - s0_hz, -rng.uniform(0.01, 0.04)),
        (morph_targets[1] - s1_hz, -rng.uniform(0.01, 0.03)),
        (morph_targets[2] - s2_hz, -rng.uniform(0.01, 0.03)),
        (morph_targets[3] - s3_hz, -rng.uniform(0.01, 0.04)),
        (rng.uniform(-800, 800), rng.uniform(-0.01, 0.005)),
        (rng.uniform(-800, 800), rng.uniform(-0.01, 0.005)),
    ]

    # Q: moderate sharpening + frequency compression toward 500 Hz
    q_sharpen = rng.uniform(0.003, 0.010)
    q_center = rng.uniform(300, 800)
    q = [
        (rng.uniform(-30, 30), q_sharpen),
        ((q_center - s1_hz) * rng.uniform(0.1, 0.3), q_sharpen),
        ((q_center - s2_hz) * rng.uniform(0.1, 0.3), q_sharpen),
        ((q_center - s3_hz) * rng.uniform(0.1, 0.3), q_sharpen),
        (rng.uniform(-500, 500), q_sharpen * 0.3),
        (rng.uniform(-500, 500), q_sharpen * 0.3),
    ]

    bp = Blueprint(
        name=f"Speaker_Knockerz_{idx:02d}",
        sentence="Bass sculpt. 4-pole sub cluster sweeps to mids. Presence anchors. Q sharpens.",
        anchor=anchor,
        morph=morph,
        q=q,
        jitter_hz=5.0,
        jitter_r=0.001,
    )

    # Compile
    body = blueprint_to_body(bp, rng, boost=4.0)
    body = solve_zero_network(body)
    body = solve_c4_surface(body, target_db=36.0)
    body = apply_gain_budget(body)

    # Gate
    audit = midpoint_audit(body)
    morph_dist = morph_trajectory_distance(body)
    issues = validate_corners(body.corners)

    # Q differentiation
    freqs = freq_points()
    enc_q0 = body.corners.interpolate(0.5, 0.0)
    enc_q1 = body.corners.interpolate(0.5, 1.0)
    db_q0 = cascade_response_db(enc_q0, freqs, SR)
    db_q1 = cascade_response_db(enc_q1, freqs, SR)
    q_diff = float(np.sqrt(np.mean((db_q0 - db_q1) ** 2)))

    # Sub anchor gate: peak in 20-60 Hz must persist
    enc_m0 = body.corners.interpolate(0.0, 0.0)
    enc_m100 = body.corners.interpolate(1.0, 0.0)
    db_m0 = cascade_response_db(enc_m0, freqs, SR)
    db_m100 = cascade_response_db(enc_m100, freqs, SR)
    sub_mask = (freqs >= 20) & (freqs <= 80)
    sub_m0_peak = float(np.max(db_m0[sub_mask])) if np.any(sub_mask) else -120
    sub_m100_peak = float(np.max(db_m100[sub_mask])) if np.any(sub_mask) else -120
    sub_drop = sub_m0_peak - sub_m100_peak

    # Peaks
    peaks = cascade_peak_db(body)

    errors = [i for i in issues if i.severity == "error"]
    warnings = [i for i in issues if i.severity == "warning"]

    passed = (
        audit["passed"]
        and morph_dist["mean_rms_db"] >= 3.0
        and q_diff >= 2.0
        and len(errors) == 0
        and sub_drop < 12.0  # sub shouldn't vanish completely
    )

    # Save body
    path = os.path.join(OUT_DIR, f"Speaker_Knockerz_{idx:02d}.json")
    with open(path, "w") as f:
        f.write(body.to_compiled_json(provenance="hostile-sk-batch"))

    return {
        "idx": idx,
        "seed": seed,
        "name": bp.name,
        "passed": passed,
        "midpoint_peak_db": audit["worst_peak_db"],
        "midpoint_passed": audit["passed"],
        "morph_rms_db": round(morph_dist["mean_rms_db"], 1),
        "q_diff_rms_db": round(q_diff, 1),
        "sub_drop_db": round(sub_drop, 1),
        "sub_m0_peak_db": round(sub_m0_peak, 1),
        "peak_db": round(max(peaks.values()), 1),
        "warnings": len(warnings),
        "errors": len(errors),
        "anchor_hz": [round(s.freq_hz, 1) for s in anchor],
        "anchor_r": [round(s.radius, 4) for s in anchor],
        "path": path,
    }


if __name__ == "__main__":
    results = []
    for i in range(10):
        seed = 1000 + i * 137
        r = make_candidate(i, seed)
        status = "PASS" if r["passed"] else "FAIL"
        print(f"[{status}] #{i:02d} seed={seed}  "
              f"mid={r['midpoint_peak_db']:+.1f}dB  "
              f"morph={r['morph_rms_db']:.1f}dB  "
              f"Q={r['q_diff_rms_db']:.1f}dB  "
              f"sub_drop={r['sub_drop_db']:.1f}dB  "
              f"peak={r['peak_db']:+.1f}dB  "
              f"warn={r['warnings']}")
        print(f"       hz={r['anchor_hz']}  r={r['anchor_r']}")
        results.append(r)

    passed = [r for r in results if r["passed"]]
    print(f"\n{'='*70}")
    print(f"{len(passed)}/10 passed all gates.")
    if passed:
        print("\nPassing candidates:")
        for r in passed:
            print(f"  #{r['idx']:02d} — morph={r['morph_rms_db']:.1f}dB  Q={r['q_diff_rms_db']:.1f}dB  sub_drop={r['sub_drop_db']:.1f}dB")
            print(f"         {r['path']}")
