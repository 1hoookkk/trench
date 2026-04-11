"""Generate Speaker Knockerz candidates via the heritage compiler.

Uses the firmware integer grid (type, freq, gain) per stage, per morph endpoint.
This is the path that produces working P2K-grade coefficients with real zeros.

Spec: Bass sculpt. Cluster sweep.
- 4 poles tight at 80-200 Hz -> sweep to mids
- 2 poles at 2-4 kHz (presence)
- Interior zeros carving 400-800 Hz mud
- Morph pulls low cluster into mids
- Q sharpens

Heritage vocabulary (DillusionMan guide):
  type 1 = bandpass resonator
  type 2 = resonant formant (fixed numerator)
  type 3 = shaped asymmetric
  freq 0-127 = packed frequency (0=72Hz, 40=320Hz, 80=1440Hz, 100=3066Hz, 127=10259Hz)
  gain 0-127 = packed gain (64=unity, <64=cut, >64=boost)
"""
import json
import math
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pyruntime.designer_compile import (
    DesignerSection, legacy_sections_to_four_corner,
    compile_four_corner_to_body,
)
from pyruntime.heritage_coeffs import fw_freq_value
from pyruntime.body import Body
from pyruntime.analysis import midpoint_audit, morph_trajectory_distance
from pyruntime.validator import validate as validate_corners
from pyruntime.freq_response import cascade_response_db, freq_points
from pyruntime.constants import SR
import numpy as np

OUT_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "vault", "sk_candidates")
os.makedirs(OUT_DIR, exist_ok=True)


def packed_to_hz(packed: int) -> float:
    """Approximate Hz from packed frequency."""
    fv = fw_freq_value(packed)
    # fv maps to roughly: theta = pi * fv / 256 -> freq = SR/2 * fv/256
    # More accurate: from the type 1 bandpass, pole angle = acos(-c3/2)
    # Rough: fv 18=72Hz, 128=1440Hz, 236=10259Hz
    return 72.0 * (2.0 ** ((packed) / 20.0))  # rough exponential mapping


# Heritage frequency reference (from Cargo.toml comment):
# freq=0 → ~72Hz, freq=20 → ~152Hz, freq=40 → ~320Hz, freq=60 → ~686Hz
# freq=80 → ~1440Hz, freq=100 → ~3066Hz, freq=120 → ~7045Hz, freq=127 → ~10259Hz

# Speaker Knockerz template designs:
# Each is 6 sections: (type, low_freq, low_gain, high_freq, high_gain)
# low = morph 0 (closed/dark), high = morph 1 (open/bright)

CANDIDATES = [
    # SK_00: Classic sub cluster sweep
    # S0-S3: type 1 resonators at sub (freq 0-15), sweep to mids (freq 40-70)
    # S4-S5: type 2 formant anchors at presence (freq 90-100)
    {
        "name": "SK_00_Classic_Sub_Sweep",
        "sections": [
            (1, 5, 70, 50, 64),    # S0: 80Hz->320Hz, boosted at rest
            (1, 8, 68, 55, 64),    # S1: 100Hz->400Hz, slight boost
            (1, 12, 66, 60, 64),   # S2: 130Hz->500Hz
            (1, 18, 64, 68, 62),   # S3: 180Hz->800Hz, slight cut at open
            (2, 90, 58, 95, 64),   # S4: 2kHz formant, cut at rest
            (2, 100, 56, 105, 66), # S5: 3kHz formant, boosted at open
        ],
    },
    # SK_01: Wider spread, more aggressive
    {
        "name": "SK_01_Wide_Spread",
        "sections": [
            (1, 3, 72, 45, 64),    # S0: 70Hz->280Hz, heavy boost at rest
            (1, 10, 70, 58, 66),   # S1: 110Hz->440Hz
            (1, 15, 68, 65, 64),   # S2: 150Hz->600Hz
            (1, 22, 64, 75, 60),   # S3: 220Hz->1000Hz, cuts at open
            (2, 85, 60, 90, 68),   # S4: 1.8kHz formant
            (2, 105, 54, 110, 64), # S5: 3.5kHz formant
        ],
    },
    # SK_02: Tight sub, type 3 presence
    {
        "name": "SK_02_Tight_Sub_T3_Presence",
        "sections": [
            (1, 2, 74, 40, 64),    # S0: 65Hz->280Hz, heavy sub
            (1, 6, 72, 48, 64),    # S1: 85Hz->330Hz
            (1, 10, 70, 55, 66),   # S2: 110Hz->400Hz
            (1, 14, 68, 62, 64),   # S3: 140Hz->530Hz
            (3, 95, 60, 100, 66),  # S4: type 3 at 2.5kHz
            (3, 108, 56, 112, 64), # S5: type 3 at 4kHz
        ],
    },
    # SK_03: DillusionMan style — shelf-like low morph, belch at high
    {
        "name": "SK_03_DillusionMan_Belch",
        "sections": [
            (1, 4, 76, 55, 68),    # S0: sub->mid, big boost both ends
            (1, 7, 74, 60, 66),    # S1: sub->mid
            (1, 11, 70, 65, 64),   # S2: sub->mid
            (1, 16, 66, 72, 62),   # S3: sub->upper-mid
            (2, 88, 52, 93, 70),   # S4: formant punch at open
            (1, 95, 50, 100, 72),  # S5: presence peak at open
        ],
    },
    # SK_04: Minimal movement, maximum weight
    {
        "name": "SK_04_Minimal_Heavy",
        "sections": [
            (1, 3, 68, 30, 64),    # S0: sub stays low, subtle sweep
            (1, 6, 68, 35, 64),    # S1: similar
            (1, 10, 66, 40, 64),   # S2: slightly wider
            (1, 15, 66, 45, 64),   # S3: still conservative
            (2, 92, 60, 95, 66),   # S4: formant
            (2, 100, 58, 102, 64), # S5: formant
        ],
    },
    # SK_05: Crossing stages — S0 goes high, S3 goes low
    {
        "name": "SK_05_Stage_Cross",
        "sections": [
            (1, 5, 70, 70, 64),    # S0: sub->upper-mid (crosses S2/S3)
            (1, 8, 68, 50, 66),    # S1: sub->mid
            (1, 12, 66, 35, 68),   # S2: mid->low (reverse!)
            (1, 20, 64, 20, 64),   # S3: stationary anchor
            (2, 88, 58, 95, 68),   # S4: formant opens
            (2, 102, 54, 108, 64), # S5: presence
        ],
    },
    # SK_06: All type 2 (formant character)
    {
        "name": "SK_06_All_Formant",
        "sections": [
            (2, 5, 68, 45, 64),    # S0: formant sub->mid
            (2, 8, 68, 52, 64),    # S1: formant sub->mid
            (2, 12, 66, 58, 66),   # S2: formant
            (2, 16, 64, 65, 64),   # S3: formant
            (2, 90, 58, 95, 66),   # S4: presence formant
            (2, 100, 56, 105, 64), # S5: hi presence
        ],
    },
    # SK_07: Asymmetric gains — big cut at rest, big boost at open
    {
        "name": "SK_07_Asymmetric_Gains",
        "sections": [
            (1, 4, 48, 50, 76),    # S0: cut at rest, boost at open
            (1, 8, 50, 55, 74),    # S1: similar
            (1, 13, 52, 62, 72),   # S2: less extreme
            (1, 18, 54, 68, 70),   # S3: still asymmetric
            (2, 90, 64, 95, 64),   # S4: flat formant
            (2, 100, 64, 105, 64), # S5: flat formant
        ],
    },
    # SK_08: Mixed types — type 1 sub + type 3 mids
    {
        "name": "SK_08_Mixed_Types",
        "sections": [
            (1, 3, 70, 42, 64),    # S0: type 1 sub->mid
            (1, 7, 70, 50, 64),    # S1: type 1 sub->mid
            (3, 12, 66, 58, 66),   # S2: type 3 shaped
            (3, 18, 64, 65, 64),   # S3: type 3 shaped
            (2, 92, 58, 97, 66),   # S4: formant
            (1, 100, 56, 105, 64), # S5: resonator presence
        ],
    },
    # SK_09: Extreme spread — sub to 1.5kHz
    {
        "name": "SK_09_Extreme_Spread",
        "sections": [
            (1, 2, 72, 60, 64),    # S0: 60Hz->530Hz
            (1, 5, 70, 68, 64),    # S1: 80Hz->800Hz
            (1, 10, 68, 75, 62),   # S2: 110Hz->1kHz
            (1, 15, 66, 82, 60),   # S3: 150Hz->1.5kHz, cuts at open
            (2, 85, 62, 92, 68),   # S4: formant punch
            (2, 100, 54, 108, 66), # S5: presence
        ],
    },
]


def run_candidate(spec: dict) -> dict:
    name = spec["name"]
    sections = spec["sections"]

    # Build section tuples for legacy compiler
    section_tuples = [(t, lf, lg, hf, hg) for t, lf, lg, hf, hg in sections]

    # Compile via heritage path
    template = legacy_sections_to_four_corner(name, section_tuples)
    body = compile_four_corner_to_body(template, boost=4.0)

    # Gate
    audit = midpoint_audit(body)
    morph = morph_trajectory_distance(body)
    issues = validate_corners(body.corners)

    # Q diff
    freqs = freq_points()
    enc_q0 = body.corners.interpolate(0.5, 0.0)
    enc_q1 = body.corners.interpolate(0.5, 1.0)
    db_q0 = cascade_response_db(enc_q0, freqs, SR)
    db_q1 = cascade_response_db(enc_q1, freqs, SR)
    q_diff = float(np.sqrt(np.mean((db_q0 - db_q1) ** 2)))

    errors = [i for i in issues if i.severity == "error"]

    passed = (
        audit["passed"]
        and morph["mean_rms_db"] >= 3.0
        and len(errors) == 0
    )

    # Save
    path = os.path.join(OUT_DIR, f"{name}.json")
    with open(path, "w") as f:
        f.write(body.to_compiled_json(provenance="heritage-sk-batch"))

    return {
        "name": name,
        "passed": passed,
        "midpoint_peak_db": audit["worst_peak_db"],
        "midpoint_passed": audit["passed"],
        "morph_rms_db": round(morph["mean_rms_db"], 1),
        "q_diff_rms_db": round(q_diff, 1),
        "errors": len(errors),
        "path": path,
    }


if __name__ == "__main__":
    results = []
    for spec in CANDIDATES:
        r = run_candidate(spec)
        status = "PASS" if r["passed"] else "FAIL"
        print(f"[{status}] {r['name']:35s}  "
              f"mid={r['midpoint_peak_db']:+.1f}dB  "
              f"morph={r['morph_rms_db']:.1f}dB  "
              f"Q={r['q_diff_rms_db']:.1f}dB  "
              f"err={r['errors']}")
        results.append(r)

    passed = [r for r in results if r["passed"]]
    print(f"\n{'='*70}")
    print(f"{len(passed)}/10 passed all gates.")
    if passed:
        print("\nPassing candidates (ready for audition):")
        for r in passed:
            print(f"  {r['name']}  morph={r['morph_rms_db']:.1f}dB  Q={r['q_diff_rms_db']:.1f}dB")
            print(f"    -> {r['path']}")
