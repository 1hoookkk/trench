"""Promote all P2K skins to generation-ready anatomy format.

Reads compiled P2K bodies (StageParams or kernel-form), reverse-solves
pole/zero positions, assigns roles, computes cascade metrics, and writes
a clean JSON per body suitable for constellation-aware generation.

Output schema (per body):
{
  "name": "P2k_012",
  "strategy": null,          # filled manually or by classifier
  "corners": {
    "M0_Q0": {
      "shared_c4": 0.56,
      "cascade_peak_db": 18.4,
      "cancellation_db": 41.3,
      "stages": [
        {
          "pole_freq_hz": 9320.9,
          "radius": 0.975,
          "role": "HF/Nyquist",
          "peak_db": 26.0,
          "zero_freq_hz": 4508.1,
          "zero_r": 0.701,
          "zero_family": "INTERIOR_ZERO"
        }, ...
      ]
    }, ...
  }
}

Usage:
    python tools/promote_p2k.py
    python tools/promote_p2k.py --source cartridges/p2k --out datasets/p2k_anatomy
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np

from pyruntime.body import Body
from pyruntime.constants import SR, TWO_PI, NUM_BODY_STAGES
from pyruntime.freq_response import cascade_response_db, freq_points
from pyruntime.corner import CornerName

FREQS = freq_points()

ROLE_BANDS = [
    (20, 200, "Anchor"),
    (200, 900, "LowMid"),
    (900, 3000, "FormantBite"),
    (3000, 6000, "Air"),
    (6000, 12000, "HF/Nyquist"),
    (12000, 20000, "Ceiling"),
]


def _freq_from_a1_r(a1: float, r: float) -> float:
    """Reverse-solve pole frequency from a1 and r."""
    if r < 0.001:
        return 0.0
    cos_theta = -a1 / (2.0 * r)
    cos_theta = max(-1.0, min(1.0, cos_theta))
    theta = math.acos(cos_theta)
    return theta * SR / TWO_PI


def _zero_from_vals(a1: float, r: float, val2: float, val3: float) -> dict:
    """Reverse-solve zero position from val2/val3.

    val2 = b1 - a1, val3 = r² - b2
    b1 = val2 + a1, b2 = r² - val3
    zero_r = sqrt(b2), zero_theta = acos(-b1 / (2*zero_r))
    """
    if r < 0.001:
        return {"zero_freq_hz": 0.0, "zero_r": 0.0, "zero_family": "NONE"}

    b1 = val2 + a1
    b2 = r * r - val3

    if b2 <= 0.001:
        return {"zero_freq_hz": 0.0, "zero_r": 0.0, "zero_family": "ALLPOLE"}

    zero_r = math.sqrt(b2)

    # Classify family
    if abs(zero_r - 1.0) < 0.02:
        family = "UNIT_CIRCLE"
    elif zero_r < r - 0.01:
        family = "INTERIOR_ZERO"
    elif zero_r > r + 0.01:
        family = "EXTERIOR_ZERO"
    else:
        family = "NEAR_POLE"

    denom = 2.0 * zero_r
    if denom < 0.001:
        return {"zero_freq_hz": 0.0, "zero_r": round(zero_r, 6), "zero_family": family}

    cos_z = -b1 / denom
    cos_z = max(-1.0, min(1.0, cos_z))
    zero_theta = math.acos(cos_z)
    zero_freq = zero_theta * SR / TWO_PI

    return {
        "zero_freq_hz": round(zero_freq, 1),
        "zero_r": round(zero_r, 6),
        "zero_family": family,
    }


def _assign_role(freq_hz: float) -> str:
    for lo, hi, role in ROLE_BANDS:
        if lo <= freq_hz < hi:
            return role
    if freq_hz < 20:
        return "DC"
    return "Ceiling"


def _stage_peak_db(a1: float, r: float, val1: float, val2: float, val3: float, c4: float) -> tuple[float, float]:
    """Compute single-stage peak dB and peak frequency.

    Uses the encoded kernel form to evaluate |H(z)|.
    """
    from pyruntime.encode import raw_to_encoded
    from pyruntime.stage_params import StageParams

    sp = StageParams(a1=a1, r=r, val1=val1, val2=val2, val3=val3)
    enc = raw_to_encoded(sp)

    # Override c0 with the shared c4 (which IS the shared b0)
    z_inv = np.exp(-2j * np.pi * FREQS / SR)
    z2 = z_inv * z_inv
    num = c4 + enc.c1 * z_inv + enc.c2 * z2
    den = 1.0 + enc.c3 * z_inv + enc.c4 * z2
    H = np.abs(num / np.where(np.abs(den) > 1e-30, den, 1e-30))
    db = 20.0 * np.log10(H + 1e-30)

    peak_idx = int(np.argmax(db))
    return round(float(db[peak_idx]), 1), round(float(FREQS[peak_idx]), 1)


def _is_passthrough(stages_data: dict) -> bool:
    """Check if a stage is passthrough from its raw values."""
    if "c0" in stages_data:
        return (abs(stages_data["c0"] - 1.0) < 0.01 and
                abs(stages_data.get("c1", 0)) < 0.01 and
                abs(stages_data.get("c2", 0)) < 0.01)
    return abs(stages_data.get("r", 0)) < 0.01 and abs(stages_data.get("a1", 0)) < 0.01


def promote_body(body: Body, name: str) -> dict:
    """Promote a Body to anatomy format."""
    corner_map = {
        CornerName.A: "M0_Q0",
        CornerName.B: "M0_Q100",
        CornerName.C: "M100_Q0",
        CornerName.D: "M100_Q100",
    }

    result = {"name": name, "strategy": None, "corners": {}}

    for cn, label in corner_map.items():
        corner = body.corners.corner(cn)
        encoded = corner.encode()

        # Find shared c4 (= c0 = b0) from first active stage
        shared_c4 = None
        for enc in encoded:
            if abs(enc.c0 - 1.0) > 0.01 and abs(enc.c0) > 0.001:
                shared_c4 = round(enc.c0, 6)
                break
        if shared_c4 is None:
            shared_c4 = 1.0

        # Cascade response for this corner
        db_cascade = cascade_response_db(encoded, FREQS, SR)
        cascade_peak = round(float(np.max(db_cascade)), 1)

        stages_out = []
        active_peaks = []

        for i, (sp, enc) in enumerate(zip(corner.stages, encoded)):
            if sp.r < 0.01 and abs(sp.a1) < 0.01:
                continue  # skip passthrough

            freq = round(_freq_from_a1_r(sp.a1, sp.r), 1)
            role = _assign_role(freq)
            zeros = _zero_from_vals(sp.a1, sp.r, sp.val2, sp.val3)

            # Single-stage peak
            z_inv = np.exp(-2j * np.pi * FREQS / SR)
            z2 = z_inv * z_inv
            num = enc.c0 + enc.c1 * z_inv + enc.c2 * z2
            den = 1.0 + enc.c3 * z_inv + enc.c4 * z2
            H = np.abs(num / np.where(np.abs(den) > 1e-30, den, 1e-30))
            stage_db = 20.0 * np.log10(H + 1e-30)
            peak_db = round(float(np.max(stage_db)), 1)
            active_peaks.append(peak_db)

            stages_out.append({
                "pole_freq_hz": freq,
                "radius": round(sp.r, 6),
                "role": role,
                "peak_db": peak_db,
                "zero_freq_hz": zeros["zero_freq_hz"],
                "zero_r": zeros["zero_r"],
                "zero_family": zeros["zero_family"],
            })

        # Cancellation = max single stage peak - cascade peak
        max_single = max(active_peaks) if active_peaks else 0.0
        cancellation = round(max_single - cascade_peak, 1)

        result["corners"][label] = {
            "shared_c4": shared_c4,
            "cascade_peak_db": cascade_peak,
            "cancellation_db": cancellation,
            "stages": stages_out,
        }

    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", default="cartridges/p2k")
    parser.add_argument("--out", default="datasets/p2k_anatomy")
    args = parser.parse_args()

    src = Path(args.source)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    files = sorted(src.glob("*.json"))
    print(f"  Promoting {len(files)} bodies from {src} → {out}\n")

    for f in files:
        name = f.stem
        try:
            body = Body.from_json(str(f))
            anatomy = promote_body(body, name)

            with open(out / f"{name}.json", "w") as fp:
                json.dump(anatomy, fp, indent=2)

            # Summary line
            m0 = anatomy["corners"]["M0_Q0"]
            n_active = len(m0["stages"])
            freqs = [s["pole_freq_hz"] for s in m0["stages"]]
            fmin = min(freqs) if freqs else 0
            fmax = max(freqs) if freqs else 0
            families = set(s["zero_family"] for s in m0["stages"])

            print(f"  {name:<12} {n_active} stages  {fmin:>6.0f}–{fmax:>6.0f} Hz  "
                  f"peak={m0['cascade_peak_db']:>5.1f}  cancel={m0['cancellation_db']:>5.1f}  "
                  f"c4={m0['shared_c4']:.4f}  zeros={','.join(sorted(families))}")

        except Exception as e:
            print(f"  {name:<12} FAILED: {e}")

    print(f"\n  Done. {len(files)} bodies written to {out}/")


if __name__ == "__main__":
    main()
