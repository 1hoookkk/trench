"""Crack P2K — take heritage bodies and push them harder.

Reads a P2K body, amplifies its existing character:
- Wider morph sweep (exaggerate freq deltas between M0 and M100)
- Tighter radii (push active poles toward 0.999)
- Deeper cancellation (push zeros closer to poles they oppose)
- More Q differentiation (exaggerate Q corner differences)

The cascade interaction is heritage-correct by construction.
We just turn up the dials on what's already there.

Usage:
    python tools/crack_p2k.py --body P2k_012 --intensity 1.5
    python tools/crack_p2k.py --all --intensity 1.3
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
from pyruntime.corner import CornerArray, CornerName, CornerState
from pyruntime.stage_params import StageParams
from pyruntime.encode import raw_to_encoded, EncodedCoeffs
from pyruntime.freq_response import cascade_response_db, freq_points
from pyruntime.forge_generator import solve_c4_surface, apply_gain_budget

FREQS = freq_points()
CARTRIDGE_DIR = Path(__file__).parent.parent / "cartridges" / "p2k"
OUT_DIR = Path(__file__).parent.parent / "vault" / "_cracked"
MAX_R = 0.9995
MIN_R = 0.005


def _freq_from_a1_r(a1: float, r: float) -> float:
    if r < 0.001:
        return 0.0
    cos_theta = max(-1.0, min(1.0, -a1 / (2.0 * r)))
    return math.acos(cos_theta) * SR / TWO_PI


def _a1_from_freq_r(freq_hz: float, r: float) -> float:
    freq_hz = max(20.0, min(SR / 2.0 - 1.0, freq_hz))
    return -2.0 * r * math.cos(TWO_PI * freq_hz / SR)


def _is_active(sp: StageParams) -> bool:
    return sp.r > 0.01 or abs(sp.a1) > 0.01


def crack_corner(
    corner: CornerState,
    ref_corner: CornerState,
    intensity: float,
    radius_push: float,
) -> CornerState:
    """Crack one corner relative to a reference (M0_Q0).

    - Exaggerate frequency deltas from reference
    - Push radii tighter
    - Preserve zero structure (val2/val3 scale with pole changes)
    """
    new_stages = []
    new_encoded = []

    for i, (sp, ref_sp) in enumerate(zip(corner.stages, ref_corner.stages)):
        if not _is_active(sp):
            new_stages.append(sp)
            new_encoded.append(raw_to_encoded(sp))
            continue

        freq = _freq_from_a1_r(sp.a1, sp.r)
        ref_freq = _freq_from_a1_r(ref_sp.a1, ref_sp.r)

        # Exaggerate freq delta from reference
        freq_delta = freq - ref_freq
        new_freq = ref_freq + freq_delta * intensity
        new_freq = max(20.0, min(SR / 2.0 - 1.0, new_freq))

        # Push radius tighter (toward MAX_R for active, preserve lossy stages)
        if sp.r > 0.5:
            headroom = MAX_R - sp.r
            new_r = min(MAX_R, sp.r + headroom * radius_push)
        else:
            new_r = sp.r  # lossy/diffuser stage — leave it

        new_a1 = _a1_from_freq_r(new_freq, new_r)

        # Scale val2/val3 proportionally to maintain zero-pole relationship
        # val2 = b1 - a1, val3 = r² - b2
        # When a1 changes, b1 should track: new_b1 = val2 + new_a1
        # When r changes, b2 should track: new_b2 = new_r² - val3
        # This preserves the zero frequency and adjusts zero radius with pole radius
        new_val2 = sp.val2  # zero angle offset stays (preserves zero freq relative to pole)
        new_val3 = new_r * new_r - (sp.r * sp.r - sp.val3)  # adjust b2 for new r

        new_sp = StageParams(
            a1=new_a1,
            r=new_r,
            val1=sp.val1,
            val2=new_val2,
            val3=new_val3,
        )
        new_stages.append(new_sp)
        new_encoded.append(raw_to_encoded(new_sp))

    return CornerState(stages=new_stages, boost=corner.boost, _pre_encoded=new_encoded)


def crack_body(body: Body, intensity: float, radius_push: float, target_db: float) -> Body:
    """Crack an entire body.

    intensity: >1 exaggerates morph/Q deltas. 1.0 = no change. 2.0 = double.
    radius_push: 0-1 fraction of headroom to add. 0 = no change. 1.0 = max.
    """
    ref = body.corners.corner(CornerName.A)  # M0_Q0 is the reference

    corners = {}
    for cn in CornerName:
        c = body.corners.corner(cn)
        corners[cn] = crack_corner(c, ref, intensity, radius_push)

    # The reference corner itself only gets radius push, no freq exaggeration
    corners[CornerName.A] = crack_corner(ref, ref, 1.0, radius_push)

    cracked = Body(
        name=body.name + "_cracked",
        corners=CornerArray(
            a=corners[CornerName.A],
            b=corners[CornerName.B],
            c=corners[CornerName.C],
            d=corners[CornerName.D],
        ),
        boost=body.boost,
    )

    cracked = solve_c4_surface(cracked, target_db=target_db)
    cracked = apply_gain_budget(cracked)
    return cracked


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--body", type=str, default=None, help="e.g. P2k_012")
    parser.add_argument("--all", action="store_true")
    parser.add_argument("--intensity", type=float, default=1.5, help="Morph/Q exaggeration (1.0=stock, 2.0=double)")
    parser.add_argument("--radius", type=float, default=0.3, help="Radius push (0=stock, 1=max)")
    parser.add_argument("--target-db", type=float, default=36.0, help="c4 solver target peak dB")
    args = parser.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    if args.all:
        files = sorted(CARTRIDGE_DIR.glob("*.json"))
    elif args.body:
        files = [CARTRIDGE_DIR / f"{args.body}.json"]
    else:
        print("  Specify --body P2k_NNN or --all")
        return

    _VOCAL_MASK = (FREQS >= 200) & (FREQS <= 5000)

    print(f"\n  Cracking {len(files)} bodies  intensity={args.intensity}  radius_push={args.radius}")
    print(f"  {'name':<20} {'orig_peak':>9} {'crack_peak':>10} {'orig_dyn':>8} {'crack_dyn':>9} {'orig_dist':>9} {'crack_dist':>10}")
    print(f"  {'─'*20} {'─'*9} {'─'*10} {'─'*8} {'─'*9} {'─'*9} {'─'*10}")

    for f in files:
        name = f.stem
        try:
            orig = Body.from_json(str(f))
            cracked = crack_body(orig, args.intensity, args.radius, args.target_db)

            # Compare
            def _metrics(b):
                db0 = cascade_response_db(b.corners.interpolate(0, 0.5), FREQS, SR)
                db1 = cascade_response_db(b.corners.interpolate(1, 0.5), FREQS, SR)
                dbm = cascade_response_db(b.corners.interpolate(0.5, 0.5), FREQS, SR)
                peak = max(float(np.max(db0)), float(np.max(db1)), float(np.max(dbm)))
                dyn = float(np.max(dbm) - np.min(dbm))
                dist = float(np.sqrt(np.mean((db0 - db1)**2)))
                return peak, dyn, dist

            op, od, odist = _metrics(orig)
            cp, cd, cdist = _metrics(cracked)

            out_name = f"{name}_x{args.intensity:.1f}"
            cracked.name = out_name
            with open(OUT_DIR / f"{out_name}.json", "w") as fp:
                fp.write(cracked.to_compiled_json(provenance="crack-p2k"))

            print(f"  {out_name:<20} {op:>9.1f} {cp:>10.1f} {od:>8.1f} {cd:>9.1f} {odist:>9.1f} {cdist:>10.1f}")

        except Exception as e:
            print(f"  {name:<20} FAILED: {e}")

    print(f"\n  Done → {OUT_DIR}/")


if __name__ == "__main__":
    main()
