"""
calibration_to_body.py — Convert calibration JSON directly to compiled-v1 body.

No heritage compiler. No back-solving packed integers. The calibration data
IS the reverse engineering — just convert StageParams → EncodedCoeffs → JSON.

Usage:
    python tools/calibration_to_body.py docs/calibration/Ear_Bender.json
    python tools/calibration_to_body.py docs/calibration/*.json --out vault/sift_queue
"""

import json
import math
import os
import sys
import argparse
from glob import glob

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pyruntime.stage_params import StageParams
from pyruntime.encode import EncodedCoeffs, raw_to_encoded
from pyruntime.corner import CornerName, CornerState, CornerArray
from pyruntime.body import Body
from pyruntime.constants import NUM_BODY_STAGES


def calibration_to_body(cal_path: str) -> Body:
    with open(cal_path) as f:
        cal = json.load(f)

    name = cal["name"]
    boost = cal.get("boost", 4.0)

    corner_map = {}
    for cn in CornerName:
        key = cn.json_key()
        raw_stages = cal["corners"][key]["stages"]

        stages = []
        encoded = []
        for s in raw_stages:
            sp = StageParams(
                a1=s["val2"] - (-2.0 * s["radius"] * math.cos(2.0 * math.pi * s["pole_freq_hz"] / 39062.5)),
                r=s["radius"],
                val1=s["val1"],
                val2=s["val2"],
                val3=s["val3"],
            )
            # Recompute a1 properly from pole position
            sp = StageParams(
                a1=-2.0 * s["radius"] * math.cos(2.0 * math.pi * s["pole_freq_hz"] / 39062.5),
                r=s["radius"],
                val1=s["val1"],
                val2=s["val2"],
                val3=s["val3"],
            )
            stages.append(sp)
            enc = raw_to_encoded(sp)
            # Override c4 with the calibration's shared b0
            enc = EncodedCoeffs(c0=enc.c0, c1=enc.c1, c2=enc.c2, c3=enc.c3, c4=s["c4_b0"])
            encoded.append(enc)

        # Pad to 12 stages
        while len(stages) < NUM_BODY_STAGES:
            stages.append(StageParams.passthrough())
            encoded.append(EncodedCoeffs(c0=1.0, c1=0.0, c2=0.0, c3=0.0, c4=0.0))

        corner_map[cn] = CornerState(stages=stages, boost=boost, _pre_encoded=encoded)

    corners = CornerArray(
        a=corner_map[CornerName.A],
        b=corner_map[CornerName.B],
        c=corner_map[CornerName.C],
        d=corner_map[CornerName.D],
    )
    return Body(name=name, corners=corners, boost=boost)


def main():
    parser = argparse.ArgumentParser(description="Convert calibration JSON to compiled-v1 body")
    parser.add_argument("inputs", nargs="+", help="Calibration JSON file(s)")
    parser.add_argument("--out", type=str, default="vault/sift_queue", help="Output directory")
    args = parser.parse_args()

    files = []
    for pattern in args.inputs:
        files.extend(glob(pattern))

    os.makedirs(args.out, exist_ok=True)
    for path in files:
        body = calibration_to_body(path)
        out_path = os.path.join(args.out, f"P2K_{body.name.replace(' ', '_')}.json")
        with open(out_path, "w") as f:
            f.write(body.to_compiled_json(provenance="calibration-direct"))
        print(f"  {body.name} → {out_path}")


if __name__ == "__main__":
    main()
