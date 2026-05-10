"""Convert compiled-v1 -> vault heritage format for trench-mcp consumption.

Inverse of compile_raw's encode_stage:
  c0 = 1 + val1     -> val1 = c0 - 1
  c1 = a1 + val2    -> val2 = c1 - c3
  c2 = r^2 - val3   -> val3 = c4 - c2
  c3 = a1
  c4 = r^2
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path


def convert_stage(s: dict) -> dict:
    c0, c1, c2, c3, c4 = (float(s[f"c{i}"]) for i in range(5))
    if c3 == 0.0 and c4 == 0.0:
        return {"a1": 0.0, "r": 0.0, "val1": 0.0, "val2": 0.0, "val3": 0.0, "flag": 1}
    return {
        "a1": c3,
        "r": math.sqrt(c4),
        "val1": c0 - 1.0,
        "val2": c1 - c3,
        "val3": c4 - c2,
        "flag": 1,
    }


def convert(in_path: Path, out_path: Path) -> None:
    cart = json.loads(in_path.read_text(encoding="utf-8"))
    out = {
        "name": cart.get("name", "unnamed"),
        "sampleRate": int(round(float(cart.get("sampleRate", 39062.5)))),
        "stages": 12,
        "keyframes": [],
    }
    for kf in cart["keyframes"]:
        out["keyframes"].append({
            "morph": kf["morph"],
            "q": kf["q"],
            "label": kf["label"],
            "boost": kf.get("boost", 1.0),
            "stages": [convert_stage(s) for s in kf["stages"]],
        })
    out_path.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(out_path)


if __name__ == "__main__":
    in_path = Path(sys.argv[1])
    out_path = Path(sys.argv[2])
    convert(in_path, out_path)
