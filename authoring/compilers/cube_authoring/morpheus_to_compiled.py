"""Morpheus pole-byte cube → trench.compiled.cube_surface.v1.

Reads `trench_re_vault/artifacts/morpheus_cubes_decoded.json` (289 cubes,
8 corners × 7 stages, freq_byte+rad_byte → freq_hz+radius), and emits one
compiled-v1 cube_surface per cube to a staging directory.

Numerator policy: passthrough zeros (b0=1, b1=0, b2=0). Pure pole audition.
Stage padding: 7 morpheus stages, padded with passthrough to 12 output stages
(matches OUTPUT_STAGE_COUNT in tools/compile_cube.py).

Sample rate: morpheus freq_hz values reported at the original 44064 Hz
Morpheus sample rate. Cube_gate runs at 44100 Hz. We treat freq_hz as
absolute Hz and use SR=44100 in the bilinear pole→coeff formula. Resulting
~0.08% pitch shift at the very top of the band is below audibility for the
admission gate.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "tools"))

from pole_math import pole_to_kernel_stage  # noqa: E402

VAULT_DEFAULT = REPO_ROOT.parent / "trench_re_vault" / "artifacts" / "morpheus_cubes_decoded.json"
OUT_DEFAULT = REPO_ROOT / "cartridges" / "cubes" / "_morpheus"

CORNER_LABELS = ("c000", "c100", "c010", "c110", "c001", "c101", "c011", "c111")
RUNTIME_SR = 44100
ACTIVE_STAGE_COUNT = 6
OUTPUT_STAGE_COUNT = 12
PASSTHROUGH_STAGE = {"c0": 1.0, "c1": 0.0, "c2": 0.0, "c3": 0.0, "c4": 0.0}


def pole_to_coeffs(freq_hz: float, radius: float, sr: int) -> dict:
    c0, c1, c2, c3, c4 = pole_to_kernel_stage(freq_hz, radius, sr)
    return {"c0": c0, "c1": c1, "c2": c2, "c3": c3, "c4": c4}


def pad_stages(stages: list[dict]) -> list[dict]:
    out = list(stages[:OUTPUT_STAGE_COUNT])
    while len(out) < OUTPUT_STAGE_COUNT:
        out.append(dict(PASSTHROUGH_STAGE))
    return out


def convert_corner(corner_stages: list[dict]) -> dict:
    coeffs = [pole_to_coeffs(s["freq_hz"], s["radius"], RUNTIME_SR) for s in corner_stages]
    return {"boost": 1.0, "stages": pad_stages(coeffs)}


def convert_cube(cube: dict) -> dict:
    corners_in = cube["corners"]
    if len(corners_in) != 8:
        raise ValueError(f"cube {cube.get('index')} has {len(corners_in)} corners, expected 8")
    corners_out = {}
    resolution = {}
    for i, label in enumerate(CORNER_LABELS):
        corners_out[label] = convert_corner(corners_in[i])
        resolution[label] = {
            "source_kind": "morpheus_pole_bytes",
            "morpheus_corner_index": i,
            "stage_count": len(corners_in[i]),
            "authoring_sample_rate_hz": RUNTIME_SR,
        }
    cube_id = f"morph_cube_{cube['index']:03d}"
    return {
        "schema": "trench.compiled.cube_surface.v1",
        "name": cube_id,
        "derived_from": {
            "schema": "morpheus_cubes_decoded.v1",
            "id": cube_id,
            "morpheus_cube_index": cube["index"],
            "source": "trench_re_vault/artifacts/morpheus_cubes_decoded.json",
            "numerator_policy": "passthrough_zeros_b0_1",
        },
        "compiler": {"name": "morpheus_to_compiled.py", "version": "1"},
        "representation": "engine_ready_coeff_packs",
        "exactness": "morpheus_native_pole_bytes",
        "control_mode": "modern_live_xyz",
        "corner_resolution": resolution,
        "corners": corners_out,
    }


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--src", type=Path, default=VAULT_DEFAULT, help="decoded morpheus cubes JSON")
    parser.add_argument("--out", type=Path, default=OUT_DEFAULT, help="output staging directory")
    parser.add_argument("--limit", type=int, default=None, help="convert only first N cubes")
    parser.add_argument("--start", type=int, default=0, help="start index")
    args = parser.parse_args(argv)

    with args.src.open(encoding="utf-8") as fh:
        data = json.load(fh)
    cubes = data["cubes"]
    args.out.mkdir(parents=True, exist_ok=True)

    end = len(cubes) if args.limit is None else min(len(cubes), args.start + args.limit)
    written = 0
    for cube in cubes[args.start:end]:
        compiled = convert_cube(cube)
        out_path = args.out / f"morph_cube_{cube['index']:03d}.cube.json"
        out_path.write_text(json.dumps(compiled, indent=2) + "\n", encoding="utf-8")
        written += 1
    print(f"wrote {written} cubes to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(__import__("sys").argv[1:]))
