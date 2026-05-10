"""Mine morpheus cubes for face-pill candidates.

Each cube has 6 axis-aligned faces; each face is a 4-corner square — i.e. a
2D filter, the natural compiled-v1 pill format. This carves out every face
of every audited cube, hard-pre-filters using the cube's native edge
verdict, and emits compiled-v1 pills for survivors.

Pre-filter (hard, per user instruction):
  - drop face if ≥3 of its 4 face-edges are *universally* dead in the cube
    native verdict (i.e. abrupt_collapse hit on every probe of that edge)
  - drop face if its corresponding face-quadrant in the cube has any corner
    whose stage 0 pole hits radius > 0.997 (near-stability-edge, blow_up
    source — usually unrecoverable on solo audition)

Each survivor → cartridges/pills/_morpheus/morph_<cube>_<face>.compiled.json
in the compiled-v1 format consumed by pill_gate.

Pill labels chosen by face: lowest-index varying axis = morph, higher = Q.
Numerator policy: passthrough zeros (b0=1, b1=0, b2=0).
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "tools" / "cube_authoring"))

DEFAULT_VAULT = REPO_ROOT.parent / "trench_re_vault" / "artifacts" / "morpheus_cubes_decoded.json"
DEFAULT_VERDICTS = REPO_ROOT / "cartridges" / "cubes" / "_morpheus" / "_verdicts_native"
DEFAULT_OUT = REPO_ROOT / "cartridges" / "pills" / "_morpheus"

RUNTIME_SR = 44100
ACTIVE_STAGES = 7   # morpheus has 7 pole stages
OUTPUT_STAGES = 12  # pad to 12 like compile_pill / compile_cube
PASSTHROUGH = {"c0": 1.0, "c1": 0.0, "c2": 0.0, "c3": 0.0, "c4": 0.0}
RADIUS_BLOWUP_THRESHOLD = 0.997

CORNER_LABELS = ("c000", "c100", "c010", "c110", "c001", "c101", "c011", "c111")
CORNER_INDEX = {label: i for i, label in enumerate(CORNER_LABELS)}

# Each face: (face_id, varying_axis_lo, varying_axis_hi, fixed_axis, fixed_value,
#             [4 cube corner labels in (M0_Q0, M100_Q0, M0_Q100, M100_Q100) order],
#             [4 face-edges as canonical EDGES tuples])
# Convention: lower-index axis = morph (X<Y<Z).
FACES = [
    # Z=0 face (vary X, Y), Z fixed at 0
    ("z0", "X", "Y", "Z", 0,
     ("c000", "c100", "c010", "c110"),
     [("c000", "c100"), ("c010", "c110"), ("c000", "c010"), ("c100", "c110")]),
    # Z=1 face
    ("z1", "X", "Y", "Z", 1,
     ("c001", "c101", "c011", "c111"),
     [("c001", "c101"), ("c011", "c111"), ("c001", "c011"), ("c101", "c111")]),
    # Y=0 face (vary X, Z)
    ("y0", "X", "Z", "Y", 0,
     ("c000", "c100", "c001", "c101"),
     [("c000", "c100"), ("c001", "c101"), ("c000", "c001"), ("c100", "c101")]),
    # Y=1 face
    ("y1", "X", "Z", "Y", 1,
     ("c010", "c110", "c011", "c111"),
     [("c010", "c110"), ("c011", "c111"), ("c010", "c011"), ("c110", "c111")]),
    # X=0 face (vary Y, Z)
    ("x0", "Y", "Z", "X", 0,
     ("c000", "c010", "c001", "c011"),
     [("c000", "c010"), ("c001", "c011"), ("c000", "c001"), ("c010", "c011")]),
    # X=1 face
    ("x1", "Y", "Z", "X", 1,
     ("c100", "c110", "c101", "c111"),
     [("c100", "c110"), ("c101", "c111"), ("c100", "c101"), ("c110", "c111")]),
]
PILL_KF_LABELS = ("M0_Q0", "M100_Q0", "M0_Q100", "M100_Q100")


def pole_to_stage(freq_hz: float, radius: float, sr: int) -> dict:
    theta = 2.0 * math.pi * (freq_hz / sr)
    return {"c0": 1.0, "c1": 0.0, "c2": 0.0,
            "c3": -2.0 * radius * math.cos(theta),
            "c4": radius * radius}


def corner_to_keyframe_stages(corner_in: list[dict], sr: int) -> list[dict]:
    stages = [pole_to_stage(s["freq_hz"], s["radius"], sr) for s in corner_in]
    while len(stages) < OUTPUT_STAGES:
        stages.append(dict(PASSTHROUGH))
    return stages[:OUTPUT_STAGES]


def corner_max_radius(corner_in: list[dict]) -> float:
    return max(float(s["radius"]) for s in corner_in)


def edge_universally_collapsed(edge_label: str, native_verdict: dict) -> bool:
    """True iff edge has abrupt_collapse fail on every probe in native verdict."""
    probes = [p for p in native_verdict.get("per_probe", [])
               if p.get("status") in ("passed", "failed")]
    if not probes:
        return False
    for p in probes:
        match = next((e for e in p.get("edges", []) if e["edge"] == edge_label), None)
        if match is None or "abrupt_collapse" not in match.get("fails", []):
            return False
    return True


def face_passes_filter(face, cube_in, native_verdict) -> tuple[bool, str]:
    _, _, _, _, _, corner_labels, edges = face
    # radius blow-up check
    for cl in corner_labels:
        idx = CORNER_INDEX[cl]
        if corner_max_radius(cube_in["corners"][idx]) > RADIUS_BLOWUP_THRESHOLD:
            return False, f"blowup_radius>{RADIUS_BLOWUP_THRESHOLD}@{cl}"
    if native_verdict is None:
        return True, "no_verdict"
    bad_edges = sum(1 for a, b in edges
                     if edge_universally_collapsed(f"{a}-{b}", native_verdict))
    if bad_edges >= 3:
        return False, f"{bad_edges}/4_edges_universally_collapsed"
    return True, "kept"


def build_face_pill(face, cube_in, sr: int) -> dict:
    face_id, lo_axis, hi_axis, fixed_axis, fixed_val, corner_labels, _ = face
    cube_index = cube_in["index"]
    keyframes = []
    for kf_label, cl in zip(PILL_KF_LABELS, corner_labels):
        morph = 0.0 if "M0" in kf_label else 1.0
        q = 0.0 if "Q0" in kf_label.split("_")[1] else 1.0
        idx = CORNER_INDEX[cl]
        stages = corner_to_keyframe_stages(cube_in["corners"][idx], sr)
        keyframes.append({"label": kf_label, "morph": morph, "q": q,
                           "boost": 1.0, "stages": stages,
                           "source_cube_corner": cl})
    pill_id = f"morph_{cube_index:03d}_face_{face_id}"
    return {
        "format": "compiled-v1",
        "name": pill_id,
        "sampleRate": sr,
        "provenance": {
            "id": pill_id,
            "schema": "morpheus_face_pill.v1",
            "source_cube_index": cube_index,
            "face": {"id": face_id, "morph_axis": lo_axis, "q_axis": hi_axis,
                       "fixed_axis": fixed_axis, "fixed_value": fixed_val},
            "numerator_policy": "passthrough_zeros_b0_1",
            "scope": "hf_resonance",
        },
        "keyframes": keyframes,
    }


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--src", type=Path, default=DEFAULT_VAULT)
    parser.add_argument("--verdicts-dir", type=Path, default=DEFAULT_VERDICTS)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--require-verdict", action="store_true",
                         help="Skip cubes that have no native verdict (default: include with no_verdict tag)")
    args = parser.parse_args(argv)

    args.out.mkdir(parents=True, exist_ok=True)
    with args.src.open(encoding="utf-8") as fh:
        cubes = json.load(fh)["cubes"]
    end = len(cubes) if args.limit is None else min(len(cubes), args.start + args.limit)
    cubes = cubes[args.start:end]

    log_path = args.out / "_mine_log.csv"
    with log_path.open("w", newline="", encoding="utf-8") as logfh:
        writer = csv.writer(logfh)
        writer.writerow(["cube_index", "face", "kept", "reason", "out_path"])

        kept = 0
        dropped = 0
        skipped_no_verdict = 0
        for cube in cubes:
            ci = cube["index"]
            vp = args.verdicts_dir / f"morph_cube_{ci:03d}.verdict_native.json"
            verdict = None
            if vp.exists():
                with vp.open(encoding="utf-8") as fh:
                    verdict = json.load(fh)
            elif args.require_verdict:
                skipped_no_verdict += 1
                continue
            for face in FACES:
                ok, reason = face_passes_filter(face, cube, verdict)
                if ok:
                    pill = build_face_pill(face, cube, RUNTIME_SR)
                    out_path = args.out / f"{pill['name']}.compiled.json"
                    out_path.write_text(json.dumps(pill, indent=2) + "\n", encoding="utf-8")
                    writer.writerow([ci, face[0], "kept", reason, str(out_path.name)])
                    kept += 1
                else:
                    writer.writerow([ci, face[0], "dropped", reason, ""])
                    dropped += 1
    print(f"audited cubes: {len(cubes)}")
    print(f"kept faces:   {kept}")
    print(f"dropped:      {dropped}")
    if skipped_no_verdict:
        print(f"skipped (no verdict, --require-verdict): {skipped_no_verdict}")
    print(f"log → {log_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
