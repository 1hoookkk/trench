#!/usr/bin/env python3
"""Compile a cube authoring path into a cube surface cartridge.

Input: `trench.authoring_path.cube.v1`
Output: `trench.compiled.cube_surface.v1`

This is a deterministic bridge. It resolves each cube corner to a concrete
heritage object and lowers that heritage object into the existing engine-ready
kernel packs already used elsewhere in TRENCH.
"""

from __future__ import annotations

import argparse
import json
import sys
from copy import deepcopy
from pathlib import Path

from bake_hedz_const import DEFAULT_INVENTORY, compile_corner, load_template
from compile_grid import load_shape
from pole_math import pole_to_kernel_stage

REPO_ROOT = Path(__file__).resolve().parent.parent
CORNER_LABELS = ("c000", "c100", "c010", "c110", "c001", "c101", "c011", "c111")
DEFAULT_EXACTNESS = "modern_cleanroom_not_native_verified"
DEFAULT_CONTROL_MODE = "modern_live_xyz"
COMPILE_SR = 44100
ACTIVE_STAGE_COUNT = 6
OUTPUT_STAGE_COUNT = 12
PASSTHROUGH_STAGE = {"c0": 1.0, "c1": 0.0, "c2": 0.0, "c3": 0.0, "c4": 0.0}
MAX_NATIVE_POLES_STAGES = OUTPUT_STAGE_COUNT


class CompileError(Exception):
    pass


def load_doc(path: Path | None) -> dict:
    try:
        if path is None:
            return json.load(sys.stdin)
        with path.open(encoding="utf-8") as handle:
            return json.load(handle)
    except json.JSONDecodeError as exc:
        raise CompileError(f"cube authoring parse error: {exc}") from exc
    except OSError as exc:
        raise CompileError(f"cube authoring read error: {exc}") from exc


def sample_from_payload(payload: dict) -> dict:
    sample = payload.get("sample", {})
    if not isinstance(sample, dict):
        raise CompileError("heritage sample must be an object")
    morph = float(sample.get("morph", 0.0))
    q = float(sample.get("q", 0.0))
    return {"morph": morph, "q": q}


def stage_list_to_objects(stages: list[list[float] | dict]) -> list[dict]:
    out = []
    for stage in stages:
        if isinstance(stage, dict):
            out.append(
                {
                    "c0": float(stage["c0"]),
                    "c1": float(stage["c1"]),
                    "c2": float(stage["c2"]),
                    "c3": float(stage["c3"]),
                    "c4": float(stage["c4"]),
                }
            )
            continue

        if not isinstance(stage, list) or len(stage) != 5:
            raise CompileError("compiled stage must be an object or 5-element list")
        out.append(
            {
                "c0": float(stage[0]),
                "c1": float(stage[1]),
                "c2": float(stage[2]),
                "c3": float(stage[3]),
                "c4": float(stage[4]),
            }
        )
    return out


def reject_runtime_fields(value: object, path: str) -> None:
    if isinstance(value, dict):
        if value.get("format") == "compiled-v1":
            raise CompileError(f"{path} must not embed compiled-v1 payloads")
        if "keyframes" in value or "representation" in value or "corners" in value:
            raise CompileError(f"{path} must not embed compiled/runtime surface payloads")
        for key, nested in value.items():
            if key in {"c0", "c1", "c2", "c3", "c4"}:
                raise CompileError(f"{path} must not embed runtime coefficients")
            reject_runtime_fields(nested, f"{path}.{key}")
    elif isinstance(value, list):
        for index, nested in enumerate(value):
            reject_runtime_fields(nested, f"{path}[{index}]")


def pad_stages(stages: list[dict]) -> list[dict]:
    out = [deepcopy(stage) for stage in stages[:ACTIVE_STAGE_COUNT]]
    while len(out) < ACTIVE_STAGE_COUNT:
        out.append(deepcopy(PASSTHROUGH_STAGE))
    while len(out) < OUTPUT_STAGE_COUNT:
        out.append(deepcopy(PASSTHROUGH_STAGE))
    return out


def interpolate_compiled_shape(doc: dict, sample: dict) -> tuple[list[dict], float]:
    if doc.get("format") != "compiled-v1":
        raise CompileError(
            f"peak/shelf heritage must resolve to compiled-v1, got {doc.get('format')!r}"
        )

    keyframes = {kf["label"]: kf for kf in doc.get("keyframes", [])}
    try:
        m0_q0 = stage_list_to_objects(keyframes["M0_Q0"]["stages"])
        m100_q0 = stage_list_to_objects(keyframes["M100_Q0"]["stages"])
        m0_q100 = stage_list_to_objects(keyframes["M0_Q100"]["stages"])
        m100_q100 = stage_list_to_objects(keyframes["M100_Q100"]["stages"])
    except KeyError as exc:
        raise CompileError(f"compiled-v1 heritage missing keyframe {exc}") from exc

    morph = sample["morph"]
    q = sample["q"]
    out = []
    for stage_index in range(ACTIVE_STAGE_COUNT):
        stage_out = {}
        for coeff in ("c0", "c1", "c2", "c3", "c4"):
            q_m0 = m0_q0[stage_index][coeff] + (m0_q100[stage_index][coeff] - m0_q0[stage_index][coeff]) * q
            q_m1 = m100_q0[stage_index][coeff] + (m100_q100[stage_index][coeff] - m100_q0[stage_index][coeff]) * q
            stage_out[coeff] = q_m0 + (q_m1 - q_m0) * morph
        out.append(stage_out)

    boosts = {
        "M0_Q0": float(keyframes["M0_Q0"].get("boost", 1.0)),
        "M100_Q0": float(keyframes["M100_Q0"].get("boost", 1.0)),
        "M0_Q100": float(keyframes["M0_Q100"].get("boost", 1.0)),
        "M100_Q100": float(keyframes["M100_Q100"].get("boost", 1.0)),
    }
    q_m0 = boosts["M0_Q0"] + (boosts["M0_Q100"] - boosts["M0_Q0"]) * q
    q_m1 = boosts["M100_Q0"] + (boosts["M100_Q100"] - boosts["M100_Q0"]) * q
    boost = q_m0 + (q_m1 - q_m0) * morph

    return pad_stages(out), boost


def compile_morph_designer_template(template: dict, sample: dict) -> tuple[list[dict], float]:
    stages = compile_corner(template["sections"], morph=sample["morph"], sr=COMPILE_SR)
    compiled = [
        {"c0": float(c0), "c1": float(c1), "c2": float(c2), "c3": float(c3), "c4": float(c4)}
        for (c0, c1, c2, c3, c4) in stages
    ]
    return pad_stages(compiled), 1.0


def compile_native_poles(payload: dict) -> tuple[list[dict], float, dict]:
    """ROM-native pole coordinates → kernel-form SOS rows, no Hz/Q reconstruction.

    Payload shape:
        {"kind": "native_poles", "sr": 44100, "boost": 1.0,
         "stages": [{"freq_hz": ..., "radius": ...}, ...]}
    """
    stages_in = payload.get("stages")
    if not isinstance(stages_in, list) or not stages_in:
        raise CompileError("native_poles requires non-empty stages list")
    if len(stages_in) > MAX_NATIVE_POLES_STAGES:
        raise CompileError(
            f"native_poles supports at most {MAX_NATIVE_POLES_STAGES} stages, got {len(stages_in)}"
        )
    try:
        sr = float(payload.get("sr", COMPILE_SR))
    except (TypeError, ValueError) as exc:
        raise CompileError(f"native_poles sr must be numeric: {exc}") from exc
    if sr <= 0.0:
        raise CompileError(f"native_poles sr must be positive, got {sr}")
    try:
        boost = float(payload.get("boost", 1.0))
    except (TypeError, ValueError) as exc:
        raise CompileError(f"native_poles boost must be numeric: {exc}") from exc

    compiled: list[dict] = []
    for index, stage in enumerate(stages_in):
        if not isinstance(stage, dict):
            raise CompileError(f"native_poles.stages[{index}] must be an object")
        if "freq_hz" not in stage or "radius" not in stage:
            raise CompileError(
                f"native_poles.stages[{index}] requires freq_hz and radius"
            )
        try:
            freq_hz = float(stage["freq_hz"])
            radius = float(stage["radius"])
        except (TypeError, ValueError) as exc:
            raise CompileError(
                f"native_poles.stages[{index}] freq_hz/radius must be numeric: {exc}"
            ) from exc
        c0, c1, c2, c3, c4 = pole_to_kernel_stage(freq_hz, radius, sr)
        compiled.append({"c0": c0, "c1": c1, "c2": c2, "c3": c3, "c4": c4})

    resolution = {
        "source_kind": "native_poles",
        "sr": sr,
        "stage_count": len(stages_in),
    }
    return pad_stages(compiled), boost, resolution


def resolve_peak_shelf_ref(path_text: str) -> dict:
    path = Path(path_text)
    if not path.is_absolute():
        path = REPO_ROOT / path
    try:
        with path.open(encoding="utf-8") as handle:
            return json.load(handle)
    except OSError as exc:
        raise CompileError(f"peak/shelf heritage read error: {exc}") from exc


def resolve_corner(payload: dict, inventory: dict) -> tuple[dict, dict]:
    kind = payload.get("kind")
    if not kind:
        raise CompileError("cube corner payload missing kind")

    if kind == "native_poles":
        stages, boost, resolution = compile_native_poles(payload)
        return {"boost": boost, "stages": stages}, resolution

    sample = sample_from_payload(payload)
    if kind == "morph_designer_ref":
        template_name = payload.get("template_name")
        if not template_name:
            raise CompileError("morph_designer_ref requires template_name")
        template = load_template(DEFAULT_INVENTORY, template_name)
        stages, boost = compile_morph_designer_template(template, sample)
        resolution = {
            "source_kind": kind,
            "template_name": template_name,
            "sample": sample,
            "authoring_sample_rate_hz": COMPILE_SR,
        }
    elif kind == "morph_designer_inline":
        template = payload.get("template")
        if not isinstance(template, dict) or "sections" not in template:
            raise CompileError("morph_designer_inline requires template.sections")
        stages, boost = compile_morph_designer_template(template, sample)
        resolution = {
            "source_kind": kind,
            "template_name": template.get("name", "inline"),
            "sample": sample,
            "authoring_sample_rate_hz": COMPILE_SR,
        }
    elif kind == "peak_shelf_ref":
        path_text = payload.get("path")
        if not path_text:
            raise CompileError("peak_shelf_ref requires path")
        doc = resolve_peak_shelf_ref(path_text)
        stages, boost = interpolate_compiled_shape(doc, sample)
        resolution = {
            "source_kind": kind,
            "path": path_text,
            "sample": sample,
        }
    elif kind == "peak_shelf_inline":
        document = payload.get("document")
        if not isinstance(document, dict):
            raise CompileError("peak_shelf_inline requires document object")
        reject_runtime_fields(document, "peak_shelf_inline.document")
        stages, boost = interpolate_compiled_shape(document, sample)
        resolution = {
            "source_kind": kind,
            "sample": sample,
        }
    else:
        raise CompileError(f"unsupported heritage object kind: {kind!r}")

    return {"boost": boost, "stages": stages}, resolution


def compile_cube(doc: dict) -> dict:
    if doc.get("schema") != "trench.authoring_path.cube.v1":
        raise CompileError(
            f"unsupported cube authoring schema: {doc.get('schema')!r}"
        )

    corners = doc.get("corners")
    if not isinstance(corners, dict):
        raise CompileError("cube authoring requires a corners object")

    compiled_corners = {}
    resolution = {}
    for label in CORNER_LABELS:
        if label not in corners:
            raise CompileError(f"cube authoring missing corner {label!r}")
        compiled_corner, corner_resolution = resolve_corner(corners[label], {})
        compiled_corners[label] = compiled_corner
        resolution[label] = corner_resolution

    return {
        "schema": "trench.compiled.cube_surface.v1",
        "name": doc.get("name", doc.get("id", "cube_surface")),
        "derived_from": {
            "schema": doc["schema"],
            "id": doc.get("id"),
            "name": doc.get("name"),
            "scope": doc.get("scope"),
            "provenance": doc.get("provenance"),
        },
        "compiler": {
            "name": "compile_cube.py",
            "version": "1",
        },
        "representation": "engine_ready_coeff_packs",
        "exactness": doc.get("exactness", DEFAULT_EXACTNESS),
        "control_mode": doc.get("control_mode_default", DEFAULT_CONTROL_MODE),
        "corner_resolution": resolution,
        "corners": compiled_corners,
    }


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(
        description="Compile trench.authoring_path.cube.v1 into trench.compiled.cube_surface.v1."
    )
    parser.add_argument(
        "cube",
        nargs="?",
        default="-",
        help="Path to cube authoring JSON, or '-' for stdin.",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=None,
        help="Output file path (stdout if omitted).",
    )
    args = parser.parse_args(argv)

    try:
        path = None if args.cube == "-" else Path(args.cube)
        authored = load_doc(path)
        compiled = compile_cube(authored)
    except CompileError as exc:
        print(f"compile error: {exc}", file=sys.stderr)
        return 2

    output = json.dumps(compiled, indent=2)
    if args.out is None:
        sys.stdout.write(output + "\n")
    else:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(output + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
