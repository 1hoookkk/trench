#!/usr/bin/env python3
"""Compile a pill authoring path into a `compiled-v1` cartridge.

Input: `trench.authoring_path.pill.v1`
Output: `compiled-v1` (4-keyframe pill), the same format the shipping phoneme
pills use (see `cartridges/engine/**/*.json`).

This is the pill analogue of `compile_cube.py`. A pill is a 2D compiled
surface (morph × Q, 4 keyframes, bilinear); a cube is 3D (X×Y×Z, 8 corners,
trilinear). If the filter's honest identity is 2D, ship as pill; see
CUBE_GATE.md §11 (demotion to pill) and the E-mu `.4` convention.

Corner resolution reuses `compile_cube.resolve_corner` — same `kind`
vocabulary (morph_designer_inline, morph_designer_ref, peak_shelf_inline,
peak_shelf_ref), same sample semantics.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from compile_cube import CompileError, load_doc, resolve_corner

KEYFRAME_LABELS = ("M0_Q0", "M100_Q0", "M0_Q100", "M100_Q100")
LABEL_TO_MQ = {
    "M0_Q0":     (0.0, 0.0),
    "M100_Q0":   (1.0, 0.0),
    "M0_Q100":   (0.0, 1.0),
    "M100_Q100": (1.0, 1.0),
}
DEFAULT_SAMPLE_RATE = 44100


def compile_pill(doc: dict) -> dict:
    if doc.get("schema") != "trench.authoring_path.pill.v1":
        raise CompileError(
            f"unsupported pill authoring schema: {doc.get('schema')!r}"
        )
    keyframes_in = doc.get("keyframes")
    if not isinstance(keyframes_in, dict):
        raise CompileError("pill authoring requires a keyframes object")

    compiled_keyframes = []
    resolution = {}
    for label in KEYFRAME_LABELS:
        if label not in keyframes_in:
            raise CompileError(f"pill authoring missing keyframe {label!r}")
        compiled_corner, corner_resolution = resolve_corner(keyframes_in[label], {})
        morph, q = LABEL_TO_MQ[label]
        compiled_keyframes.append({
            "label": label,
            "morph": morph,
            "q": q,
            "boost": compiled_corner["boost"],
            "stages": compiled_corner["stages"],
        })
        resolution[label] = corner_resolution

    return {
        "format": "compiled-v1",
        "name": doc.get("name", doc.get("id", "pill")),
        "sampleRate": doc.get("sampleRate", DEFAULT_SAMPLE_RATE),
        "provenance": {
            "schema": doc["schema"],
            "id": doc.get("id"),
            "name": doc.get("name"),
            "scope": doc.get("scope"),
            "dimensions": 2,
            "compiler": {"name": "compile_pill.py", "version": "1"},
            "exactness": doc.get("exactness", "modern_cleanroom_not_native_verified"),
            "keyframe_resolution": resolution,
            "source_provenance": doc.get("provenance"),
        },
        "keyframes": compiled_keyframes,
    }


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(
        description="Compile trench.authoring_path.pill.v1 into a compiled-v1 pill.",
    )
    parser.add_argument(
        "pill", nargs="?", default="-",
        help="Path to pill authoring JSON, or '-' for stdin.",
    )
    parser.add_argument("--out", type=Path, default=None,
                         help="Output file path (stdout if omitted).")
    args = parser.parse_args(argv)

    try:
        path = None if args.pill == "-" else Path(args.pill)
        authored = load_doc(path)
        compiled = compile_pill(authored)
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
