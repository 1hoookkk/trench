#!/usr/bin/env python3
"""Compile an authored 2-token grid into a compiled-v1 cartridge.

Input: JSON grid with exactly 2 tokens (start + end).
Output: compiled-v1 cartridge JSON that loads via trench-core::Cartridge::from_json.

Spec: docs/superpowers/specs/compiler_grid_to_cartridge.md
Decision log: §5 Path A, strict 2-token limit.

Stdlib only. No RBJ cookbook derivation — shape coefficients are taken
verbatim from `cartridges/engine/_source/shapes/<category>/<key>.json`.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
INVENTORY_PATH = REPO_ROOT / "cartridges" / "engine" / "_source" / "token_inventory_unified_v2.json"
SHAPES_ROOT = REPO_ROOT / "cartridges" / "engine" / "_source" / "shapes"

CORNER_LABELS = ("M0_Q0", "M0_Q100", "M100_Q0", "M100_Q100")
FIRST_TOKEN_LABELS = ("M0_Q0", "M0_Q100")
LAST_TOKEN_LABELS = ("M100_Q0", "M100_Q100")


class CompileError(Exception):
    pass


def load_inventory() -> dict:
    if not INVENTORY_PATH.exists():
        raise CompileError(f"token inventory missing: {INVENTORY_PATH}")
    with INVENTORY_PATH.open() as f:
        inv = json.load(f)
    tokens = inv.get("tokens")
    if not isinstance(tokens, dict):
        raise CompileError("inventory has no 'tokens' dict")
    return tokens


def resolve_shape_path(token_entry: dict) -> Path:
    """Compute shape file path from category + key. Ignores stale `shape_path`."""
    category = token_entry.get("category")
    key = token_entry.get("key")
    if not category or not key:
        raise CompileError(f"token entry missing category/key: {token_entry}")
    return SHAPES_ROOT / category / f"{key}.json"


def load_shape(token_name: str, inventory: dict) -> dict:
    entry = inventory.get(token_name)
    if entry is None:
        raise CompileError(f"unknown token: {token_name!r}")
    path = resolve_shape_path(entry)
    if not path.exists():
        raise CompileError(f"shape file missing for token {token_name!r}: {path}")
    with path.open() as f:
        shape = json.load(f)
    if shape.get("format") != "compiled-v1":
        raise CompileError(
            f"shape file {path} is not compiled-v1 (got {shape.get('format')!r})"
        )
    return shape


def extract_keyframe(shape: dict, label: str, source_token: str) -> dict:
    for kf in shape.get("keyframes", []):
        if kf.get("label") == label:
            return kf
    raise CompileError(
        f"shape for token {source_token!r} missing required keyframe {label!r}"
    )


def compile_grid(grid_doc: dict) -> dict:
    name = grid_doc.get("name", "unnamed_body")
    grid = grid_doc.get("grid")
    if not isinstance(grid, list):
        raise CompileError("grid document must have a 'grid' list")
    if len(grid) != 2:
        raise CompileError(
            f"v1 compiler requires exactly 2 tokens (start + end); got {len(grid)}. "
            f"Multi-waypoint compilation is gated on compiled-v2 schema evolution "
            f"(see docs/superpowers/specs/compiler_grid_to_cartridge.md §5)."
        )

    first_entry, last_entry = grid[0], grid[1]
    first_token = first_entry.get("token")
    last_token = last_entry.get("token")
    if not first_token or not last_token:
        raise CompileError("each grid entry must have a 'token' field")

    inventory = load_inventory()
    first_shape = load_shape(first_token, inventory)
    last_shape = load_shape(last_token, inventory)

    out_keyframes = []
    for label in FIRST_TOKEN_LABELS:
        kf = extract_keyframe(first_shape, label, first_token)
        out_keyframes.append(
            {
                "label": label,
                "boost": kf.get("boost", 1.0),
                "stages": kf["stages"],
            }
        )
    for label in LAST_TOKEN_LABELS:
        kf = extract_keyframe(last_shape, label, last_token)
        out_keyframes.append(
            {
                "label": label,
                "boost": kf.get("boost", 1.0),
                "stages": kf["stages"],
            }
        )

    out_labels = {kf["label"] for kf in out_keyframes}
    missing = set(CORNER_LABELS) - out_labels
    if missing:
        raise CompileError(f"output cartridge missing corners: {sorted(missing)}")

    return {
        "format": "compiled-v1",
        "name": name,
        "provenance": "compile_grid",
        "grid_source": [first_token, last_token],
        "keyframes": out_keyframes,
    }


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(
        description="Compile a 2-token grid into a compiled-v1 cartridge."
    )
    parser.add_argument(
        "grid",
        type=Path,
        help="Path to grid JSON file, or '-' for stdin.",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=None,
        help="Output cartridge path (stdout if omitted).",
    )
    args = parser.parse_args(argv)

    try:
        if str(args.grid) == "-":
            grid_doc = json.load(sys.stdin)
        else:
            with args.grid.open() as f:
                grid_doc = json.load(f)
    except json.JSONDecodeError as e:
        print(f"grid parse error: {e}", file=sys.stderr)
        return 1
    except OSError as e:
        print(f"grid read error: {e}", file=sys.stderr)
        return 1

    try:
        cart = compile_grid(grid_doc)
    except CompileError as e:
        print(f"compile error: {e}", file=sys.stderr)
        return 2

    out_json = json.dumps(cart, indent=2)
    if args.out is None:
        sys.stdout.write(out_json + "\n")
    else:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(out_json + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
