"""Bake a user-authored Morph Designer body → compiled-v1 JSON for the live plugin hot-swap path.

Reads a sections JSON file (the same shape as the inner payload of
`heritage_designer_sections.json` templates) and writes a 4-corner
compiled-v1 cartridge JSON that the engine's `loadLiveOverrideFromJson`
can consume via the `trench_live.json` poller.

This is the UI-side authoring bake: the MorphDesignerPanel in JUCE edits
sections in memory, drops them to a temp JSON, and invokes this script
to produce the compiled cartridge. Q axis is degenerate for 2-frame
Morph Designer bodies — corners 0/2 and corners 1/3 share coefficients
(`build_cartridge` in `bake_hedz_const.py` already does this duplication
as `[m0, m1, m0, m1]`).

The compile math comes from `tools.bake_hedz_const` — `compile_corner` is
the proven path (gated by `trench-core/tests/hedz_cascade.rs` and
`talking_hedz_parity.rs`). Any drift from that path will break heritage
character, so this script imports rather than inlines.

Usage:

    python tools/bake_md_live.py <sections_in.json> <compiled_v1_out.json> [--name NAME]

Input JSON schema (flat, no wrapping template envelope):

    {
      "name": "my-body",
      "sections": [
        {"index": 1, "type": 3, "low_freq": 0,  "low_gain": 127,
                                "high_freq": 116, "high_gain": 87},
        ... up to 6 sections
      ]
    }
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from tools.bake_hedz_const import compile_corner, NUM_STAGES, SR  # noqa: E402

PASSTHROUGH_STAGE = (1.0, 0.0, 0.0, 0.0, 0.0)
RUNTIME_STAGES = 12      # engine topology: 6 active + 6 passthrough
RUNTIME_BOOST = 4.0
RUNTIME_SR = 39062.5     # native cascade SR — see project_native_sr_resampling memory


def pad_to_runtime(stages: list[tuple[float, ...]]) -> list[tuple[float, ...]]:
    """Pad the 6 MD stages with passthrough to fill the 12-stage cascade."""
    padded = list(stages)
    while len(padded) < RUNTIME_STAGES:
        padded.append(PASSTHROUGH_STAGE)
    return padded[:RUNTIME_STAGES]


def stages_to_keyframe(label: str, stages: list[tuple[float, ...]]) -> dict:
    return {
        "label": label,
        "boost": RUNTIME_BOOST,
        "stages": [
            {"c0": s[0], "c1": s[1], "c2": s[2], "c3": s[3], "c4": s[4]}
            for s in pad_to_runtime(stages)
        ],
    }


def build_compiled_v1(name: str, sections: list[dict]) -> dict:
    """Compile 2-frame MD sections into 4-corner compiled-v1 JSON.

    The Q axis is degenerate: keyframes M0_Q0 == M0_Q100 and
    M100_Q0 == M100_Q100. This matches the `[m0, m1, m0, m1]` shape
    produced by `build_cartridge` in `bake_hedz_const.py`.
    """
    m0 = compile_corner(sections, morph=0.0, sr=SR)
    m1 = compile_corner(sections, morph=1.0, sr=SR)
    return {
        "format": "compiled-v1",
        "name": name,
        "provenance": "md-live",
        "sampleRate": RUNTIME_SR,
        "stages": RUNTIME_STAGES,
        "keyframes": [
            stages_to_keyframe("M0_Q0", m0),
            stages_to_keyframe("M100_Q0", m1),
            stages_to_keyframe("M0_Q100", m0),
            stages_to_keyframe("M100_Q100", m1),
        ],
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Bake MD sections → compiled-v1 JSON")
    ap.add_argument("input", type=Path, help="sections input JSON")
    ap.add_argument("output", type=Path, help="compiled-v1 output JSON")
    ap.add_argument("--name", default=None, help="override body name")
    args = ap.parse_args(argv)

    payload = json.loads(args.input.read_text(encoding="utf-8"))
    sections = payload.get("sections", [])
    if not sections:
        print("no sections in input", file=sys.stderr)
        return 2
    name = args.name or payload.get("name", "md-live")

    cartridge = build_compiled_v1(name, sections)
    args.output.write_text(
        json.dumps(cartridge, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
