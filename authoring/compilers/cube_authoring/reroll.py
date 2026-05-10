"""Reroll generator — sample a single cube corner (or a full cube) from a scope.

One corner is a 6-section Morph Designer template. The reroll picks fresh
section freq/gain values inside the scope's per-corner envelope, respecting
the layout rule (cluster vs distributed) implied by the corner's 3-bit
address.

Inputs
  - scope name (hardcoded envelopes below; move to scopes.json later)
  - corner label (c000..c111) for --corner mode, or --cube for all 8
  - optional RNG seed for reproducibility

Output
  - stdout: one corner template (JSON) or a full cube authoring doc (JSON)

This generator honors the narrow envelope of the declared scope. It does not
invent new scopes. Adding a scope means adding an entry in SCOPE_ENVELOPES.
"""
from __future__ import annotations

import argparse
import json
import random
import sys
from typing import Iterable

CORNER_LABELS = ("c000", "c100", "c010", "c110", "c001", "c101", "c011", "c111")

# Pill keyframes: 4 edges of the morph×Q plane. Maps name → (corner_label, morph, q)
PILL_KEYFRAMES = {
    "M0_Q0":    ("c000", 0.0, 0.0),
    "M100_Q0":  ("c100", 0.0, 0.0),
    "M0_Q100":  ("c010", 0.0, 0.0),
    "M100_Q100": ("c110", 0.0, 0.0),
}

# -----------------------------------------------------------------------------
# Scope envelopes
#
# freq_center / gain_center: Morph Designer packed-byte targets (0..127).
# freq_spread / gain_spread: half-width of the sampling window around center.
# layout: "cluster" = all 6 sections near one freq; "distributed" = 6 sections
#         evenly spread across [center - spread .. center + spread].
# -----------------------------------------------------------------------------

SCOPE_ENVELOPES = {
    "hf_resonance": {
        "section_type": 1,          # E-mu firmware type=1 (symmetric Peak)
        "freq_min": 75,             # HF floor — below this is not hf_resonance
        "freq_max": 127,            # packed maximum
        "gain_min": 64,
        "gain_max": 127,
        "in_scope_probes": ["sine_sweep", "reese", "drum_loop", "noise"],
        "corners": {
            "c000": {"freq_center":  95, "freq_spread":  4, "gain_center":  80, "gain_spread": 6,  "layout": "cluster"},
            "c100": {"freq_center": 120, "freq_spread":  4, "gain_center":  80, "gain_spread": 6,  "layout": "cluster"},
            "c010": {"freq_center":  95, "freq_spread":  4, "gain_center": 120, "gain_spread": 5,  "layout": "cluster"},
            "c110": {"freq_center": 120, "freq_spread":  4, "gain_center": 120, "gain_spread": 5,  "layout": "cluster"},
            "c001": {"freq_center":  95, "freq_spread": 12, "gain_center":  80, "gain_spread": 6,  "layout": "distributed"},
            "c101": {"freq_center": 118, "freq_spread":  9, "gain_center":  80, "gain_spread": 6,  "layout": "distributed"},
            "c011": {"freq_center":  95, "freq_spread": 12, "gain_center": 120, "gain_spread": 5,  "layout": "distributed"},
            "c111": {"freq_center": 118, "freq_spread":  9, "gain_center": 120, "gain_spread": 5,  "layout": "distributed"},
        },
    },
}

# Small jitter applied to every sampled value regardless of layout, to avoid
# two reroll calls producing identical sections when the envelope is tight.
JITTER = 2


# -----------------------------------------------------------------------------
# Sampling
# -----------------------------------------------------------------------------

def clamp(v: int, lo: int, hi: int) -> int:
    return max(lo, min(hi, v))


def sample_cluster_freqs(rng: random.Random, env: dict, scope: dict) -> list[int]:
    center = env["freq_center"]
    spread = env["freq_spread"]
    lo = clamp(center - spread, scope["freq_min"], scope["freq_max"])
    hi = clamp(center + spread, scope["freq_min"], scope["freq_max"])
    return [rng.randint(lo, hi) for _ in range(6)]


def sample_distributed_freqs(rng: random.Random, env: dict, scope: dict) -> list[int]:
    center = env["freq_center"]
    spread = env["freq_spread"]
    lo = clamp(center - spread, scope["freq_min"], scope["freq_max"])
    hi = clamp(center + spread, scope["freq_min"], scope["freq_max"])
    if hi <= lo:
        return [lo] * 6
    step = (hi - lo) / 5.0
    freqs = [round(lo + i * step) for i in range(6)]
    return [clamp(f + rng.randint(-JITTER, JITTER),
                   scope["freq_min"], scope["freq_max"]) for f in freqs]


def sample_gains(rng: random.Random, env: dict, scope: dict) -> list[int]:
    center = env["gain_center"]
    spread = env["gain_spread"]
    lo = clamp(center - spread, scope["gain_min"], scope["gain_max"])
    hi = clamp(center + spread, scope["gain_min"], scope["gain_max"])
    return [rng.randint(lo, hi) for _ in range(6)]


def sample_corner(
    label: str, scope_name: str, rng: random.Random, scope_env: dict | None = None
) -> dict:
    if scope_env is None:
        if scope_name not in SCOPE_ENVELOPES:
            raise ValueError(f"unknown scope {scope_name!r}")
        scope_env = SCOPE_ENVELOPES[scope_name]
    scope = scope_env
    if label not in scope["corners"]:
        raise ValueError(f"scope {scope_name!r} has no envelope for corner {label!r}")
    env = scope["corners"][label]

    if env["layout"] == "cluster":
        freqs = sample_cluster_freqs(rng, env, scope)
    elif env["layout"] == "distributed":
        freqs = sample_distributed_freqs(rng, env, scope)
    else:
        raise ValueError(f"unknown layout {env['layout']!r}")

    gains = sample_gains(rng, env, scope)
    section_type = scope["section_type"]

    sections = []
    for i in range(6):
        sections.append({
            "index": i + 1,
            "type": section_type,
            "low_freq": int(freqs[i]),
            "low_gain": int(gains[i]),
            "high_freq": int(freqs[i]),
            "high_gain": int(gains[i]),
        })
    return {
        "name": f"{scope_name}.{label}.reroll",
        "sections": sections,
    }


def corner_authoring_block(label: str, scope_name: str, rng: random.Random) -> dict:
    template = sample_corner(label, scope_name, rng)
    return {
        "kind": "morph_designer_inline",
        "sample": {"morph": 0.0, "q": 0.0},
        "template": template,
    }


def sample_pill(
    scope_name: str,
    rng: random.Random,
    *,
    pill_id: str,
    name: str,
    scope_env: dict | None = None,
) -> dict:
    scope = scope_env if scope_env is not None else SCOPE_ENVELOPES[scope_name]
    keyframes = {}
    for kf_label, (corner_label, morph, q) in PILL_KEYFRAMES.items():
        template = sample_corner(corner_label, scope_name, rng, scope_env=scope_env)
        template["name"] = f"{scope_name}.{kf_label}.reroll"
        keyframes[kf_label] = {
            "kind": "morph_designer_inline",
            "sample": {"morph": morph, "q": q},
            "template": template,
        }
    return {
        "schema": "trench.authoring_path.pill.v1",
        "id": pill_id,
        "name": name,
        "scope": scope_name,
        "exactness": "modern_cleanroom_not_native_verified",
        "sampleRate": 44100,
        "provenance": {
            "generator": "tools/cube_authoring/reroll.py",
            "derived_from": "trench_re_vault/analysis/cube_laws/cube_laws.json",
            "rulebook_version": "cube_laws_v0",
            "scope_envelope_version": "v0",
            "in_scope_probes": scope.get("in_scope_probes", []),
        },
        "keyframes": keyframes,
    }


def sample_cube(scope_name: str, rng: random.Random, *, cube_id: str, name: str) -> dict:
    scope = SCOPE_ENVELOPES[scope_name]
    corners = {label: corner_authoring_block(label, scope_name, rng)
                for label in CORNER_LABELS}
    return {
        "schema": "trench.authoring_path.cube.v1",
        "id": cube_id,
        "name": name,
        "scope": scope_name,
        "exactness": "modern_cleanroom_not_native_verified",
        "control_mode_default": "modern_live_xyz",
        "provenance": {
            "generator": "tools/cube_authoring/reroll.py",
            "derived_from": "trench_re_vault/analysis/cube_laws/cube_laws.json",
            "rulebook_version": "cube_laws_v0",
            "scope_envelope_version": "v0",
        },
        "corners": corners,
    }


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------

def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(
        description="Reroll a single cube corner or a full cube from a scope envelope.",
    )
    parser.add_argument("--scope", default="hf_resonance",
                         help="Scope name (default: hf_resonance).")
    parser.add_argument("--seed", type=int, default=None,
                         help="RNG seed for reproducibility (default: random).")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--corner", choices=CORNER_LABELS,
                        help="Emit a single corner template.")
    group.add_argument("--cube", action="store_true",
                        help="Emit a full 8-corner cube authoring document.")
    group.add_argument("--pill", action="store_true",
                        help="Emit a 4-keyframe pill authoring document (trench.authoring_path.pill.v1).")
    parser.add_argument("--cube-id", default="cube.rolled",
                         help="id field for --cube output (default: cube.rolled).")
    parser.add_argument("--cube-name", default=None,
                         help="name field for --cube output (default: derived from id).")
    parser.add_argument("--pill-id", default="pill.rolled",
                         help="id field for --pill output (default: pill.rolled).")
    parser.add_argument("--pill-name", default=None,
                         help="name field for --pill output (default: derived from id).")
    parser.add_argument("--out", default=None,
                         help="Output path (stdout if omitted).")
    args = parser.parse_args(argv)

    rng = random.Random(args.seed)

    if args.corner:
        out = corner_authoring_block(args.corner, args.scope, rng)
    elif args.pill:
        pill_name = args.pill_name or args.pill_id
        out = sample_pill(args.scope, rng, pill_id=args.pill_id, name=pill_name)
    else:
        cube_name = args.cube_name or args.cube_id
        out = sample_cube(args.scope, rng, cube_id=args.cube_id, name=cube_name)

    text = json.dumps(out, indent=2)
    if args.out:
        from pathlib import Path
        p = Path(args.out)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(text + "\n", encoding="utf-8")
    else:
        sys.stdout.write(text + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
