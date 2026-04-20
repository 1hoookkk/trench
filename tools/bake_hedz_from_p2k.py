"""Bake `hedz_rom.rs` + `hedz_golden.rs` directly from a compiled-v1 P2K
cartridge JSON.

Registry identities (see `FILTER_ARTIFACTS.json`):
    source  = compiled_p2k_013   (default)
    output  = rom_hedz           (trench-core/src/hedz_rom.rs)

Replaces the heritage-template compile path in `bake_hedz_const.py` (which
emitted `md_hedz`-sourced coefficients — stale, morph-only) with a
straight ingest of the compiled-v1 runtime truth. The committed Rust rom
and the JSON the plugin loads for P2K index 013 are now the same bytes.

Doctrinally still COMPILER-side: reads a repo JSON artifact, emits Rust
source. No pyruntime imports, stdlib only.

Usage::

    python tools/bake_hedz_from_p2k.py
    # override source:
    python tools/bake_hedz_from_p2k.py --source trench-juce/cartridges/p2k/P2k_013.json
"""
from __future__ import annotations

import argparse
import json
import struct
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
DEFAULT_SOURCE = REPO / "trench-juce" / "cartridges" / "p2k" / "P2k_013.json"
DEFAULT_ROM = REPO / "trench-core" / "src" / "hedz_rom.rs"
DEFAULT_GOLDEN = REPO / "trench-core" / "src" / "hedz_golden.rs"

# Runtime cascade shape (must match trench-core/src/cartridge.rs).
NUM_CORNERS = 4
NUM_STAGES = 6
NUM_COEFFS = 5
PASSTHROUGH = [1.0, 0.0, 0.0, 0.0, 0.0]
IMPULSE_LEN = 256

# JSON → rom corner ordering.
JSON_LABELS = ["M0_Q0", "M0_Q100", "M100_Q0", "M100_Q100"]
ROM_LABELS = ["M0_Q0", "M100_Q0", "M0_Q100", "M100_Q100"]


def f32(x: float) -> float:
    """Round-trip through IEEE-754 single precision (matches np.float32)."""
    return struct.unpack("<f", struct.pack("<f", x))[0]


# ---------------------------------------------------------------------------
# Compiled-v1 JSON ingest
# ---------------------------------------------------------------------------

def load_compiled_v1(path: Path) -> tuple[list[list[tuple[float, ...]]], float, str]:
    """Return (corners_in_rom_order, boost, source_name).

    `corners_in_rom_order` is a list of 4 corners, each a list of 6 stages,
    each a 5-tuple of kernel coefficients. The Q-axis collapse claim of
    the heritage template is validated here: if a corner past stage 6 is
    not a passthrough, this raises.
    """
    data = json.loads(path.read_text(encoding="utf-8"))
    fmt = data.get("format")
    if fmt != "compiled-v1":
        raise SystemExit(f"{path}: expected format 'compiled-v1', got {fmt!r}")

    keyframes = {kf["label"]: kf for kf in data["keyframes"]}
    missing = [l for l in JSON_LABELS if l not in keyframes]
    if missing:
        raise SystemExit(f"{path}: missing corner labels {missing}")

    boosts = {label: keyframes[label]["boost"] for label in JSON_LABELS}
    boost_vals = set(boosts.values())
    if len(boost_vals) != 1:
        raise SystemExit(
            f"{path}: non-uniform boost across corners ({boosts}); "
            "rom layout assumes one scalar boost."
        )
    boost = next(iter(boost_vals))

    # Extract the first 6 stages per corner and verify 7..12 are passthrough.
    per_label: dict[str, list[tuple[float, ...]]] = {}
    for label in JSON_LABELS:
        stages_json = keyframes[label]["stages"]
        if len(stages_json) != NUM_STAGES * 2:
            raise SystemExit(
                f"{path} {label}: expected {NUM_STAGES * 2} stages, "
                f"got {len(stages_json)}"
            )
        head = []
        for s in stages_json[:NUM_STAGES]:
            head.append(tuple(f32(s[k]) for k in ("c0", "c1", "c2", "c3", "c4")))
        for i, s in enumerate(stages_json[NUM_STAGES:], start=NUM_STAGES + 1):
            vals = [s["c0"], s["c1"], s["c2"], s["c3"], s["c4"]]
            if vals != PASSTHROUGH:
                raise SystemExit(
                    f"{path} {label} stage {i}: expected passthrough "
                    f"{PASSTHROUGH}, got {vals}. The rom only stores 6 stages."
                )
        per_label[label] = head

    corners_rom = [per_label[label] for label in ROM_LABELS]
    return corners_rom, float(boost), data.get("name", path.stem)


# ---------------------------------------------------------------------------
# DF2T cascade emulator (mirrors trench-core/src/cascade.rs math exactly —
# static coefficients, sample-accurate, no ramping).
# ---------------------------------------------------------------------------

def run_impulse(
    corner_stages: list[tuple[float, ...]],
    length: int,
    boost: float,
) -> list[float]:
    state = [(0.0, 0.0) for _ in range(NUM_STAGES)]
    out = []
    for n in range(length):
        x = 1.0 if n == 0 else 0.0
        sig = x
        for i, (c0, c1, c2, c3, c4) in enumerate(corner_stages):
            w1, w2 = state[i]
            y = c0 * sig + w1
            new_w1 = c1 * sig - c3 * y + w2
            new_w2 = c2 * sig - c4 * y
            state[i] = (new_w1, new_w2)
            sig = y
        sig *= boost
        out.append(f32(sig))
    return out


# ---------------------------------------------------------------------------
# Rust codegen
# ---------------------------------------------------------------------------

_ROM_HEADER = '''//! Talking Hedz cartridge — baked from `trench-juce/cartridges/p2k/P2k_013.json`.
//!
//! Registry ID: `rom_hedz` (derived_from `compiled_p2k_013`).
//! See `FILTER_ARTIFACTS.json` — do not compare this rom against any
//! canonical render outside its listed `null_test_pairs`.
//!
//! DO NOT EDIT BY HAND. Regenerate with `python tools/bake_hedz_from_p2k.py`.
//!
//! Source: the compiled-v1 P2K cartridge the plugin runtime loads for the
//! Talking Hedz skin (P2K index 013). This rom is a byte-compatible
//! restatement of the JSON so the plugin default and the Rust cascade
//! tests consume identical coefficients — no heritage XML compile step
//! in between.
//!
//! Corner order matches `Cartridge::corners`:
//!     0: M0_Q0, 1: M100_Q0, 2: M0_Q100, 3: M100_Q100
//!
//! The source JSON ships 12 stages per corner; stages 7–12 are all
//! passthrough (`[1,0,0,0,0]`). The bake script asserts that before
//! truncating to the active 6 so this rom can never silently drop a
//! non-trivial stage.

use crate::cartridge::{{CornerData, NUM_CORNERS}};

/// Authoring name — one heap allocation at plugin init, never on the
/// audio thread.
pub const HEDZ_NAME: &str = "Talking Hedz";

/// Per-corner post-cascade boost. The source JSON carries a single
/// scalar boost at every corner; the bake script refuses to emit if
/// that invariant breaks.
pub const HEDZ_BOOSTS: [f64; NUM_CORNERS] = [{boost_lit}, {boost_lit}, {boost_lit}, {boost_lit}];

/// Full cartridge — 4 corners × 6 stages × 5 coefficients.
pub const HEDZ_CORNERS: [CornerData; NUM_CORNERS] = [
'''

_ROM_FOOTER = "];\n"


def emit_rom(
    corners: list[list[tuple[float, ...]]],
    boost: float,
    out_path: Path,
) -> None:
    boost_lit = f"{boost!r}"
    lines = [_ROM_HEADER.format(boost_lit=boost_lit)]
    for corner_idx, corner in enumerate(corners):
        lines.append(f"    // {ROM_LABELS[corner_idx]}\n")
        lines.append("    [\n")
        for stage_idx, stage in enumerate(corner):
            c0, c1, c2, c3, c4 = stage
            lines.append(
                f"        [{c0!r}, {c1!r}, {c2!r}, {c3!r}, {c4!r}], // stage {stage_idx + 1}\n"
            )
        lines.append("    ],\n")
    lines.append(_ROM_FOOTER)
    out_path.write_text("".join(lines), encoding="utf-8")


_GOLDEN_HEADER = '''//! Talking Hedz cascade golden impulse response.
//!
//! DO NOT EDIT BY HAND. Regenerate with `python tools/bake_hedz_from_p2k.py`.
//!
//! Each `HEDZ_GOLDEN_<LABEL>` array is the first {length} samples of a
//! unit-impulse response through the static hedz cascade at the named
//! corner, computed by a stdlib DF2T emulator in Python. The Rust
//! cross-language parity test (`trench-core/tests/hedz_cascade.rs`)
//! feeds the same impulse through `trench_core::Cascade` with
//! `set_targets(chunk_size=1)` so coefficients snap immediately and no
//! ramping interpolation runs, and asserts the two outputs match
//! byte-close (< 1e-6 peak error). If this test fails, either the
//! bake script or the Rust cascade math has drifted — fix the root
//! cause, do not relax the tolerance.

pub const HEDZ_GOLDEN_LEN: usize = {length};

'''


def emit_golden(
    impulses: dict[str, list[float]],
    out_path: Path,
    length: int,
) -> None:
    lines = [_GOLDEN_HEADER.format(length=length)]
    for label, samples in impulses.items():
        lines.append(f"pub const HEDZ_GOLDEN_{label}: [f32; HEDZ_GOLDEN_LEN] = [\n")
        for i in range(0, length, 4):
            chunk = samples[i : i + 4]
            formatted = ", ".join(f"{s!r}f32" for s in chunk)
            lines.append(f"    {formatted},\n")
        lines.append("];\n\n")
    out_path.write_text("".join(lines), encoding="utf-8")


# ---------------------------------------------------------------------------
# Entrypoint
# ---------------------------------------------------------------------------

def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--rom-out", type=Path, default=DEFAULT_ROM)
    parser.add_argument("--golden-out", type=Path, default=DEFAULT_GOLDEN)
    parser.add_argument("--impulse-len", type=int, default=IMPULSE_LEN)
    args = parser.parse_args(argv)

    corners, boost, source_name = load_compiled_v1(args.source)

    impulses: dict[str, list[float]] = {}
    for label, stages in zip(ROM_LABELS, corners):
        impulses[label] = run_impulse(stages, args.impulse_len, boost)

    args.rom_out.parent.mkdir(parents=True, exist_ok=True)
    emit_rom(corners, boost, args.rom_out)
    emit_golden(impulses, args.golden_out, args.impulse_len)

    print(f"source: {args.source}  ({source_name})")
    print(f"boost:  {boost}")
    print(f"wrote {args.rom_out}")
    print(f"wrote {args.golden_out}")
    for label, samples in impulses.items():
        peak = max(abs(s) for s in samples)
        energy = sum(s * s for s in samples)
        print(f"  {label}: peak={peak:.6f}  energy={energy:.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
