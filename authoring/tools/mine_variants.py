"""forge/mine_variants.py — Path 2: Generative Variant Mining.

Bake N variants of an arc, randomising within the doctrinal Q-law envelope.
Every variant gate-passes by construction. Saves a contact-sheet plot grid
plus the individual .json + .magnitude.png files so you can scroll and
audition. Comparable to what the heritage E-MU sound designers were doing
when they threw a thousand variants and kept the winners.

What's randomised:
  * uniform_radius      ∈ [0.940, 0.999]   — pole strength
  * m0_zero_radius      ∈ [0.70, 1.00]      — M0 carve depth
  * m100_zero_radius    ∈ [0.30, 0.85]      — M100 dissolve point
  * per-stage zero offsets (treble-basket jitter ± 800 Hz)
  * eq_zero_offset_semitones (cascade only) ∈ [3.0, 9.0]
  * anchor pole jitter (vocal arcs) ± 30%

What's locked:
  * stage_gain ≡ 1.0 (Rule 5)
  * shared c4 across active stages (cascade) / lattice law (lattice)
  * topology (per recipe — vocal = cascade_6, sweep = lattice_3)
  * gate must pass

Usage:
    python forge/mine_variants.py ah_to_ee --count 16 --seed 42
"""
from __future__ import annotations

import argparse
import copy
import json
import math
import random
import struct
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ── reuse the forge core ──────────────────────────────────────────────────
sys.path.insert(0, str(Path(__file__).resolve().parent))
import forge_cli  # noqa: E402


REPO = Path(__file__).resolve().parent.parent
OUT_DIR = REPO / "cartridges" / "factory" / "generated" / "qlaw" / "variants"


@dataclass
class VariantParams:
    radius: float
    m0_zero_radius: float
    m100_zero_radius: float
    eq_zero_offset_semitones: float
    treble_basket_jitter_hz: float


def random_params(rng: random.Random) -> VariantParams:
    return VariantParams(
        radius=rng.uniform(0.940, 0.999),
        m0_zero_radius=rng.uniform(0.70, 1.00),
        m100_zero_radius=rng.uniform(0.30, 0.85),
        eq_zero_offset_semitones=rng.uniform(3.0, 9.0),
        treble_basket_jitter_hz=rng.uniform(-800.0, 800.0),
    )


def apply_params(p: VariantParams) -> dict[str, Any]:
    """Patch forge_cli module-level constants with this variant's settings.
    Snapshots originals so we can restore between variants."""
    snap = {
        "UNIFORM_RADIUS": forge_cli.UNIFORM_RADIUS,
        "M0_ZERO_RADIUS": forge_cli.M0_ZERO_RADIUS,
        "M100_ZERO_RADIUS": forge_cli.M100_ZERO_RADIUS,
        "EQ_ZERO_OFFSET_SEMITONES": forge_cli.EQ_ZERO_OFFSET_SEMITONES,
        "ROLE_ZERO_PLACEMENT": copy.deepcopy(forge_cli.ROLE_ZERO_PLACEMENT),
    }
    forge_cli.UNIFORM_RADIUS = p.radius
    forge_cli.M0_ZERO_RADIUS = p.m0_zero_radius
    forge_cli.M100_ZERO_RADIUS = p.m100_zero_radius
    forge_cli.EQ_ZERO_OFFSET_SEMITONES = p.eq_zero_offset_semitones
    # jitter the treble-basket zero positions
    for role, spec in forge_cli.ROLE_ZERO_PLACEMENT.items():
        if spec.get("mode") == "fixed_hz" and spec.get("value") is not None:
            base = float(spec["value"])
            jittered = max(20.0, min(forge_cli.NYQUIST_LIMIT_HZ * 0.985,
                                     base + p.treble_basket_jitter_hz))
            spec["value"] = jittered
    return snap


def restore_params(snap: dict[str, Any]) -> None:
    forge_cli.UNIFORM_RADIUS = snap["UNIFORM_RADIUS"]
    forge_cli.M0_ZERO_RADIUS = snap["M0_ZERO_RADIUS"]
    forge_cli.M100_ZERO_RADIUS = snap["M100_ZERO_RADIUS"]
    forge_cli.EQ_ZERO_OFFSET_SEMITONES = snap["EQ_ZERO_OFFSET_SEMITONES"]
    forge_cli.ROLE_ZERO_PLACEMENT = snap["ROLE_ZERO_PLACEMENT"]


def composite_db(cart: dict, label: str, freqs: np.ndarray) -> np.ndarray:
    kf = next(k for k in cart["keyframes"] if k["label"] == label)
    sr = float(cart.get("sampleRate", 39062.5))
    boost = float(kf.get("boost", 1.0))
    w = 2.0 * np.pi * freqs / sr
    e1 = np.exp(-1j * w)
    e2 = np.exp(-2j * w)
    mag = np.ones_like(freqs)
    for s in kf["stages"]:
        c0 = float(s["c0"]); c1 = float(s["c1"]); c2 = float(s["c2"])
        c3 = float(s["c3"]); c4 = float(s["c4"])
        if c3 == 0.0 and c4 == 0.0:
            continue
        num = c0 + c1 * e1 + c2 * e2
        den = 1.0 + c3 * e1 + c4 * e2
        mag *= np.abs(num / den)
    return 20.0 * np.log10(np.maximum(mag * boost, 1e-9))


def trajectory_db(cart: dict, freqs: np.ndarray, n_steps: int = 9) -> list[np.ndarray]:
    """Sample the morph trajectory at N evenly-spaced positions between M0_Q0
    and M100_Q0. Mirrors the runtime's linear coefficient lerp."""
    m0_kf = next(k for k in cart["keyframes"] if k["label"] == "M0_Q0")
    m1_kf = next(k for k in cart["keyframes"] if k["label"] == "M100_Q0")
    sr = float(cart.get("sampleRate", 39062.5))
    boost_a = float(m0_kf.get("boost", 1.0))
    boost_b = float(m1_kf.get("boost", 1.0))
    w = 2.0 * np.pi * freqs / sr
    e1 = np.exp(-1j * w)
    e2 = np.exp(-2j * w)
    out = []
    for step in range(n_steps):
        t = step / (n_steps - 1)
        boost = (1.0 - t) * boost_a + t * boost_b
        mag = np.ones_like(freqs)
        for sa, sb in zip(m0_kf["stages"], m1_kf["stages"]):
            c0 = (1.0 - t) * float(sa["c0"]) + t * float(sb["c0"])
            c1 = (1.0 - t) * float(sa["c1"]) + t * float(sb["c1"])
            c2 = (1.0 - t) * float(sa["c2"]) + t * float(sb["c2"])
            c3 = (1.0 - t) * float(sa["c3"]) + t * float(sb["c3"])
            c4 = (1.0 - t) * float(sa["c4"]) + t * float(sb["c4"])
            if c3 == 0.0 and c4 == 0.0:
                continue
            num = c0 + c1 * e1 + c2 * e2
            den = 1.0 + c3 * e1 + c4 * e2
            mag *= np.abs(num / den)
        out.append(20.0 * np.log10(np.maximum(mag * boost, 1e-9)))
    return out


def trajectory_score(traces: list[np.ndarray]) -> float:
    """Higher = more progressive / coherent morph motion.

    Combines:
      * spectral motion magnitude (M0 must differ from M100)
      * trajectory smoothness (consecutive frames close to each other)
    """
    if len(traces) < 3:
        return 0.0
    arr = np.stack(traces)              # shape (n_steps, n_freqs)
    motion = np.linalg.norm(arr[-1] - arr[0]) / arr.shape[1]
    # smoothness: low avg jump between adjacent frames (relative to total motion)
    jumps = np.linalg.norm(arr[1:] - arr[:-1], axis=1) / arr.shape[1]
    smoothness = 1.0 / (1.0 + jumps.std())
    return float(motion * smoothness)


def render_contact_sheet(
    arc: forge_cli.Arc,
    variants: list[tuple[int, VariantParams, dict, float]],
    out_path: Path,
) -> None:
    """Each cell shows the FULL morph trajectory — N spectra at N evenly-
    spaced morph positions, color-graded cool→warm. A coherent variant
    sweeps smoothly; a broken one shows messy crossings or stays static."""
    n = len(variants)
    cols = 4
    rows = (n + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(16.0, 2.7 * rows), dpi=140, sharex=True)
    fig.patch.set_facecolor("#16191a")
    fig.suptitle(
        f"FORGE :: variant trajectories — {arc.name}  ({arc.recipe})\n"
        f"each line = one morph position (cool=M0, warm=M100). "
        f"smooth sweep = whole. messy = broken.",
        color="#ffba00", fontweight="bold", fontsize=12,
    )
    freqs = np.geomspace(40.0, 18000.0, 768)
    flat = axes.flatten() if hasattr(axes, "flatten") else [axes]

    n_steps = 9
    cmap = plt.get_cmap("plasma")  # cool→warm
    line_colors = [cmap(i / (n_steps - 1)) for i in range(n_steps)]

    for ax, (idx, p, cart, score) in zip(flat, variants):
        ax.set_facecolor("#0d0f10")
        traces = trajectory_db(cart, freqs, n_steps=n_steps)
        for col, trace in zip(line_colors, traces):
            ax.semilogx(freqs, trace, color=col, lw=0.9, alpha=0.85)
        # bold the endpoints so M0 and M100 still pop
        ax.semilogx(freqs, traces[0],  color="#22ddff", lw=1.6, alpha=0.95)
        ax.semilogx(freqs, traces[-1], color="#ffba00", lw=1.6, alpha=0.95)
        ax.set_title(
            f"#{idx:02d}  score={score:.2f}  r={p.radius:.3f}",
            color="white", fontsize=9, fontweight="bold",
        )
        ax.grid(True, which="both", alpha=0.18, color="white")
        ax.tick_params(colors="white", labelsize=7)
        for spine in ax.spines.values():
            spine.set_color("#3a3e42")

    for j in range(n, len(flat)):
        flat[j].axis("off")

    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, facecolor=fig.get_facecolor())
    plt.close(fig)


def main(argv: list[str]) -> int:
    p = argparse.ArgumentParser()
    p.add_argument("arc", help="arc name (run `forge_cli list` to see them)")
    p.add_argument("--count", type=int, default=16,
                   help="number of variants to mine (default 16)")
    p.add_argument("--seed", type=int, default=None,
                   help="reproducibility seed (default: random)")
    p.add_argument("--policy", default="A", choices=["A", "B"])
    args = p.parse_args(argv)

    rng = random.Random(args.seed)
    arc = forge_cli.find_arc(args.arc)
    arc_dir = OUT_DIR / arc.name
    arc_dir.mkdir(parents=True, exist_ok=True)

    variants: list[tuple[int, VariantParams, dict, float]] = []
    rejects = 0
    attempts = 0
    score_freqs = np.geomspace(40.0, 18000.0, 256)

    while len(variants) < args.count and attempts < args.count * 5:
        attempts += 1
        params = random_params(rng)
        snap = apply_params(params)
        try:
            frame_map = forge_cli.arc_to_frame_map(arc)
            stages, meta = forge_cli.plan(frame_map, arc.recipe)
            solved = forge_cli.solve(stages, args.policy)
            cart = forge_cli.build_cartridge(meta, solved, args.policy)
            gate_result = forge_cli.gate(cart)
        finally:
            restore_params(snap)

        if not gate_result.passed:
            rejects += 1
            continue

        traces = trajectory_db(cart, score_freqs, n_steps=9)
        score = trajectory_score(traces)
        variants.append((len(variants), params, cart, score))

    # Sort by score (best wholeness first), then re-write JSON files in score order
    variants.sort(key=lambda v: v[3], reverse=True)
    sorted_variants: list[tuple[int, VariantParams, dict, float]] = []
    for new_idx, (_old_idx, p, cart, score) in enumerate(variants):
        cart_path = arc_dir / f"{arc.name}_v{new_idx:02d}_score{score:.2f}.json"
        cart_path.write_text(json.dumps(cart, indent=2) + "\n", encoding="utf-8")
        sorted_variants.append((new_idx, p, cart, score))
    variants = sorted_variants

    if not variants:
        print(f"no variants survived the gate after {attempts} attempts", file=sys.stderr)
        return 1

    sheet_path = arc_dir / f"_contact_sheet.png"
    render_contact_sheet(arc, variants, sheet_path)

    print(f"mined {len(variants)} variants ({rejects} rejected) for {arc.name}")
    print(f"contact sheet → {sheet_path}")
    print(f"individual cartridges → {arc_dir}/")

    if sys.platform == "win32":
        try:
            import os
            os.startfile(str(sheet_path))  # type: ignore[attr-defined]
        except Exception:
            pass
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
