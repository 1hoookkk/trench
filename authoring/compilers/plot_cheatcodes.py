"""Magnitude-response plots for cheatcode cartridges.

One PNG per cartridge: 2×2 grid of corners (M0_Q0, M100_Q0, M0_Q100, M100_Q100).
Each panel shows the cascade magnitude response plus each active stage's
individual contribution, so you can see which stage does what.

Matches the shipping runtime's math exactly — plain-DF2T kernel form where
the numerator is [c0, c1, c2] and the denominator is [1, c3, c4]. Evaluates
at host SR (44100 Hz) per memory:compiled_v1_sample_rate.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

SR = 44100.0
NFFT = 1024
F_MIN = 20.0
F_MAX = SR / 2.0
DB_MIN, DB_MAX = -40.0, 40.0

REPO = Path(__file__).resolve().parents[1]
IN_DIR = REPO / "looperator" / "cheatcode_presets"
OUT_DIR = IN_DIR / "plots"

CORNER_LAYOUT = [
    ("M0_Q0",     0, 0),
    ("M100_Q0",   0, 1),
    ("M0_Q100",   1, 0),
    ("M100_Q100", 1, 1),
]


def stage_response(stage: dict, freqs_hz: np.ndarray) -> np.ndarray:
    """|H(e^jω)| for one DF2T biquad stage on the given Hz grid."""
    c0, c1, c2, c3, c4 = stage["c0"], stage["c1"], stage["c2"], stage["c3"], stage["c4"]
    w = 2.0 * np.pi * freqs_hz / SR
    e1 = np.exp(-1j * w)
    e2 = np.exp(-2j * w)
    num = c0 + c1 * e1 + c2 * e2
    den = 1.0 + c3 * e1 + c4 * e2
    return np.abs(num / den)


def is_passthrough(stage: dict) -> bool:
    return (
        abs(stage["c0"] - 1.0) < 1e-9
        and abs(stage["c1"]) < 1e-9
        and abs(stage["c2"]) < 1e-9
        and abs(stage["c3"]) < 1e-9
        and abs(stage["c4"]) < 1e-9
    )


def describe_stage(stage: dict) -> str:
    """Decode pole/zero Hz for a legend label."""
    if is_passthrough(stage):
        return "passthrough"
    c0, c1, c2, c3, c4 = stage["c0"], stage["c1"], stage["c2"], stage["c3"], stage["c4"]
    # Denominator [1, c3, c4] → pole pair at r=√c4, θ=acos(-c3/(2r))
    pole_r = math.sqrt(max(0.0, c4))
    if pole_r > 1e-6:
        arg = max(-1.0, min(1.0, -c3 / (2.0 * pole_r)))
        pole_hz = math.acos(arg) * SR / (2.0 * math.pi)
    else:
        pole_hz = 0.0
    # Numerator [c0, c1, c2] — only meaningful if non-degenerate
    has_zero = abs(c1) > 1e-9 or abs(c2) > 1e-9
    if has_zero and abs(c0) > 1e-9:
        # Roots of c2 + c1·z + c0·z² = 0 (z-domain, after multiplying by z²)
        # → z = (-c1 ± √(c1² - 4·c0·c2)) / (2·c0)
        disc = c1 * c1 - 4.0 * c0 * c2
        if disc < 0:
            # Complex conjugate pair: z = (-c1/(2c0)) ± j·√(-disc)/(2c0)
            re = -c1 / (2.0 * c0)
            im = math.sqrt(-disc) / (2.0 * c0)
            zero_r = math.sqrt(re * re + im * im)
            zero_theta = math.atan2(abs(im), re)
            zero_hz = zero_theta * SR / (2.0 * math.pi)
            return f"pole {pole_hz:.0f}Hz r={pole_r:.2f}  zero {zero_hz:.0f}Hz r={zero_r:.2f}"
        else:
            return f"pole {pole_hz:.0f}Hz r={pole_r:.2f}  zero(real)"
    return f"pole {pole_hz:.0f}Hz r={pole_r:.2f}  gain={c0:.2f}"


def plot_cartridge(path: Path, out_path: Path) -> None:
    cart = json.loads(path.read_text(encoding="utf-8"))
    freqs = np.logspace(math.log10(F_MIN), math.log10(F_MAX), NFFT)

    by_label = {kf["label"]: kf for kf in cart["keyframes"]}

    fig, axes = plt.subplots(2, 2, figsize=(14, 8), sharex=True, sharey=True)
    fig.suptitle(
        f"{cart['name']}  ·  SR={SR:.0f} Hz  ·  plain-DF2T kernel form",
        fontsize=12, family="monospace",
    )

    stage_colors = plt.cm.viridis(np.linspace(0.15, 0.85, 6))

    for label, row, col in CORNER_LAYOUT:
        ax = axes[row][col]
        kf = by_label.get(label)
        if kf is None:
            ax.set_title(f"{label} — MISSING", color="red")
            continue

        # Overall cascade = product of stage responses
        cascade = np.ones_like(freqs)
        per_stage = []
        for si, stage in enumerate(kf["stages"]):
            r = stage_response(stage, freqs)
            cascade = cascade * r
            per_stage.append((si, stage, r))

        # Plot per-stage (thin, faint)
        for si, stage, r in per_stage:
            if is_passthrough(stage):
                continue
            ax.plot(
                freqs, 20.0 * np.log10(np.maximum(r, 1e-9)),
                color=stage_colors[si], linewidth=0.9, alpha=0.75,
                label=f"S{si+1}: {describe_stage(stage)}",
            )

        # Overall cascade (thick, white with dark outline for contrast)
        ax.plot(
            freqs, 20.0 * np.log10(np.maximum(cascade, 1e-9)),
            color="black", linewidth=2.8, alpha=0.9, label="cascade",
        )

        boost = kf.get("boost", 1.0)
        ax.set_title(f"{label}  ·  boost={boost:g}", family="monospace")
        ax.set_xscale("log")
        ax.set_xlim(F_MIN, F_MAX)
        ax.set_ylim(DB_MIN, DB_MAX)
        ax.grid(True, which="both", alpha=0.25)
        ax.axhline(0, color="gray", lw=0.5, alpha=0.5)
        if row == 1:
            ax.set_xlabel("Hz")
        if col == 0:
            ax.set_ylabel("dB")
        ax.legend(loc="lower left", fontsize=7, framealpha=0.9)

    plt.tight_layout()
    fig.savefig(out_path, dpi=110, bbox_inches="tight")
    plt.close(fig)


def plot_per_stage_strip(path: Path, out_path: Path) -> None:
    """One tall strip per cartridge: rows = stages, cols = corners.

    Easier than the overlay for seeing which stage is doing what at each
    corner — one panel per stage-corner cell, 4 corners across, N stages down.
    """
    cart = json.loads(path.read_text(encoding="utf-8"))
    freqs = np.logspace(math.log10(F_MIN), math.log10(F_MAX), NFFT)

    by_label = {kf["label"]: kf for kf in cart["keyframes"]}
    num_stages = len(next(iter(by_label.values()))["stages"])

    fig, axes = plt.subplots(
        num_stages, 4, figsize=(14, 1.3 * num_stages + 1), sharex=True, sharey=True,
    )
    fig.suptitle(f"{cart['name']}  ·  per-stage × per-corner", fontsize=11, family="monospace")

    corner_order = [lbl for lbl, _, _ in CORNER_LAYOUT]

    for si in range(num_stages):
        for ci, label in enumerate(corner_order):
            ax = axes[si, ci] if num_stages > 1 else axes[ci]
            stage = by_label[label]["stages"][si]
            r = stage_response(stage, freqs)
            passthrough = is_passthrough(stage)
            color = "lightgray" if passthrough else plt.cm.plasma(si / max(1, num_stages - 1))
            ax.plot(freqs, 20.0 * np.log10(np.maximum(r, 1e-9)),
                    color=color, linewidth=1.4)
            ax.set_xscale("log")
            ax.set_xlim(F_MIN, F_MAX)
            ax.set_ylim(DB_MIN, DB_MAX)
            ax.grid(True, which="both", alpha=0.2)
            ax.axhline(0, color="gray", lw=0.4, alpha=0.5)
            desc = describe_stage(stage)
            if si == 0:
                ax.set_title(label, family="monospace", fontsize=10)
            if ci == 0:
                ax.set_ylabel(f"S{si+1}\n{desc}", fontsize=7, family="monospace", rotation=0,
                              ha="right", va="center")
            if si == num_stages - 1:
                ax.set_xlabel("Hz", fontsize=8)

    plt.tight_layout()
    fig.savefig(out_path, dpi=110, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    carts = sorted(IN_DIR.glob("Cheat_*.json"))
    if not carts:
        print(f"No cartridges at {IN_DIR}")
        return 1

    print(f"Plotting {len(carts)} cartridges → {OUT_DIR}")
    for path in carts:
        overlay = OUT_DIR / f"{path.stem}_overlay.png"
        strip = OUT_DIR / f"{path.stem}_stages.png"
        plot_cartridge(path, overlay)
        plot_per_stage_strip(path, strip)
        print(f"  {path.stem}: {overlay.name} + {strip.name}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
