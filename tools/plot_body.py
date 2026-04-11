"""Plot frequency response of a TRENCH body at multiple morph/Q positions.

Usage:
    python tools/plot_body.py cartridges/Acid_Squelch_77.json
"""

import json
import math
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams.update({
    "figure.facecolor": "#0a0a14",
    "axes.facecolor": "#0a0a14",
    "axes.edgecolor": "#333355",
    "axes.labelcolor": "#8888aa",
    "xtick.color": "#666688",
    "ytick.color": "#666688",
    "text.color": "#ccccee",
    "grid.color": "#1a1a2e",
    "grid.alpha": 0.6,
    "font.family": "monospace",
    "font.size": 9,
})

SR = 39062.5
NUM_FREQS = 2048


def load_body(path: str) -> dict:
    with open(path) as f:
        return json.load(f)


def get_corner_coeffs(body: dict, label: str) -> list[dict]:
    """Extract stage coefficients for a corner label (e.g. 'M0_Q0')."""
    for kf in body.get("keyframes", []):
        if kf["label"] == label:
            return kf["stages"]
    # Try corners format
    if "corners" in body:
        corner = body["corners"].get(label, [])
        if isinstance(corner, list) and len(corner) > 0:
            if isinstance(corner[0], list):
                return [{"c0": s[0], "c1": s[1], "c2": s[2], "c3": s[3], "c4": s[4]}
                        for s in corner]
            return corner
    return []


def interpolate_coeffs(body: dict, morph: float, q: float) -> list[dict]:
    """Bilinear interpolation: Q first, then morph (per spec)."""
    corners = {
        "M0_Q0": get_corner_coeffs(body, "M0_Q0"),
        "M100_Q0": get_corner_coeffs(body, "M100_Q0"),
        "M0_Q100": get_corner_coeffs(body, "M0_Q100"),
        "M100_Q100": get_corner_coeffs(body, "M100_Q100"),
    }

    num_stages = len(corners["M0_Q0"])
    result = []
    for i in range(num_stages):
        stage = {}
        for key in ("c0", "c1", "c2", "c3", "c4"):
            m0_q0 = corners["M0_Q0"][i][key]
            m100_q0 = corners["M100_Q0"][i][key]
            m0_q100 = corners["M0_Q100"][i][key]
            m100_q100 = corners["M100_Q100"][i][key]

            # Q first
            q_m0 = m0_q0 + (m0_q100 - m0_q0) * q
            q_m1 = m100_q0 + (m100_q100 - m100_q0) * q
            # Then morph
            stage[key] = q_m0 + (q_m1 - q_m0) * morph
        result.append(stage)
    return result


def cascade_response(stages: list[dict], freqs: np.ndarray) -> np.ndarray:
    """Compute cascade magnitude response in dB.

    Kernel-form DF2T transfer function:
    H(z) = (c0 + c1*z^-1 + c2*z^-2) / (1 + c3*z^-1 + c4*z^-2)
    """
    z_inv = np.exp(-2j * np.pi * freqs / SR)
    z_inv2 = z_inv * z_inv

    total = np.ones(len(freqs), dtype=complex)
    for s in stages:
        c0, c1, c2, c3, c4 = s["c0"], s["c1"], s["c2"], s["c3"], s["c4"]
        # Skip passthrough stages: c0=1, c1-c4=0
        if (abs(c0 - 1.0) < 1e-9 and abs(c1) < 1e-9 and
                abs(c2) < 1e-9 and abs(c3) < 1e-9 and abs(c4) < 1e-9):
            continue

        num = c0 + c1 * z_inv + c2 * z_inv2
        den = 1.0 + c3 * z_inv + c4 * z_inv2
        H = num / np.where(np.abs(den) > 1e-30, den, 1e-30)
        total *= H

    mag_db = 20.0 * np.log10(np.maximum(np.abs(total), 1e-30))
    return mag_db


def plot_body(path: str):
    body = load_body(path)
    name = body.get("name", os.path.basename(path))

    freqs = np.geomspace(20.0, SR / 2.0 - 1.0, NUM_FREQS)

    # Morph sweep at Q=0 (7 positions)
    morph_positions = [0.0, 0.17, 0.33, 0.5, 0.67, 0.83, 1.0]
    cmap = plt.cm.cool

    fig, axes = plt.subplots(2, 2, figsize=(14, 9))
    fig.suptitle(f"TRENCH · {name}", fontsize=14, fontweight="bold", color="#3a7aff")

    # --- Panel 1: Morph sweep at Q=0 ---
    ax = axes[0, 0]
    for i, m in enumerate(morph_positions):
        stages = interpolate_coeffs(body, m, 0.0)
        mag = cascade_response(stages, freqs)
        color = cmap(i / (len(morph_positions) - 1))
        ax.semilogx(freqs, mag, color=color, linewidth=1.2,
                     label=f"M={m:.0%}", alpha=0.85)
    ax.set_title("Morph sweep · Q=0%", fontsize=10)
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Magnitude (dB)")
    ax.set_xlim(20, SR / 2)
    ax.set_ylim(-60, 80)
    ax.legend(fontsize=7, loc="upper right", framealpha=0.3)
    ax.grid(True, which="both", linewidth=0.5)

    # --- Panel 2: Morph sweep at Q=100 ---
    ax = axes[0, 1]
    for i, m in enumerate(morph_positions):
        stages = interpolate_coeffs(body, m, 1.0)
        mag = cascade_response(stages, freqs)
        color = cmap(i / (len(morph_positions) - 1))
        ax.semilogx(freqs, mag, color=color, linewidth=1.2,
                     label=f"M={m:.0%}", alpha=0.85)
    ax.set_title("Morph sweep · Q=100%", fontsize=10)
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Magnitude (dB)")
    ax.set_xlim(20, SR / 2)
    ax.set_ylim(-60, 80)
    ax.legend(fontsize=7, loc="upper right", framealpha=0.3)
    ax.grid(True, which="both", linewidth=0.5)

    # --- Panel 3: Q sweep at M=0 ---
    ax = axes[1, 0]
    q_positions = [0.0, 0.25, 0.5, 0.75, 1.0]
    cmap_q = plt.cm.autumn
    for i, qv in enumerate(q_positions):
        stages = interpolate_coeffs(body, 0.0, qv)
        mag = cascade_response(stages, freqs)
        color = cmap_q(i / (len(q_positions) - 1))
        ax.semilogx(freqs, mag, color=color, linewidth=1.2,
                     label=f"Q={qv:.0%}", alpha=0.85)
    ax.set_title("Q sweep · M=0%", fontsize=10)
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Magnitude (dB)")
    ax.set_xlim(20, SR / 2)
    ax.set_ylim(-60, 80)
    ax.legend(fontsize=7, loc="upper right", framealpha=0.3)
    ax.grid(True, which="both", linewidth=0.5)

    # --- Panel 4: Diagonal sweep (morph=Q) ---
    ax = axes[1, 1]
    diag_positions = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
    cmap_d = plt.cm.spring
    for i, d in enumerate(diag_positions):
        stages = interpolate_coeffs(body, d, d)
        mag = cascade_response(stages, freqs)
        color = cmap_d(i / (len(diag_positions) - 1))
        ax.semilogx(freqs, mag, color=color, linewidth=1.2,
                     label=f"M=Q={d:.0%}", alpha=0.85)
    ax.set_title("Diagonal sweep · M=Q", fontsize=10)
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Magnitude (dB)")
    ax.set_xlim(20, SR / 2)
    ax.set_ylim(-60, 80)
    ax.legend(fontsize=7, loc="upper right", framealpha=0.3)
    ax.grid(True, which="both", linewidth=0.5)

    plt.tight_layout()
    out_path = path.rsplit(".", 1)[0] + "_plot.png"
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Saved: {out_path}")
    plt.show()


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python tools/plot_body.py <body.json>")
        sys.exit(1)
    plot_body(sys.argv[1])
