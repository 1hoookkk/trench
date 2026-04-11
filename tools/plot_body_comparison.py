"""Plot frequency responses comparing two bodies at key morph/Q positions.

Usage: PYTHONPATH=. python tools/plot_body_comparison.py \
         cartridges/p2k/P2k_006.json \
         cartridges/candidates/Speaker_Knockerz_v3.keyframe.json
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pyruntime.body import Body
from pyruntime.constants import SR
from pyruntime.freq_response import cascade_response_db, freq_points


POSITIONS = [
    (0.0, 0.0, "M0 Q0 (Vault)"),
    (0.0, 1.0, "M0 Q1 (Vault+Pressure)"),
    (0.5, 0.0, "M0.5 Q0 (Midpoint)"),
    (0.5, 0.5, "M0.5 Q0.5 (Center)"),
    (0.75, 0.0, "M0.75 Q0"),
    (1.0, 0.0, "M1 Q0 (Open)"),
    (1.0, 1.0, "M1 Q1 (Open+Pressure)"),
]


def plot_comparison(path_a: str, path_b: str, out_path: str) -> None:
    body_a = Body.from_json(path_a)
    body_b = Body.from_json(path_b)
    freqs = freq_points(n=512, sr=SR)

    fig, axes = plt.subplots(4, 2, figsize=(16, 20))
    fig.suptitle(f"{body_a.name}  vs  {body_b.name}", fontsize=14, fontweight="bold")
    axes_flat = axes.flatten()

    for i, (m, q, label) in enumerate(POSITIONS):
        ax = axes_flat[i]
        enc_a = body_a.corners.interpolate(m, q)
        enc_b = body_b.corners.interpolate(m, q)
        db_a = cascade_response_db(enc_a, freqs, SR)
        db_b = cascade_response_db(enc_b, freqs, SR)

        ax.semilogx(freqs, db_a, color="#666666", linewidth=1.0, alpha=0.7, label=body_a.name)
        ax.semilogx(freqs, db_b, color="#e04040", linewidth=1.5, label=body_b.name)
        ax.set_title(label, fontsize=11)
        ax.set_xlim(20, 19000)
        ax.set_ylim(-60, 50)
        ax.set_xlabel("Hz")
        ax.set_ylabel("dB")
        ax.axhline(0, color="#444", linewidth=0.5, linestyle="--")
        ax.axhline(35, color="#cc0000", linewidth=0.5, linestyle=":", alpha=0.5)
        ax.legend(fontsize=8, loc="upper right")
        ax.grid(True, alpha=0.3)

    # Hide unused subplot
    axes_flat[-1].set_visible(False)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"wrote {out_path}")
    plt.close()


def main() -> None:
    if len(sys.argv) < 3:
        print("Usage: plot_body_comparison.py <body_a.json> <body_b.json>")
        sys.exit(1)

    path_a, path_b = sys.argv[1], sys.argv[2]
    out = str(ROOT / "cartridges" / "candidates" / "Speaker_Knockerz_v3_comparison.png")
    plot_comparison(path_a, path_b, out)


if __name__ == "__main__":
    main()
