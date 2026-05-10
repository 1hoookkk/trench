"""Magnitude audit for the Q-law compiled AH->EE cartridge.

Plots per-corner composite magnitude + per-stage thin lines + role labels,
so we can read whether the Q-law (uniform c4, log-modulated zeros) still
delivers six legible stage roles after gate passage.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


SCRATCH = Path(__file__).resolve().parent
COMPILED = SCRATCH / "ah_to_ee_qlaw.compiled.json"
PLAN = SCRATCH / "ah_to_ee_planner_v2.md.json"
F_MIN, F_MAX, NPTS = 40.0, 18000.0, 2048


def stage_mag(stage: dict, w: np.ndarray) -> np.ndarray:
    c0, c1, c2, c3, c4 = (float(stage[f"c{i}"]) for i in range(5))
    e1 = np.exp(-1j * w)
    e2 = np.exp(-2j * w)
    return np.abs((c0 + c1 * e1 + c2 * e2) / (1.0 + c3 * e1 + c4 * e2))


def main() -> int:
    cart = json.loads(COMPILED.read_text(encoding="utf-8"))
    plan = json.loads(PLAN.read_text(encoding="utf-8"))
    sr = float(cart["sampleRate"])
    freqs = np.geomspace(F_MIN, F_MAX, NPTS)
    w = 2.0 * np.pi * freqs / sr

    landmarks = [s["planner_landmark"] for s in plan["stages"]]
    shapes = [s["shape"] for s in plan["stages"]]

    corners = ["M0_Q0", "M100_Q0"]
    fig, axes = plt.subplots(len(corners), 1, figsize=(11.0, 8.0), dpi=150, sharex=True)
    fig.suptitle(f"{cart['name']} — Q-law: uniform c4, log-modulated zeros")

    for ax, corner_label in zip(axes, corners):
        kf = next(k for k in cart["keyframes"] if k["label"] == corner_label)
        boost = float(kf["boost"])
        composite = np.ones_like(freqs)
        for idx, st in enumerate(kf["stages"]):
            if st["c3"] == 0.0 and st["c4"] == 0.0:
                continue
            mag = stage_mag(st, w)
            composite *= mag
            if idx < len(landmarks):
                label = f"S{idx} {shapes[idx]:5s} {landmarks[idx]}"
            else:
                label = f"S{idx}"
            ax.plot(freqs, 20 * np.log10(np.maximum(mag, 1e-9)), lw=0.9, alpha=0.55, label=label)
        composite_db = 20 * np.log10(np.maximum(composite * boost, 1e-9))
        ax.plot(freqs, composite_db, color="#111", lw=2.2, label=f"composite (boost={boost})")
        ax.set_xscale("log")
        ax.set_title(f"{corner_label}  (c4={kf['stages'][0]['c4']:.6f} uniform across all 6 active stages)")
        ax.set_ylabel("dB")
        ax.set_ylim(-60, 100)
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(fontsize=7, loc="lower left", ncol=2)

    axes[-1].set_xlabel("Hz")
    out = SCRATCH / "ah_to_ee_qlaw.magnitude.png"
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)
    print(out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
