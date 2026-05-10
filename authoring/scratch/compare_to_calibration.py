#!/usr/bin/env python3
"""Render two cartridges (any of compiled-v1, hardware-reference, raw) into a
single 4-row x 2-col grid: rows are corners (M0_Q0, M0_Q100, M100_Q0,
M100_Q100), columns are [candidate, reference]. Same kernel math, same
frequency grid, same y-range per row -- so visual differences are real, not
stylistic.

Usage:
    python compare_to_calibration.py <candidate.json> <reference.json> --out <png>
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Reuse the parser/math in plot_magnitude.
TOOLS = Path(__file__).resolve().parents[1]
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))
import plot_magnitude as pm  # noqa: E402

CORNERS = ("M0_Q0", "M0_Q100", "M100_Q0", "M100_Q100")
F_MIN = 40.0
F_MAX = 20000.0
NPTS = 2048

BG, FG, GRID, AX = "#0b1020", "#cbd5f5", "#1e293b", "#64748b"
STAGE_PALETTE = ("#f97316", "#a855f7", "#22d3ee", "#84cc16", "#ec4899", "#facc15")
CASCADE_LEFT = "#7dd3fc"
CASCADE_RIGHT = "#f59e0b"
PASSTHROUGH_C = (1.0, 0.0, 0.0, 0.0, 0.0)


def cascade_db_per_stage(stages, boost, freqs, sr):
    w = 2.0 * math.pi * freqs / sr
    cascade = np.ones_like(freqs)
    stage_curves = []
    for i, coefs in enumerate(stages):
        if coefs == PASSTHROUGH_C:
            continue
        m = pm._stage_mag(coefs, w)
        cascade = cascade * m
        stage_curves.append((i, m))
    cascade = cascade * max(boost, 1e-9)
    return cascade, stage_curves


def panel(ax, name, label, stages, boost, sr, freqs, cascade_color):
    cascade, stage_curves = cascade_db_per_stage(stages, boost, freqs, sr)

    for i, m in stage_curves:
        db = 20.0 * np.log10(np.maximum(m, 1e-12))
        c = STAGE_PALETTE[i % len(STAGE_PALETTE)]
        ax.plot(freqs, db, color=c, linewidth=0.8, alpha=0.55,
                linestyle="--", label=f"S{i+1}")

    cascade_db = 20.0 * np.log10(np.maximum(cascade, 1e-12))
    ax.plot(freqs, cascade_db, color=cascade_color, linewidth=2.0,
            label="cascade")

    y_lo = float(np.min(cascade_db))
    y_hi = float(np.max(cascade_db))
    for _, m in stage_curves:
        sdb = 20.0 * np.log10(np.maximum(m, 1e-12))
        y_lo = min(y_lo, float(np.min(sdb)))
        y_hi = max(y_hi, float(np.max(sdb)))
    y_lo = max(y_lo - 4.0, -100.0)
    y_hi = min(y_hi + 4.0, 80.0)

    ax.set_xscale("log")
    ax.set_xlim(F_MIN, F_MAX)
    ax.set_ylim(y_lo, y_hi)
    ax.set_facecolor(BG)
    ax.grid(True, which="both", color=GRID, linewidth=0.4, alpha=0.65)
    ax.axhline(0.0, color=AX, linewidth=0.5, alpha=0.55)
    ax.tick_params(colors=AX, labelsize=8)
    for sp in ax.spines.values():
        sp.set_color(GRID)
    ax.set_title(f"{name} | {label}", color=FG, fontsize=10, pad=4)
    ax.set_xlabel("Hz", color=AX, fontsize=9)
    ax.set_ylabel("dB", color=AX, fontsize=9)
    ax.legend(loc="upper left", fontsize=7, facecolor=BG,
              edgecolor=GRID, labelcolor=FG, framealpha=0.8, ncol=2)
    return cascade_db, y_lo, y_hi


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("candidate", type=Path)
    ap.add_argument("reference", type=Path)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--candidate-name", default=None)
    ap.add_argument("--reference-name", default=None)
    args = ap.parse_args()

    import json
    cand_doc = json.loads(args.candidate.read_text(encoding="utf-8"))
    ref_doc = json.loads(args.reference.read_text(encoding="utf-8"))

    cand_sr, cand_map = pm._normalize(cand_doc)
    ref_sr, ref_map = pm._normalize(ref_doc)

    cand_name = args.candidate_name or cand_doc.get("name") or args.candidate.stem
    ref_name = args.reference_name or ref_doc.get("name") or args.reference.stem

    freqs = np.logspace(math.log10(F_MIN), math.log10(F_MAX), NPTS)

    fig, axes = plt.subplots(4, 2, figsize=(13.5, 14.0), dpi=140, facecolor=BG)
    fig.suptitle(f"{cand_name}  vs  {ref_name}",
                 color=FG, fontsize=14, y=0.995)

    for row, label in enumerate(CORNERS):
        if label in cand_map:
            d = cand_map[label]
            panel(axes[row][0], cand_name, label,
                  d["stages"], d["boost"], cand_sr, freqs, CASCADE_LEFT)
        else:
            axes[row][0].set_visible(False)
        if label in ref_map:
            d = ref_map[label]
            panel(axes[row][1], ref_name, label,
                  d["stages"], d["boost"], ref_sr, freqs, CASCADE_RIGHT)
        else:
            axes[row][1].set_visible(False)

    fig.tight_layout(rect=(0, 0, 1, 0.985))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, facecolor=BG)
    plt.close(fig)
    print(f"OUT -> {args.out}")


if __name__ == "__main__":
    main()
