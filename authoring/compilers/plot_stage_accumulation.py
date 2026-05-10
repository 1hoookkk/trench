#!/usr/bin/env python3
"""Plot cumulative stage buildup for compiled-v1 cartridge endpoints."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


F_MIN = 40.0
F_MAX = 20000.0
NPTS = 2048
DEFAULT_STAGE_COUNTS = (1, 2, 3, 6)


def _stage_response(stage: dict, freqs_hz: np.ndarray, sample_rate: float) -> np.ndarray:
    c0 = float(stage["c0"])
    c1 = float(stage["c1"])
    c2 = float(stage["c2"])
    c3 = float(stage["c3"])
    c4 = float(stage["c4"])
    w = 2.0 * np.pi * freqs_hz / sample_rate
    e1 = np.exp(-1j * w)
    e2 = np.exp(-2j * w)
    num = c0 + c1 * e1 + c2 * e2
    den = 1.0 + c3 * e1 + c4 * e2
    return np.abs(num / den)


def _cumulative_curves(keyframe: dict, sample_rate: float, stage_counts: tuple[int, ...]) -> tuple[np.ndarray, dict[int, np.ndarray]]:
    freqs = np.logspace(math.log10(F_MIN), math.log10(F_MAX), NPTS)
    active_stages = []
    for stage in keyframe["stages"]:
        vals = tuple(float(stage[k]) for k in ("c0", "c1", "c2", "c3", "c4"))
        if vals == (1.0, 0.0, 0.0, 0.0, 0.0):
            continue
        active_stages.append(stage)

    cascade = np.ones_like(freqs)
    out: dict[int, np.ndarray] = {}
    stage_limit_set = set(stage_counts)
    for idx, stage in enumerate(active_stages, 1):
        cascade = cascade * _stage_response(stage, freqs, sample_rate)
        if idx in stage_limit_set:
            out[idx] = cascade.copy()
    return freqs, out


def _probe_db(freqs: np.ndarray, mag: np.ndarray, target_hz: float) -> float:
    idx = int(np.argmin(np.abs(freqs - target_hz)))
    return 20.0 * math.log10(max(float(mag[idx]), 1e-12))


def main() -> int:
    ap = argparse.ArgumentParser(description="Plot cumulative response after selected stage counts for M0/M100 endpoints.")
    ap.add_argument("compiled", type=Path, help="compiled-v1 cartridge JSON")
    ap.add_argument("--out", type=Path, required=True, help="output PNG")
    ap.add_argument("--labels", default="M0_Q0,M100_Q0", help="comma-separated keyframe labels to plot")
    ap.add_argument("--stage-counts", default="1,2,3,6", help="comma-separated cumulative stage counts")
    ap.add_argument("--split-stages", action="store_true", help="render each stage-count in its own subplot instead of overlaying")
    args = ap.parse_args()

    cart = json.loads(args.compiled.read_text(encoding="utf-8"))
    keyframes = {kf["label"]: kf for kf in cart["keyframes"]}
    labels = tuple(part.strip() for part in args.labels.split(",") if part.strip())
    stage_counts = tuple(sorted({int(part.strip()) for part in args.stage_counts.split(",") if part.strip()}))
    sample_rate = float(cart.get("sampleRate", 44100.0))

    colors = {
        1: "#f59e0b",
        2: "#ef4444",
        3: "#8b5cf6",
        6: "#7dd3fc",
    }

    summaries: list[str] = []
    overall_floor = 1e9
    overall_ceil = -1e9

    if args.split_stages:
        rows = len(labels)
        cols = len(stage_counts)
        fig, axes = plt.subplots(rows, cols, figsize=(4.2 * cols, 3.2 * rows), dpi=140, facecolor="#0b1020")
        axes_grid = np.atleast_2d(axes).reshape(rows, cols)
        all_axes = list(axes_grid.flat)
    else:
        fig, axes = plt.subplots(1, len(labels), figsize=(6.2 * len(labels), 4.4), dpi=140, facecolor="#0b1020")
        axes_grid = np.atleast_1d(axes).reshape(1, -1)
        all_axes = list(axes_grid.flat)

    for row_idx, label in enumerate(labels):
        freqs, curves = _cumulative_curves(keyframes[label], sample_rate, stage_counts)

        if args.split_stages:
            for col_idx, stage_count in enumerate(stage_counts):
                ax = axes_grid[row_idx, col_idx]
                if stage_count in curves:
                    db = 20.0 * np.log10(np.maximum(curves[stage_count], 1e-12))
                    overall_floor = min(overall_floor, float(np.min(db)))
                    overall_ceil = max(overall_ceil, float(np.max(db)))
                    color = colors.get(stage_count, "#cbd5f5")
                    ax.plot(freqs, db, color=color, linewidth=1.8, label=f"after stage {stage_count}")
                ax.set_title(f"{label} — after stage {stage_count}", color="#cbd5f5", fontsize=9, pad=4)
                ax.set_xscale("log")
                ax.set_xlim(F_MIN, F_MAX)
                ax.set_facecolor("#0b1020")
                ax.grid(True, which="both", color="#1e293b", linewidth=0.4, alpha=0.7)
                ax.tick_params(colors="#64748b", labelsize=7)
                for sp in ax.spines.values():
                    sp.set_color("#1e293b")
                ax.set_xlabel("Hz", color="#64748b", fontsize=8)
                ax.set_ylabel("dB", color="#64748b", fontsize=8)
        else:
            ax = axes_grid[0, row_idx]
            for stage_count in stage_counts:
                if stage_count not in curves:
                    continue
                db = 20.0 * np.log10(np.maximum(curves[stage_count], 1e-12))
                overall_floor = min(overall_floor, float(np.min(db)))
                overall_ceil = max(overall_ceil, float(np.max(db)))
                color = colors.get(stage_count, "#cbd5f5")
                ax.plot(freqs, db, color=color, linewidth=1.8, label=f"after stage {stage_count}")
            ax.set_xscale("log")
            ax.set_xlim(F_MIN, F_MAX)
            ax.set_facecolor("#0b1020")
            ax.grid(True, which="both", color="#1e293b", linewidth=0.4, alpha=0.7)
            ax.tick_params(colors="#64748b", labelsize=7)
            for sp in ax.spines.values():
                sp.set_color("#1e293b")
            ax.set_title(label, color="#cbd5f5", fontsize=10, pad=5)
            ax.set_xlabel("Hz", color="#64748b", fontsize=8)
            ax.set_ylabel("dB", color="#64748b", fontsize=8)
            ax.legend(loc="upper left", fontsize=7, facecolor="#0b1020", edgecolor="#1e293b", labelcolor="#cbd5f5")

        final_curve = curves[max(curves)]
        low_80 = _probe_db(freqs, final_curve, 80.0)
        low_120 = _probe_db(freqs, final_curve, 120.0)
        two_stage = curves.get(2, final_curve)
        two_80 = _probe_db(freqs, two_stage, 80.0)
        summaries.append(f"{label}: stage2@80Hz={two_80:.1f} dB  full@80Hz={low_80:.1f} dB  full@120Hz={low_120:.1f} dB")

    y_lo = max(overall_floor - 6.0, -80.0)
    y_hi = min(overall_ceil + 4.0, 100.0)
    for ax in all_axes:
        ax.set_ylim(y_lo, y_hi)

    summary_text = "\n".join(summaries)
    fig.suptitle(f"{cart['name']} — cumulative stage buildup\n{summary_text}", color="#cbd5f5", fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, facecolor="#0b1020", dpi=140)
    print(args.out)
    for line in summaries:
        print(line)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
