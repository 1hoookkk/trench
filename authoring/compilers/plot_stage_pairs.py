#!/usr/bin/env python3
"""Plot each stage's M0 vs M100 pair response individually for a compiled-v1 cartridge.

One subplot per active stage. Each subplot overlays the stage's M0_Q0 response
and M100_Q0 response, so you can read how the Lo/Hi pair moves per slot.
"""

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
PASSTHROUGH = (1.0, 0.0, 0.0, 0.0, 0.0)
M0_COLOR = "#7dd3fc"
M100_COLOR = "#86efac"


def _stage_mag(stage: dict, freqs_hz: np.ndarray, sample_rate: float) -> np.ndarray:
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


def _active_stages(keyframe: dict) -> list[dict]:
    out = []
    for stage in keyframe["stages"]:
        vals = tuple(float(stage[k]) for k in ("c0", "c1", "c2", "c3", "c4"))
        if vals == PASSTHROUGH:
            continue
        out.append(stage)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description="Plot M0 vs M100 per-stage pair responses for a compiled-v1 cartridge.")
    ap.add_argument("compiled", type=Path, help="compiled-v1 cartridge JSON")
    ap.add_argument("--out", type=Path, required=True, help="output PNG")
    ap.add_argument("--m0-label", default="M0_Q0")
    ap.add_argument("--m100-label", default="M100_Q0")
    args = ap.parse_args()

    cart = json.loads(args.compiled.read_text(encoding="utf-8"))
    keyframes = {kf["label"]: kf for kf in cart["keyframes"]}
    sample_rate = float(cart.get("sampleRate", 44100.0))

    m0_stages = _active_stages(keyframes[args.m0_label])
    m100_stages = _active_stages(keyframes[args.m100_label])
    n = min(len(m0_stages), len(m100_stages))
    if n == 0:
        raise SystemExit("no active stages to plot")

    freqs = np.logspace(math.log10(F_MIN), math.log10(F_MAX), NPTS)

    raw_doc_hint = cart.get("name", "")
    roles: list[str] = []
    src_candidates = [
        Path("cartridges/factory/_source/bench_raw_bank") / f"{raw_doc_hint}.raw.json",
    ]
    for candidate in src_candidates:
        if candidate.exists():
            try:
                raw = json.loads(candidate.read_text(encoding="utf-8"))
                m0_raw = raw.get("frames", {}).get("M0", {}).get("stages", [])
                roles = [str(st.get("role", "")) for st in m0_raw]
                break
            except Exception:
                pass

    cols = min(n, 3)
    rows = math.ceil(n / cols)
    fig, axes = plt.subplots(rows, cols, figsize=(4.4 * cols, 3.2 * rows), dpi=140, facecolor="#0b1020")
    axes_arr = np.atleast_1d(axes).reshape(rows, cols)

    floor, ceil = 1e9, -1e9
    for idx in range(n):
        ax = axes_arr.flat[idx]
        m0_mag = _stage_mag(m0_stages[idx], freqs, sample_rate)
        m100_mag = _stage_mag(m100_stages[idx], freqs, sample_rate)
        m0_db = 20.0 * np.log10(np.maximum(m0_mag, 1e-12))
        m100_db = 20.0 * np.log10(np.maximum(m100_mag, 1e-12))
        floor = min(floor, float(np.min(m0_db)), float(np.min(m100_db)))
        ceil = max(ceil, float(np.max(m0_db)), float(np.max(m100_db)))

        ax.plot(freqs, m0_db, color=M0_COLOR, linewidth=1.6, label="M0 (Lo)")
        ax.plot(freqs, m100_db, color=M100_COLOR, linewidth=1.6, label="M100 (Hi)")
        role = roles[idx] if idx < len(roles) and roles[idx] else f"stage {idx+1}"
        ax.set_title(f"stage {idx+1} — {role}", color="#cbd5f5", fontsize=10, pad=4)
        ax.set_xscale("log")
        ax.set_xlim(F_MIN, F_MAX)
        ax.set_facecolor("#0b1020")
        ax.grid(True, which="both", color="#1e293b", linewidth=0.4, alpha=0.7)
        ax.tick_params(colors="#64748b", labelsize=7)
        for sp in ax.spines.values():
            sp.set_color("#1e293b")
        ax.set_xlabel("Hz", color="#64748b", fontsize=8)
        ax.set_ylabel("dB", color="#64748b", fontsize=8)
        ax.legend(loc="upper left", fontsize=7, facecolor="#0b1020", edgecolor="#1e293b", labelcolor="#cbd5f5")

    for extra_ax in axes_arr.flat[n:]:
        extra_ax.axis("off")

    y_lo = max(floor - 4.0, -80.0)
    y_hi = min(ceil + 3.0, 90.0)
    for ax in axes_arr.flat[:n]:
        ax.set_ylim(y_lo, y_hi)

    fig.suptitle(f"{cart.get('name','')} — per-stage Lo/Hi pair", color="#cbd5f5", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, facecolor="#0b1020", dpi=140)
    print(args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
