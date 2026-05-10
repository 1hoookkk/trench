#!/usr/bin/env python3
"""Trusted oscilloscope for TRENCH filter cartridges.

Accepts two JSON shapes:
  1. compiled-v1 cartridge          -- `format: "compiled-v1"` + `keyframes`
                                       (each stage: c0..c4 direct)
  2. hardware-reference JSON        -- no `format` field; `corners` with
                                       val1/val2/val3/pole_freq_hz/radius
                                       (used for Ear_Bender, Talking_Hedz etc.)

Magnitude math matches tools/compile_raw.py kernel form:
  H(z) = (c0 + c1 z^-1 + c2 z^-2) / (1 + c3 z^-1 + c4 z^-2)

Per-corner panel draws:
  - thin dashed per-stage magnitudes
  - bold composite cascade (boost applied here)

Usage:
  python tools/plot_magnitude.py path/to/cartridge.json
  python tools/plot_magnitude.py path/to/cartridge.json --corner M0_Q0
  python tools/plot_magnitude.py path/to/cartridge.json --split
  python tools/plot_magnitude.py authoring/plots/Ear_Bender.json
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

_THIS = Path(__file__).resolve().parent
if str(_THIS) not in sys.path:
    sys.path.insert(0, str(_THIS))
import bark  # noqa: E402

F_MIN = 40.0
F_MAX = 20000.0
NPTS = 2048
LANDMARK_COLOR = "#94a3b8"
LANDMARK_LINE = "#475569"
CORNERS = ("M0_Q0", "M0_Q100", "M100_Q0", "M100_Q100")
CORNER_COLORS = {
    "M0_Q0":     "#7dd3fc",
    "M0_Q100":   "#38bdf8",
    "M100_Q0":   "#fb7185",
    "M100_Q100": "#f59e0b",
}
STAGE_PALETTE = ("#f97316", "#a855f7", "#22d3ee", "#84cc16", "#ec4899", "#facc15")
BG, FG, GRID, AX = "#0b1020", "#cbd5f5", "#1e293b", "#64748b"
PASSTHROUGH_C = (1.0, 0.0, 0.0, 0.0, 0.0)


def _stage_from_hwref(stage: dict, sr: float) -> tuple[float, float, float, float, float]:
    f = float(stage["pole_freq_hz"])
    r = float(stage["radius"])
    v1 = float(stage.get("val1", 0.0))
    v2 = float(stage.get("val2", 0.0))
    v3 = float(stage.get("val3", 0.0))
    theta = 2.0 * math.pi * f / sr
    a1 = -2.0 * r * math.cos(theta)
    a2 = r * r
    return (1.0 + v1, a1 + v2, a2 - v3, a1, a2)


def _stage_from_compiled(stage: dict) -> tuple[float, float, float, float, float]:
    return (
        float(stage["c0"]), float(stage["c1"]), float(stage["c2"]),
        float(stage["c3"]), float(stage["c4"]),
    )

def _stage_from_raw(stage: dict) -> tuple[float, float, float, float, float]:
    # raw stages use a1, r, val1, val2, val3
    # mapping to c0..c4:
    # c0 = val1, c1 = val2, c2 = val3, c3 = a1, c4 = r
    return (
        float(stage.get("val1", 0.0)),
        float(stage.get("val2", 0.0)),
        float(stage.get("val3", 0.0)),
        float(stage.get("a1", 0.0)),
        float(stage.get("r", 0.0))
    )


def _stage_mag(coefs: tuple, w: np.ndarray) -> np.ndarray:
    c0, c1, c2, c3, c4 = coefs
    e1 = np.exp(-1j * w)
    e2 = np.exp(-2j * w)
    num = c0 + c1 * e1 + c2 * e2
    den = 1.0 + c3 * e1 + c4 * e2
    return np.abs(num / den)


def _is_passthrough(coefs: tuple) -> bool:
    return coefs == PASSTHROUGH_C


def _normalize(doc: dict) -> tuple[float, dict[str, dict]]:
    """Return (sr_hz, {label: {"stages": [(c0..c4), ...], "boost": float}})."""
    fmt = doc.get("format")
    sr = float(doc.get("sampleRate", doc.get("sample_rate", 39062.5)))
    out: dict[str, dict] = {}
    
    # Handle keyframes (compiled-v1 or raw survivors)
    if "keyframes" in doc:
        for kf in doc["keyframes"]:
            label = kf.get("label")
            if label not in CORNERS:
                continue
            
            # Detect stage format
            raw_stages = kf.get("stages", [])
            if not raw_stages:
                continue
                
            if "c0" in raw_stages[0]:
                stages = [_stage_from_compiled(s) for s in raw_stages]
            else:
                stages = [_stage_from_raw(s) for s in raw_stages]
                
            boost = float(kf.get("boost", 1.0))
            out[label] = {"stages": stages, "boost": boost}
            
        if out:
            return sr, out

    # hardware-reference
    global_boost = float(doc.get("boost", 1.0))
    for label, corner in doc.get("corners", {}).items():
        if label not in CORNERS:
            continue
        stages = [_stage_from_hwref(s, sr) for s in corner.get("stages", [])]
        out[label] = {"stages": stages, "boost": global_boost}
    return sr, out


def _draw_landmark_grid(ax, x_vals: np.ndarray, y_lo: float, y_hi: float,
                        xaxis: str) -> None:
    grid = bark.plot_grid()
    for g in grid:
        if xaxis == "bark":
            x = g["bark"]
        else:
            x = g["freq_hz"]
        if x < x_vals[0] or x > x_vals[-1]:
            continue
        ax.axvline(x, color=LANDMARK_LINE, linewidth=0.55, alpha=0.55,
                   linestyle=":")
        ax.text(x, y_hi - (y_hi - y_lo) * 0.04, g["label"],
                color=LANDMARK_COLOR, fontsize=7, ha="center",
                va="top", alpha=0.85)


def _render(labels: list[str], corner_map: dict, sr: float, name: str,
            out_path: Path, xaxis: str = "hz") -> Path:
    n = len(labels)
    if n == 1:
        rows, cols, fig_w, fig_h = 1, 1, 7.2, 4.6
    elif n == 2:
        rows, cols, fig_w, fig_h = 1, 2, 12.0, 4.4
    else:
        rows, cols, fig_w, fig_h = 2, 2, 12.0, 8.6

    fig, axes = plt.subplots(rows, cols, figsize=(fig_w, fig_h), dpi=140, facecolor=BG)
    axes_arr = np.atleast_1d(axes).flatten()

    freqs = np.logspace(math.log10(F_MIN), math.log10(F_MAX), NPTS)
    w = 2.0 * math.pi * freqs / sr
    x_vals = np.asarray(bark.hz_to_bark(freqs)) if xaxis == "bark" else freqs
    x_lo = float(x_vals[0])
    x_hi = float(x_vals[-1])

    for ax, label in zip(axes_arr[:n], labels):
        data = corner_map[label]
        stages = data["stages"]
        boost = data["boost"]

        cascade = np.ones_like(freqs)
        stage_curves: list[tuple[int, np.ndarray]] = []
        for i, coefs in enumerate(stages):
            if _is_passthrough(coefs):
                continue
            m = _stage_mag(coefs, w)
            cascade = cascade * m
            stage_curves.append((i, m))
        cascade = cascade * max(boost, 1e-9)

        for i, m in stage_curves:
            db = 20.0 * np.log10(np.maximum(m, 1e-12))
            color = STAGE_PALETTE[i % len(STAGE_PALETTE)]
            ax.plot(x_vals, db, color=color, linewidth=0.9, alpha=0.7,
                    linestyle="--", label=f"S{i+1}")

        cascade_db = 20.0 * np.log10(np.maximum(cascade, 1e-12))
        corner_color = CORNER_COLORS.get(label, "#cbd5f5")
        ax.plot(x_vals, cascade_db, color=corner_color, linewidth=2.0,
                label="cascade")

        y_lo = float(np.min(cascade_db))
        y_hi = float(np.max(cascade_db))
        for _, m in stage_curves:
            sdb = 20.0 * np.log10(np.maximum(m, 1e-12))
            y_lo = min(y_lo, float(np.min(sdb)))
            y_hi = max(y_hi, float(np.max(sdb)))
        y_lo = max(y_lo - 4.0, -100.0)
        y_hi = min(y_hi + 4.0, 80.0)

        if xaxis == "bark":
            ax.set_xscale("linear")
        else:
            ax.set_xscale("log")
        ax.set_xlim(x_lo, x_hi)
        ax.set_ylim(y_lo, y_hi)
        ax.set_facecolor(BG)
        ax.grid(True, which="both", color=GRID, linewidth=0.4, alpha=0.65)
        ax.axhline(0.0, color=AX, linewidth=0.55, alpha=0.65)
        ax.tick_params(colors=AX, labelsize=8)
        for sp in ax.spines.values():
            sp.set_color(GRID)
        _draw_landmark_grid(ax, x_vals, y_lo, y_hi, xaxis)
        ax.set_title(f"{name} - {label}", color=FG, fontsize=10, pad=6)
        ax.set_xlabel("Bark" if xaxis == "bark" else "Hz",
                      color=AX, fontsize=9)
        ax.set_ylabel("dB", color=AX, fontsize=9)
        ax.legend(loc="upper left", fontsize=7, facecolor=BG, edgecolor=GRID,
                  labelcolor=FG, framealpha=0.8, ncol=2)

    for ax in axes_arr[n:]:
        ax.set_visible(False)

    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, facecolor=BG)
    plt.close(fig)
    return out_path


def plot(doc: dict, out_path: Path, selected: list[str] | None = None,
         split: bool = False, xaxis: str = "hz") -> list[Path]:
    sr, corner_map = _normalize(doc)
    name = doc.get("name", out_path.stem)

    labels = [l for l in CORNERS if l in corner_map]
    if selected:
        labels = [l for l in selected if l in corner_map]
    if not labels:
        raise SystemExit("no matching corners found")

    if split:
        outs: list[Path] = []
        stem = out_path.stem if out_path.suffix == ".png" else out_path.name
        base_dir = out_path.parent if out_path.suffix == ".png" else out_path
        base_dir.mkdir(parents=True, exist_ok=True)
        for label in labels:
            path = base_dir / f"{stem}_{label}.png"
            _render([label], corner_map, sr, name, path, xaxis=xaxis)
            outs.append(path)
        return outs

    _render(labels, corner_map, sr, name, out_path, xaxis=xaxis)
    return [out_path]


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Plot magnitude response for compiled-v1 or hardware-reference JSON."
    )
    ap.add_argument("json_path", type=Path,
                    help="compiled-v1 cartridge or hardware-reference JSON")
    ap.add_argument("--out", type=Path, default=None,
                    help="Output PNG path (or directory when --split is given).")
    ap.add_argument("--corner", action="append", default=None,
                    help="Restrict to one corner label. Repeat for multiple.")
    ap.add_argument("--split", action="store_true",
                    help="Write one PNG per corner instead of a grid.")
    ap.add_argument("--xaxis", choices=("hz", "bark"), default="hz",
                    help="Frequency axis: 'hz' log-frequency (default) "
                         "or 'bark' Traunmueller critical-band scale with "
                         "TRENCH landmark overlays.")
    args = ap.parse_args()

    doc = json.loads(args.json_path.read_text(encoding="utf-8"))
    default_out = args.json_path.with_name(f"{args.json_path.stem}__mag.png")
    out = args.out or default_out
    paths = plot(doc, out, selected=args.corner, split=args.split,
                 xaxis=args.xaxis)
    for p in paths:
        print(p)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
