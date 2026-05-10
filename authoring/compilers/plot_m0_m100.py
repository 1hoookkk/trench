#!/usr/bin/env python3
"""Plot full-cascade M0 and M100 responses at a chosen Q slice."""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

_TOOLS = Path(__file__).resolve().parent
if str(_TOOLS) not in sys.path:
    sys.path.insert(0, str(_TOOLS))

from cube_authoring.thumbnail import magnitude_spectrum_db, run_impulse  # noqa: E402


def _interp_stage(a: dict, b: dict, t: float) -> tuple[float, float, float, float, float]:
    return tuple(float(a[k]) * (1.0 - t) + float(b[k]) * t for k in ("c0", "c1", "c2", "c3", "c4"))


def _curve(cart: dict, morph_label: str, q: float) -> tuple[np.ndarray, np.ndarray]:
    keyframes = {kf["label"]: kf for kf in cart["keyframes"]}
    q0 = keyframes[f"{morph_label}_Q0"]
    q1 = keyframes[f"{morph_label}_Q100"]
    stages = [_interp_stage(a, b, q) for a, b in zip(q0["stages"], q1["stages"])]
    boost = float(q0.get("boost", 1.0)) * (1.0 - q) + float(q1.get("boost", 1.0)) * q
    sample_rate = int(round(float(cart.get("sampleRate", 44100.0))))
    ir = run_impulse(stages, 4096, boost)
    f, db = magnitude_spectrum_db(ir, sample_rate)
    mask = (f >= 40.0) & (f <= 20000.0)
    return f[mask], db[mask]


def main() -> int:
    ap = argparse.ArgumentParser(description="Plot M0/M100 full-cascade responses for compiled-v1 cartridges")
    ap.add_argument("--manifest", type=Path, required=True, help="manifest.json with entries[].file")
    ap.add_argument("--source-dir", type=Path, required=True, help="directory holding cartridge json files")
    ap.add_argument("--out", type=Path, required=True, help="output PNG")
    ap.add_argument("--q", type=float, default=0.5, help="Q slice in [0,1]")
    args = ap.parse_args()

    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    entries = list(manifest.get("entries", []))
    if not entries:
        raise SystemExit("manifest has no entries")
    q = max(0.0, min(1.0, float(args.q)))

    cols = 2
    rows = math.ceil(len(entries) / cols)
    fig, axes = plt.subplots(rows, cols, figsize=(cols * 5.4, rows * 4.0), dpi=140, facecolor="#0b1020")
    axes_arr = np.atleast_1d(axes).reshape(rows, cols)

    for ax, entry in zip(axes_arr.flat, entries):
        cart = json.loads((args.source_dir / str(entry["file"])).read_text(encoding="utf-8"))
        f0, db0 = _curve(cart, "M0", q)
        f1, db1 = _curve(cart, "M100", q)
        ax.plot(f0, db0, color="#7dd3fc", linewidth=1.6, label="M0")
        ax.plot(f1, db1, color="#86efac", linewidth=1.6, label="M100")
        peak = max(float(np.max(db0)), float(np.max(db1)))
        floor = min(float(np.min(db0)), float(np.min(db1)))
        top = min(peak + 4.0, 90.0)
        bottom = max(floor - 6.0, -90.0)
        ax.set_xscale("log")
        ax.set_xlim(40.0, 20000.0)
        ax.set_ylim(bottom, top)
        ax.set_facecolor("#0b1020")
        ax.grid(True, which="both", color="#1e293b", linewidth=0.4, alpha=0.7)
        ax.tick_params(colors="#64748b", labelsize=7)
        for sp in ax.spines.values():
            sp.set_color("#1e293b")
        ax.set_title(str(entry["file"]).replace(".json", ""), color="#cbd5f5", fontsize=10, pad=4)
        ax.set_xlabel("Hz", color="#64748b", fontsize=8)
        ax.set_ylabel("dB", color="#64748b", fontsize=8)
        ax.legend(loc="upper left", fontsize=7, facecolor="#0b1020", edgecolor="#1e293b", labelcolor="#cbd5f5")

    for ax in list(axes_arr.flat)[len(entries):]:
        ax.axis("off")

    fig.suptitle(f"{args.manifest.parent.name} — full cascade endpoints @ Q={q:.2f}", color="#cbd5f5", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, facecolor="#0b1020", dpi=140)
    print(args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
