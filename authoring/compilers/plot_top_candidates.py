#!/usr/bin/env python3
"""Plot top generated candidates as 4-corner overlay contact sheets."""

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


SR = 44100
IR_LEN = 4096
FREQ_MIN = 40.0
FREQ_MAX = 20000.0
DB_MIN = -72.0
DB_MAX = 18.0

BG = "#0b1020"
GRID = "#1e293b"
TICK = "#64748b"
TITLE = "#cbd5f5"

LAYOUT = (
    ("M0_Q0", "#7dd3fc"),
    ("M100_Q0", "#fbbf24"),
    ("M0_Q100", "#f472b6"),
    ("M100_Q100", "#86efac"),
)


def _load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _curve_from_keyframe(kf: dict) -> tuple[np.ndarray, np.ndarray]:
    stages = [
        (
            float(s["c0"]),
            float(s["c1"]),
            float(s["c2"]),
            float(s["c3"]),
            float(s["c4"]),
        )
        for s in kf["stages"]
    ]
    ir = run_impulse(stages, IR_LEN, float(kf.get("boost", 1.0)))
    return magnitude_spectrum_db(ir, SR)


def _style(ax, title: str) -> None:
    ax.set_xscale("log")
    ax.set_xlim(FREQ_MIN, FREQ_MAX)
    ax.set_ylim(DB_MIN, DB_MAX)
    ax.axhline(0, color="#334155", linewidth=0.6, alpha=0.7)
    ax.set_facecolor(BG)
    ax.grid(True, which="both", color=GRID, linewidth=0.4, alpha=0.7)
    ax.tick_params(colors=TICK, labelsize=7)
    for sp in ax.spines.values():
        sp.set_color(GRID)
    ax.set_title(title, color=TITLE, fontsize=9, pad=3)
    ax.set_xlabel("Hz", color=TICK, fontsize=8)
    ax.set_ylabel("dB", color=TICK, fontsize=8)


def _top_entries(entries: list[dict], limit: int, intensity: float | None) -> list[dict]:
    filtered = entries
    if intensity is not None:
        filtered = [e for e in entries if abs(float(e["intensity"]) - intensity) < 1e-9]

    filtered.sort(
        key=lambda e: (
            -float(e.get("selection_score", 0.0)),
            str(e.get("category", "")),
            str(e.get("target", "")),
            str(e.get("recipe", "")),
            str(e.get("variant", "")),
            str(e.get("file", "")),
        )
    )
    return filtered[:limit]


def main() -> int:
    ap = argparse.ArgumentParser(description="Plot top generated candidates from a manifest")
    ap.add_argument("--manifest", type=Path, required=True, help="manifest.json or manifest.shortlist.json")
    ap.add_argument("--source-dir", type=Path, required=True, help="directory containing the cartridge JSON files")
    ap.add_argument("--out", type=Path, required=True, help="output PNG path")
    ap.add_argument("--top", type=int, default=12, help="number of candidates to plot")
    ap.add_argument("--intensity", type=float, default=None, help="optional exact intensity filter")
    args = ap.parse_args()

    manifest = _load_json(args.manifest)
    entries = list(manifest.get("entries", []))
    if not entries:
        raise SystemExit("manifest has no entries")
    if args.top <= 0:
        raise SystemExit("--top must be > 0")

    chosen = _top_entries(entries, args.top, args.intensity)
    if not chosen:
        raise SystemExit("no entries matched requested filters")

    cols = 3
    rows = math.ceil(len(chosen) / cols)
    fig, axes = plt.subplots(rows, cols, figsize=(cols * 5.0, rows * 3.7), dpi=130, facecolor=BG)
    axes_arr = np.atleast_1d(axes).reshape(rows, cols)

    for ax, entry in zip(axes_arr.flat, chosen):
        cart = _load_json(args.source_dir / str(entry["file"]))
        keyframes = {kf["label"]: kf for kf in cart["keyframes"]}
        for label, color in LAYOUT:
            freqs, db = _curve_from_keyframe(keyframes[label])
            mask = freqs >= FREQ_MIN
            ax.plot(freqs[mask], db[mask], color=color, linewidth=1.25, alpha=0.95)
        title = (
            f"{entry['recipe']}:{entry.get('variant', 'core')} s{entry['skin']} "
            f"i{int(round(float(entry['intensity']) * 100)):03d}\n"
            f"score {float(entry.get('selection_score', 0.0)):.3f}"
        )
        _style(ax, title)

    for ax in list(axes_arr.flat)[len(chosen):]:
        ax.axis("off")

    handles = [plt.Line2D([0], [0], color=color, lw=2) for _, color in LAYOUT]
    labels = [label for label, _ in LAYOUT]
    fig.legend(
        handles,
        labels,
        loc="upper center",
        ncol=4,
        fontsize=8,
        facecolor=BG,
        edgecolor=GRID,
        labelcolor=TITLE,
    )
    fig.suptitle("Top Candidate Overlays", color=TITLE, fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, facecolor=BG, dpi=130)
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
