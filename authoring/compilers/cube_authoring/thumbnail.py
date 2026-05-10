"""Spectrum thumbnailer — render a compiled cube corner to a small PNG.

Runs a unit impulse through the 6-stage DF2T cascade (mirrors the shipping
runtime: trench-core/src/cascade.rs and tools/bake_hedz_const.py), FFTs the
output, and saves a dark log-frequency magnitude-response thumbnail.

Used by the Cube Net UI as the visual currency per corner. Also usable as a
CLI to dump a full cube or a single corner.
"""
from __future__ import annotations

import argparse
import json
import struct
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

CORNER_LABELS = ("c000", "c100", "c010", "c110", "c001", "c101", "c011", "c111")
SR = 44100
IMPULSE_LEN = 512
THUMB_WIDTH_PX = 256
THUMB_HEIGHT_PX = 128
DPI = 80
FREQ_MIN_HZ = 40.0
FREQ_MAX_HZ = 22000.0
DB_MIN = -60.0
DB_MAX = 40.0


# -----------------------------------------------------------------------------
# DF2T cascade — matches trench-core/src/cascade.rs:38-68 and
# tools/bake_hedz_const.py:run_impulse.
# -----------------------------------------------------------------------------

def f32(x: float) -> float:
    return struct.unpack("<f", struct.pack("<f", x))[0]


def run_impulse(stages: list[tuple[float, float, float, float, float]],
                 length: int, boost: float) -> np.ndarray:
    state = [(0.0, 0.0) for _ in range(len(stages))]
    out = np.zeros(length, dtype=np.float32)
    for n in range(length):
        sig = 1.0 if n == 0 else 0.0
        for i, (c0, c1, c2, c3, c4) in enumerate(stages):
            w1, w2 = state[i]
            y = c0 * sig + w1
            new_w1 = c1 * sig - c3 * y + w2
            new_w2 = c2 * sig - c4 * y
            state[i] = (new_w1, new_w2)
            sig = y
        out[n] = f32(sig * boost)
    return out


def corner_stage_tuples(corner: dict) -> list[tuple[float, float, float, float, float]]:
    out = []
    for s in corner["stages"]:
        out.append((float(s["c0"]), float(s["c1"]), float(s["c2"]),
                    float(s["c3"]), float(s["c4"])))
    return out


def magnitude_spectrum_db(impulse_response: np.ndarray, sr: int
                           ) -> tuple[np.ndarray, np.ndarray]:
    n = impulse_response.size
    spec = np.fft.rfft(impulse_response)
    mag = np.abs(spec)
    # avoid log(0)
    mag = np.maximum(mag, 1e-12)
    db = 20.0 * np.log10(mag)
    freqs = np.fft.rfftfreq(n, d=1.0 / sr)
    return freqs, db


# -----------------------------------------------------------------------------
# Render
# -----------------------------------------------------------------------------

def plot_response(ax, freqs: np.ndarray, db: np.ndarray, *, title: str | None = None):
    mask = freqs >= FREQ_MIN_HZ
    ax.plot(freqs[mask], db[mask], color="#7dd3fc", linewidth=1.4)
    ax.set_xscale("log")
    ax.set_xlim(FREQ_MIN_HZ, FREQ_MAX_HZ)
    ax.set_ylim(DB_MIN, DB_MAX)
    ax.axhline(0, color="#334155", linewidth=0.6, alpha=0.5)
    ax.set_facecolor("#0b1020")
    ax.grid(True, which="both", color="#1e293b", linewidth=0.4, alpha=0.6)
    ax.tick_params(colors="#64748b", labelsize=6)
    for spine in ax.spines.values():
        spine.set_color("#1e293b")
    if title:
        ax.set_title(title, color="#cbd5f5", fontsize=8, pad=2)


def render_corner_thumbnail(corner: dict, out_path: Path,
                              *, sr: int = SR, title: str | None = None) -> None:
    stages = corner_stage_tuples(corner)
    boost = float(corner.get("boost", 1.0))
    ir = run_impulse(stages, IMPULSE_LEN, boost)
    freqs, db = magnitude_spectrum_db(ir, sr)

    fig = plt.figure(figsize=(THUMB_WIDTH_PX / DPI, THUMB_HEIGHT_PX / DPI),
                      dpi=DPI, facecolor="#0b1020")
    ax = fig.add_axes([0.02, 0.05, 0.96, 0.90])
    plot_response(ax, freqs, db, title=title)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=DPI, facecolor=fig.get_facecolor())
    plt.close(fig)


def render_cube_composite(compiled: dict, out_path: Path, *, sr: int = SR) -> None:
    """2x4 grid — row 1: cluster corners (c000/c100/c010/c110),
                     row 2: distributed corners (c001/c101/c011/c111)."""
    rows = (("c000", "c100", "c010", "c110"),
            ("c001", "c101", "c011", "c111"))
    fig = plt.figure(figsize=(12, 4), dpi=DPI, facecolor="#0b1020")
    for r, row in enumerate(rows):
        for c, label in enumerate(row):
            ax = fig.add_subplot(2, 4, r * 4 + c + 1)
            corner = compiled["corners"][label]
            stages = corner_stage_tuples(corner)
            boost = float(corner.get("boost", 1.0))
            ir = run_impulse(stages, IMPULSE_LEN, boost)
            freqs, db = magnitude_spectrum_db(ir, sr)
            plot_response(ax, freqs, db, title=label)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=DPI, facecolor=fig.get_facecolor())
    plt.close(fig)


def thumbnail_cube(compiled: dict, out_dir: Path, *, sr: int = SR) -> list[Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = []
    for label in CORNER_LABELS:
        p = out_dir / f"{label}.png"
        render_corner_thumbnail(compiled["corners"][label], p, sr=sr, title=label)
        paths.append(p)
    return paths


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------

def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(
        description="Render spectrum thumbnails for a compiled cube surface.",
    )
    parser.add_argument("compiled", help="Path to a compiled cube surface JSON "
                                             "(trench.compiled.cube_surface.v1).")
    parser.add_argument("--out-dir", default=None,
                         help="Directory to write per-corner thumbnails "
                              "(default: alongside the compiled file, in a "
                              "'thumbs_<stem>/' subdir).")
    parser.add_argument("--composite", action="store_true",
                         help="Also emit a single 2x4 composite PNG.")
    parser.add_argument("--sr", type=int, default=SR,
                         help=f"Sample rate used for the impulse render (default: {SR}).")
    args = parser.parse_args(argv)

    compiled_path = Path(args.compiled)
    with compiled_path.open() as f:
        compiled = json.load(f)
    if compiled.get("schema") != "trench.compiled.cube_surface.v1":
        print(f"error: expected trench.compiled.cube_surface.v1, got "
              f"{compiled.get('schema')!r}", file=sys.stderr)
        return 2

    out_dir = Path(args.out_dir) if args.out_dir else (
        compiled_path.parent / f"thumbs_{compiled_path.stem}"
    )
    paths = thumbnail_cube(compiled, out_dir, sr=args.sr)
    print(f"wrote {len(paths)} corner thumbnails to {out_dir}")

    if args.composite:
        comp_path = out_dir / "_composite.png"
        render_cube_composite(compiled, comp_path, sr=args.sr)
        print(f"wrote composite {comp_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
