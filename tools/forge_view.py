#!/usr/bin/env python3
"""Live filter-response viewer for the trench-forge-dev watch slot.

Renders what makes a preset *a preset*: the four corner magnitude responses,
peaks (resonances / formants), notches (carved zeros), per-corner boost,
and the active-stage count. Polls `forge-watch/active.json` on a 200 ms
cadence and redraws on every write — so the loop from `forge_write.py`
to ear-by-eye is sub-half-second.

Usage:
    python tools/forge_view.py

Env overrides:
    TRENCH_FORGE_WATCH   path to compiled-v1 JSON (default forge-watch/active.json)
    TRENCH_FORGE_POLL_MS poll interval in ms (default 200)
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

CORNERS = ["M0_Q0", "M0_Q100", "M100_Q0", "M100_Q100"]
NFFT = 4096
DB_FLOOR = 1e-12

# A stage is "passthrough" when c0=1, c1..c4=0. Compile output uses
# exact binary zeros for the inactive tail, so an exact compare is
# sufficient — no epsilon needed.
PASSTHROUGH = (1.0, 0.0, 0.0, 0.0, 0.0)


def watch_path() -> Path:
    return Path(os.environ.get("TRENCH_FORGE_WATCH", "forge-watch/active.json"))


def poll_ms() -> int:
    return int(os.environ.get("TRENCH_FORGE_POLL_MS", "200"))


def cascade_response(stages: list[dict], sr: float) -> tuple[np.ndarray, np.ndarray]:
    """Cascade the biquads (DF2T → rational transfer function)
    and return (freqs_hz, magnitude_db) from 0 to Nyquist."""
    w = np.linspace(0.0, np.pi, NFFT)
    z_inv = np.exp(-1j * w)
    h = np.ones_like(w, dtype=complex)
    for s in stages:
        c0, c1, c2, c3, c4 = s["c0"], s["c1"], s["c2"], s["c3"], s["c4"]
        if (c0, c1, c2, c3, c4) == PASSTHROUGH:
            continue
        num = c0 + c1 * z_inv + c2 * (z_inv * z_inv)
        den = 1.0 + c3 * z_inv + c4 * (z_inv * z_inv)
        h *= num / den
    freqs = w * sr / (2 * np.pi)
    mag_db = 20.0 * np.log10(np.maximum(np.abs(h), DB_FLOOR))
    return freqs, mag_db


def find_extrema(
    freqs: np.ndarray, mag_db: np.ndarray, min_hz: float = 40.0
) -> tuple[list[tuple[float, float]], list[tuple[float, float]]]:
    """Local maxima (peaks) and minima (notches) above min_hz, ranked
    by prominence against the cascade's median level."""
    mask = freqs >= min_hz
    idx = np.where(mask)[0]
    if idx.size < 3:
        return [], []
    start = idx[0]
    sub = mag_db[start:]
    f_sub = freqs[start:]
    diff = np.diff(sub)
    sign = np.sign(diff)
    # zero-crossings of the derivative
    maxima = []
    minima = []
    for i in range(1, len(sign)):
        if sign[i - 1] > 0 and sign[i] <= 0:
            maxima.append(i)
        elif sign[i - 1] < 0 and sign[i] >= 0:
            minima.append(i)
    median = float(np.median(sub))
    # prominence = distance from median; filter thin wobbles
    peaks = [(f_sub[i], sub[i]) for i in maxima if sub[i] - median > 3.0]
    notches = [(f_sub[i], sub[i]) for i in minima if median - sub[i] > 3.0]
    peaks.sort(key=lambda x: -x[1])
    notches.sort(key=lambda x: x[1])
    return peaks[:3], notches[:2]


def active_stage_count(stages: list[dict]) -> int:
    n = 0
    for s in stages:
        key = (s["c0"], s["c1"], s["c2"], s["c3"], s["c4"])
        if key != PASSTHROUGH:
            n += 1
    return n


def render_corner(ax, cart: dict, label: str, sr: float) -> None:
    ax.clear()
    keyframe = next((k for k in cart["keyframes"] if k["label"] == label), None)
    if keyframe is None:
        ax.set_title(f"{label}  (missing)")
        return
    stages = keyframe["stages"]
    boost = keyframe.get("boost", 1.0)
    active = active_stage_count(stages)
    freqs, mag = cascade_response(stages, sr)
    # skip DC for log-x
    ax.semilogx(freqs[1:], mag[1:], color="#d14", linewidth=1.6)
    peaks, notches = find_extrema(freqs, mag)
    for f, m in peaks:
        ax.plot(f, m, marker="^", color="#d14", markersize=7)
        ax.annotate(
            f"{int(round(f))} Hz\n{m:+.1f} dB",
            xy=(f, m),
            xytext=(6, 4),
            textcoords="offset points",
            fontsize=7,
            color="#d14",
        )
    for f, m in notches:
        ax.plot(f, m, marker="v", color="#147", markersize=6)
        ax.annotate(
            f"{int(round(f))} Hz\n{m:+.1f} dB",
            xy=(f, m),
            xytext=(6, -14),
            textcoords="offset points",
            fontsize=7,
            color="#147",
        )
    ax.set_title(f"{label}   boost={boost:.2f}x   active={active}/6", fontsize=10)
    ax.set_xlim(20.0, sr / 2.0)
    ax.set_ylim(-60.0, 40.0)
    ax.set_xlabel("Hz")
    ax.set_ylabel("dB")
    ax.grid(True, which="both", alpha=0.25)


def render_error(fig, axes, message: str) -> None:
    for ax in axes:
        ax.clear()
        ax.axis("off")
    fig.suptitle(message, color="#b00", fontsize=11)


def main() -> int:
    path = watch_path()
    interval = poll_ms() / 1000.0

    plt.rcParams["figure.constrained_layout.use"] = True
    fig, axes = plt.subplots(2, 2, figsize=(13, 8))
    flat = list(axes.flatten())
    fig.canvas.manager.set_window_title("TRENCH forge-view")
    plt.ion()
    plt.show(block=False)

    last_mtime = -1.0
    while plt.fignum_exists(fig.number):
        try:
            mtime = path.stat().st_mtime if path.exists() else -1.0
        except OSError:
            mtime = -1.0

        if mtime != last_mtime:
            last_mtime = mtime
            if mtime < 0:
                render_error(fig, flat, f"waiting for {path}")
            else:
                try:
                    cart = json.loads(path.read_text())
                    name = cart.get("name", "?")
                    fmt = cart.get("format", "?")
                    sr = float(cart.get("sampleRate", 48000.0))
                    fig.suptitle(
                        f"{name}   ·   {fmt}   ·   {sr:.1f} Hz",
                        fontsize=12,
                    )
                    for ax, label in zip(flat, CORNERS):
                        render_corner(ax, cart, label, sr)
                except (OSError, json.JSONDecodeError, KeyError, ValueError) as exc:
                    render_error(fig, flat, f"{type(exc).__name__}: {exc}")
            fig.canvas.draw_idle()

        plt.pause(interval)

    return 0


if __name__ == "__main__":
    sys.exit(main())
