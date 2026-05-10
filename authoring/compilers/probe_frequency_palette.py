"""Probe reachable peak/notch frequencies across MorphLPX + Morph2 params.

For each (filter_type, skin, q_byte) combo, compile a cartridge, render the
M100_Q100 corner impulse response, find the peak location and the minimum
(notch) location in 80 Hz–20 kHz. Build a map so we can pick candidates
targeting specific Aluminum Siding frequency anchors (1 kHz scoop, 12 kHz
ring, 18 kHz peak).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

_TOOLS = Path(__file__).resolve().parent
sys.path.insert(0, str(_TOOLS))

from compile_morphlp import compile_morphlp_filter  # noqa: E402
from compile_morph2 import compile_morph2_filter    # noqa: E402
from cube_authoring.thumbnail import (               # noqa: E402
    magnitude_spectrum_db,
    run_impulse,
)

SR = 44100
IR_LEN = 8192
FMIN, FMAX = 80.0, 20000.0


def corner_curve(kf: dict) -> tuple[np.ndarray, np.ndarray]:
    stages = [(float(s["c0"]), float(s["c1"]), float(s["c2"]),
               float(s["c3"]), float(s["c4"])) for s in kf["stages"]]
    ir = run_impulse(stages, IR_LEN, float(kf.get("boost", 1.0)))
    return magnitude_spectrum_db(ir, SR)


def feature_summary(kf: dict) -> dict:
    f, db = corner_curve(kf)
    mask = (f >= FMIN) & (f <= FMAX)
    f, db = f[mask], db[mask]
    peak_i = int(np.argmax(db))
    notch_i = int(np.argmin(db))
    # Kilohertz band averages
    def band_db(lo, hi):
        m = (f >= lo) & (f <= hi)
        return float(np.mean(db[m])) if m.any() else float("nan")
    return {
        "peak_hz":   float(f[peak_i]),
        "peak_db":   float(db[peak_i]),
        "notch_hz":  float(f[notch_i]),
        "notch_db":  float(db[notch_i]),
        "db_1k":     band_db(800, 1200),
        "db_5k":     band_db(4500, 5500),
        "db_10k":    band_db(9000, 11000),
        "db_12k":    band_db(11000, 13000),
        "db_18k":    band_db(16000, 19500),
    }


def probe_morphlpx():
    print("=== MorphLPX M100_Q100 feature map ===")
    print(f"{'skin':>4} {'q0':>3} {'q100':>4} | "
          f"{'peak_hz':>8} {'peak_dB':>7} | {'1k':>7} {'5k':>7} {'10k':>7} {'12k':>7} {'18k':>7}")
    for skin in range(4):
        for q100 in (16, 32, 63, 96, 127, 160, 200):
            cart = compile_morphlp_filter("morphlpx", skin, q0_byte=0, q100_byte=q100)
            kf = next(k for k in cart["keyframes"] if k["label"] == "M100_Q100")
            s = feature_summary(kf)
            print(f"{skin:>4} {0:>3} {q100:>4} | "
                  f"{s['peak_hz']:>8.0f} {s['peak_db']:>7.1f} | "
                  f"{s['db_1k']:>7.1f} {s['db_5k']:>7.1f} {s['db_10k']:>7.1f} "
                  f"{s['db_12k']:>7.1f} {s['db_18k']:>7.1f}")
        print()


def probe_morph2():
    print("=== Morph2 M100_Q0 feature map (Q axis unused, morph=1) ===")
    print(f"{'skin':>4} {'cq':>3} {'qs':>3} {'pit':>4} | "
          f"{'peak_hz':>8} {'peak_dB':>7} | {'1k':>7} {'5k':>7} {'10k':>7} {'12k':>7} {'18k':>7}")
    for skin in range(4):
        for cq in (64, 128, 192, 224):
            for qs in (16, 48, 96):
                cart = compile_morph2_filter(
                    skin, m0_center_q=64, m0_q_spread=32, m0_pitch=0x40,
                    m100_center_q=cq, m100_q_spread=qs, m100_pitch=0x40,
                )
                kf = next(k for k in cart["keyframes"] if k["label"] == "M100_Q0")
                s = feature_summary(kf)
                print(f"{skin:>4} {cq:>3} {qs:>3} {0x40:>4} | "
                      f"{s['peak_hz']:>8.0f} {s['peak_db']:>7.1f} | "
                      f"{s['db_1k']:>7.1f} {s['db_5k']:>7.1f} {s['db_10k']:>7.1f} "
                      f"{s['db_12k']:>7.1f} {s['db_18k']:>7.1f}")
        print()


if __name__ == "__main__":
    probe_morphlpx()
    probe_morph2()
