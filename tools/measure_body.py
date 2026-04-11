"""Measure a body against the Speaker Knockerz shipping gates.

Reports: body_profile, peak/trough frequencies at key morph points.
Usage:  PYTHONPATH=. python tools/measure_body.py cartridges/p2k/P2k_006.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pyruntime.analysis import body_profile
from pyruntime.body import Body
from pyruntime.constants import SR
from pyruntime.freq_response import cascade_response_db, freq_points


def peak_trough_report(body: Body, morph: float, q: float) -> dict:
    freqs = freq_points(n=512, sr=SR)
    enc = body.corners.interpolate(morph, q)
    db = cascade_response_db(enc, freqs, SR)

    peak_idx = int(np.argmax(db))
    trough_idx = int(np.argmin(db))

    return {
        "morph": morph,
        "q": q,
        "peak_hz": round(float(freqs[peak_idx]), 1),
        "peak_db": round(float(db[peak_idx]), 1),
        "trough_hz": round(float(freqs[trough_idx]), 1),
        "trough_db": round(float(db[trough_idx]), 1),
        "range_db": round(float(db[peak_idx] - db[trough_idx]), 1),
    }


def pole_freqs(body: Body, corner_label: str) -> list[dict]:
    from pyruntime.corner import CornerName
    cmap = {"M0_Q0": CornerName.A, "M0_Q100": CornerName.B,
            "M100_Q0": CornerName.C, "M100_Q100": CornerName.D}
    cn = cmap[corner_label]
    cs = body.corners.corner(cn)
    result = []
    for i, sp in enumerate(cs.stages):
        if sp.r < 0.01:
            continue
        result.append({"stage": i, "freq_hz": round(sp.freq_hz(), 1), "r": round(sp.r, 4)})
    return result


def main() -> None:
    if len(sys.argv) < 2:
        print("Usage: measure_body.py <body.json> [<body2.json>]")
        sys.exit(1)

    paths = sys.argv[1:]
    for p in paths:
        body = Body.from_json(p)
        print(f"\n{'='*60}")
        print(f"  {body.name}  ({p})")
        print(f"{'='*60}")

        # body_profile
        prof = body_profile(body)
        md = prof["morph_distance"]
        print(f"\n  morph_distance:   q0={md['q0_rms_db']:.1f}  q1={md['q1_rms_db']:.1f}  mean={md['mean_rms_db']:.1f}")
        cx = prof["crossings"]
        print(f"  crossings:        m0q0={cx['m0_q0']}  m0q100={cx['m0_q100']}  m100q0={cx['m100_q0']}  m100q100={cx['m100_q100']}")
        print(f"  spectral_tilt_db: {prof['spectral_tilt_db']:.1f}")
        print(f"  active_stages:    {prof['active_stages']}")
        audit = prof["midpoint_audit"]
        status = "PASS" if audit["passed"] else "FAIL"
        print(f"  midpoint_audit:   {status}  worst={audit['worst_peak_db']:.1f} dB @ morph={audit['worst_point']['morph']} q={audit['worst_point']['q']}")
        for pt in audit["points"]:
            print(f"    m={pt['morph']:.2f} q={pt['q']:.1f} peak={pt['peak_db']:.1f} dB")

        # Peak / trough at key positions
        print(f"\n  Peak / trough:")
        for m, q in [(0, 0), (0, 1), (0.5, 0), (0.5, 0.5), (0.5, 1), (0.75, 0), (1, 0), (1, 1)]:
            r = peak_trough_report(body, m, q)
            print(f"    M{m:.2f}_Q{q:.1f}:  peak={r['peak_hz']:>8.1f} Hz / {r['peak_db']:>6.1f} dB   "
                  f"trough={r['trough_hz']:>8.1f} Hz / {r['trough_db']:>7.1f} dB   range={r['range_db']:.1f}")

        # Pole skeleton at key corners
        print(f"\n  Pole skeleton:")
        for label in ["M0_Q0", "M0_Q100", "M100_Q0", "M100_Q100"]:
            poles = pole_freqs(body, label)
            freqs_str = ", ".join(f"S{p['stage']}:{p['freq_hz']}Hz(r={p['r']})" for p in poles)
            print(f"    {label}: {freqs_str}")


if __name__ == "__main__":
    main()
