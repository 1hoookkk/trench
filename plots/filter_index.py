"""Generate plots/p2k_filter_index.json — LLM-readable summary of every
P2K filter's frequency response at all 4 corners.

Computes for each filter × corner:
  - shape classification (lowpass / highpass / bandpass / bandstop /
    peaking / notched / complex)
  - global peak (hz, dB) and global trough (hz, dB)
  - formants: local peaks with prominence >= 3 dB, sorted by frequency
  - notches: local troughs with depth >= 6 dB, sorted by frequency
  - anchor readings at 100 Hz, 1 kHz, 8 kHz
  - energy-weighted spectral centroid

Stdlib only. Reads filter data from git history (same source as
png_plotter.py — all 33 P2k JSON files deleted in f146fb3).
"""
from __future__ import annotations

import json
import math
import subprocess
from pathlib import Path

SR = 39062.5


# ---------------------------------------------------------------------------
# DSP (shared with png_plotter.py)
# ---------------------------------------------------------------------------

def _biquad_mag(c, omega):
    cos1, sin1 = math.cos(omega), math.sin(omega)
    cos2, sin2 = math.cos(2 * omega), math.sin(2 * omega)
    num_re = c[0] + c[1] * cos1 + c[2] * cos2
    num_im = -c[1] * sin1 - c[2] * sin2
    den_re = 1.0 + c[3] * cos1 + c[4] * cos2
    den_im = -c[3] * sin1 - c[4] * sin2
    num = math.sqrt(num_re ** 2 + num_im ** 2)
    den = math.sqrt(den_re ** 2 + den_im ** 2)
    return num / den if den > 1e-18 else 1e6


def _cascade_db(stages, boost, f):
    omega = 2 * math.pi * f / SR
    m = boost
    for s in stages:
        m *= _biquad_mag(s, omega)
    return 20 * math.log10(m) if m > 1e-10 else -200.0


# ---------------------------------------------------------------------------
# Analysis
# ---------------------------------------------------------------------------

F_MIN, F_MAX = 50.0, 12000.0
N = 2048  # enough resolution for peak/notch detection; not a visual render


def _sample(stages, boost):
    """Return (freqs, dbs) arrays of length N."""
    log_lo = math.log10(F_MIN)
    log_hi = math.log10(F_MAX)
    freqs = [10 ** (log_lo + (log_hi - log_lo) * i / (N - 1)) for i in range(N)]
    dbs = [_cascade_db(stages, boost, f) for f in freqs]
    return freqs, dbs


def _prominence(dbs, i):
    """Topographic prominence of a local maximum at index i."""
    # Left key col: minimum between i and the nearest higher point to the left
    left_min = dbs[i]
    for j in range(i - 1, -1, -1):
        if dbs[j] > dbs[i]:
            break
        if dbs[j] < left_min:
            left_min = dbs[j]
    # Right key col
    right_min = dbs[i]
    for j in range(i + 1, len(dbs)):
        if dbs[j] > dbs[i]:
            break
        if dbs[j] < right_min:
            right_min = dbs[j]
    return dbs[i] - max(left_min, right_min)


def _depth(dbs, i):
    """Topographic depth of a local minimum at index i (inverse prominence)."""
    left_max = dbs[i]
    for j in range(i - 1, -1, -1):
        if dbs[j] < dbs[i]:
            break
        if dbs[j] > left_max:
            left_max = dbs[j]
    right_max = dbs[i]
    for j in range(i + 1, len(dbs)):
        if dbs[j] < dbs[i]:
            break
        if dbs[j] > right_max:
            right_max = dbs[j]
    return min(left_max, right_max) - dbs[i]


def analyze_corner(stages, boost):
    freqs, dbs = _sample(stages, boost)

    # --- global peak / trough (unclamped) ---
    peak_i = max(range(N), key=lambda i: dbs[i])
    trough_i = min(range(N), key=lambda i: dbs[i])

    # --- formants: local maxima with prominence >= 3 dB, top 8 by prominence ---
    formants = []
    for i in range(1, N - 1):
        if dbs[i] > dbs[i - 1] and dbs[i] > dbs[i + 1]:
            p = _prominence(dbs, i)
            if p >= 3.0:
                formants.append({"hz": round(freqs[i]), "db": round(dbs[i], 1),
                                  "prominence_db": round(p, 1)})
    formants.sort(key=lambda x: x["prominence_db"], reverse=True)
    formants = formants[:8]
    formants.sort(key=lambda x: x["hz"])

    # --- notches: local minima with depth >= 6 dB, top 4 by depth ---
    notches = []
    for i in range(1, N - 1):
        if dbs[i] < dbs[i - 1] and dbs[i] < dbs[i + 1]:
            d = _depth(dbs, i)
            if d >= 6.0:
                notches.append({"hz": round(freqs[i]), "db": round(dbs[i], 1),
                                 "depth_db": round(d, 1)})
    notches.sort(key=lambda x: x["depth_db"], reverse=True)
    notches = notches[:4]
    notches.sort(key=lambda x: x["hz"])

    # --- anchor readings ---
    d100  = _cascade_db(stages, boost, 100.0)
    d1k   = _cascade_db(stages, boost, 1000.0)
    d8k   = _cascade_db(stages, boost, 8000.0)

    # --- shape classification ---
    THRESH = 6.0
    if d100 > d8k + THRESH and d100 > d1k + THRESH:
        shape = "lowpass"
    elif d8k > d100 + THRESH and d8k > d1k + THRESH:
        shape = "highpass"
    elif d1k > d100 + THRESH and d1k > d8k + THRESH:
        shape = "bandpass"
    elif d100 > d1k + THRESH and d8k > d1k + THRESH:
        shape = "bandstop"
    elif dbs[peak_i] > 6.0 and len(formants) <= 2:
        shape = "peaking"
    elif notches and dbs[trough_i] < -12.0:
        shape = "notched"
    else:
        shape = "complex"

    # --- spectral centroid (energy-weighted) ---
    total_w = 0.0
    weighted_f = 0.0
    for i in range(N):
        lin = 10 ** (dbs[i] / 20.0)
        e = lin * lin
        weighted_f += freqs[i] * e
        total_w += e
    centroid_hz = round(weighted_f / total_w) if total_w > 0 else 0

    return {
        "shape": shape,
        "peak_hz": round(freqs[peak_i]),
        "peak_db": round(dbs[peak_i], 1),
        "trough_hz": round(freqs[trough_i]),
        "trough_db": round(dbs[trough_i], 1),
        "centroid_hz": centroid_hz,
        "db_at_100hz": round(d100, 1),
        "db_at_1khz": round(d1k, 1),
        "db_at_8khz": round(d8k, 1),
        "formants": formants,
        "notches": notches,
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    repo = Path(__file__).resolve().parents[1]
    corners = ["M0_Q0", "M100_Q0", "M0_Q100", "M100_Q100"]

    all_filters = []
    for n in range(33):
        num = f"{n:03d}"
        result = subprocess.run(
            ["git", "-C", str(repo), "show",
             f"f146fb3^:cartridges/p2k/P2k_{num}.json"],
            capture_output=True, text=True, check=True,
        )
        p2k = json.loads(result.stdout)

        entry = {"id": f"p2k_{num}", "corners": {}}
        for corner in corners:
            kf = next(k for k in p2k["keyframes"] if k["label"] == corner)
            stages = [[s["c0"], s["c1"], s["c2"], s["c3"], s["c4"]]
                      for s in kf["stages"]]
            boost = kf.get("boost", 1.0)
            entry["corners"][corner] = analyze_corner(stages, boost)

        all_filters.append(entry)
        print(f"indexed p2k_{num}")

    index = {
        "schema_version": "1",
        "description": (
            "LLM-readable frequency-response index for all 33 E-mu Proteus 2000 "
            "filters (P2k_000–P2k_032). Each filter has 4 corners: "
            "M=Morph (0=start, 100=end), Q=resonance (0=flat, 100=max). "
            "Frequencies in Hz, levels in dBFS. "
            "formants: local peaks with topographic prominence >= 3 dB. "
            "notches: local troughs with depth >= 6 dB. "
            "shape: lowpass|highpass|bandpass|bandstop|peaking|notched|complex. "
            "centroid_hz: energy-weighted spectral center of mass. "
            "Source: design-time 6-stage biquad model at SR=39062.5 Hz "
            "(not the runtime 5-stage plugin model)."
        ),
        "sample_rate_hz": SR,
        "freq_range_hz": [F_MIN, F_MAX],
        "filters": all_filters,
    }

    out_path = repo / "plots" / "p2k_filter_index.json"
    out_path.write_text(json.dumps(index, indent=2))
    size = out_path.stat().st_size
    print(f"wrote {out_path.relative_to(repo)}  ({size:,} bytes, {len(all_filters)} filters)")


if __name__ == "__main__":
    main()
