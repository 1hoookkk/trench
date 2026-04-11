"""Profile a single P2K skin — extract pole/zero anatomy for all 4 corners.

Usage:
    python tools/profile_p2k_skin.py --skin P2k_013
    python tools/profile_p2k_skin.py --skin P2k_026 --out vault/_profiles

Outputs a structured JSON profile with per-stage pole/zero frequencies,
radii, and corner deltas. No classification — that's for the reader.
"""
from __future__ import annotations
import argparse, json, math, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from pyruntime.constants import SR, TWO_PI
from pyruntime.body import Body
from pyruntime.freq_response import cascade_response_db, freq_points

FREQS = freq_points(sr=SR)
CORNER_LABELS = ["M0_Q0", "M0_Q100", "M100_Q0", "M100_Q100"]
CORNER_COORDS = [(0.0, 0.0), (0.0, 1.0), (1.0, 0.0), (1.0, 1.0)]


def pole_zero_from_kernel(c0, c1, c2, c3, c4):
    """Extract pole and zero freq/radius from kernel coefficients."""
    r = math.sqrt(abs(c4))
    if r > 0.001:
        cos_t = max(-1, min(1, -c3 / (2 * r)))
        pole_hz = math.acos(cos_t) * SR / TWO_PI
    else:
        pole_hz = 0.0

    if abs(c0) > 0.001:
        zr = math.sqrt(abs(c2 / c0))
        if zr > 0.001:
            zcos = max(-1, min(1, -(c1 / c0) / (2 * zr)))
            zero_hz = math.acos(zcos) * SR / TWO_PI
        else:
            zero_hz = 0.0
    else:
        zero_hz = 0.0
        zr = 0.0

    return {
        "pole_hz": round(pole_hz, 1),
        "pole_r": round(r, 6),
        "zero_hz": round(zero_hz, 1),
        "zero_r": round(zr, 6),
    }


def response_stats(body: Body, morph: float, q: float):
    """Get peak/notch dB and spectral centroid at a morph/q point."""
    coeffs = body.corners.interpolate(morph, q)
    db = cascade_response_db(coeffs, FREQS, SR)
    peak_db = float(db.max())
    notch_db = float(db.min())
    dyn_range = peak_db - notch_db

    # spectral centroid (linear weighting)
    import numpy as np
    linear = 10 ** (db / 20)
    total = linear.sum()
    if total > 0:
        centroid_hz = float((FREQS * linear).sum() / total)
    else:
        centroid_hz = 0.0

    # count peaks and notches (local extrema > 3dB prominence)
    from scipy.signal import find_peaks
    peaks, _ = find_peaks(db, prominence=3)
    notches, _ = find_peaks(-db, prominence=3)

    return {
        "peak_db": round(peak_db, 1),
        "notch_db": round(notch_db, 1),
        "dynamic_range_db": round(dyn_range, 1),
        "centroid_hz": round(centroid_hz, 1),
        "num_peaks": len(peaks),
        "num_notches": len(notches),
    }


def profile_skin(skin_path: Path):
    body = Body.from_json(str(skin_path))
    result = {
        "name": body.name,
        "source": str(skin_path),
        "sample_rate": SR,
        "boost": body.boost,
        "corners": {},
    }

    for label, (m, q) in zip(CORNER_LABELS, CORNER_COORDS):
        coeffs = body.corners.interpolate(m, q)
        stages = []
        for i, enc in enumerate(coeffs[:6]):
            pz = pole_zero_from_kernel(enc.c0, enc.c1, enc.c2, enc.c3, enc.c4)
            pz["stage"] = i + 1
            stages.append(pz)

        stats = response_stats(body, m, q)
        result["corners"][label] = {"stages": stages, "response": stats}

    # compute motion deltas
    def _delta(a, b):
        return {
            k: round(b["response"][k] - a["response"][k], 2)
            for k in ["peak_db", "dynamic_range_db", "centroid_hz"]
        }

    c = result["corners"]
    result["motion"] = {
        "morph_at_Q0": _delta(c["M0_Q0"], c["M100_Q0"]),
        "morph_at_Q100": _delta(c["M0_Q100"], c["M100_Q100"]),
        "Q_at_M0": _delta(c["M0_Q0"], c["M0_Q100"]),
        "Q_at_M100": _delta(c["M100_Q0"], c["M100_Q100"]),
        "diagonal": _delta(c["M0_Q0"], c["M100_Q100"]),
    }

    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--skin", required=True, help="Skin name e.g. P2k_013")
    parser.add_argument("--out", default="vault/_profiles", help="Output dir")
    args = parser.parse_args()

    cartridge = Path("cartridges/p2k") / f"{args.skin}.json"
    if not cartridge.exists():
        print(f"ERROR: {cartridge} not found")
        sys.exit(1)

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    profile = profile_skin(cartridge)
    out_path = out_dir / f"{args.skin}_profile.json"
    with open(out_path, "w") as f:
        json.dump(profile, f, indent=2)
    print(f"Wrote {out_path}")

    # Print summary for the reader
    print(f"\n=== {profile['name']} ===")
    for label in CORNER_LABELS:
        corner = profile["corners"][label]
        r = corner["response"]
        print(f"\n  {label}: {r['num_peaks']}pk/{r['num_notches']}notch  "
              f"range={r['dynamic_range_db']}dB  centroid={r['centroid_hz']}Hz")
        for s in corner["stages"]:
            print(f"    Stage {s['stage']}: pole={s['pole_hz']:>7.1f}Hz r={s['pole_r']:.4f}"
                  f"  |  zero={s['zero_hz']:>7.1f}Hz zr={s['zero_r']:.4f}")

    print(f"\n  Motion:")
    for axis, d in profile["motion"].items():
        print(f"    {axis}: peak={d['peak_db']:+.1f}dB  "
              f"range={d['dynamic_range_db']:+.1f}dB  centroid={d['centroid_hz']:+.1f}Hz")


if __name__ == "__main__":
    main()
