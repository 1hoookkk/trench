#!/usr/bin/env python3
"""DF-II direct-coefficient demo authoring.

Writes a compiled-v1 cartridge from a single RBJ-cookbook biquad. The DF-II
runtime cascade in `runtime/trench-core/src/cascade.rs` accepts
(c0, c1, c2, c3, c4) per stage directly, where:

    c0 = b0/a0   c1 = b1/a0   c2 = b2/a0   c3 = a1/a0   c4 = a2/a0

so any standard biquad designed via the RBJ cookbook drops in cleanly.

Scope and doctrine
------------------
This is the **demo / content path** for showing what DF-II does on camera.
DF-II is a Musical Filter. The faceplate reads "Musical Filter". This script
demonstrates that the runtime accepts any biquad you hand it, by going from
a one-line description ("warm 200 Hz lowpass with a little resonance") to a
compiled-v1 cartridge in one CLI call.

Per `authoring/FILTER_WORK.md`, shipping body cartridges (Speaker Knockerz,
Aluminum Siding, Small Talk Ah-Ee, Cul-De-Sac) require **body-first**
authoring against the DF-II cascade — that compiler is not yet built. The
project doctrine bans the RBJ cookbook for character filters; this script
intentionally stays clear of that ban by emitting **generic** filters
(lp/hp/peak/shelf/bp/notch), not body cartridges. Output goes to
`cartridges/factory/my_designs/`, not into the shipping body tree.

Usage
-----
    python authoring/tools/design_df2.py NAME TYPE FREQ Q [GAIN_DB]
        [--sr SR] [--out PATH]

Filter types:
    lp         lowpass         args: FREQ, Q
    hp         highpass        args: FREQ, Q
    peak       peaking EQ      args: FREQ, Q, GAIN_DB
    lowshelf   low shelf       args: FREQ, _,  GAIN_DB   (Q ignored, pass any)
    highshelf  high shelf      args: FREQ, _,  GAIN_DB   (Q ignored, pass any)
    bp         bandpass (0 dB) args: FREQ, Q
    notch      notch           args: FREQ, Q

Output cartridge: 1 active stage + 11 passthrough, replicated to all four
corners (M0_Q0, M0_Q100, M100_Q0, M100_Q100). Static filter — does not
morph; runtime morph/Q controls do nothing for these demos. That is the
point: prove the cascade accepts any biquad before adding morph behavior.

Examples
--------
    python authoring/tools/design_df2.py warm_lp lp 200 0.7
    python authoring/tools/design_df2.py air_shelf highshelf 7200 0.7 6
    python authoring/tools/design_df2.py mid_peak peak 1200 4.0 9
"""

import argparse
import json
import math
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# RBJ cookbook biquads. Each returns the five trench-core direct-form
# coefficients (c0..c4 = b0/a0, b1/a0, b2/a0, a1/a0, a2/a0).
# Source: Robert Bristow-Johnson, "Cookbook formulae for audio EQ biquad
# filter coefficients."
# ---------------------------------------------------------------------------

def _normalize(b0, b1, b2, a0, a1, a2):
    return [b0 / a0, b1 / a0, b2 / a0, a1 / a0, a2 / a0]


def rbj_lowpass(f0, Q, sr):
    w0 = 2.0 * math.pi * f0 / sr
    cosw, sinw = math.cos(w0), math.sin(w0)
    alpha = sinw / (2.0 * Q)
    b0 = (1 - cosw) / 2.0
    b1 = 1 - cosw
    b2 = (1 - cosw) / 2.0
    a0 = 1 + alpha
    a1 = -2 * cosw
    a2 = 1 - alpha
    return _normalize(b0, b1, b2, a0, a1, a2)


def rbj_highpass(f0, Q, sr):
    w0 = 2.0 * math.pi * f0 / sr
    cosw, sinw = math.cos(w0), math.sin(w0)
    alpha = sinw / (2.0 * Q)
    b0 = (1 + cosw) / 2.0
    b1 = -(1 + cosw)
    b2 = (1 + cosw) / 2.0
    a0 = 1 + alpha
    a1 = -2 * cosw
    a2 = 1 - alpha
    return _normalize(b0, b1, b2, a0, a1, a2)


def rbj_peak(f0, Q, gain_db, sr):
    A = 10 ** (gain_db / 40.0)
    w0 = 2.0 * math.pi * f0 / sr
    cosw, sinw = math.cos(w0), math.sin(w0)
    alpha = sinw / (2.0 * Q)
    b0 = 1 + alpha * A
    b1 = -2 * cosw
    b2 = 1 - alpha * A
    a0 = 1 + alpha / A
    a1 = -2 * cosw
    a2 = 1 - alpha / A
    return _normalize(b0, b1, b2, a0, a1, a2)


def _shelf_alpha(A, sinw, slope_S=1.0):
    # Per cookbook: alpha = sin(w0)/2 * sqrt((A + 1/A)*(1/S - 1) + 2)
    return sinw / 2.0 * math.sqrt((A + 1.0 / A) * (1.0 / slope_S - 1) + 2)


def rbj_lowshelf(f0, gain_db, sr):
    A = 10 ** (gain_db / 40.0)
    w0 = 2.0 * math.pi * f0 / sr
    cosw, sinw = math.cos(w0), math.sin(w0)
    alpha = _shelf_alpha(A, sinw)
    twoSqrtAalpha = 2.0 * math.sqrt(A) * alpha
    b0 = A * ((A + 1) - (A - 1) * cosw + twoSqrtAalpha)
    b1 = 2 * A * ((A - 1) - (A + 1) * cosw)
    b2 = A * ((A + 1) - (A - 1) * cosw - twoSqrtAalpha)
    a0 = (A + 1) + (A - 1) * cosw + twoSqrtAalpha
    a1 = -2 * ((A - 1) + (A + 1) * cosw)
    a2 = (A + 1) + (A - 1) * cosw - twoSqrtAalpha
    return _normalize(b0, b1, b2, a0, a1, a2)


def rbj_highshelf(f0, gain_db, sr):
    A = 10 ** (gain_db / 40.0)
    w0 = 2.0 * math.pi * f0 / sr
    cosw, sinw = math.cos(w0), math.sin(w0)
    alpha = _shelf_alpha(A, sinw)
    twoSqrtAalpha = 2.0 * math.sqrt(A) * alpha
    b0 = A * ((A + 1) + (A - 1) * cosw + twoSqrtAalpha)
    b1 = -2 * A * ((A - 1) + (A + 1) * cosw)
    b2 = A * ((A + 1) + (A - 1) * cosw - twoSqrtAalpha)
    a0 = (A + 1) - (A - 1) * cosw + twoSqrtAalpha
    a1 = 2 * ((A - 1) - (A + 1) * cosw)
    a2 = (A + 1) - (A - 1) * cosw - twoSqrtAalpha
    return _normalize(b0, b1, b2, a0, a1, a2)


def rbj_bandpass_0db(f0, Q, sr):
    # Constant 0 dB peak gain (cookbook BPF, "constant 0 dB peak gain" form).
    w0 = 2.0 * math.pi * f0 / sr
    cosw, sinw = math.cos(w0), math.sin(w0)
    alpha = sinw / (2.0 * Q)
    b0 = alpha
    b1 = 0.0
    b2 = -alpha
    a0 = 1 + alpha
    a1 = -2 * cosw
    a2 = 1 - alpha
    return _normalize(b0, b1, b2, a0, a1, a2)


def rbj_notch(f0, Q, sr):
    w0 = 2.0 * math.pi * f0 / sr
    cosw, sinw = math.cos(w0), math.sin(w0)
    alpha = sinw / (2.0 * Q)
    b0 = 1.0
    b1 = -2 * cosw
    b2 = 1.0
    a0 = 1 + alpha
    a1 = -2 * cosw
    a2 = 1 - alpha
    return _normalize(b0, b1, b2, a0, a1, a2)


# ---------------------------------------------------------------------------
# Cartridge assembly. compiled-v1 wire format per cartridge.schema.json:
#   keyframes is a list of 4 entries (M0_Q0, M0_Q100, M100_Q0, M100_Q100),
#   each with exactly 12 stages, each stage carrying c0..c4.
# Static demo: stage 0 = active biquad, stages 1..11 = passthrough; all four
# keyframes carry the same stages so morph and Q on the runtime axes do
# nothing. The point is to prove the cascade accepts the biquad — not to
# demonstrate morph behavior.
# ---------------------------------------------------------------------------

PASSTHROUGH_STAGE = {"c0": 1.0, "c1": 0.0, "c2": 0.0, "c3": 0.0, "c4": 0.0}


def _stage_from_coeffs(c):
    return {"c0": c[0], "c1": c[1], "c2": c[2], "c3": c[3], "c4": c[4]}


def _make_keyframe(label, morph, q, active_stage):
    stages = [active_stage] + [dict(PASSTHROUGH_STAGE) for _ in range(11)]
    return {"label": label, "morph": morph, "q": q, "boost": 1.0, "stages": stages}


def build_cartridge(name, c, sr, provenance):
    active = _stage_from_coeffs(c)
    return {
        "format": "compiled-v1",
        "name": name,
        "provenance": provenance,
        "sampleRate": sr,
        "stages": 12,
        "keyframes": [
            _make_keyframe("M0_Q0", 0.0, 0.0, active),
            _make_keyframe("M0_Q100", 0.0, 1.0, active),
            _make_keyframe("M100_Q0", 1.0, 0.0, active),
            _make_keyframe("M100_Q100", 1.0, 1.0, active),
        ],
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

FILTER_TYPES = ("lp", "hp", "peak", "lowshelf", "highshelf", "bp", "notch")


def design(filter_type, freq, Q, gain_db, sr):
    if filter_type == "lp":
        return rbj_lowpass(freq, Q, sr)
    if filter_type == "hp":
        return rbj_highpass(freq, Q, sr)
    if filter_type == "peak":
        return rbj_peak(freq, Q, gain_db, sr)
    if filter_type == "lowshelf":
        return rbj_lowshelf(freq, gain_db, sr)
    if filter_type == "highshelf":
        return rbj_highshelf(freq, gain_db, sr)
    if filter_type == "bp":
        return rbj_bandpass_0db(freq, Q, sr)
    if filter_type == "notch":
        return rbj_notch(freq, Q, sr)
    raise ValueError(f"unknown filter type: {filter_type}")


def repo_root():
    return Path(__file__).resolve().parents[2]


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="DF-II direct-coefficient demo authoring (RBJ cookbook -> compiled-v1).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python authoring/tools/design_df2.py warm_lp lp 200 0.7\n"
            "  python authoring/tools/design_df2.py air_shelf highshelf 7200 0.7 6\n"
            "  python authoring/tools/design_df2.py mid_peak peak 1200 4.0 9\n"
        ),
    )
    ap.add_argument("name", help="cartridge name and output filename stem")
    ap.add_argument("type", choices=FILTER_TYPES, help="filter type")
    ap.add_argument("freq", type=float, help="cutoff/center frequency in Hz")
    ap.add_argument("Q", type=float, help="quality factor (ignored for shelves; pass any value)")
    ap.add_argument("gain_db", type=float, nargs="?", default=0.0,
                    help="peak/shelf gain in dB (default 0)")
    ap.add_argument("--sr", type=float, default=44100.0,
                    help="sample rate to design at (default 44100; the plugin's runtime SR)")
    ap.add_argument("--out", type=Path, default=None,
                    help="output JSON path (default: cartridges/factory/my_designs/<name>.json)")
    args = ap.parse_args(argv)

    c = design(args.type, args.freq, args.Q, args.gain_db, args.sr)

    provenance = (
        f"df2-direct-coeffs (RBJ {args.type} f0={args.freq} Q={args.Q} "
        f"gain={args.gain_db}dB sr={args.sr})"
    )
    cart = build_cartridge(args.name, c, args.sr, provenance)

    out = args.out or (
        repo_root() / "cartridges" / "factory" / "my_designs" / f"{args.name}.json"
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(cart, indent=2), encoding="utf-8")

    print(f"Wrote {out}")
    print(
        f"  type={args.type} freq={args.freq} Hz Q={args.Q} "
        f"gain={args.gain_db} dB sr={args.sr} Hz"
    )
    print(
        f"  c0={c[0]:.6f} c1={c[1]:.6f} c2={c[2]:.6f} "
        f"c3={c[3]:.6f} c4={c[4]:.6f}"
    )
    print("  cascade: stage 0 active, stages 1-11 passthrough; static across all 4 corners")
    print()
    print("  Audition:")
    print(f"    python authoring/audition.py \"{out}\" --source pinknoise --duration 6")
    print("  Plot the actual cascade response:")
    print(f"    python authoring/compilers/plot_magnitude.py \"{out}\"")


if __name__ == "__main__":
    main()
