"""Genetic splice — cross-breed heritage XML donor DNA into 4 mutant bodies.

Donors:
  Bass Boost LP 1  — sub pressure, type 2 formants, type 1 bandpass
  Complex BP 1     — alternating kill/boost, type 3 asymmetric
  Steep 8 Pole     — cascaded type 2 frequency wall
  Crazy Rez        — triple max-gain resonance
  Six Pole Extreme Q — stacked identical resonance
  On the Edge      — varied type 2/3 with moderate gains

StepSweeper64FUN — 64-step alternating decay used as morph corruption source.
"""

import sys, os, json, math
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from pyruntime.designer_compile import (
    make_four_corner_template,
    compile_four_corner_to_body,
    DesignerCell,
)
from pyruntime.analysis import body_profile, shipping_gate, midpoint_audit
from pyruntime.freq_response import freq_points, cascade_response_db
from pyruntime.constants import SR
import numpy as np

# ---------------------------------------------------------------------------
# StepSweeper64FUN — raw 64-step array
# ---------------------------------------------------------------------------

STEP_SWEEPER = [
    1.0, -1.0, 0.984375, -0.9375, 0.9375, -0.875, 0.859375, -0.8125,
    0.796875, -0.734375, 0.6875, -0.65625, 0.640625, -0.578125, 0.5625,
    -0.484375, 0.46875, -0.390625, 0.40625, -0.328125, 0.34375, -0.265625,
    0.296875, -0.21875, 0.234375, -0.15625, 0.140625, -0.109375, 0.171875,
    -0.140625, 0.203125, -0.21875, 0.25, -0.171875, 0.1875, -0.140625,
    0.109375, -0.171875, 0.171875, -0.234375, 0.234375, -0.296875, 0.3125,
    -0.359375, 0.390625, -0.4375, 0.46875, -0.5, 0.53125, -0.5625,
    0.59375, -0.625, 0.65625, -0.6875, 0.734375, -0.734375, 0.828125,
    -0.796875, 0.90625, -0.890625, 0.96875, -0.953125, 1.0, -1.0,
]


def corrupt_value(base: int, step_idx: int, intensity: float = 0.5) -> int:
    """Corrupt a 0-127 value using StepSweeper as mutation source."""
    step_val = STEP_SWEEPER[step_idx % 64]
    mutation = int(step_val * intensity * 64)
    return max(0, min(127, base + mutation))


def corrupt_freq(base: int, step_idx: int, intensity: float = 0.7) -> int:
    """Corrupt frequency — allow wider swings."""
    step_val = STEP_SWEEPER[step_idx % 64]
    mutation = int(step_val * intensity * 80)
    return max(0, min(127, base + mutation))


# ---------------------------------------------------------------------------
# Donor DNA — extracted raw sections from XML
# ---------------------------------------------------------------------------

# Bass Boost LP 1 sections: (type, low_freq, low_gain, high_freq, high_gain)
BASS_BOOST = [
    (2, 25, 92, 103, 0),     # sec 1: type2 formant, sub boost
    (0, 19, 0, 34, 0),       # sec 2: bypass morph
    (1, 0, 107, 0, 109),     # sec 3: type1 bandpass, locked freq=0, high gain
    (0, 12, 127, 12, 100),   # sec 4: bypass but max gain DNA
    (1, 5, 101, 12, 101),    # sec 5: type1 bandpass, sub-freq sweep
    (2, 19, 0, 103, 127),    # sec 6: type2 formant, darkness to brightness
]

# Complex BP 1 sections
COMPLEX_BP = [
    (3, 0, 127, 116, 87),    # sec 1: type3 asymmetric, max gain
    (1, 15, 0, 33, 0),       # sec 2: type1, killed gain
    (1, 33, 127, 51, 127),   # sec 3: type1, max gain both endpoints
    (1, 51, 0, 70, 0),       # sec 4: type1, killed gain (void)
    (1, 70, 127, 89, 127),   # sec 5: type1, max gain sweep
    (2, 0, 127, 127, 127),   # sec 6: type2, full range max gain
]

# Steep 8 Pole sections
STEEP_8 = [
    (2, 0, 0, 127, 0),       # sec 1: type2, full sweep, zero gain
    (2, 110, 0, 127, 0),     # sec 2: type2, near-nyquist
    (2, 94, 0, 127, 0),      # sec 3: type2, high freq
    (2, 82, 0, 127, 0),      # sec 4: type2, descending
    (0, 0, 0, 0, 0),         # sec 5: off
    (0, 0, 0, 0, 0),         # sec 6: off
]

# Crazy Rez sections
CRAZY_REZ = [
    (1, 21, 127, 68, 127),   # sec 1: max gain, wide sweep
    (1, 23, 127, 69, 114),   # sec 2: max gain, wide sweep
    (1, 31, 127, 102, 127),  # sec 3: max gain, massive sweep
    (0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0),
]

# Six Pole Extreme Q sections
SIXPOLE_Q = [
    (2, 26, 127, 127, 76),   # sec 1: type2, sub to nyquist
    (2, 26, 127, 127, 76),   # sec 2: identical
    (2, 26, 127, 127, 76),   # sec 3: identical
    (0, 21, 127, 127, 0),    # sec 4: bypass DNA
    (0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0),
]

# On the Edge sections
ON_EDGE = [
    (2, 0, 45, 54, 48),      # sec 1: type2, moderate
    (3, 14, 27, 69, 75),     # sec 2: type3, climbing
    (2, 41, 69, 82, 88),     # sec 3: type2, mid push
    (2, 68, 109, 85, 102),   # sec 4: type2, high pressure
    (0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0),
]


def cell(t: int, f: int, g: int) -> dict:
    return {"type": max(0, min(3, t)), "freq": max(0, min(127, f)), "gain": max(0, min(127, g))}


# ---------------------------------------------------------------------------
# MUTANT 1: "TubeJam Splice"
#
# Concept: Bass Boost's sub-pressure foundation (type2 formant at freq 25)
# crossed with Complex BP's type3 asymmetric aggression and Steep 8 Pole's
# cascaded frequency wall. StepSweeper corrupts the Q axis.
#
# Strategy: Fixed Anchor (#9) — stage 1 is Bass Boost's anchored sub,
# stages 2-6 are Complex BP / Steep hybrids that pile up at morph=1.
# ---------------------------------------------------------------------------

def build_tubejam():
    # M0_Q0: dark anchor state — Bass Boost sub + Complex BP voids
    m0_q0 = [
        cell(2, 25, 92),           # Bass Boost sec1: sub formant anchor
        cell(3, 0, 127),           # Complex BP sec1: type3 asymmetric max
        cell(1, 0, 107),           # Bass Boost sec3: locked-freq bandpass
        cell(1, 33, 127),          # Complex BP sec3: mid bandpass max gain
        cell(2, 110, 0),           # Steep sec2: near-nyquist formant (dead)
        cell(2, 94, 0),            # Steep sec3: high formant (dead)
    ]

    # M100_Q0: morph open — Steep wall descends, Complex BP peaks activate
    m100_q0 = [
        cell(2, 25, 92),           # Sub anchor HOLDS (invariant)
        cell(3, 116, 87),          # Complex BP sec1 morphed: type3 pushed high
        cell(1, 0, 109),           # Bass Boost sec3 morphed: gain push
        cell(1, 51, 127),          # Complex BP sec3 morphed: freq climbs
        cell(2, 127, 76),          # Steep wall: top of range, Six Pole gain
        cell(2, 127, 76),          # Steep wall: stacked at Nyquist
    ]

    # M0_Q100: Q corrupted by StepSweeper — tighten sub, break highs
    m0_q100 = [
        cell(2, corrupt_freq(25, 3), corrupt_value(120, 7)),   # sub tightened
        cell(3, corrupt_freq(0, 11), corrupt_value(127, 15)),  # type3 mutated
        cell(1, corrupt_freq(0, 19), corrupt_value(107, 23)),  # bandpass warped
        cell(1, corrupt_freq(33, 27), corrupt_value(127, 31)), # mid peak shifted
        cell(2, corrupt_freq(110, 35), corrupt_value(0, 39)),  # near-nyquist (still dead)
        cell(2, corrupt_freq(94, 43), corrupt_value(0, 47)),   # high (still dead)
    ]

    # M100_Q100: full chaos — everything pushed, StepSweeper corruption max
    m100_q100 = [
        cell(2, corrupt_freq(25, 5, 0.3), corrupt_value(92, 9, 0.8)),  # sub anchor survives
        cell(3, corrupt_freq(116, 13, 0.9), corrupt_value(127, 17, 0.9)),  # type3 peak chaos
        cell(1, corrupt_freq(0, 21, 0.5), corrupt_value(127, 25, 0.9)),    # bandpass maxed
        cell(1, corrupt_freq(70, 29, 0.8), corrupt_value(127, 33, 0.5)),   # freq from Complex BP sec5
        cell(2, corrupt_freq(127, 37, 0.6), corrupt_value(127, 41, 0.9)),  # Nyquist alive at max gain
        cell(2, corrupt_freq(127, 45, 0.7), corrupt_value(127, 49, 0.9)),  # stacked Nyquist screamer
    ]

    return make_four_corner_template("TubeJam Splice", {
        "M0_Q0": m0_q0, "M0_Q100": m0_q100,
        "M100_Q0": m100_q0, "M100_Q100": m100_q100,
    })


# ---------------------------------------------------------------------------
# MUTANT 2: "PoleCross Fracture"
#
# Concept: Crazy Rez's triple max-gain resonances in the low state, morphing
# into On the Edge's type2/3 mix. Q axis crosses Complex BP's alternating
# kill/boost pattern onto the Crazy Rez fundamentals.
#
# Strategy: Resonance Phalanx (#5) at morph=0, Frequency Swap (#1) at morph=1
# ---------------------------------------------------------------------------

def build_polecross():
    # M0_Q0: Crazy Rez phalanx — 3 max-gain resonances, close-packed
    m0_q0 = [
        cell(1, 21, 127),          # Crazy Rez sec1: low peak
        cell(1, 23, 127),          # Crazy Rez sec2: adjacent peak (2 semitones!)
        cell(1, 31, 127),          # Crazy Rez sec3: third peak
        cell(2, 0, 127),           # Complex BP sec6 low: type2 at bottom
        cell(1, 70, 127),          # Complex BP sec5 low: upper peak
        cell(3, 0, 127),           # Complex BP sec1 low: type3 anchor
    ]

    # M100_Q0: frequency swap — peaks rearrange violently
    m100_q0 = [
        cell(1, 68, 127),          # Crazy sec1 morphed
        cell(1, 69, 114),          # Crazy sec2 morphed (gain drops)
        cell(1, 102, 127),         # Crazy sec3 morphed: massive freq jump
        cell(2, 127, 127),         # Complex BP sec6 high: max range
        cell(1, 89, 127),          # Complex BP sec5 high
        cell(3, 116, 87),          # Complex BP sec1 high: type3 at top
    ]

    # M0_Q100: On the Edge DNA crosses onto Crazy Rez frequencies
    m0_q100 = [
        cell(2, 21, 45),           # Crazy freq, On the Edge type2/gain
        cell(3, 23, 27),           # Crazy freq, On the Edge type3/gain
        cell(2, 31, 69),           # Crazy freq, On the Edge type2/gain
        cell(2, 0, 109),           # On the Edge sec4 gain onto Complex freq
        cell(1, 70, 0),            # Complex BP sec4 kill pattern
        cell(3, 14, 75),           # On the Edge sec2 high values
    ]

    # M100_Q100: total crossing — frequencies from morph=1, gains from Q=1
    m100_q100 = [
        cell(2, 68, corrupt_value(48, 0, 0.9)),     # Edge gain + Crazy morphed freq
        cell(3, 69, corrupt_value(75, 8, 0.8)),      # type3 + Crazy morphed freq
        cell(2, 102, corrupt_value(88, 16, 0.7)),    # mid push at Crazy's swept freq
        cell(2, 127, corrupt_value(102, 24, 0.6)),   # Complex high + Edge pressure
        cell(1, 89, corrupt_value(0, 32, 0.5)),      # kill zone at Complex freq
        cell(3, 116, corrupt_value(75, 40, 0.9)),    # type3 screamer
    ]

    return make_four_corner_template("PoleCross Fracture", {
        "M0_Q0": m0_q0, "M0_Q100": m0_q100,
        "M100_Q0": m100_q0, "M100_Q100": m100_q100,
    })


# ---------------------------------------------------------------------------
# MUTANT 3: "Nyquist Shatter"
#
# Concept: Steep 8 Pole's descending frequency wall (110, 94, 82) stacked
# with Six Pole Extreme Q's identical resonance and Bass Boost's sub anchor.
# At morph=1 the wall folds DOWN into the sub range — total spectral inversion.
#
# Strategy: HF Wall Fold (#3) — 4 HF stages folding to mids/subs
# ---------------------------------------------------------------------------

def build_nyquist_shatter():
    # M0_Q0: the wall — 4 type2 stages stacked at HF + sub anchor
    m0_q0 = [
        cell(2, 12, 127),          # Bass Boost sec4: sub anchor, max gain forced
        cell(2, 110, 76),          # Steep sec2 + Six Pole gain
        cell(2, 94, 76),           # Steep sec3 + Six Pole gain
        cell(2, 82, 76),           # Steep sec4 + Six Pole gain
        cell(2, 26, 127),          # Six Pole sec1: low anchor duplicate
        cell(1, 5, 101),           # Bass Boost sec5: sub bandpass
    ]

    # M100_Q0: the fold — HF wall descends into low-mids
    m100_q0 = [
        cell(2, 12, 100),          # Sub anchor HOLDS (slight gain drop)
        cell(2, 33, 127),          # was 110 → folds to Complex BP sec3 freq
        cell(2, 25, 127),          # was 94 → folds to Bass Boost sec1 freq
        cell(2, 15, 127),          # was 82 → folds below, gain maxed
        cell(2, 127, 76),          # Six Pole morphed: swaps to Nyquist
        cell(1, 12, 101),          # Bass Boost sec5 morphed: slight climb
    ]

    # M0_Q100: identical wall but Six Pole triplication kicks in
    m0_q100 = [
        cell(2, 12, 127),          # Sub anchor
        cell(2, 110, 127),         # Steep + max gain (Six Pole x3 DNA)
        cell(2, 110, 127),         # DUPLICATE — identical resonance phalanx
        cell(2, 110, 127),         # TRIPLICATE — constructive interference
        cell(2, 26, 127),          # Six Pole anchor
        cell(1, 5, 127),           # Sub bandpass maxed
    ]

    # M100_Q100: wall fold + triplication = total spectral collapse
    m100_q100 = [
        cell(2, corrupt_freq(12, 2, 0.2), 127),   # Sub barely moves
        cell(2, 33, corrupt_value(127, 10, 0.9)),  # folded + corrupted
        cell(2, 33, corrupt_value(127, 18, 0.8)),  # folded + stacked
        cell(2, 33, corrupt_value(127, 26, 0.7)),  # folded + tripled at same freq
        cell(2, corrupt_freq(127, 34, 0.9), corrupt_value(127, 42, 0.9)),  # Nyquist chaos
        cell(1, corrupt_freq(12, 50, 0.5), corrupt_value(127, 58, 0.8)),   # sub bandpass mutant
    ]

    return make_four_corner_template("Nyquist Shatter", {
        "M0_Q0": m0_q0, "M0_Q100": m0_q100,
        "M100_Q0": m100_q0, "M100_Q100": m100_q100,
    })


# ---------------------------------------------------------------------------
# MUTANT 4: "Comb Destroyer"
#
# Concept: Complex BP's alternating kill/boost bandpass isolation spliced
# with StepSweeper-corrupted frequency placement. Every stage frequency is
# derived from the 64-step array mapped onto 0-127, creating non-harmonic
# comb spacing that shifts across morph.
#
# Strategy: Broadband Diffuser (#10) into Cluster Sweep (#6)
# ---------------------------------------------------------------------------

def _step_to_freq(step_idx: int) -> int:
    """Map StepSweeper value to 0-127 frequency."""
    val = STEP_SWEEPER[step_idx % 64]
    return max(0, min(127, int((val + 1.0) * 63.5)))


def build_comb_destroyer():
    # Build frequency arrays from StepSweeper — non-harmonic spacing
    # Pick 6 evenly-spaced but phase-corrupted step positions
    freq_steps_low  = [0, 10, 20, 30, 40, 50]
    freq_steps_high = [5, 15, 25, 35, 45, 55]
    freq_steps_q    = [3, 13, 23, 33, 43, 53]
    freq_steps_qh   = [7, 17, 27, 37, 47, 57]

    # Alternating type pattern from Complex BP: kill, peak, kill, peak, kill, peak
    # But use type2 for peaks and type1 for kills (inverted)
    type_pattern = [1, 2, 1, 2, 1, 2]
    # Gain pattern: Complex BP alternating 0/127
    gain_lo = [127, 0, 127, 0, 127, 0]
    gain_hi = [0, 127, 0, 127, 0, 127]

    m0_q0 = [
        cell(type_pattern[i], _step_to_freq(freq_steps_low[i]), gain_lo[i])
        for i in range(6)
    ]

    m100_q0 = [
        cell(type_pattern[i], _step_to_freq(freq_steps_high[i]), gain_hi[i])
        for i in range(6)
    ]

    m0_q100 = [
        cell(type_pattern[i], _step_to_freq(freq_steps_q[i]), gain_lo[5-i])
        for i in range(6)
    ]

    m100_q100 = [
        cell(type_pattern[i], _step_to_freq(freq_steps_qh[i]), gain_hi[5-i])
        for i in range(6)
    ]

    return make_four_corner_template("Comb Destroyer", {
        "M0_Q0": m0_q0, "M0_Q100": m0_q100,
        "M100_Q0": m100_q0, "M100_Q100": m100_q100,
    })


# ---------------------------------------------------------------------------
# Build, measure, export
# ---------------------------------------------------------------------------

def measure_body(body):
    """Run full measurement suite."""
    freqs = freq_points()
    sr = SR
    results = {}

    # Corner responses
    for m_val, m_label in [(0.0, "M0"), (0.5, "M50"), (1.0, "M100")]:
        for q_val, q_label in [(0.0, "Q0"), (0.5, "Q50"), (1.0, "Q100")]:
            enc = body.corners.interpolate(m_val, q_val)
            db = cascade_response_db(enc, freqs, sr)
            peak = float(np.max(db))
            valley = float(np.min(db))
            results[f"{m_label}_{q_label}_peak"] = round(peak, 1)
            results[f"{m_label}_{q_label}_valley"] = round(valley, 1)

    # Profile
    profile = body_profile(body)
    results["morph_distance"] = profile["morph_distance"]
    results["active_stages"] = profile["active_stages"]
    results["spectral_tilt_db"] = profile["spectral_tilt_db"]
    results["midpoint_audit"] = profile["midpoint_audit"]

    return results


def main():
    builders = [
        ("TubeJam Splice", build_tubejam),
        ("PoleCross Fracture", build_polecross),
        ("Nyquist Shatter", build_nyquist_shatter),
        ("Comb Destroyer", build_comb_destroyer),
    ]

    cartridge_dir = os.path.join(os.path.dirname(__file__), "..", "cartridges", "mutants")
    os.makedirs(cartridge_dir, exist_ok=True)

    for name, builder in builders:
        print(f"\n{'='*70}")
        print(f"  MUTANT: {name}")
        print(f"{'='*70}")

        template = builder()
        body = compile_four_corner_to_body(template, boost=4.0)

        # Measure
        measurements = measure_body(body)

        print(f"\n  Morph distance: Q0={measurements['morph_distance']['q0_rms_db']:.1f} dB, "
              f"Q1={measurements['morph_distance']['q1_rms_db']:.1f} dB")
        print(f"  Active stages: {measurements['active_stages']}")
        print(f"  Spectral tilt: {measurements['spectral_tilt_db']:.1f} dB")
        print(f"  Midpoint audit: {'PASS' if measurements['midpoint_audit']['passed'] else 'FAIL'} "
              f"(worst peak: {measurements['midpoint_audit']['worst_peak_db']:.1f} dB)")

        # Corner peaks
        print(f"\n  Corner peaks (dB):")
        for key, val in sorted(measurements.items()):
            if key.endswith("_peak") and "midpoint" not in key:
                print(f"    {key}: {val}")

        # Export compiled JSON
        compiled_json = body.to_compiled_json(provenance="genetic_splice")
        filename = name.replace(" ", "_") + ".compiled.json"
        filepath = os.path.join(cartridge_dir, filename)
        with open(filepath, "w") as f:
            f.write(compiled_json)
        print(f"\n  Exported: {filepath}")

        # Also export keyframe JSON
        kf_json = body.to_json()
        kf_filename = name.replace(" ", "_") + ".keyframe.json"
        kf_filepath = os.path.join(cartridge_dir, kf_filename)
        with open(kf_filepath, "w") as f:
            f.write(kf_json)
        print(f"  Exported: {kf_filepath}")

    print(f"\n{'='*70}")
    print(f"  4 mutants compiled and exported to {cartridge_dir}")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
