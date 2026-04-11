"""
forge_generator.py — Anchor + Vector body generator for TRENCH Sift workflow.

Authoring paradigm: define a resting anatomy (Anchor), then apply two
directional forces (Morph Vector, Q Vector). The compiler resolves these
into 4 corners for the Rust engine and solves shared c4 per corner to
hit a target cascade peak.

Corners are a playback format, not an authoring concept.

Pipeline:
    Blueprint (anchor + vectors)
    → 4-corner Designer Grid
    → heritage compiler (pole/zero structure)
    → c4 solver (shared b0 per corner, target cascade peak)
    → compiled-v1 JSON

Usage:
    python tools/forge_generator.py
    python tools/forge_generator.py --seed 99
    python tools/forge_generator.py --archetype acid_squelch
    python tools/forge_generator.py --target-db 40
"""

import json
import math
import os
import random
import sys
import argparse
from dataclasses import dataclass

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np

from pyruntime.corner import CornerName, CornerState, CornerArray
from pyruntime.encode import EncodedCoeffs, raw_to_encoded
from pyruntime.stage_params import StageParams
from pyruntime.body import Body
from pyruntime.constants import SR, NUM_BODY_STAGES


# ---------------------------------------------------------------------------
# Blueprint: Anchor + Morph Vector + Q Vector — direct Hz/radius space
#
# No heritage compiler. Stages are authored as (freq_hz, radius).
# val1 is derived from the c4 solver (shared b0). val2/val3 default to the
# P2K interior-zero regime unless overridden.
# ---------------------------------------------------------------------------

@dataclass
class Stage:
    freq_hz: float    # pole frequency in Hz
    radius: float     # pole radius (0.0 - 0.999)
    val2: float = 0.0 # zero angle offset (default: interior zero)
    val3: float = 0.0 # zero radius offset (default: interior zero)


@dataclass
class Blueprint:
    name: str
    sentence: str
    anchor: list[Stage]           # 6 stages at M0_Q0
    morph: list[tuple[float, float]]  # per-stage (freq_hz_delta, radius_delta)
    q: list[tuple[float, float]]      # per-stage (freq_hz_delta, radius_delta)
    jitter_hz: float = 20.0      # seed jitter in Hz
    jitter_r: float = 0.002      # seed jitter in radius


def _stage_to_params(s: Stage) -> StageParams:
    """Convert a Stage to StageParams.

    val1 is placeholder (c4 solver overrides via shared b0).
    If val2/val3 are 0, default to "zeros at origin" — pure all-pole resonator.
    """
    theta = 2.0 * math.pi * max(20.0, min(SR / 2.0 - 1.0, s.freq_hz)) / SR
    r = max(0.0, min(0.999, s.radius))
    a1 = -2.0 * r * math.cos(theta)

    val2 = s.val2 if s.val2 != 0.0 else -a1
    val3 = s.val3 if s.val3 != 0.0 else r * r

    return StageParams(a1=a1, r=r, val1=0.0, val2=val2, val3=val3)


def _apply_stage(base: Stage, df: float, dr: float, jf: float, jr: float) -> Stage:
    """Apply morph/Q delta + jitter to a base stage."""
    return Stage(
        freq_hz=max(20.0, min(SR / 2.0 - 1.0, base.freq_hz + df + jf)),
        radius=max(0.0, min(0.999, base.radius + dr + jr)),
        val2=base.val2,
        val3=base.val3,
    )


def blueprint_to_body(bp: Blueprint, rng: random.Random, boost: float = 4.0) -> Body:
    """Resolve a Blueprint directly to a Body in Hz/radius space.

    Bypasses the heritage compiler entirely. StageParams → EncodedCoeffs.
    c4 is set to 1.0 here (placeholder); the c4 solver overrides it later.
    """
    corner_map = {}
    for cn in CornerName:
        use_morph = cn in (CornerName.C, CornerName.D)
        use_q = cn in (CornerName.B, CornerName.D)

        stages = []
        encoded = []
        for i, base in enumerate(bp.anchor):
            df = (bp.morph[i][0] if use_morph else 0.0) + (bp.q[i][0] if use_q else 0.0)
            dr = (bp.morph[i][1] if use_morph else 0.0) + (bp.q[i][1] if use_q else 0.0)
            jf = rng.uniform(-bp.jitter_hz, bp.jitter_hz)
            jr = rng.uniform(-bp.jitter_r, bp.jitter_r)

            s = _apply_stage(base, df, dr, jf, jr)
            sp = _stage_to_params(s)
            stages.append(sp)
            encoded.append(raw_to_encoded(sp))

        # Pad to 12 with passthrough
        while len(stages) < NUM_BODY_STAGES:
            stages.append(StageParams.passthrough())
            encoded.append(EncodedCoeffs(c0=1.0, c1=0.0, c2=0.0, c3=0.0, c4=0.0))

        corner_map[cn] = CornerState(stages=stages, boost=boost, _pre_encoded=encoded)

    corners = CornerArray(
        a=corner_map[CornerName.A],
        b=corner_map[CornerName.B],
        c=corner_map[CornerName.C],
        d=corner_map[CornerName.D],
    )
    return Body(name=bp.name, corners=corners, boost=boost)


# ---------------------------------------------------------------------------
# Archetypes
# ---------------------------------------------------------------------------

BLUEPRINTS: dict[str, Blueprint] = {}


def _reg(bp: Blueprint):
    key = bp.name.lower().replace(" ", "_").replace("-", "_")
    BLUEPRINTS[key] = bp

# =========================================================================
# Packed freq → Hz reference (Type 1/2, gain=64):
#   0→72  10→105  20→152  30→220  40→320  50→463  60→686  70→993
#   80→1440  90→2095  100→3066  105→3800  110→4632  120→7045  127→10259
#
# Each blueprint implements a specific P2K-proven topological strategy.
# The strategy name is the first line of the sentence field.
# =========================================================================


# ── STRATEGY: FREQUENCY SWAP (BassBox 303) ──────────────────────────────
# Lows go high, highs go low. Total stage reordering during morph.
_reg(Blueprint(
    name="Sub Physicalizer",
    sentence="FREQUENCY SWAP. Sub anchor + HF wall. Morph swaps them — bass goes high, highs go low.",
    # BassBox 303: S0=79 Hz r=0.962, S1-S4 HF (10-17 kHz) r=0.964-0.990, S5 DC r=0.126
    anchor=[
        Stage(79, 0.997),       # sub anchor — P2K-grade radius
        Stage(7000, 0.990),     # HF wall 1
        Stage(5500, 0.990),     # HF wall 2
        Stage(9000, 0.985),     # HF wall 3
        Stage(8000, 0.990),     # HF wall 4
        Stage(40, 0.126),       # DC/sub — lossy anchor (BassBox S5)
    ],
    morph=[                     # FREQUENCY SWAP: sub goes high, HF drops to mids
        (5000.0, -0.03),        # sub: 79→5079 Hz
        (-5600.0, 0.005),       # HF1: 7000→1400 Hz — drops to formant
        (-4500.0, 0.005),       # HF2: 5500→1000 Hz
        (-7500.0, 0.005),       # HF3: 9000→1500 Hz — crosses HF1
        (-6500.0, 0.005),       # HF4: 8000→1500 Hz
        (460.0, 0.87),          # DC: 40→500 Hz, radius rockets 0.126→0.996
    ],
    q=[                         # Q scatters HF wall downward
        (4000.0, 0.0),          # sub: 79→4079 Hz
        (-6600.0, 0.005),       # HF1: 7000→400 Hz — crashes
        (-4700.0, 0.005),       # HF2: 5500→800 Hz
        (-8600.0, 0.005),       # HF3: 9000→400 Hz — crashes
        (-7500.0, 0.005),       # HF4: 8000→500 Hz
        (960.0, 0.50),          # DC: 40→1000 Hz, radius to 0.626
    ],
))

# ── RADIUS-ONLY MORPH (Talking Hedz) ─────────────────────────────────────
_reg(Blueprint(
    name="Alien Formant",
    sentence="Morph sharpens radii to 0.999. Frequencies locked. Sub anchor launches at Q.",
    anchor=[Stage(730, 0.975), Stage(1720, 0.977), Stage(2657, 0.970),
            Stage(3500, 0.948), Stage(180, 0.991), Stage(5500, 0.950)],
    morph=[(20.0, 0.024), (-10.0, 0.022), (30.0, 0.029),
           (-20.0, 0.049), (-20.0, 0.008), (20.0, 0.040)],
    q=[(-80.0, 0.005), (50.0, 0.005), (-50.0, 0.005),
       (30.0, 0.005), (1600.0, 0.0), (-3500.0, 0.005)],
))

# ── HF WALL FOLD (Ear Bender) ───────────────────────────────────────────
_reg(Blueprint(
    name="Pad Breather",
    sentence="Four HF stages fold to mids. Bass anchor drops to sub-bass monster.",
    anchor=[Stage(302, 0.997), Stage(3723, 0.997), Stage(10283, 0.998),
            Stage(13571, 0.997), Stage(14394, 0.997), Stage(15837, 0.997)],
    morph=[(-257.0, 0.002), (6467.0, -0.003), (2370.0, -0.008),
           (-25.0, -0.003), (166.0, -0.003), (294.0, -0.003)],
    q=[(5970.0, -0.012), (-2985.0, 0.005), (-7899.0, 0.005),
       (-7661.0, 0.005), (-8873.0, 0.005), (-15023.0, 0.005)],
))

# ── ONE-VS-FIVE (Meaty Gizmo) ───────────────────────────────────────────
_reg(Blueprint(
    name="Circuit Breaker",
    sentence="59 Hz monster at r=0.999. Five HF suppressors cancel 58 dB of it.",
    anchor=[Stage(59, 0.999), Stage(10329, 0.994), Stage(12893, 0.997),
            Stage(13615, 0.998), Stage(14670, 0.998), Stage(16277, 0.997)],
    morph=[(13420.0, 0.0), (-7952.0, -0.640), (-2544.0, 0.003),
           (3207.0, 0.002), (-1188.0, 0.002), (-13645.0, 0.0)],
    q=[(5075.0, 0.0), (-2112.0, -0.005), (-12022.0, 0.005),
       (-10128.0, 0.005), (-12306.0, 0.005), (-15821.0, 0.005)],
))

# ── RESONANCE PHALANX (Fuzzi Face) ──────────────────────────────────────
_reg(Blueprint(
    name="World Splitter",
    sentence="Five poles at 400 Hz constructive pile-up. Morph fans them outward.",
    anchor=[Stage(4036, 0.998), Stage(420, 0.997), Stage(430, 0.997),
            Stage(380, 0.997), Stage(350, 0.997), Stage(440, 0.997)],
    morph=[(0.0, 0.0), (1580.0, -0.01), (1170.0, -0.01),
           (-160.0, -0.01), (-200.0, -0.01), (2060.0, -0.01)],
    q=[(-2000.0, 0.0), (-100.0, 0.002), (50.0, 0.002),
       (80.0, 0.002), (100.0, 0.002), (-150.0, 0.002)],
))

# ── CLUSTER SWEEP (acid) ────────────────────────────────────────────────
_reg(Blueprint(
    name="Acid Squelch",
    sentence="Three poles at 400 Hz sweep to 3500 Hz. DC anchor wakes up.",
    anchor=[Stage(400, 0.997), Stage(420, 0.997), Stage(370, 0.997),
            Stage(460, 0.990), Stage(685, 0.985), Stage(40, 0.126)],
    morph=[(3100.0, -0.005), (2780.0, -0.005), (3630.0, -0.005),
           (2540.0, 0.005), (1415.0, -0.005), (460.0, 0.870)],
    q=[(-80.0, 0.002), (50.0, 0.002), (-50.0, 0.002),
       (500.0, 0.009), (-200.0, 0.005), (200.0, 0.10)],
))

# ── NYQUIST SHARPENING (Razor Blades) ───────────────────────────────────
_reg(Blueprint(
    name="Asymmetric Notch",
    sentence="HF cluster pushes toward 16 kHz. Q pulls everything down to mids.",
    anchor=[Stage(12648, 0.993), Stage(443, 0.980), Stage(1383, 0.987),
            Stage(13518, 0.976), Stage(6183, 0.995), Stage(10096, 0.948)],
    morph=[(1722.0, 0.0), (318.0, -0.074), (221.0, -0.056),
           (2571.0, 0.0), (5478.0, 0.0), (5899.0, 0.0)],
    q=[(-8777.0, 0.006), (-360.0, 0.0), (-248.0, 0.011),
       (-10374.0, 0.022), (-738.0, 0.004), (-7556.0, 0.050)],
))

# ── CONSTRUCTIVE RADIO (Radio Craze) ────────────────────────────────────
_reg(Blueprint(
    name="Lo-Fi Telephone",
    sentence="Three stages at 5 kHz reinforce. Cascade louder than any single stage.",
    anchor=[Stage(4948, 0.988), Stage(199, 0.946), Stage(1044, 0.966),
            Stage(9965, 0.718), Stage(5274, 0.910), Stage(4836, 0.976)],
    morph=[(2740.0, -0.248), (1388.0, -0.098), (1548.0, 0.0),
           (608.0, 0.0), (-1336.0, -0.260), (228.0, 0.0)],
    q=[(-2290.0, 0.010), (1198.0, 0.0), (-693.0, 0.017),
       (-5542.0, 0.280), (-3790.0, 0.088), (-1584.0, 0.016)],
))

# ── FIXED ANCHOR (Early Rizer) ──────────────────────────────────────────
_reg(Blueprint(
    name="Liquid Sharpener",
    sentence="S1 locked at 553 Hz. HF compresses toward it. Ceiling crashes at Q.",
    anchor=[Stage(8802, 0.997), Stage(553, 0.976), Stage(3787, 0.998),
            Stage(11354, 0.997), Stage(6177, 0.997), Stage(16995, 0.910)],
    morph=[(-3883.0, 0.002), (0.0, 0.023), (-2278.0, 0.001),
           (-7484.0, 0.002), (-3702.0, 0.002), (0.0, 0.0)],
    q=[(3021.0, 0.0), (1401.0, 0.0), (99.0, 0.0),
       (3836.0, 0.0), (365.0, 0.0), (-15494.0, -0.035)],
))

# ── VOWEL FORMANT (Ooh to Eee) ──────────────────────────────────────────
_reg(Blueprint(
    name="Vocal Fold",
    sentence="/oo/ to /ee/. Three formants, parallel upward shift.",
    anchor=[Stage(2483, 0.992), Stage(567, 0.996), Stage(222, 0.999),
            Stage(8000, 0.950), Stage(12000, 0.940), Stage(15000, 0.930)],
    morph=[(2455.0, 0.0), (535.0, 0.0), (212.0, 0.0),
           (0.0, 0.0), (0.0, 0.0), (0.0, 0.0)],
    q=[(-486.0, 0.005), (869.0, 0.0), (345.0, 0.0),
       (-3000.0, 0.010), (-5000.0, 0.010), (-7000.0, 0.010)],
))

# ── BROADBAND DIFFUSER (Freak Shifta) ───────────────────────────────────
_reg(Blueprint(
    name="Reece Growl",
    sentence="S5 at r=0.559 diffuser. S1 locked at 858 Hz. Highest c4 body.",
    anchor=[Stage(15645, 0.983), Stage(858, 0.880), Stage(1860, 0.901),
            Stage(12668, 0.901), Stage(10324, 0.974), Stage(12641, 0.559)],
    morph=[(-5872.0, 0.015), (12.0, 0.0), (-615.0, -0.044),
           (384.0, 0.096), (-4481.0, 0.0), (-168.0, 0.014)],
    q=[(-5814.0, 0.005), (-104.0, 0.109), (423.0, 0.094),
       (-2161.0, 0.005), (-5894.0, 0.021), (-4623.0, 0.435)],
))

# ── FREQUENCY SWAP (inverted) ───────────────────────────────────────────
_reg(Blueprint(
    name="Hollow Bone",
    sentence="All HF at rest. Morph brings everything down into mids/bass.",
    anchor=[Stage(10000, 0.997), Stage(8000, 0.997), Stage(6500, 0.997),
            Stage(5500, 0.997), Stage(5000, 0.997), Stage(9000, 0.990)],
    morph=[(-9800.0, 0.0), (-7600.0, 0.0), (-6000.0, 0.0),
           (-4900.0, 0.0), (-4000.0, 0.0), (-8650.0, 0.0)],
    q=[(0.0, 0.0), (-6500.0, 0.0), (0.0, 0.0),
       (-5000.0, 0.0), (1000.0, 0.0), (-6500.0, 0.0)],
))

# ── PHALANX at 1500 Hz ──────────────────────────────────────────────────
_reg(Blueprint(
    name="Swarm Gate",
    sentence="Six poles at 1450 Hz reinforce. Morph fans from 350 Hz to 5 kHz.",
    anchor=[Stage(1440, 0.997), Stage(1450, 0.997), Stage(1430, 0.997),
            Stage(1460, 0.997), Stage(1420, 0.997), Stage(1470, 0.997)],
    morph=[(0.0, 0.0), (960.0, -0.01), (-730.0, -0.01),
           (2060.0, -0.01), (-1070.0, -0.01), (3560.0, -0.01)],
    q=[(0.0, 0.002), (-200.0, 0.002), (200.0, 0.002),
       (-400.0, 0.002), (400.0, 0.002), (-500.0, 0.002)],
))

# ── WALL-TO-BANDPASS ────────────────────────────────────────────────────
_reg(Blueprint(
    name="Iron Curtain",
    sentence="Brickwall LP morphs to bandpass. HP shelf + resonance wake up.",
    anchor=[Stage(1440, 0.997), Stage(1500, 0.997), Stage(1560, 0.997),
            Stage(1620, 0.997), Stage(150, 0.500), Stage(3065, 0.500)],
    morph=[(-1240.0, 0.0), (-1250.0, 0.0), (-1210.0, 0.0),
           (-1170.0, 0.0), (1290.0, 0.490), (3935.0, 0.490)],
    q=[(-100.0, 0.002), (50.0, 0.002), (-50.0, 0.002),
       (100.0, 0.002), (-50.0, -0.10), (200.0, 0.10)],
))

# ── NASAL BODY ──────────────────────────────────────────────────────────
_reg(Blueprint(
    name="Chest Cavity",
    sentence="Two formants + two anti-formants. Morph shrinks the body.",
    anchor=[Stage(250, 0.997), Stage(1100, 0.997), Stage(1000, 0.990),
            Stage(3000, 0.985), Stage(80, 0.999), Stage(5000, 0.950)],
    morph=[(650.0, 0.0), (1400.0, 0.0), (-650.0, 0.005),
           (-1500.0, 0.005), (-30.0, 0.0), (-4100.0, 0.020)],
    q=[(-50.0, 0.002), (200.0, 0.002), (200.0, -0.05),
       (100.0, -0.05), (300.0, 0.0), (-1500.0, 0.010)],
))

# ── HF CLUSTER DROP ─────────────────────────────────────────────────────
_reg(Blueprint(
    name="Wire Brush",
    sentence="Four poles at 5-8 kHz drop to mids. Stages cross during sweep.",
    anchor=[Stage(180, 0.985), Stage(4600, 0.998), Stage(5100, 0.998),
            Stage(5600, 0.998), Stage(6200, 0.998), Stage(9000, 0.990)],
    morph=[(100.0, 0.005), (-4140.0, 0.0), (-4500.0, 0.0),
           (-5200.0, 0.0), (-5400.0, 0.0), (-7000.0, 0.005)],
    q=[(-50.0, 0.005), (200.0, 0.001), (-200.0, 0.001),
       (200.0, 0.001), (-200.0, 0.001), (100.0, 0.005)],
))

# ── FORMANT CLAMP ───────────────────────────────────────────────────────
_reg(Blueprint(
    name="Throat Clamp",
    sentence="F1 rises, F3 drops, they CROSS through F2. Vocal constriction.",
    anchor=[Stage(250, 0.997), Stage(1000, 0.997), Stage(2343, 0.997),
            Stage(3800, 0.995), Stage(686, 0.990), Stage(7000, 0.960)],
    morph=[(1190.0, 0.0), (200.0, 0.0), (-1743.0, 0.0),
           (-2800.0, 0.002), (1314.0, -0.020), (-6540.0, 0.030)],
    q=[(500.0, 0.002), (0.0, 0.002), (500.0, 0.002),
       (300.0, 0.002), (-200.0, 0.005), (300.0, -0.010)],
))

# ── HARMONIC OVERBLOW ───────────────────────────────────────────────────
_reg(Blueprint(
    name="Rust Flute",
    sentence="Bell ratios at rest. Morph kills fundamental, partials take over.",
    anchor=[Stage(220, 0.999), Stage(330, 0.995), Stage(550, 0.990),
            Stage(993, 0.985), Stage(463, 0.990), Stage(2343, 0.960)],
    morph=[(50.0, -0.20), (30.0, -0.15), (-50.0, 0.009),
           (-100.0, 0.014), (500.0, -0.10), (-1993.0, 0.039)],
    q=[(-30.0, 0.0), (30.0, 0.0), (-50.0, 0.005),
       (50.0, 0.005), (-80.0, 0.005), (100.0, 0.005)],
))

# ── ONE HOT POLE (formant range) ────────────────────────────────────────
_reg(Blueprint(
    name="Splinter",
    sentence="730 Hz monster at r=0.999. Five HF suppressors. Morph scatters.",
    anchor=[Stage(730, 0.999), Stage(4600, 0.990), Stage(5500, 0.985),
            Stage(7000, 0.990), Stage(9000, 0.985), Stage(6500, 0.960)],
    morph=[(6270.0, -0.010), (-4140.0, -0.20), (-4300.0, 0.010),
           (-5000.0, 0.005), (-7500.0, 0.005), (-4900.0, 0.030)],
    q=[(2000.0, 0.0), (-4400.0, 0.005), (500.0, 0.005),
       (-6000.0, 0.005), (500.0, 0.005), (-4500.0, 0.005)],
))

# ---------------------------------------------------------------------------
# Shared-b0 c4 solver
#
# Calibration data (Millennium, ft=26) proves:
#   - All stages within a corner share the same c4 (= b0)
#   - c4 is a consequence of pole structure, not a design parameter
#   - Given poles, c4 is whatever keeps cascade peak at target level
# ---------------------------------------------------------------------------

DEFAULT_TARGET_DB = 36.0
GAIN_BUDGET_DB = 48.0


def _is_passthrough(enc: EncodedCoeffs) -> bool:
    return (abs(enc.c0 - 1.0) < 0.01 and abs(enc.c1) < 0.01 and
            abs(enc.c2) < 0.01 and abs(enc.c3) < 0.01 and
            abs(enc.c4) < 0.01)


def _solve_b0_for_corner(
    encoded: list[EncodedCoeffs],
    target_db: float = DEFAULT_TARGET_DB,
    num_freqs: int = 1024,
) -> float:
    """Find the shared b0 (= c0) that produces target_db cascade peak.

    DF2T: H_i(z) = (c0 + c1*z^-1 + c2*z^-2) / (1 + c3*z^-1 + c4*z^-2)

    For all-pole stages (c1≈0, c2≈0), c0 factors out cleanly:
        cascade_peak = c0^N * shape_peak
        c0 = (target / shape_peak)^(1/N)

    For stages with zeros, the numerator shape at c0=1 is used as the
    approximation — accurate when c1, c2 are small relative to c0.
    """
    freqs = np.linspace(20.0, SR / 2.0 - 1.0, num_freqs)
    z_inv = np.exp(-2j * np.pi * freqs / SR)
    z_inv2 = z_inv * z_inv

    # Compute shape peak with c0=1: ∏ |(1 + c1*z^-1 + c2*z^-2) / D_i(z)|
    shape = np.ones(num_freqs)
    n = 0
    for enc in encoded:
        if _is_passthrough(enc):
            continue
        num = 1.0 + enc.c1 * z_inv + enc.c2 * z_inv2
        den = 1.0 + enc.c3 * z_inv + enc.c4 * z_inv2
        H_shape = np.abs(num / np.where(np.abs(den) > 1e-30, den, 1e-30))
        shape *= H_shape
        n += 1

    if n == 0:
        return 1.0

    shape_peak = float(np.max(shape))
    if shape_peak < 1e-30:
        return 0.5

    target_lin = 10.0 ** (target_db / 20.0)
    b0 = (target_lin / shape_peak) ** (1.0 / n)
    return max(0.01, min(4.0, b0))


def _apply_shared_b0(enc: EncodedCoeffs, target_b0: float) -> EncodedCoeffs:
    """Set c0 to target b0. c1..c4 stay fixed — pole/zero positions don't move."""
    if _is_passthrough(enc):
        return enc
    return EncodedCoeffs(c0=target_b0, c1=enc.c1, c2=enc.c2, c3=enc.c3, c4=enc.c4)


# ---------------------------------------------------------------------------
# Body-level zero solver
#
# Constraints:
#   1. Interior zeros only (zero_r < pole_r, NOT unit circle)
#   2. Zero cluster centroid in 4500-8500 Hz band
#   3. Minimum 2-octave pole-zero spread per stage
#   4. Zero positions (frequencies) are morph-invariant — only zero radii
#      change under morph (matching P2K Talking Hedz behavior)
#
# The solver computes a single set of zero frequencies for the M0_Q0 corner,
# then applies them to all 4 corners with per-corner zero radius scaling.
# ---------------------------------------------------------------------------

ZERO_BAND_LO = 4500.0
ZERO_BAND_HI = 8500.0
MIN_POLE_ZERO_SPREAD_OCTAVES = 2.0
ZERO_RADIUS_RATIO = 0.72  # default zero_r / pole_r — interior, not on unit circle
# From Talking Hedz: the target pole (where zeros converge) has the lowest
# ratio (0.53-0.59), other stages have higher ratios (0.70-0.76).
# And zero radii are approximately constant across morph — they don't track
# pole radius changes. This means zero_r is an absolute value, not a ratio.
ZERO_R_ABSOLUTE = 0.71  # mean zero_r from TH M0_Q0 (excluding S4 target pole)
ZERO_R_TARGET_POLE = 0.55  # zero_r for the stage whose pole IS the target


def _compute_zero_freq(pole_hz: float, band_lo: float, band_hi: float) -> float:
    """Place zero frequency inside the [band_lo, band_hi] cluster band.

    Clustering is the primary constraint (from Talking Hedz analysis: all 6 zeros
    land in 4500-8500 Hz regardless of pole position). The 2-octave spread is a
    soft preference, not a hard constraint — Talking Hedz has spreads from 0.8
    to 5.5 octaves, but clustering always wins.

    Within the band, zero placement is proportional to pole position:
    lower poles → lower end of zero band, higher poles → higher end.
    This distributes zeros evenly across the band instead of piling them
    all at the centroid.
    """
    if pole_hz < 1.0:
        return (band_lo + band_hi) / 2.0

    # Map pole frequency to a position within the zero band.
    # Use log-frequency to spread evenly in perceptual space.
    # Pole range for typical bodies: ~100 Hz to ~16000 Hz.
    pole_log = math.log2(max(20.0, pole_hz))
    pole_log_min = math.log2(20.0)
    pole_log_max = math.log2(SR / 2.0)
    t = (pole_log - pole_log_min) / (pole_log_max - pole_log_min)
    t = max(0.0, min(1.0, t))

    # Map to zero band (inverted: low poles → high zeros for maximum spread)
    zero_f = band_hi - t * (band_hi - band_lo)

    # If the pole is inside the zero band, push the zero to the opposite edge
    # to avoid destructive near-zero spread (pole ≈ zero = cancellation, not character)
    if band_lo <= pole_hz <= band_hi:
        mid = (band_lo + band_hi) / 2.0
        if pole_hz < mid:
            zero_f = band_hi  # pole in lower half → zero goes high
        else:
            zero_f = band_lo  # pole in upper half → zero goes low

    return max(band_lo, min(band_hi, zero_f))


def _inject_zeros_into_stage(
    sp: StageParams,
    zero_freq_hz: float,
    zero_r: float = ZERO_R_ABSOLUTE,
) -> StageParams:
    """Inject an interior zero into a StageParams by setting val2/val3.

    val2 = b1 - a1, val3 = r² - b2
    where b1 = -2 * zero_r * cos(zero_theta), b2 = zero_r²

    zero_r is absolute (not a ratio of pole_r). From Talking Hedz: zero radii
    are approximately constant (~0.71) across all corners regardless of pole
    radius changes under morph. This creates the zero dissolution effect —
    as poles sharpen (r→0.999), zeros stay at fixed radius, weakening relative
    to the poles.
    """
    pole_r = sp.r
    if pole_r < 0.01:
        return sp

    # Clamp zero_r below pole_r (interior zero requirement)
    zr = min(zero_r, pole_r - 0.001)
    if zr < 0.01:
        return sp

    zero_theta = 2.0 * math.pi * max(20.0, min(SR / 2.0 - 1.0, zero_freq_hz)) / SR

    b1 = -2.0 * zr * math.cos(zero_theta)
    b2 = zr * zr

    val2 = b1 - sp.a1
    val3 = pole_r * pole_r - b2

    return StageParams(a1=sp.a1, r=sp.r, val1=sp.val1, val2=val2, val3=val3)


def _compute_foreign_pole_targets(pole_freqs: list[float]) -> list[float]:
    """Compute zero target frequencies using foreign-pole alignment.

    From Talking Hedz RE analysis:
    - The PRIMARY target is the pole sitting inside the zero band (4500-8500 Hz).
      All other stages aim their zeros at this pole. (5-to-1 convergence)
    - The primary target stage's own zero targets the HIGHEST-frequency pole
      (cross-targeting: the two anchor stages cancel each other).

    If no pole sits inside the zero band, fall back to highest-frequency targeting.
    """
    # Find the primary target: the pole inside the zero band
    primary_idx = -1
    primary_freq = 0.0
    for i, f in enumerate(pole_freqs):
        if ZERO_BAND_LO <= f <= ZERO_BAND_HI:
            if f > primary_freq:  # pick highest in-band pole
                primary_freq = f
                primary_idx = i

    # Find the highest-frequency pole (for cross-targeting)
    hi_idx = -1
    hi_freq = 0.0
    for i, f in enumerate(pole_freqs):
        if f > hi_freq:
            hi_freq = f
            hi_idx = i

    # Fallback: if no in-band pole, use second-highest as primary
    if primary_idx < 0:
        indexed = sorted([(f, i) for i, f in enumerate(pole_freqs) if f > 1.0],
                         key=lambda x: -x[0])
        if len(indexed) >= 2:
            primary_idx = indexed[1][1]
            primary_freq = indexed[1][0]
        elif len(indexed) >= 1:
            primary_idx = indexed[0][1]
            primary_freq = indexed[0][0]
        else:
            centroid = (ZERO_BAND_LO + ZERO_BAND_HI) / 2.0
            return [centroid] * len(pole_freqs)

    # If primary and highest are the same stage, pick second-highest for cross
    if primary_idx == hi_idx:
        indexed = sorted([(f, i) for i, f in enumerate(pole_freqs) if f > 1.0],
                         key=lambda x: -x[0])
        if len(indexed) >= 2:
            hi_idx = indexed[1][1]
            hi_freq = indexed[1][0]

    zero_targets = []
    for i, pf in enumerate(pole_freqs):
        if pf < 1.0:
            zero_targets.append(0.0)
        elif i == primary_idx:
            # Primary target stage → its zero targets the highest pole
            zero_targets.append(hi_freq)
        else:
            # Everyone else → zero targets the primary (in-band) pole
            zero_targets.append(primary_freq)

    return zero_targets


def solve_zero_network(body: Body) -> Body:
    """Body-level zero solver using foreign-pole alignment.

    From Talking Hedz RE analysis: stages aim their zeros at other stages' poles.
    The two highest-frequency poles serve as targets; all other stages' zeros
    converge on them. This creates 41+ dB of cancellation via constructive
    zero pileup at specific frequencies.

    Zero frequencies are computed once from M0_Q0 (morph-invariant).
    Zero radii are absolute (not ratios of pole_r) — from TH data, zero_r ≈ 0.71
    for non-target stages and ≈ 0.55 for the target-pole stage. This is constant
    across corners, creating the zero dissolution effect as poles sharpen under morph.
    """
    # Step 1: Read pole frequencies from M0_Q0 anchor
    anchor_stages = body.corners.corner(CornerName.A).stages
    pole_freqs = [sp.freq_hz() if sp.r > 0.01 else 0.0 for sp in anchor_stages[:6]]

    # Step 2: Compute zero targets via foreign-pole alignment
    zero_freqs = _compute_foreign_pole_targets(pole_freqs)

    # Step 3: Identify target pole indices (stages whose poles are targeted by others)
    indexed = [(f, i) for i, f in enumerate(pole_freqs) if f > 1.0]
    indexed.sort(key=lambda x: -x[0])
    target_indices = set()
    if len(indexed) >= 2:
        target_indices = {indexed[0][1], indexed[1][1]}

    # Step 4: Apply morph-invariant zero frequencies to all corners
    # Zero_r is absolute: ~0.55 for target pole stages, ~0.71 for others
    new_corners = {}
    for cn in CornerName:
        corner = body.corners.corner(cn)
        new_stages = []
        new_encoded = []
        for i, sp in enumerate(corner.stages):
            if i < 6 and zero_freqs[i] > 0:
                zr = ZERO_R_TARGET_POLE if i in target_indices else ZERO_R_ABSOLUTE
                sp_z = _inject_zeros_into_stage(sp, zero_freqs[i], zr)
            else:
                sp_z = sp
            new_stages.append(sp_z)
            new_encoded.append(raw_to_encoded(sp_z))

        new_corners[cn] = CornerState(
            stages=new_stages, boost=corner.boost, _pre_encoded=new_encoded,
        )

    new_ca = CornerArray(
        a=new_corners[CornerName.A],
        b=new_corners[CornerName.B],
        c=new_corners[CornerName.C],
        d=new_corners[CornerName.D],
    )
    return Body(name=body.name, corners=new_ca, boost=body.boost)


def solve_c4_surface(body: Body, target_db: float = DEFAULT_TARGET_DB) -> Body:
    """Solve shared b0 (= c0) per corner to hit target cascade peak."""
    new_corners = {}
    for cn in CornerName:
        corner = body.corners.corner(cn)
        encoded = corner.encode()
        b0 = _solve_b0_for_corner(encoded, target_db=target_db)
        new_encoded = [_apply_shared_b0(enc, b0) for enc in encoded]
        new_corners[cn] = CornerState(
            stages=corner.stages, boost=corner.boost, _pre_encoded=new_encoded,
        )
    new_ca = CornerArray(
        a=new_corners[CornerName.A],
        b=new_corners[CornerName.B],
        c=new_corners[CornerName.C],
        d=new_corners[CornerName.D],
    )
    return Body(name=body.name, corners=new_ca, boost=body.boost)


def apply_constant_peak_compensation(body: Body) -> Body:
    """Apply (1-r) per-stage gain compensation to tame resonance explosions.

    For each active stage, scales the feedforward coefficients (c0, c1, c2)
    by (1-r) where r is the stage's pole radius (sqrt(c4)).
    This pins the peak gain at the resonant frequency while narrowing the
    bandwidth — high-Q stages carve inward instead of exploding upward.

    Applied AFTER solve_c4_surface so the shared b0 is already set.
    """
    new_corners = {}
    for cn in CornerName:
        corner = body.corners.corner(cn)
        encoded = corner.encode()
        compensated = []
        for enc in encoded:
            if _is_passthrough(enc):
                compensated.append(enc)
                continue
            r = math.sqrt(max(0.0, enc.c4))
            comp = 1.0 - r
            compensated.append(EncodedCoeffs(
                c0=enc.c0 * comp,
                c1=enc.c1 * comp,
                c2=enc.c2 * comp,
                c3=enc.c3,
                c4=enc.c4,
            ))
        new_corners[cn] = CornerState(
            stages=corner.stages, boost=corner.boost, _pre_encoded=compensated,
        )
    new_ca = CornerArray(
        a=new_corners[CornerName.A],
        b=new_corners[CornerName.B],
        c=new_corners[CornerName.C],
        d=new_corners[CornerName.D],
    )
    return Body(name=body.name, corners=new_ca, boost=body.boost)


# ---------------------------------------------------------------------------
# Cascade gain measurement
# ---------------------------------------------------------------------------

def cascade_peak_db(body: Body, num_freqs: int = 1024) -> dict:
    """Compute peak cascade gain (dB) for each corner.

    Evaluates the actual DF2T transfer function:
    H(z) = (c0 + c1*z^-1 + c2*z^-2) / (1 + c3*z^-1 + c4*z^-2)
    """
    freqs = np.linspace(20.0, SR / 2.0 - 1.0, num_freqs)
    z_inv = np.exp(-2j * np.pi * freqs / SR)
    z_inv2 = z_inv * z_inv

    results = {}
    for cn in CornerName:
        corner = body.corners.corner(cn)
        encoded = corner.encode()
        total = np.ones(num_freqs)
        for enc in encoded:
            if _is_passthrough(enc):
                continue
            num = enc.c0 + enc.c1 * z_inv + enc.c2 * z_inv2
            den = 1.0 + enc.c3 * z_inv + enc.c4 * z_inv2
            H = np.abs(num / np.where(np.abs(den) > 1e-30, den, 1e-30))
            total *= H
        peak = float(np.max(total))
        results[cn.json_key()] = 20.0 * math.log10(max(peak, 1e-30))
    return results


def apply_gain_budget(body: Body, budget_db: float = GAIN_BUDGET_DB) -> Body:
    """If cascade peak exceeds budget, reduce boost."""
    peaks = cascade_peak_db(body)
    worst_db = max(peaks.values())
    if worst_db <= budget_db:
        return body
    excess_linear = 10.0 ** ((worst_db - budget_db) / 20.0)
    return Body(name=body.name, corners=body.corners, boost=body.boost / excess_linear)


# ---------------------------------------------------------------------------
# Grid → compiled Body
# ---------------------------------------------------------------------------

def generate_body(
    archetype: str,
    seed: int = 42,
    boost: float = 4.0,
    target_db: float = DEFAULT_TARGET_DB,
    solve: bool = True,
) -> Body:
    """Generate a Body from a named archetype."""
    bp = BLUEPRINTS[archetype]
    rng = random.Random(seed)
    body = blueprint_to_body(bp, rng, boost=boost)
    if solve:
        body = solve_c4_surface(body, target_db=target_db)
        body = apply_gain_budget(body)
    return body


def generate_batch(
    seed: int = 42,
    boost: float = 4.0,
    target_db: float = DEFAULT_TARGET_DB,
    solve: bool = True,
) -> list[tuple[str, Body]]:
    """Generate one body per archetype."""
    bodies = []
    for i, (key, bp) in enumerate(BLUEPRINTS.items()):
        rng = random.Random(seed + i)
        body = blueprint_to_body(bp, rng, boost=boost)
        if solve:
            body = solve_c4_surface(body, target_db=target_db)
            body = apply_gain_budget(body)
        bodies.append((key, body))
    return bodies


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

def _report_body(body: Body):
    peaks = cascade_peak_db(body)
    worst = max(peaks.values())
    # Solved c4 per corner
    c4s = {}
    for cn in CornerName:
        encoded = body.corners.corner(cn).encode()
        for enc in encoded:
            if not _is_passthrough(enc):
                c4s[cn.json_key()] = enc.c4
                break
    c4_str = "  ".join(f"{k}={v:.4f}" for k, v in c4s.items())
    peak_str = "  ".join(f"{k}={v:+.1f}" for k, v in peaks.items())
    print(f"    c4:   {c4_str}")
    print(f"    peak: {peak_str}  worst={worst:+.1f} dB")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="TRENCH Anchor+Vector body generator"
    )
    parser.add_argument(
        "--archetype", type=str, choices=list(BLUEPRINTS.keys()),
        help="Generate a single named archetype",
    )
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--boost", type=float, default=4.0)
    parser.add_argument("--target-db", type=float, default=DEFAULT_TARGET_DB,
                        help=f"Target cascade peak dB (default {DEFAULT_TARGET_DB})")
    parser.add_argument("--out", type=str, default="vault/sift_queue")
    parser.add_argument("--no-solve", action="store_true",
                        help="Skip c4 solver (raw heritage compiler output)")
    args = parser.parse_args()

    do_solve = not args.no_solve

    if args.archetype:
        body = generate_body(args.archetype, seed=args.seed,
                             boost=args.boost, target_db=args.target_db,
                             solve=do_solve)
        os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
        with open(args.out, "w") as f:
            f.write(body.to_compiled_json(provenance="forge-generator"))
        print(f"  {body.name}")
        _report_body(body)
    else:
        os.makedirs(args.out, exist_ok=True)
        bodies = generate_batch(seed=args.seed, boost=args.boost,
                                target_db=args.target_db, solve=do_solve)
        for i, (key, body) in enumerate(bodies):
            fname = f"{i:03d}_{body.name.replace(' ', '_')}.json"
            path = os.path.join(args.out, fname)
            with open(path, "w") as f:
                f.write(body.to_compiled_json(provenance="forge-generator"))
            print(f"  {i:03d}. {body.name}")
            _report_body(body)

        print(f"\nGenerated {len(bodies)} bodies in {args.out}/")


if __name__ == "__main__":
    main()
