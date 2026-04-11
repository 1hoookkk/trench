"""Forge shipping — P2K direct-synthesis optimizer for the 4 flagship bodies.

Bypasses the heritage 0-127 integer compiler entirely. Searches continuous
Hz/radius/zero space with pymoo NSGA-II, seeded from acoustic formant targets.

Each body has a 6-stage allocation with roles:
  - Anchor: locked frequency, immune to morph
  - Character: sweeps aggressively with morph
  - Boundary: tilt/shelf shaping

Usage:
    python pyruntime/forge_shipping.py --body small_talk
    python pyruntime/forge_shipping.py --body speaker_knockerz
    python pyruntime/forge_shipping.py --body all
    python pyruntime/forge_shipping.py --body small_talk --generations 100 --pop 60
"""
from __future__ import annotations

import argparse
import json
import math
import os
from datetime import datetime
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.core.problem import Problem
from pymoo.optimize import minimize
from pymoo.operators.crossover.sbx import SBX
from pymoo.operators.mutation.pm import PM
from pymoo.operators.sampling.rnd import FloatRandomSampling

from pyruntime.body import Body
from pyruntime.constants import SR, TWO_PI, NUM_BODY_STAGES
from pyruntime.corner import CornerArray, CornerName, CornerState
from pyruntime.stage_params import StageParams
from pyruntime.encode import raw_to_encoded, EncodedCoeffs
from pyruntime.freq_response import cascade_response_db, freq_points
from pyruntime.forge_generator import solve_c4_surface, apply_gain_budget, cascade_peak_db

VAULT_DIR = Path(__file__).parent.parent / "vault"
FREQS = freq_points(sr=SR)
_VOCAL_MASK = (FREQS >= 200) & (FREQS <= 5000)
_VOCAL_FREQS = FREQS[_VOCAL_MASK]


# =============================================================================
# Stage role constants
# =============================================================================

ROLE_ANCHOR = "anchor"       # frequency locked across morph
ROLE_CHARACTER = "character"  # sweeps with morph
ROLE_BOUNDARY = "boundary"   # tilt/shelf, lower radius


# =============================================================================
# Shipping body definitions — acoustic targets + stage allocation
# =============================================================================

from dataclasses import dataclass, field


@dataclass(frozen=True)
class StageTarget:
    """Per-stage acoustic target for one corner."""
    freq_hz: float
    radius: float
    role: str
    zero_type: str = "allpole"  # allpole, unit_circle, interior


@dataclass(frozen=True)
class ShippingBodyDef:
    """Complete shipping body definition with 4-corner acoustic targets."""
    name: str
    key: str
    m0_q0: list[StageTarget]     # Corner A: morph=0, q=0
    m100_q0: list[StageTarget]   # Corner C: morph=1, q=0
    q_compression: float         # 0-1: how much Q compresses freqs toward center
    q_center_hz: float           # compression target frequency
    q_radius_mod: float          # radius delta for Q corners
    min_radius: float = 0.990    # minimum pole radius for active stages
    target_db: float = 36.0      # cascade peak target for c4 solver


# =============================================================================
# FROZEN SPECS — do not modify during optimization runs
# =============================================================================
# All radii chosen so per-stage gain comp (b0=1-r²) keeps individual
# stage peaks near 0 dB. Cascade peak determined by stage clustering,
# not by which stage has the highest radius.
# =============================================================================

# ---------------------------------------------------------------------------
# BODY 1: Small Talk Ah-Ee
# ---------------------------------------------------------------------------
# A throat transitioning from open "Ah" to tight "Ee."
# Two dominant formant peaks always present, always moving.
# F1 drops, F2 rises. Q compresses toward nasal center.
# Stages 0-1: F1/F2 (the vocal identity — these carry the vowel)
# Stages 2-3: F3/F4 (brightness/nasality — supporting color)
# Stages 4-5: Tilt boundaries (low warmth + HF rolloff)

SMALL_TALK = ShippingBodyDef(
    name="Small Talk Ah-Ee",
    key="small_talk",
    m0_q0=[
        StageTarget(768,  0.985, ROLE_CHARACTER),   # F1 "Ah"
        StageTarget(1189, 0.985, ROLE_CHARACTER),   # F2 "Ah"
        StageTarget(2555, 0.983, ROLE_CHARACTER),   # F3
        StageTarget(3508, 0.980, ROLE_CHARACTER),   # F4
        StageTarget(150,  0.970, ROLE_BOUNDARY),    # low warmth
        StageTarget(8000, 0.950, ROLE_BOUNDARY),    # HF rolloff
    ],
    m100_q0=[
        StageTarget(342,  0.985, ROLE_CHARACTER),   # F1 "Ee"
        StageTarget(2322, 0.985, ROLE_CHARACTER),   # F2 "Ee"
        StageTarget(3000, 0.983, ROLE_CHARACTER),   # F3
        StageTarget(3657, 0.980, ROLE_CHARACTER),   # F4
        StageTarget(150,  0.970, ROLE_BOUNDARY),    # low warmth (stationary)
        StageTarget(8000, 0.950, ROLE_BOUNDARY),    # HF rolloff (stationary)
    ],
    q_compression=0.5,
    q_center_hz=1200.0,
    q_radius_mod=0.003,
    min_radius=0.970,
)


# ---------------------------------------------------------------------------
# BODY 2: Speaker Knockerz
# ---------------------------------------------------------------------------
# Sub stays. Mids fan upward and sharpen. r ≤ 0.985 everywhere.

SPEAKER_KNOCKERZ = ShippingBodyDef(
    name="Speaker Knockerz",
    key="speaker_knockerz",
    m0_q0=[
        StageTarget(45,   0.975, ROLE_ANCHOR),       # sub anchor
        StageTarget(120,  0.965, ROLE_CHARACTER),     # chest bloom
        StageTarget(250,  0.963, ROLE_CHARACTER),     # choke seed
        StageTarget(500,  0.960, ROLE_CHARACTER),     # rip seed
        StageTarget(1000, 0.958, ROLE_CHARACTER),     # rattle seed
        StageTarget(8000, 0.940, ROLE_BOUNDARY),      # HF ceiling
    ],
    m100_q0=[
        StageTarget(45,   0.975, ROLE_ANCHOR),        # sub — locked
        StageTarget(400,  0.985, ROLE_CHARACTER),      # cardboard rip
        StageTarget(800,  0.985, ROLE_CHARACTER),      # rattle
        StageTarget(1200, 0.985, ROLE_CHARACTER),      # cone cry
        StageTarget(2500, 0.983, ROLE_CHARACTER),      # upper fracture
        StageTarget(8000, 0.940, ROLE_BOUNDARY),       # HF ceiling
    ],
    q_compression=0.5,
    q_center_hz=400.0,
    q_radius_mod=0.003,
    min_radius=0.955,
)


# ---------------------------------------------------------------------------
# BODY 3: Aluminum Siding
# ---------------------------------------------------------------------------
# 1kHz void locked. HF peaks sharpen with morph. r ≤ 0.985 everywhere.

ALUMINUM_SIDING = ShippingBodyDef(
    name="Aluminum Siding",
    key="aluminum_siding",
    # No active notch — the 1kHz void is structural: nothing resonates there.
    # All character energy lives above 5kHz, boundary at 300Hz.
    m0_q0=[
        StageTarget(300,  0.960, ROLE_BOUNDARY),      # low anchor — everything below is dark
        StageTarget(5000, 0.975, ROLE_CHARACTER),     # dull sheen
        StageTarget(7000, 0.973, ROLE_CHARACTER),     # glass seed
        StageTarget(9000, 0.970, ROLE_CHARACTER),     # sibilant seed
        StageTarget(11000, 0.968, ROLE_CHARACTER),    # tear seed
        StageTarget(14000, 0.960, ROLE_BOUNDARY),     # air ceiling
    ],
    m100_q0=[
        StageTarget(300,  0.960, ROLE_BOUNDARY),      # low anchor (stationary)
        StageTarget(7000, 0.975, ROLE_CHARACTER),     # glass stress
        StageTarget(10000, 0.975, ROLE_CHARACTER),    # aluminium tear
        StageTarget(12000, 0.973, ROLE_CHARACTER),    # dog whistle
        StageTarget(15000, 0.970, ROLE_CHARACTER),    # shatter point
        StageTarget(17000, 0.950, ROLE_BOUNDARY),     # air ceiling
    ],
    q_compression=0.4,
    q_center_hz=9000.0,
    q_radius_mod=0.003,
    min_radius=0.958,
)


# ---------------------------------------------------------------------------
# BODY 4: Cul-De-Sac
# ---------------------------------------------------------------------------
# Clustered pipe → scattered comb. r ≤ 0.985 everywhere.

CUL_DE_SAC = ShippingBodyDef(
    name="Cul-De-Sac",
    key="cul_de_sac",
    m0_q0=[
        StageTarget(100,  0.973, ROLE_ANCHOR),        # hum anchor
        StageTarget(300,  0.965, ROLE_CHARACTER),     # cluster 1
        StageTarget(400,  0.963, ROLE_CHARACTER),     # cluster 2
        StageTarget(500,  0.963, ROLE_CHARACTER),     # cluster 3
        StageTarget(650,  0.960, ROLE_CHARACTER),     # cluster 4
        StageTarget(800,  0.958, ROLE_CHARACTER),     # cluster 5
    ],
    m100_q0=[
        StageTarget(100,  0.973, ROLE_ANCHOR),         # hum — LOCKED
        StageTarget(1500, 0.985, ROLE_CHARACTER),     # shard 1
        StageTarget(3500, 0.985, ROLE_CHARACTER),     # shard 2
        StageTarget(6000, 0.983, ROLE_CHARACTER),     # shard 3
        StageTarget(10000, 0.980, ROLE_CHARACTER),    # shard 4
        StageTarget(14000, 0.978, ROLE_CHARACTER),    # shard 5
    ],
    q_compression=0.6,
    q_center_hz=2000.0,
    q_radius_mod=0.004,
    min_radius=0.955,
)


ALL_BODIES = {
    "small_talk": SMALL_TALK,
    "speaker_knockerz": SPEAKER_KNOCKERZ,
    "aluminum_siding": ALUMINUM_SIDING,
    "cul_de_sac": CUL_DE_SAC,
}


# =============================================================================
# Builder: definition → Body
# =============================================================================

def _hz_to_params(freq_hz: float, radius: float, zero_type: str = "allpole") -> StageParams:
    """Convert Hz/radius to StageParams with specified zero placement."""
    freq_hz = max(20.0, min(SR / 2.0 - 1.0, freq_hz))
    radius = max(0.0, min(0.999, radius))
    theta = TWO_PI * freq_hz / SR
    a1 = -2.0 * radius * math.cos(theta)

    if zero_type == "unit_circle":
        # Zero on unit circle at pole frequency — deep notch
        val2 = -2.0 * math.cos(theta) - a1
        val3 = radius * radius - 1.0
    elif zero_type == "interior":
        # Interior zero at 60% of pole radius
        zr = radius * 0.6
        val2 = -2.0 * zr * math.cos(theta) - a1
        val3 = radius * radius - zr * zr
    elif zero_type == "bandpass":
        # Zero at DC (z=1) — kills DC buildup for low-frequency poles.
        # b0 = gain_comp, b1 = -2*b0, b2 = b0. H(DC) = 0.
        b0 = 1.0 - radius * radius
        val2 = -2.0 * b0 - a1   # b1 = a1 + val2 → val2 = b1 - a1 = -2*b0 - a1
        val3 = radius * radius - b0  # b2 = r² - val3 → val3 = r² - b0
        # val1 is set below via gain comp, which gives b0 = 1-r² = b0. Consistent.
    else:
        # All-pole: zeros at origin — maximum resonance
        val2 = -a1
        val3 = radius * radius

    # Per-stage gain compensation: b0 = (1 - r²) so peak ≈ 0 dB regardless of Q.
    # Skip for unit_circle (notch) — gain comp inverts the notch behavior.
    if zero_type != "unit_circle":
        val1 = (1.0 - radius * radius) - 1.0  # val1 = -r², so b0 = 1 + val1 = 1 - r²
    else:
        val1 = 0.0  # notch stage: no gain comp, b0 = 1.0
    return StageParams(a1=a1, r=radius, val1=val1, val2=val2, val3=val3)


def _build_corner(targets: list[StageTarget]) -> CornerState:
    """Build a corner from stage targets."""
    stages = []
    encoded = []
    for t in targets:
        sp = _hz_to_params(t.freq_hz, t.radius, t.zero_type)
        stages.append(sp)
        encoded.append(raw_to_encoded(sp))
    # Pad to 12 with passthrough
    while len(stages) < NUM_BODY_STAGES:
        stages.append(StageParams.passthrough())
        encoded.append(EncodedCoeffs(c0=1.0, c1=0.0, c2=0.0, c3=0.0, c4=0.0))
    return CornerState(stages=stages, boost=4.0, _pre_encoded=encoded)


def _q_corner(base: list[StageTarget], comp: float, center: float, r_mod: float) -> list[StageTarget]:
    """Derive a Q corner by compressing frequencies toward center and modifying radii."""
    result = []
    for t in base:
        if t.role == ROLE_ANCHOR:
            # Anchors don't move with Q either
            result.append(t)
        else:
            hz_q = t.freq_hz + (center - t.freq_hz) * comp
            hz_q = max(20.0, min(SR / 2.0 - 1.0, hz_q))
            r_q = min(0.999, t.radius + r_mod)
            result.append(StageTarget(hz_q, r_q, t.role, t.zero_type))
    return result


def build_body_from_def(bdef: ShippingBodyDef) -> Body:
    """Build a complete 4-corner Body from a shipping body definition."""
    # Corner A: M0_Q0
    corner_a = _build_corner(bdef.m0_q0)
    # Corner C: M100_Q0
    corner_c = _build_corner(bdef.m100_q0)
    # Corner B: M0_Q100 (Q-compressed version of A)
    m0_q100 = _q_corner(bdef.m0_q0, bdef.q_compression, bdef.q_center_hz, bdef.q_radius_mod)
    corner_b = _build_corner(m0_q100)
    # Corner D: M100_Q100 (Q-compressed version of C)
    m100_q100 = _q_corner(bdef.m100_q0, bdef.q_compression, bdef.q_center_hz, bdef.q_radius_mod)
    corner_d = _build_corner(m100_q100)

    corners = CornerArray(a=corner_a, b=corner_b, c=corner_c, d=corner_d)
    return Body(name=bdef.name, corners=corners, boost=4.0)


# =============================================================================
# Fitness evaluation
# =============================================================================

_OCC_STEPS = [i / 10.0 for i in range(11)]  # 11x11 = 121 points across morph surface


def _extract_formants(db: np.ndarray, n: int = 3) -> list[float]:
    """Return up to n formant frequencies (Hz) from the vocal band."""
    vocal_db = db[_VOCAL_MASK]
    if len(vocal_db) < 3:
        return []
    mean_db = float(np.mean(db))
    peaks: list[tuple[float, float]] = []
    for i in range(1, len(vocal_db) - 1):
        if vocal_db[i] > vocal_db[i - 1] and vocal_db[i] > vocal_db[i + 1]:
            if vocal_db[i] > mean_db + 3.0:
                peaks.append((float(_VOCAL_FREQS[i]), float(vocal_db[i])))
    peaks.sort(key=lambda p: -p[1])
    return sorted(p[0] for p in peaks[:n])


def _trajectory_score(tracks: list[list[float]]) -> float:
    """Score vowel-path quality from per-sweep-point formant lists."""
    valid = [ft for ft in tracks if len(ft) >= 2]
    if len(valid) < 3:
        return 0.0
    f1 = [ft[0] for ft in valid]
    f2 = [ft[1] for ft in valid]
    range_score = min(1.0, (max(f1) - min(f1) + max(f2) - min(f2)) / 1500.0)
    hits = sum(1 for a, b in zip(f1, f2) if 200 < a < 1000 and 500 < b < 3000)
    vowel_score = hits / len(valid)
    arc = 0.0
    for i in range(1, len(valid)):
        d1 = math.log(f1[i] / f1[i - 1]) if f1[i] > 0 and f1[i - 1] > 0 else 0.0
        d2 = math.log(f2[i] / f2[i - 1]) if f2[i] > 0 and f2[i - 1] > 0 else 0.0
        arc += math.sqrt(d1 * d1 + d2 * d2)
    arc_score = min(1.0, arc / 2.0)
    return 0.4 * range_score + 0.3 * vowel_score + 0.3 * arc_score


def _continuity_score(tracks: list[list[float]]) -> float:
    """Reward smooth formant migration, penalise jumps."""
    valid = [ft for ft in tracks if len(ft) >= 2]
    if len(valid) < 2:
        return 0.0
    n_f = min(2, min(len(ft) for ft in valid))
    penalty = 0.0
    for fi in range(n_f):
        trk = [ft[fi] for ft in valid]
        for i in range(1, len(trk)):
            lo, hi = min(trk[i], trk[i - 1]), max(trk[i], trk[i - 1])
            ratio = hi / max(1.0, lo)
            if ratio > 1.5:
                penalty += ratio - 1.0
    return max(0.0, min(1.0, 1.0 - penalty / 3.0))


def _talkingness(db: np.ndarray) -> float:
    """Vocal peak count in 200-5000 Hz, normalised to [0, 1]."""
    vocal_db = db[_VOCAL_MASK]
    if len(vocal_db) < 3:
        return 0.0
    mean_db = float(np.mean(db))
    peaks = 0
    for i in range(1, len(vocal_db) - 1):
        if vocal_db[i] > vocal_db[i - 1] and vocal_db[i] > vocal_db[i + 1]:
            if vocal_db[i] > mean_db + 3.0:
                peaks += 1
    return min(1.0, peaks / 4.0)


def eval_shipping_body(body: Body, category: str = "vocal") -> dict:
    """Layered fitness evaluation for shipping bodies.

    L1: hard gate (peak > 50 dB, dead midpoint)
    L2: occupancy (fraction of 5x5 bins alive)
    L3: trajectory (formant path quality)
    L4: continuity (smooth migration)
    L5: morph distance (RMS spectral change)
    L6: category-specific metrics
    """
    responses: dict[tuple[float, float], np.ndarray] = {}
    for mi in _OCC_STEPS:
        for qi in _OCC_STEPS:
            enc = body.corners.interpolate(mi, qi)
            responses[(mi, qi)] = cascade_response_db(enc, FREQS, SR)

    # L1: hard gate — 40 dB ceiling across ENTIRE 11x11 surface (121 evals)
    max_peak = max(float(np.max(db)) for db in responses.values())
    db_mid = responses[(0.5, 0.5)]
    mid_talk = _talkingness(db_mid)
    gate_fail = max_peak > 40.0

    # L2: occupancy
    n_bins = len(responses)
    talk_grid = {k: _talkingness(db) for k, db in responses.items()}
    dead = sum(1 for t in talk_grid.values() if t == 0.0)
    occupancy = (n_bins - dead) / float(n_bins)

    # L3/L4: trajectory + continuity (Q=0.5 morph sweep)
    sweep_keys = [(m, 0.5) for m in _OCC_STEPS]
    formant_tracks = [_extract_formants(responses[k]) for k in sweep_keys]
    trajectory = _trajectory_score(formant_tracks)
    continuity = _continuity_score(formant_tracks)

    # L5: morph distance
    morph_dist = float(0.5 * (
        np.sqrt(np.mean((responses[(0.0, 0.0)] - responses[(1.0, 0.0)]) ** 2))
        + np.sqrt(np.mean((responses[(0.0, 1.0)] - responses[(1.0, 1.0)]) ** 2))
    ))

    # L6: category-specific
    ridge_vals = []
    for db in responses.values():
        pk = float(np.max(db))
        mn = float(np.mean(db))
        ridge_vals.append(max(0.0, min(1.0, (pk - mn) / 40.0)))
    ridge = float(np.mean(ridge_vals))

    dynamic_range = float(np.max(db_mid) - np.min(db_mid))
    ruggedness = min(1.0, float(np.std(db_mid)) / 30.0)

    tvs = list(talk_grid.values())
    alive = [t for t in tvs if t > 0]
    talkingness_mean = float(np.mean(alive)) if alive else 0.0

    return {
        "gate_fail": gate_fail,
        "occupancy": occupancy,
        "trajectory": trajectory,
        "continuity": continuity,
        "morph_distance": morph_dist,
        "ridge": ridge,
        "talkingness": talkingness_mean,
        "dynamic_range": dynamic_range,
        "ruggedness": ruggedness,
        "peak_db": max_peak,
    }


# =============================================================================
# Pymoo optimizer — searches Hz/radius/zero space around seeded targets
# =============================================================================

def _unpack_body(params: np.ndarray, bdef: ShippingBodyDef) -> Body:
    """Decode optimizer parameter vector into a Body.

    Parameter layout (6 active stages):
      [0:6]   anchor freq offsets (Hz delta from seed)
      [6:12]  anchor radius offsets (delta from seed)
      [12:18] morph target freq offsets (Hz delta from seed)
      [18:24] morph target radius offsets (delta from seed)
      [24:30] zero placement offsets per stage (val2 scale factor)
      [30:36] zero depth per stage (val3 scale factor)
      [36]    Q compression factor
      [37]    Q center offset (Hz)
      [38]    Q radius mod
    Total: 39 params
    """
    n = len(bdef.m0_q0)

    # Decode M0_Q0 corner
    m0_targets = []
    for i in range(n):
        seed = bdef.m0_q0[i]
        freq = seed.freq_hz + params[i]
        freq = max(20.0, min(SR / 2.0 - 1.0, freq))
        radius = seed.radius + params[6 + i]
        radius = max(0.5, min(0.999, radius))
        # Enforce minimum radius for active stages
        if seed.role != ROLE_BOUNDARY:
            radius = max(bdef.min_radius, radius)
        m0_targets.append(StageTarget(freq, radius, seed.role, seed.zero_type))

    # Decode M100_Q0 corner
    m100_targets = []
    for i in range(n):
        seed = bdef.m100_q0[i]
        if seed.role == ROLE_ANCHOR:
            # Anchor stages are locked — use exact seed values
            m100_targets.append(seed)
        else:
            freq = seed.freq_hz + params[12 + i]
            freq = max(20.0, min(SR / 2.0 - 1.0, freq))
            radius = seed.radius + params[18 + i]
            radius = max(bdef.min_radius, min(0.999, radius))
            m100_targets.append(StageTarget(freq, radius, seed.role, seed.zero_type))

    # Build corners
    corner_a = _build_corner(m0_targets)
    corner_c = _build_corner(m100_targets)

    # Q corners
    q_comp = max(0.0, min(1.0, params[36]))
    q_center = bdef.q_center_hz + params[37]
    q_center = max(100.0, min(15000.0, q_center))
    q_r_mod = max(-0.01, min(0.01, params[38]))

    m0_q100 = _q_corner(m0_targets, q_comp, q_center, q_r_mod)
    m100_q100 = _q_corner(m100_targets, q_comp, q_center, q_r_mod)
    corner_b = _build_corner(m0_q100)
    corner_d = _build_corner(m100_q100)

    corners = CornerArray(a=corner_a, b=corner_b, c=corner_c, d=corner_d)
    return Body(name=bdef.name, corners=corners, boost=4.0)


# Category-specific objective sets
BODY_OBJECTIVES = {
    "small_talk":       ["trajectory", "continuity", "morph_distance"],
    "speaker_knockerz": ["dynamic_range", "morph_distance", "ridge"],
    "aluminum_siding":  ["ridge", "morph_distance", "ruggedness"],
    "cul_de_sac":       ["morph_distance", "ridge", "dynamic_range"],
}


class ShippingProblem(Problem):
    """Pymoo problem for shipping body optimization.

    Searches around seeded acoustic targets with bounded perturbations.
    Anchor stages are locked. Character stages have wider search bounds.
    """

    def __init__(self, bdef: ShippingBodyDef):
        self.bdef = bdef
        n = len(bdef.m0_q0)
        objectives = BODY_OBJECTIVES[bdef.key]
        self.objectives = objectives

        # Build bounds: [freq_offsets, radius_offsets, morph_freq, morph_r, zero, Q]
        xl, xu = [], []

        # Anchor freq offsets [0:6]
        for i in range(n):
            if bdef.m0_q0[i].role == ROLE_ANCHOR:
                xl.append(0.0); xu.append(0.0)  # locked
            elif bdef.m0_q0[i].role == ROLE_BOUNDARY:
                xl.append(-200.0); xu.append(200.0)
            else:
                xl.append(-150.0); xu.append(150.0)

        # Anchor radius offsets [6:12]
        for i in range(n):
            if bdef.m0_q0[i].role == ROLE_ANCHOR:
                xl.append(0.0); xu.append(0.0)
            else:
                xl.append(-0.008); xu.append(0.008)

        # Morph target freq offsets [12:18]
        for i in range(n):
            if bdef.m100_q0[i].role == ROLE_ANCHOR:
                xl.append(0.0); xu.append(0.0)
            elif bdef.m100_q0[i].role == ROLE_BOUNDARY:
                xl.append(-300.0); xu.append(300.0)
            else:
                xl.append(-250.0); xu.append(250.0)

        # Morph target radius offsets [18:24]
        for i in range(n):
            if bdef.m100_q0[i].role == ROLE_ANCHOR:
                xl.append(0.0); xu.append(0.0)
            else:
                xl.append(-0.008); xu.append(0.008)

        # Zero placement scale [24:30] (unused for now, reserved)
        for _ in range(n):
            xl.append(-0.5); xu.append(0.5)

        # Zero depth scale [30:36] (unused for now, reserved)
        for _ in range(n):
            xl.append(-0.5); xu.append(0.5)

        # Q params [36:39]
        xl.append(max(0.0, bdef.q_compression - 0.3))
        xu.append(min(1.0, bdef.q_compression + 0.3))
        xl.append(-500.0); xu.append(500.0)  # Q center offset
        xl.append(-0.005); xu.append(0.008)   # Q radius mod

        super().__init__(
            n_var=39,
            n_obj=len(objectives),
            n_ieq_constr=2,
            xl=np.array(xl),
            xu=np.array(xu),
        )

    def _evaluate(self, X, out, *args, **kwargs):
        F, G = [], []
        for x in X:
            try:
                body = _unpack_body(x, self.bdef)
                body = solve_c4_surface(body, target_db=self.bdef.target_db)
                body = apply_gain_budget(body)
                metrics = eval_shipping_body(body, self.bdef.key)
                # Negate for minimization
                row = [-metrics[obj] for obj in self.objectives]
                F.append(row)
                # Constraints: gate fail, minimum occupancy
                g1 = 1.0 if metrics["gate_fail"] else -1.0
                g2 = 0.4 - metrics["occupancy"]
                G.append([g1, g2])
            except Exception:
                F.append([0.0] * len(self.objectives))
                G.append([1.0, 1.0])
        out["F"] = np.array(F)
        out["G"] = np.array(G)


# =============================================================================
# Plotting
# =============================================================================

def _plot_morph_sweep(body: Body, out_path: Path, n_steps: int = 11):
    """Plot frequency response at multiple morph positions (Q=0.5)."""
    fig, ax = plt.subplots(figsize=(12, 6))
    cmap = plt.cm.coolwarm
    for i in range(n_steps):
        morph = i / (n_steps - 1)
        enc = body.corners.interpolate(morph, 0.5)
        db = cascade_response_db(enc, FREQS, SR)
        color = cmap(morph)
        label = f"M={morph:.1f}" if i % 2 == 0 else None
        ax.semilogx(FREQS, db, color=color, alpha=0.7, linewidth=1.5, label=label)
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Magnitude (dB)")
    ax.set_title(f"{body.name} — Morph Sweep (Q=0.5)")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(20, SR / 2)
    ax.set_ylim(-40, 50)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()


def _plot_q_sweep(body: Body, out_path: Path, n_steps: int = 7):
    """Plot frequency response at multiple Q positions (Morph=0.5)."""
    fig, ax = plt.subplots(figsize=(12, 6))
    cmap = plt.cm.viridis
    for i in range(n_steps):
        q = i / (n_steps - 1)
        enc = body.corners.interpolate(0.5, q)
        db = cascade_response_db(enc, FREQS, SR)
        color = cmap(q)
        ax.semilogx(FREQS, db, color=color, alpha=0.7, linewidth=1.5, label=f"Q={q:.2f}")
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Magnitude (dB)")
    ax.set_title(f"{body.name} — Q Sweep (Morph=0.5)")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(20, SR / 2)
    ax.set_ylim(-40, 50)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()


def _plot_corners(body: Body, out_path: Path):
    """Plot all 4 corners overlaid."""
    fig, ax = plt.subplots(figsize=(12, 6))
    colors = {"M0_Q0": "blue", "M0_Q100": "cyan", "M100_Q0": "red", "M100_Q100": "orange"}
    for cn in CornerName:
        enc = body.corners.corner(cn).encode()
        db = cascade_response_db(enc, FREQS, SR)
        ax.semilogx(FREQS, db, color=colors[cn.json_key()], linewidth=2,
                     label=cn.json_key(), alpha=0.8)
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Magnitude (dB)")
    ax.set_title(f"{body.name} — 4 Corners")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(20, SR / 2)
    ax.set_ylim(-40, 50)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()


def _plot_pareto(result, objectives, out_path):
    """Plot Pareto front."""
    F = -result.F
    n_obj = len(objectives)
    n_plots = min(3, n_obj - 1)
    fig, axes = plt.subplots(1, n_plots, figsize=(5 * n_plots, 4))
    if n_plots == 1:
        axes = [axes]
    for i, ax in enumerate(axes):
        if i + 1 < F.shape[1]:
            ax.scatter(F[:, 0], F[:, i + 1], c="steelblue", s=20, alpha=0.7)
            ax.set_xlabel(objectives[0])
            ax.set_ylabel(objectives[i + 1])
            ax.set_title("Pareto front")
            ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()


# =============================================================================
# Main optimizer entry point
# =============================================================================

def optimize_shipping_body(
    key: str,
    generations: int = 50,
    pop_size: int = 40,
    seed: int = 42,
) -> Path:
    """Run pymoo optimization for one shipping body."""
    bdef = ALL_BODIES[key]
    objectives = BODY_OBJECTIVES[key]

    print(f"\n{'='*60}")
    print(f"  SHIPPING BODY: {bdef.name}")
    print(f"  Pipeline: P2K Direct StageParams (no heritage compiler)")
    print(f"  Objectives: {objectives}")
    print(f"  Min radius: {bdef.min_radius}")
    print(f"  Generations: {generations}, Population: {pop_size}")
    print(f"{'='*60}")

    # First: build and evaluate the seed body (acoustic targets, no optimization)
    seed_body = build_body_from_def(bdef)
    seed_body = solve_c4_surface(seed_body, target_db=bdef.target_db)
    seed_body = apply_gain_budget(seed_body)
    seed_metrics = eval_shipping_body(seed_body, key)
    seed_peaks = cascade_peak_db(seed_body)

    print(f"\n  Seed body metrics:")
    for k, v in seed_metrics.items():
        if isinstance(v, float):
            print(f"    {k}: {v:.3f}")
        else:
            print(f"    {k}: {v}")
    print(f"  Seed peaks: {', '.join(f'{k}={v:+.1f}dB' for k, v in seed_peaks.items())}")

    # Run optimization
    problem = ShippingProblem(bdef)
    algorithm = NSGA2(
        pop_size=pop_size,
        sampling=FloatRandomSampling(),
        crossover=SBX(prob=0.9, eta=15),
        mutation=PM(eta=20),
    )
    result = minimize(problem, algorithm, ("n_gen", generations),
                      seed=seed, verbose=False)

    # Build and save Pareto-optimal bodies
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = VAULT_DIR / f"shipping_{key}_{ts}"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Save seed body first
    seed_path = out_dir / f"SEED_{bdef.name.replace(' ', '_')}.json"
    seed_path.write_text(seed_body.to_compiled_json(provenance="forge-shipping-seed"))
    print(f"\n  Seed body: {seed_path}")

    bodies = []
    names = []
    metrics_list = []

    if result.X is None:
        print("\n  WARNING: No feasible Pareto-optimal solutions found.")
        print("  All candidates violated constraints (gate_fail or low occupancy).")
        print(f"  Seed body saved to: {seed_path}")
        return out_dir

    for i, x in enumerate(result.X):
        body = _unpack_body(x, bdef)
        body = solve_c4_surface(body, target_db=bdef.target_db)
        body = apply_gain_budget(body)
        name = f"{key}_{i:03d}"
        body_named = Body(name=name, corners=body.corners, boost=body.boost)
        metrics = eval_shipping_body(body_named, key)

        path = out_dir / f"{name}.json"
        path.write_text(body_named.to_compiled_json(provenance="forge-shipping-opt"))

        bodies.append(body_named)
        names.append(name)
        metrics_list.append(metrics)

    # Rank by composite of all objectives
    def _rank_score(m):
        return sum(m.get(obj, 0.0) for obj in objectives)

    ranked = sorted(zip(names, metrics_list, bodies), key=lambda x: -_rank_score(x[1]))

    print(f"\n  {len(bodies)} Pareto-optimal bodies saved to {out_dir}")
    print(f"\n  Top 10 by composite score:")
    for name, m, _ in ranked[:10]:
        vals = "  ".join(f"{obj}={m[obj]:.3f}" for obj in objectives)
        print(f"    {name}: {vals}")

    # Generate plots for top candidate
    if ranked:
        best_name, best_metrics, best_body = ranked[0]
        _plot_morph_sweep(best_body, out_dir / "best_morph_sweep.png")
        _plot_q_sweep(best_body, out_dir / "best_q_sweep.png")
        _plot_corners(best_body, out_dir / "best_corners.png")
        _plot_pareto(result, objectives, out_dir / "pareto.png")
        print(f"\n  Plots saved to {out_dir}")

        # Save best as the promoted candidate
        best_path = out_dir / f"BEST_{bdef.name.replace(' ', '_')}.json"
        best_body_final = Body(name=bdef.name, corners=best_body.corners, boost=best_body.boost)
        best_path.write_text(best_body_final.to_compiled_json(provenance="forge-shipping-best"))
        print(f"  Best candidate: {best_path}")
        print(f"  Best metrics: {', '.join(f'{k}={v:.3f}' for k, v in best_metrics.items() if isinstance(v, float))}")

    return out_dir


# =============================================================================
# CLI
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description="Forge shipping body optimizer")
    parser.add_argument("--body", default="small_talk",
                        help="Body to optimize: small_talk, speaker_knockerz, aluminum_siding, cul_de_sac, all")
    parser.add_argument("--generations", type=int, default=50)
    parser.add_argument("--pop", type=int, default=40)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    keys = list(ALL_BODIES.keys()) if args.body == "all" else [args.body]

    for key in keys:
        if key not in ALL_BODIES:
            print(f"Unknown body: {key}. Available: {list(ALL_BODIES.keys())}")
            continue
        optimize_shipping_body(key, args.generations, args.pop, args.seed)


if __name__ == "__main__":
    main()
