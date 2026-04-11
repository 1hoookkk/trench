"""Trajectory-first forge — author stress narratives, solve for corners.

The human authors landmarks (what audible events happen where in the sweep)
and invariants (what must survive). The solver finds the 4-corner StageParams
that produce that trajectory in the actual runtime cascade.

Corners are an output format, not the authoring interface.

Usage:
    python pyruntime/forge_trajectory.py --body speaker_knockerz
    python pyruntime/forge_trajectory.py --body all --generations 80 --pop 60
"""
from __future__ import annotations

import argparse
import json
import math
import os
from datetime import datetime
from pathlib import Path
from dataclasses import dataclass, field

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
from pyruntime.forge_generator import solve_c4_surface, apply_gain_budget

VAULT_DIR = Path(__file__).parent.parent / "vault"
FREQS = freq_points(sr=SR)

# Hard limits
MAX_RADIUS = 0.985
GAIN_CEILING_DB = 40.0
B0_FLOOR = 0.01


# =============================================================================
# Trajectory spec — the authoring interface
# =============================================================================

@dataclass(frozen=True)
class Event:
    """An audible event at a specific morph position."""
    type: str        # "peak", "notch", "void", "shelf"
    freq_hz: float   # target center frequency
    bandwidth: float = 0.5  # octaves — how wide the match window is
    weight: float = 1.0     # importance in the fitness

@dataclass(frozen=True)
class Landmark:
    """What must happen at a specific morph position."""
    morph: float         # 0.0 - 1.0
    events: tuple[Event, ...]

@dataclass(frozen=True)
class Invariant:
    """A condition that must hold across the ENTIRE morph sweep."""
    type: str            # "sub_floor", "gain_ceiling", "midpoint_alive", "no_hf_takeover"
    freq_hz: float = 0.0
    threshold_db: float = 0.0

@dataclass(frozen=True)
class QBehavior:
    """How Q affects the trajectory."""
    compression_center_hz: float  # Q compresses formants toward this frequency
    compression_amount: float     # 0.0 = no effect, 1.0 = full collapse
    tightening: float             # how much Q narrows peaks (radius increase)

@dataclass(frozen=True)
class TrajectorySpec:
    """Complete body specification as a stress narrative."""
    name: str
    key: str
    landmarks: tuple[Landmark, ...]
    invariants: tuple[Invariant, ...]
    q_behavior: QBehavior
    target_db: float = 36.0


# =============================================================================
# The 4 shipping bodies — authored as trajectories
# =============================================================================

SPEAKER_KNOCKERZ = TrajectorySpec(
    name="Speaker Knockerz",
    key="speaker_knockerz",
    landmarks=(
        Landmark(0.0, (
            Event("peak", 45, 1.0, 1.5),      # sub: the vault — heavy, dark
            Event("shelf", 120, 2.0, 0.5),     # warm low-mid wash
        )),
        Landmark(0.2, (
            Event("peak", 45, 1.0, 1.0),       # sub stays
            Event("peak", 150, 0.8, 1.0),      # chest resonance emerges
        )),
        Landmark(0.4, (
            Event("peak", 45, 1.0, 0.8),       # sub stays
            Event("notch", 200, 0.5, 1.2),     # the choke — steep notch separates sub from overtones
            Event("peak", 350, 0.6, 1.0),      # pressure building above the choke
        )),
        Landmark(0.6, (
            Event("peak", 45, 1.0, 0.6),       # sub stays (quieter in the mix)
            Event("peak", 400, 0.5, 1.5),      # cardboard rip — peaky, unstable
        )),
        Landmark(0.8, (
            Event("peak", 45, 1.0, 0.5),       # sub stays
            Event("peak", 800, 0.4, 1.3),      # the rattle — comb-like flutter
            Event("peak", 1200, 0.5, 1.0),     # cone cry starting
        )),
        Landmark(1.0, (
            Event("peak", 45, 1.0, 0.5),       # sub stays
            Event("void", 200, 2.0, 1.0),      # low-mids completely gone
            Event("peak", 1200, 0.3, 1.5),     # cone cry — narrow scream
            Event("peak", 2500, 0.5, 0.8),     # upper fracture
        )),
    ),
    invariants=(
        Invariant("sub_floor", 45.0, 5.0),         # sub must stay above +5 dB at all morph
        Invariant("gain_ceiling", 0.0, GAIN_CEILING_DB),  # no point above 40 dB
        Invariant("midpoint_alive", 0.0, 0.0),     # m=0.5 must have spectral activity
    ),
    q_behavior=QBehavior(400.0, 0.5, 0.003),
)

SMALL_TALK = TrajectorySpec(
    name="Small Talk Ah-Ee",
    key="small_talk",
    landmarks=(
        Landmark(0.0, (
            Event("peak", 768, 0.4, 1.5),     # F1 "Ah" — open throat
            Event("peak", 1189, 0.4, 1.5),    # F2 "Ah"
            Event("peak", 2555, 0.5, 0.8),    # F3 brightness
        )),
        Landmark(0.2, (
            Event("peak", 660, 0.5, 1.2),     # F1 shifting toward "Eh"
            Event("peak", 1400, 0.5, 1.2),    # F2 rising
        )),
        Landmark(0.4, (
            Event("peak", 550, 0.5, 1.0),     # F1 dropping — "Eh" to "Ee" transition
            Event("peak", 1800, 0.5, 1.0),    # F2 climbing
        )),
        Landmark(0.6, (
            Event("peak", 450, 0.5, 1.2),     # F1 approaching "Ee"
            Event("peak", 2100, 0.4, 1.2),    # F2 nearing target
        )),
        Landmark(0.8, (
            Event("peak", 380, 0.5, 1.0),     # F1 nearly at "Ee"
            Event("peak", 2250, 0.4, 1.0),    # F2 nearly at "Ee"
        )),
        Landmark(1.0, (
            Event("peak", 342, 0.4, 1.5),     # F1 "Ee" — tight throat
            Event("peak", 2322, 0.4, 1.5),    # F2 "Ee"
            Event("peak", 3000, 0.5, 0.8),    # F3 brightness
        )),
    ),
    invariants=(
        Invariant("gain_ceiling", 0.0, GAIN_CEILING_DB),
        Invariant("midpoint_alive", 0.0, 0.0),
    ),
    q_behavior=QBehavior(1200.0, 0.5, 0.003),
)

ALUMINUM_SIDING = TrajectorySpec(
    name="Aluminum Siding",
    key="aluminum_siding",
    landmarks=(
        Landmark(0.0, (
            Event("void", 1000, 2.0, 1.5),    # permanent midrange void
            Event("shelf", 5000, 1.0, 0.8),   # dull, phasey HF sheen
        )),
        Landmark(0.2, (
            Event("void", 1000, 2.0, 1.2),    # void stays
            Event("peak", 7000, 0.5, 1.0),    # first sign of shimmer
        )),
        Landmark(0.4, (
            Event("void", 1000, 2.0, 1.0),    # void stays
            Event("peak", 7000, 0.4, 1.2),    # glass stress — sharpening
            Event("peak", 10000, 0.5, 0.8),   # second peak emerging
        )),
        Landmark(0.6, (
            Event("void", 1000, 2.0, 1.0),    # void stays
            Event("peak", 8000, 0.4, 1.3),    # sibilant fold
            Event("peak", 11000, 0.4, 1.0),   # aluminium tear beginning
        )),
        Landmark(0.8, (
            Event("void", 1000, 2.0, 0.8),    # void stays
            Event("peak", 10000, 0.3, 1.5),   # harsh metallic ringing
            Event("peak", 13000, 0.4, 1.0),   # dog whistle territory
        )),
        Landmark(1.0, (
            Event("void", 1000, 2.0, 1.5),    # void stays
            Event("peak", 12000, 0.3, 1.5),   # aluminium tear — maximum
            Event("peak", 15000, 0.4, 1.0),   # shatter point
        )),
    ),
    invariants=(
        Invariant("gain_ceiling", 0.0, GAIN_CEILING_DB),
        Invariant("midpoint_alive", 0.0, 0.0),
        Invariant("no_hf_takeover", 0.0, 0.0),     # HF must not exceed low-end by > 20 dB
    ),
    q_behavior=QBehavior(9000.0, 0.4, 0.003),
)

CUL_DE_SAC = TrajectorySpec(
    name="Cul-De-Sac",
    key="cul_de_sac",
    landmarks=(
        Landmark(0.0, (
            Event("peak", 100, 1.0, 0.8),     # hum anchor — the tether
            Event("peak", 450, 1.5, 1.5),      # thick clustered pipe resonance
        )),
        Landmark(0.2, (
            Event("peak", 100, 1.0, 0.6),     # hum stays
            Event("peak", 500, 1.2, 1.2),      # pipe starting to thin
            Event("notch", 300, 0.5, 0.8),     # notches appearing — the rust
        )),
        Landmark(0.4, (
            Event("peak", 100, 1.0, 0.5),     # hum stays
            Event("peak", 3000, 0.8, 1.3),     # the bulge — unnatural swell
        )),
        Landmark(0.5, (
            Event("peak", 100, 1.0, 0.5),     # hum barely there
            Event("void", 500, 3.0, 2.0),      # STATE BOUNDARY — the null. phase cancellation.
        )),
        Landmark(0.7, (
            Event("peak", 100, 1.0, 0.5),     # hum re-emerges
            Event("peak", 1500, 0.5, 1.0),     # fracture — peaks reappear scattered
            Event("peak", 4000, 0.5, 1.0),     # second comb tooth
            Event("peak", 8000, 0.5, 0.8),     # third comb tooth
        )),
        Landmark(1.0, (
            Event("peak", 100, 1.0, 0.5),     # hum anchor
            Event("peak", 1500, 0.3, 1.0),     # comb tooth 1
            Event("peak", 3500, 0.3, 1.0),     # comb tooth 2
            Event("peak", 6000, 0.3, 0.8),     # comb tooth 3
            Event("peak", 10000, 0.4, 0.6),    # comb tooth 4
        )),
    ),
    invariants=(
        Invariant("sub_floor", 100.0, 0.0),        # hum must stay present
        Invariant("gain_ceiling", 0.0, GAIN_CEILING_DB),
        Invariant("midpoint_alive", 0.0, 0.0),
    ),
    q_behavior=QBehavior(2000.0, 0.6, 0.004),
)

ALL_SPECS = {
    "speaker_knockerz": SPEAKER_KNOCKERZ,
    "small_talk": SMALL_TALK,
    "aluminum_siding": ALUMINUM_SIDING,
    "cul_de_sac": CUL_DE_SAC,
}


# =============================================================================
# Stage builder — Hz/radius to StageParams with gain comp
# =============================================================================

def _make_stage(freq_hz: float, radius: float) -> StageParams:
    """Build an all-pole resonator stage. r capped at MAX_RADIUS.

    No per-stage gain comp — the c4 solver normalizes the cascade globally.
    Resonance character comes from the pole radius directly.
    """
    freq_hz = max(20.0, min(SR / 2.0 - 1.0, freq_hz))
    radius = max(0.0, min(MAX_RADIUS, radius))
    theta = TWO_PI * freq_hz / SR
    a1 = -2.0 * radius * math.cos(theta)
    val1 = 0.0       # b0 = 1.0 — no gain comp
    val2 = -a1        # zeros at origin
    val3 = radius * radius
    return StageParams(a1=a1, r=radius, val1=val1, val2=val2, val3=val3)


def _build_corner_from_stages(stage_list: list[tuple[float, float]]) -> CornerState:
    """Build a corner from (freq_hz, radius) pairs."""
    stages = []
    encoded = []
    for freq, r in stage_list:
        sp = _make_stage(freq, r)
        stages.append(sp)
        encoded.append(raw_to_encoded(sp))
    while len(stages) < NUM_BODY_STAGES:
        stages.append(StageParams.passthrough())
        encoded.append(EncodedCoeffs(c0=1.0, c1=0.0, c2=0.0, c3=0.0, c4=0.0))
    return CornerState(stages=stages, boost=4.0, _pre_encoded=encoded)


# =============================================================================
# Event scoring — how well does a response match a landmark?
# =============================================================================

def _find_peak_near(db: np.ndarray, freqs: np.ndarray, target_hz: float, bw_oct: float) -> float:
    """Find the maximum dB within bw_oct octaves of target_hz. Returns dB."""
    lo = target_hz * 2.0 ** (-bw_oct)
    hi = target_hz * 2.0 ** (bw_oct)
    mask = (freqs >= lo) & (freqs <= hi)
    if not np.any(mask):
        return -100.0
    return float(np.max(db[mask]))


def _find_valley_near(db: np.ndarray, freqs: np.ndarray, target_hz: float, bw_oct: float) -> float:
    """Find the minimum dB within bw_oct octaves of target_hz. Returns dB."""
    lo = target_hz * 2.0 ** (-bw_oct)
    hi = target_hz * 2.0 ** (bw_oct)
    mask = (freqs >= lo) & (freqs <= hi)
    if not np.any(mask):
        return 0.0
    return float(np.min(db[mask]))


def _score_event(db: np.ndarray, freqs: np.ndarray, event: Event, mean_db: float) -> float:
    """Score how well a response matches an event. Returns 0.0 (miss) to 1.0 (hit)."""
    if event.type == "peak":
        peak_db = _find_peak_near(db, freqs, event.freq_hz, event.bandwidth)
        # A peak should be above the mean by at least 3 dB
        prominence = peak_db - mean_db
        return min(1.0, max(0.0, prominence / 10.0))  # 10 dB prominence = perfect

    elif event.type == "notch":
        valley_db = _find_valley_near(db, freqs, event.freq_hz, event.bandwidth)
        # A notch should be below the mean by at least 3 dB
        depth = mean_db - valley_db
        return min(1.0, max(0.0, depth / 10.0))

    elif event.type == "void":
        # A void: the entire band should be below the mean
        lo = event.freq_hz * 2.0 ** (-event.bandwidth)
        hi = event.freq_hz * 2.0 ** (event.bandwidth)
        mask = (freqs >= lo) & (freqs <= hi)
        if not np.any(mask):
            return 0.0
        band_db = float(np.mean(db[mask]))
        deficit = mean_db - band_db
        return min(1.0, max(0.0, deficit / 15.0))  # 15 dB deficit = perfect void

    elif event.type == "shelf":
        # A shelf: the band should be close to the mean (no sharp features)
        lo = event.freq_hz * 2.0 ** (-event.bandwidth)
        hi = event.freq_hz * 2.0 ** (event.bandwidth)
        mask = (freqs >= lo) & (freqs <= hi)
        if not np.any(mask):
            return 0.0
        band_std = float(np.std(db[mask]))
        return max(0.0, 1.0 - band_std / 10.0)  # low std = good shelf

    return 0.0


def _score_landmark(body: Body, landmark: Landmark) -> float:
    """Score how well a body matches a landmark at a specific morph position."""
    enc = body.corners.interpolate(landmark.morph, 0.5)
    db = cascade_response_db(enc, FREQS, SR)
    mean_db = float(np.mean(db))

    total_weight = sum(e.weight for e in landmark.events)
    if total_weight == 0:
        return 0.0

    score = 0.0
    for event in landmark.events:
        s = _score_event(db, FREQS, event, mean_db)
        score += s * event.weight
    return score / total_weight


def _check_invariants(body: Body, spec: TrajectorySpec) -> tuple[bool, list[str]]:
    """Check all invariants across the morph sweep. Returns (pass, violations)."""
    violations = []
    n_steps = 21

    for inv in spec.invariants:
        if inv.type == "gain_ceiling":
            for i in range(n_steps):
                m = i / (n_steps - 1)
                for q in [0.0, 0.25, 0.5, 0.75, 1.0]:
                    enc = body.corners.interpolate(m, q)
                    db = cascade_response_db(enc, FREQS, SR)
                    peak = float(np.max(db))
                    if peak > inv.threshold_db:
                        violations.append(f"gain_ceiling: {peak:.1f}dB > {inv.threshold_db:.1f}dB at m={m:.2f} q={q:.2f}")
                        return False, violations

        elif inv.type == "sub_floor":
            sub_mask = (FREQS >= 20) & (FREQS <= inv.freq_hz * 1.5)
            for i in range(n_steps):
                m = i / (n_steps - 1)
                enc = body.corners.interpolate(m, 0.5)
                db = cascade_response_db(enc, FREQS, SR)
                sub_db = float(np.max(db[sub_mask])) if np.any(sub_mask) else -100
                if sub_db < inv.threshold_db:
                    violations.append(f"sub_floor: {sub_db:.1f}dB < {inv.threshold_db:.1f}dB at m={m:.2f}")
                    return False, violations

        elif inv.type == "midpoint_alive":
            enc = body.corners.interpolate(0.5, 0.5)
            db = cascade_response_db(enc, FREQS, SR)
            dyn = float(np.max(db) - np.min(db))
            if dyn < 5.0:
                violations.append(f"midpoint_alive: dynamic range {dyn:.1f}dB < 5dB")
                return False, violations

        elif inv.type == "no_hf_takeover":
            for i in range(n_steps):
                m = i / (n_steps - 1)
                enc = body.corners.interpolate(m, 0.5)
                db = cascade_response_db(enc, FREQS, SR)
                lo_mask = FREQS < 1000
                hi_mask = FREQS > 5000
                lo_mean = float(np.mean(db[lo_mask])) if np.any(lo_mask) else 0
                hi_mean = float(np.mean(db[hi_mask])) if np.any(hi_mask) else 0
                if hi_mean - lo_mean > 20.0:
                    violations.append(f"no_hf_takeover: HF {hi_mean:.1f} > LF {lo_mean:.1f} + 20dB at m={m:.2f}")
                    return False, violations

    return True, []


# =============================================================================
# Global solver — find corners from trajectory spec
# =============================================================================

def _unpack_body(params: np.ndarray, spec: TrajectorySpec) -> Body:
    """Decode parameter vector into a Body.

    6 stages parameterized as:
      params[0:6]   = log(freq_hz) at M0 for each stage
      params[6:12]  = radius at M0 for each stage
      params[12:18] = log(freq_hz) at M100 for each stage
      params[18:24] = radius at M100 for each stage
      params[24]    = Q compression amount
      params[25]    = Q radius mod
    Total: 26 params.

    Q corners derived from Q behavior + these params.
    """
    n = 6

    m0_stages = []
    m100_stages = []
    for i in range(n):
        f0 = np.exp(params[i])
        r0 = min(MAX_RADIUS, max(0.5, params[6 + i]))
        f1 = np.exp(params[12 + i])
        r1 = min(MAX_RADIUS, max(0.5, params[18 + i]))
        m0_stages.append((f0, r0))
        m100_stages.append((f1, r1))

    q_comp = max(0.0, min(1.0, params[24]))
    q_r_mod = max(-0.005, min(0.005, params[25]))

    corner_a = _build_corner_from_stages(m0_stages)
    corner_c = _build_corner_from_stages(m100_stages)

    # Q corners: compress frequencies toward center, adjust radii
    center = spec.q_behavior.compression_center_hz
    m0_q = [(f + (center - f) * q_comp, min(MAX_RADIUS, r + q_r_mod)) for f, r in m0_stages]
    m100_q = [(f + (center - f) * q_comp, min(MAX_RADIUS, r + q_r_mod)) for f, r in m100_stages]
    corner_b = _build_corner_from_stages(m0_q)
    corner_d = _build_corner_from_stages(m100_q)

    corners = CornerArray(a=corner_a, b=corner_b, c=corner_c, d=corner_d)
    body = Body(name=spec.name, corners=corners, boost=4.0)
    body = solve_c4_surface(body, target_db=spec.target_db)
    body = apply_gain_budget(body)
    return body


def _eval_trajectory(body: Body, spec: TrajectorySpec) -> dict:
    """Score a body against its trajectory spec.

    Returns dict with:
      trajectory_fit: weighted mean of landmark scores (0-1)
      invariant_pass: bool
      per_landmark: list of per-landmark scores
      violations: list of invariant violation messages
    """
    # Score each landmark
    scores = []
    for lm in spec.landmarks:
        s = _score_landmark(body, lm)
        scores.append(s)

    trajectory_fit = float(np.mean(scores)) if scores else 0.0

    # Check invariants
    inv_pass, violations = _check_invariants(body, spec)

    return {
        "trajectory_fit": trajectory_fit,
        "invariant_pass": inv_pass,
        "per_landmark": scores,
        "violations": violations,
    }


class TrajectoryProblem(Problem):
    """Pymoo problem: maximize trajectory fit, respect invariants."""

    def __init__(self, spec: TrajectorySpec):
        self.spec = spec

        # Bounds for 6 stages: log(freq) + radius at M0 and M100, plus Q params
        # Frequency range: 20 Hz to 18000 Hz in log space
        log_lo = math.log(20.0)
        log_hi = math.log(18000.0)
        r_lo = 0.920
        r_hi = MAX_RADIUS

        xl = ([log_lo] * 6 + [r_lo] * 6 +    # M0: freq + radius
              [log_lo] * 6 + [r_lo] * 6 +    # M100: freq + radius
              [0.0, -0.005])                   # Q params
        xu = ([log_hi] * 6 + [r_hi] * 6 +
              [log_hi] * 6 + [r_hi] * 6 +
              [1.0, 0.005])

        super().__init__(
            n_var=26,
            n_obj=1,       # single objective: maximize trajectory fit
            n_ieq_constr=1,  # invariant pass
            xl=np.array(xl),
            xu=np.array(xu),
        )

    def _evaluate(self, X, out, *args, **kwargs):
        F, G = [], []
        for x in X:
            try:
                body = _unpack_body(x, self.spec)
                result = _eval_trajectory(body, self.spec)
                F.append([-result["trajectory_fit"]])  # negate for minimization
                g = 0.0 if result["invariant_pass"] else 1.0
                G.append([g])
            except Exception:
                F.append([0.0])
                G.append([1.0])
        out["F"] = np.array(F)
        out["G"] = np.array(G)


# =============================================================================
# Seed from spec — initial population biased toward landmark frequencies
# =============================================================================

def _seed_from_spec(spec: TrajectorySpec, n_individuals: int = 10) -> np.ndarray:
    """Create seeded initial population from the landmark frequencies.

    Extracts all peak/notch frequencies from landmarks and uses them
    as stage frequency seeds. Better starting point than random.
    """
    rng = np.random.default_rng(42)

    # Collect all event frequencies from landmarks
    all_freqs = []
    for lm in spec.landmarks:
        for ev in lm.events:
            all_freqs.append(ev.freq_hz)
    all_freqs = sorted(set(all_freqs))

    # Pick 6 representative frequencies (spread across the range)
    if len(all_freqs) >= 6:
        indices = np.linspace(0, len(all_freqs) - 1, 6, dtype=int)
        seed_freqs = [all_freqs[i] for i in indices]
    else:
        seed_freqs = (all_freqs + [1000.0] * 6)[:6]

    # M0 and M100 freqs from first and last landmark
    m0_freqs = seed_freqs[:]
    m100_freqs = seed_freqs[:]

    # Try to use first landmark's peaks for M0, last for M100
    first_peaks = [e.freq_hz for e in spec.landmarks[0].events if e.type == "peak"]
    last_peaks = [e.freq_hz for e in spec.landmarks[-1].events if e.type == "peak"]

    for i, fp in enumerate(first_peaks[:6]):
        m0_freqs[i] = fp
    for i, lp in enumerate(last_peaks[:6]):
        m100_freqs[i] = lp

    # Build seeded population
    seeds = []
    for _ in range(n_individuals):
        params = np.zeros(26)
        for i in range(6):
            jitter = rng.uniform(0.8, 1.2)
            params[i] = math.log(max(20, min(18000, m0_freqs[i] * jitter)))
            params[6 + i] = min(MAX_RADIUS, 0.960 + rng.uniform(0.0, 0.025))
            params[12 + i] = math.log(max(20, min(18000, m100_freqs[i] * jitter)))
            params[18 + i] = min(MAX_RADIUS, 0.960 + rng.uniform(0.0, 0.025))
        params[24] = spec.q_behavior.compression_amount + rng.uniform(-0.1, 0.1)
        params[25] = spec.q_behavior.tightening + rng.uniform(-0.002, 0.002)
        seeds.append(params)

    return np.array(seeds)


# =============================================================================
# Plotting
# =============================================================================

def _plot_trajectory(body: Body, spec: TrajectorySpec, out_path: Path):
    """Plot morph sweep with landmark annotations."""
    fig, ax = plt.subplots(figsize=(14, 7))
    cmap = plt.cm.coolwarm

    # Plot sweep
    n_steps = 21
    for i in range(n_steps):
        m = i / (n_steps - 1)
        enc = body.corners.interpolate(m, 0.5)
        db = cascade_response_db(enc, FREQS, SR)
        color = cmap(m)
        label = f"M={m:.1f}" if i % 4 == 0 else None
        ax.semilogx(FREQS, db, color=color, alpha=0.6, linewidth=1.0, label=label)

    # Annotate landmarks
    for lm in spec.landmarks:
        for ev in lm.events:
            if ev.type == "peak":
                ax.axvline(ev.freq_hz, color=cmap(lm.morph), alpha=0.3, linestyle="--", linewidth=0.8)
            elif ev.type in ("notch", "void"):
                ax.axvline(ev.freq_hz, color=cmap(lm.morph), alpha=0.3, linestyle=":", linewidth=0.8)

    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Magnitude (dB)")
    ax.set_title(f"{spec.name} — Trajectory Sweep")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(20, 20000)
    ax.set_ylim(-40, 50)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()


def _plot_landmark_scores(scores: list[float], spec: TrajectorySpec, out_path: Path):
    """Bar chart of per-landmark trajectory scores."""
    fig, ax = plt.subplots(figsize=(10, 4))
    morph_positions = [lm.morph for lm in spec.landmarks]
    colors = [plt.cm.coolwarm(m) for m in morph_positions]
    bars = ax.bar(range(len(scores)), scores, color=colors, edgecolor="black", linewidth=0.5)
    ax.set_xticks(range(len(scores)))
    ax.set_xticklabels([f"m={m:.1f}" for m in morph_positions], fontsize=9)
    ax.set_ylabel("Landmark Score")
    ax.set_title(f"{spec.name} — Landmark Fit")
    ax.set_ylim(0, 1.1)
    ax.axhline(0.5, color="red", alpha=0.3, linestyle="--", label="50% threshold")
    ax.legend()
    ax.grid(axis="y", alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()


# =============================================================================
# Main solver
# =============================================================================

def solve_trajectory(
    key: str,
    generations: int = 80,
    pop_size: int = 60,
    seed: int = 42,
) -> Path:
    """Solve a trajectory spec into a 4-corner body."""
    spec = ALL_SPECS[key]

    print(f"\n{'='*60}")
    print(f"  TRAJECTORY FORGE: {spec.name}")
    print(f"  Landmarks: {len(spec.landmarks)}")
    print(f"  Invariants: {len(spec.invariants)}")
    print(f"  Generations: {generations}, Population: {pop_size}")
    print(f"  r cap: {MAX_RADIUS}, gain ceiling: {GAIN_CEILING_DB} dB")
    print(f"{'='*60}")

    # Seed population
    seed_pop = _seed_from_spec(spec, n_individuals=pop_size // 3)

    problem = TrajectoryProblem(spec)
    algorithm = NSGA2(
        pop_size=pop_size,
        sampling=FloatRandomSampling(),
        crossover=SBX(prob=0.9, eta=15),
        mutation=PM(eta=20),
    )

    result = minimize(problem, algorithm, ("n_gen", generations),
                      seed=seed, verbose=False)

    # Output
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = VAULT_DIR / f"traj_{key}_{ts}"
    out_dir.mkdir(parents=True, exist_ok=True)

    if result.X is None:
        print(f"\n  NO FEASIBLE SOLUTIONS FOUND")
        print(f"  Output: {out_dir}")
        return out_dir

    # Handle single vs multiple solutions
    solutions = result.X if result.X.ndim == 2 else result.X.reshape(1, -1)
    fitnesses = -result.F if result.F.ndim == 2 else (-result.F).reshape(1, -1)

    bodies = []
    results_list = []
    for i, x in enumerate(solutions):
        body = _unpack_body(x, spec)
        ev = _eval_trajectory(body, spec)
        bodies.append(body)
        results_list.append(ev)

        path = out_dir / f"{key}_{i:03d}.json"
        named = Body(name=f"{spec.name}_{i:03d}", corners=body.corners, boost=body.boost)
        path.write_text(named.to_compiled_json(provenance="forge-trajectory"))

    # Rank by trajectory fit
    ranked = sorted(zip(bodies, results_list, range(len(bodies))),
                    key=lambda x: -x[1]["trajectory_fit"])

    print(f"\n  {len(bodies)} solutions saved to {out_dir}")
    print(f"\n  Top candidates:")
    for body, ev, idx in ranked[:10]:
        lm_str = " ".join(f"{s:.2f}" for s in ev["per_landmark"])
        inv = "PASS" if ev["invariant_pass"] else "FAIL"
        print(f"    #{idx:03d}: fit={ev['trajectory_fit']:.3f}  inv={inv}  landmarks=[{lm_str}]")

    # Save best
    if ranked:
        best_body, best_ev, best_idx = ranked[0]
        best_path = out_dir / f"BEST_{spec.name.replace(' ', '_')}.json"
        best_named = Body(name=spec.name, corners=best_body.corners, boost=best_body.boost)
        best_path.write_text(best_named.to_compiled_json(provenance="forge-trajectory-best"))
        print(f"\n  Best: #{best_idx:03d} fit={best_ev['trajectory_fit']:.3f}")

        # Plots
        _plot_trajectory(best_body, spec, out_dir / "trajectory_sweep.png")
        _plot_landmark_scores(best_ev["per_landmark"], spec, out_dir / "landmark_scores.png")
        print(f"  Plots saved to {out_dir}")

    return out_dir


# =============================================================================
# CLI
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description="Trajectory-first forge")
    parser.add_argument("--body", default="speaker_knockerz")
    parser.add_argument("--generations", type=int, default=80)
    parser.add_argument("--pop", type=int, default=60)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    keys = list(ALL_SPECS.keys()) if args.body == "all" else [args.body]
    for key in keys:
        if key not in ALL_SPECS:
            print(f"Unknown body: {key}. Available: {list(ALL_SPECS.keys())}")
            continue
        solve_trajectory(key, args.generations, args.pop, args.seed)


if __name__ == "__main__":
    main()
