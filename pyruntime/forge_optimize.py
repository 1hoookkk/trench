"""Forge optimizer — pymoo multi-objective body search.

Per-category fitness functions, fast in-process evaluation,
Pareto-optimal output to vault with matplotlib plots.
Post-optimization taste model scoring for automatic ranking.

Usage:
    python pyruntime/forge_optimize.py --category vocal --generations 50 --pop 40
    python pyruntime/forge_optimize.py --category all
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
from pyruntime.macro_compile import (
    Actor, BodySpec, CompileMacro, PressureBehavior,
    SlotSpec, SlotState, compile_body, freq_to_place,
)
from pyruntime.stage_params import StageParams
from pyruntime.stage_math import resonator
from pyruntime.encode import raw_to_encoded
from pyruntime.freq_response import cascade_response_db, freq_points
from pyruntime.analysis import morph_trajectory_distance, shipping_gate, dense_midpoint_audit
from pyruntime.target import _empty_slot
from pyruntime.zero_law import ContourFamily, ContourKind

VAULT_DIR = Path(__file__).parent.parent / "vault"
FREQS = freq_points(sr=SR)
TASTE_MODEL_PATH = Path(__file__).parent.parent / "taste_model.json"
_VOCAL_MASK = (FREQS >= 200) & (FREQS <= 5000)
_VOCAL_FREQS = FREQS[_VOCAL_MASK]


# ---------------------------------------------------------------------------
# Taste model scorer (Python port of trench-foundry taste.rs)
# ---------------------------------------------------------------------------

class TasteScorer:
    """Calibrated logistic regression scorer matching the Rust TasteModel.

    Loads taste_model.json and scores bodies using the model-only path
    (no corpus kinship — that requires the full P2K embedding which lives
    in the Rust binary). This gives the calibrated probability from the
    logistic regression + Platt scaling, which is the 30% blend component
    in vault_seeds. Good enough for ranking optimizer output.
    """

    def __init__(self, model_path: Path):
        with open(model_path) as f:
            m = json.load(f)
        self.active_features = m["active_features"]
        self.means = np.array(m["means"])
        self.stds = np.array(m["stds"])
        self.weights = np.array(m["weights"])
        self.intercept = m["intercept"]
        self.cal_weight = m.get("calibrator_weight", 1.0)
        self.cal_bias = m.get("calibrator_bias", 0.0)
        self.calibrated = m.get("calibration_enabled", False)

    def score(self, features: dict) -> float:
        """Score a body's features dict → calibrated probability [0, 1]."""
        # Feature order: fitness, trust_score, morph_sweep, q_cliff, sci,
        #                emergence, movement, novelty
        raw = np.array([
            features.get("fitness", 0.4),
            features.get("trust_score", 0.4),
            features.get("morph_sweep", 0.3),
            features.get("q_cliff", 0.3),
            features.get("sci", 0.5),
            features.get("emergence", 0.5),
            features.get("movement", 0.3),
            features.get("novelty", 0.9),   # default to corpus mean
        ])
        # Select active features and standardize
        active = raw[self.active_features]
        z = (active - self.means) / np.maximum(self.stds, 1e-9)
        logit = self.intercept + float(self.weights @ z)
        if self.calibrated:
            logit = self.cal_bias + self.cal_weight * logit
        return 1.0 / (1.0 + math.exp(-logit))


def _sigmoid(x: float) -> float:
    return 1.0 / (1.0 + math.exp(-max(-30.0, min(30.0, x))))


def _percentile(vals: list[float], p: float) -> float:
    if not vals:
        return 0.0
    s = sorted(vals)
    idx = p * (len(s) - 1)
    lo = int(idx)
    hi = min(lo + 1, len(s) - 1)
    frac = idx - lo
    return s[lo] * (1 - frac) + s[hi] * frac


def _response_distance(db_a: np.ndarray, db_b: np.ndarray) -> float:
    return float(np.sqrt(np.mean((db_a - db_b) ** 2)))


def _extract_taste_features(body: Body) -> dict:
    """Extract the 8 taste features from a body using 4x4 grid sampling.

    Mirrors the Rust extract_base_signals → to_features pipeline.
    """
    grid = 4
    snapshots_db = []
    snapshots_meta = []  # (peak, valley, dyn_range, n_peaks, n_notches, centroid)

    for mi in range(grid):
        morph = mi / max(1, grid - 1)
        for qi in range(grid):
            q = qi / max(1, grid - 1)
            enc = body.corners.interpolate(morph, q)
            db = cascade_response_db(enc, FREQS, SR)
            snapshots_db.append(db)

            peak = float(np.max(db))
            valley = float(np.min(db))
            mean_db = float(np.mean(db))

            # Count peaks and notches
            n_peaks = n_notches = 0
            for i in range(1, len(db) - 1):
                if db[i] > db[i-1] and db[i] > db[i+1] and db[i] > valley + 3:
                    n_peaks += 1
                if db[i] < db[i-1] and db[i] < db[i+1] and db[i] < peak - 3:
                    n_notches += 1

            # Spectral centroid
            lin = 10.0 ** (db / 20.0)
            denom = float(np.sum(lin))
            centroid = float(np.sum(FREQS * lin) / denom) if denom > 0 else 1000.0

            # Ridge prominence: max peak above mean, normalized
            ridge_prom = min(1.0, max(0.0, (peak - mean_db) / 40.0))

            # Talkingness proxy: peaks in 200-5000 Hz above mean+3
            mask = (FREQS >= 200) & (FREQS <= 5000)
            vocal_db = db[mask]
            vocal_peaks = 0
            if len(vocal_db) > 2:
                for i in range(1, len(vocal_db) - 1):
                    if vocal_db[i] > vocal_db[i-1] and vocal_db[i] > vocal_db[i+1]:
                        if vocal_db[i] > mean_db + 3:
                            vocal_peaks += 1
            talkingness = min(1.0, vocal_peaks / 4.0)

            # Fragmentation: spectral variance in octave bands
            ruggedness = min(1.0, float(np.std(db)) / 30.0)

            snapshots_meta.append({
                "peak": peak, "valley": valley,
                "dyn_range": peak - valley,
                "n_peaks": n_peaks, "n_notches": n_notches,
                "centroid": centroid, "ridge_prom": ridge_prom,
                "talkingness": talkingness, "ruggedness": ruggedness,
            })

    # --- Profile (worst-corner blend) ---
    def _profile_score(m):
        dyn_score = max(0, min(1, 1.0 - math.exp(-m["dyn_range"] / 18.0)))
        shape_score = min(1.0, (m["n_peaks"] + m["n_notches"]) / 5.0)
        centroid_oct = abs(math.log2(max(40, m["centroid"]) / 1400.0))
        centroid_score = max(0, min(1, 1.0 - centroid_oct / 4.0))
        peak_score = _sigmoid((m["peak"] - 1.0) / 3.0)
        return max(0, min(1, 0.45*dyn_score + 0.20*shape_score + 0.20*centroid_score + 0.15*peak_score))

    corner_profiles = [_profile_score(snapshots_meta[i]) for i in [0, 3, 12, 15]]
    p0 = corner_profiles[0]
    worst = min(corner_profiles)
    profile = 0.4 * p0 + 0.6 * worst

    # --- Emergence (ridge prominence) ---
    ridges = [m["ridge_prom"] for m in snapshots_meta]
    emergence = max(0, min(1, 0.60 * np.mean(ridges) + 0.40 * _percentile(ridges, 0.80)))

    # --- Movement (neighbor response distances) ---
    distances = []
    for mi in range(grid):
        for qi in range(grid):
            idx = mi * grid + qi
            if qi + 1 < grid:
                distances.append(_response_distance(snapshots_db[idx], snapshots_db[idx + 1]))
            if mi + 1 < grid:
                distances.append(_response_distance(snapshots_db[idx], snapshots_db[idx + grid]))
    mean_d = float(np.mean(distances)) if distances else 0.0
    p75 = _percentile(distances, 0.75)
    movement = max(0, min(1, 1.0 - math.exp(-(0.6*mean_d + 0.4*p75) / 8.0)))

    # --- Morph sweep (M-axis distances at each Q) ---
    morph_deltas = []
    for qi in range(grid):
        for mi in range(grid - 1):
            morph_deltas.append(_response_distance(
                snapshots_db[mi * grid + qi],
                snapshots_db[(mi + 1) * grid + qi]))
        morph_deltas.append(_response_distance(
            snapshots_db[0 * grid + qi],
            snapshots_db[(grid - 1) * grid + qi]))
    morph_sweep = max(0, min(1,
        1.0 - math.exp(-(0.5 * np.mean(morph_deltas) + 0.5 * _percentile(morph_deltas, 0.80)) / 9.0)
    )) if morph_deltas else 0.0

    # --- Q cliff (Q-axis distances at each M) ---
    q_deltas = []
    for mi in range(grid):
        for qi in range(grid - 1):
            q_deltas.append(_response_distance(
                snapshots_db[mi * grid + qi],
                snapshots_db[mi * grid + qi + 1]))
    edge = (0.7 * _percentile(q_deltas, 0.90) + 0.3 * max(q_deltas)) if q_deltas else 0.0
    q_cliff = max(0, min(1, 1.0 - math.exp(-edge / 10.0)))

    # --- SCI (spectral complexity) ---
    peak_density = [min(1, (m["n_peaks"] + m["n_notches"]) / 12.0) for m in snapshots_meta]
    ruggedness_vals = [m["ruggedness"] for m in snapshots_meta]
    sci = max(0, min(1, 0.30*np.mean(peak_density) + 0.25*np.mean(ruggedness_vals)
                       + 0.25*np.mean(ruggedness_vals) + 0.20*0.5))

    # --- Content (talkingness composite) ---
    talk_vals = [m["talkingness"] for m in snapshots_meta]
    content = max(0, min(1, 0.55*np.mean(talk_vals) + 0.25*0.5 + 0.20*0.8))

    # --- Coherence ---
    coherence = max(0, min(1, 0.65*0.6 + 0.20*(1.0 - np.mean(ruggedness_vals)) + 0.15*0.8))

    # --- Low end ---
    low_cut = max(1, len(FREQS) // 8)
    mid_cut = max(low_cut + 1, len(FREQS) // 2)
    low_end_vals = []
    for db in snapshots_db:
        low_band = float(np.mean(db[:low_cut]))
        mid_band = float(np.mean(db[low_cut:mid_cut]))
        low_end_vals.append(low_band - mid_band)
    low_end = max(0, min(1, _sigmoid(np.mean(low_end_vals) / 6.0)))

    # --- Compose taste features (matches Rust BaseSignals::to_features) ---
    expressiveness = max(0, min(1, 0.45*morph_sweep + 0.35*q_cliff + 0.20*sci))
    trust_score = max(0, min(1, 0.40*profile + 0.35*coherence + 0.20*low_end + 0.05*content))
    fitness = max(0, min(1, 0.30*profile + 0.22*expressiveness + 0.16*content
                           + 0.12*sci + 0.12*coherence + 0.08*low_end))

    return {
        "fitness": fitness,
        "trust_score": trust_score,
        "morph_sweep": morph_sweep,
        "q_cliff": q_cliff,
        "sci": sci,
        "emergence": emergence,
        "movement": movement,
        "novelty": 0.9,  # no corpus in Python; use training mean
    }


# Try to load taste model at import time
_taste_scorer = None
if TASTE_MODEL_PATH.is_file():
    try:
        _taste_scorer = TasteScorer(TASTE_MODEL_PATH)
    except Exception as e:
        print(f"  Warning: could not load taste model: {e}")


VOWELS = {
    "oo": (378, 997, 2343, 3357),
    "ee": (342, 2322, 3000, 3657),
    "ah": (768, 1189, 2555, 3508),
    "eh": (660, 1720, 2410, 3500),
    "oh": (450, 1030, 2500, 3500),
    "ae": (730, 1090, 2440, 3500),
    "ih": (446, 1993, 2657, 3599),
    "er": (490, 1350, 1690, 3500),
    "aw": (570, 840, 2410, 3500),
    "uh": (640, 1190, 2390, 3500),
}
VOWEL_KEYS = list(VOWELS.keys())

ACTOR_CONTOUR = {
    Actor.FOUNDATION: (ContourKind.PURE, 0.0),
    Actor.MASS: (ContourKind.NEAR_ALLPASS, 0.3),
    Actor.THROAT: (ContourKind.INTERIOR_ZERO, 0.5),
    Actor.BITE: (ContourKind.UNIT_CIRCLE, 0.6),
    Actor.AIR: (ContourKind.INTERIOR_ZERO, 0.4),
    Actor.SCAR: (ContourKind.UNIT_CIRCLE, 0.0),
}


# ---------------------------------------------------------------------------
# Fast evaluation helpers (no file I/O)
# ---------------------------------------------------------------------------

def _talkingness_single(db: np.ndarray) -> float:
    """Vocal peak count in 200–5000 Hz, normalised to [0, 1].

    Threshold is full-spectrum mean + 3 dB (matches original metric).
    """
    vocal_db = db[_VOCAL_MASK]
    if len(vocal_db) < 3:
        return 0.0
    mean_db = float(np.mean(db))           # full spectrum baseline
    peaks = 0
    for i in range(1, len(vocal_db) - 1):
        if vocal_db[i] > vocal_db[i - 1] and vocal_db[i] > vocal_db[i + 1]:
            if vocal_db[i] > mean_db + 3.0:
                peaks += 1
    return min(1.0, peaks / 4.0)


def _extract_formants(db: np.ndarray, n: int = 3) -> list[float]:
    """Return up to *n* formant frequencies (Hz), sorted low → high.

    Uses full-spectrum mean as the threshold baseline (not vocal-band
    mean) so formant peaks in a broadband response are still detected.
    """
    vocal_db = db[_VOCAL_MASK]
    if len(vocal_db) < 3:
        return []
    mean_db = float(np.mean(db))          # full spectrum baseline
    peaks: list[tuple[float, float]] = []
    for i in range(1, len(vocal_db) - 1):
        if vocal_db[i] > vocal_db[i - 1] and vocal_db[i] > vocal_db[i + 1]:
            if vocal_db[i] > mean_db + 3.0:
                peaks.append((float(_VOCAL_FREQS[i]), float(vocal_db[i])))
    peaks.sort(key=lambda p: -p[1])
    return sorted(p[0] for p in peaks[:n])


def _trajectory_score(tracks: list[list[float]]) -> float:
    """Score vowel-path quality from per-sweep-point formant lists.

    High score = F1/F2 move through vowel space with range, stay in
    the vowel triangle, and traverse a distinct arc in log-freq space.
    """
    valid = [ft for ft in tracks if len(ft) >= 2]
    if len(valid) < 3:
        return 0.0
    f1 = [ft[0] for ft in valid]
    f2 = [ft[1] for ft in valid]

    # Range: total Hz movement of F1 + F2
    range_score = min(1.0, (max(f1) - min(f1) + max(f2) - min(f2)) / 1500.0)

    # Vowel proximity: fraction of points in the vowel triangle
    hits = sum(1 for a, b in zip(f1, f2) if 200 < a < 1000 and 500 < b < 3000)
    vowel_score = hits / len(valid)

    # Log-space arc length (longer path ≈ more distinctive transition)
    arc = 0.0
    for i in range(1, len(valid)):
        d1 = math.log(f1[i] / f1[i - 1]) if f1[i] > 0 and f1[i - 1] > 0 else 0.0
        d2 = math.log(f2[i] / f2[i - 1]) if f2[i] > 0 and f2[i - 1] > 0 else 0.0
        arc += math.sqrt(d1 * d1 + d2 * d2)
    arc_score = min(1.0, arc / 2.0)

    return 0.4 * range_score + 0.3 * vowel_score + 0.3 * arc_score


def _continuity_score(tracks: list[list[float]]) -> float:
    """Reward smooth formant migration, penalise jumps / crossovers."""
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
            if ratio > 1.5:          # > 50 % jump ≈ crossover
                penalty += ratio - 1.0
    return max(0.0, min(1.0, 1.0 - penalty / 3.0))


# ---- 5×5 morph/Q grid positions (reused across all eval calls) ----------
_OCC_STEPS = [i / 4.0 for i in range(5)]       # 0, .25, .5, .75, 1


def _eval_body(body: Body) -> dict:
    """Layered fitness stack.

    L1  hard gate   — peak ceiling, dead midpoint, invalid open state
    L2  occupancy   — fraction of 5×5 bins with talk > 0
    L3  trajectory  — F1/F2/F3 path quality over the Q = 0.5 morph sweep
    L4  continuity  — smooth formant migration across sweep
    L5  flavor      — talk · ridge · coherence (secondary composite)
    """
    # -- 5×5 grid: 25 cascade evals, cached in dict -----------------------
    responses: dict[tuple[float, float], np.ndarray] = {}
    for mi in _OCC_STEPS:
        for qi in _OCC_STEPS:
            enc = body.corners.interpolate(mi, qi)
            responses[(mi, qi)] = cascade_response_db(enc, FREQS, SR)

    # -- L1: hard gates ----------------------------------------------------
    max_peak = max(float(np.max(db)) for db in responses.values())
    db_mid = responses[(0.5, 0.5)]
    gate_fail = max_peak > 50.0

    # -- L2: occupancy -----------------------------------------------------
    talk_grid = {k: _talkingness_single(db) for k, db in responses.items()}
    dead = sum(1 for t in talk_grid.values() if t == 0.0)
    occupancy = (25 - dead) / 25.0

    # -- L3/L4: trajectory + continuity (Q = 0.5 sweep, 5 points) ---------
    sweep_keys = [(m, 0.5) for m in _OCC_STEPS]
    formant_tracks = [_extract_formants(responses[k]) for k in sweep_keys]
    trajectory = _trajectory_score(formant_tracks)
    continuity = _continuity_score(formant_tracks)

    # -- L5: secondary flavour metrics -------------------------------------
    ridge_vals = []
    for db in responses.values():
        pk = float(np.max(db))
        mn = float(np.mean(db))
        ridge_vals.append(max(0.0, min(1.0, (pk - mn) / 40.0)))
    ridge = float(np.mean(ridge_vals))

    tvs = list(talk_grid.values())
    talk_max = max(tvs)
    coherence = (1.0 - float(np.std(tvs) / talk_max)) if talk_max > 0 else 0.0

    alive = [t for t in tvs if t > 0]
    talkingness = float(np.mean(alive)) if alive else 0.0

    flavor = 0.4 * talkingness + 0.3 * ridge + 0.3 * coherence

    # -- secondary metrics for non-voice categories ------------------------
    peak = float(np.max(db_mid))
    valley = float(np.min(db_mid))
    dynamic_range = peak - valley
    ruggedness = min(1.0, float(np.std(db_mid)) / 30.0)
    tilt = float(db_mid[-1] - db_mid[0])
    morph_dist = float(0.5 * (
        np.sqrt(np.mean((responses[(0.0, 0.0)] - responses[(1.0, 0.0)]) ** 2))
        + np.sqrt(np.mean((responses[(0.0, 1.0)] - responses[(1.0, 1.0)]) ** 2))
    ))

    # -- normalized morph distance (shape change, not level) ---------------
    def _norm(db):
        return db - float(np.mean(db))

    norm_morph_dist = float(0.5 * (
        np.sqrt(np.mean((_norm(responses[(0.0, 0.0)]) - _norm(responses[(1.0, 0.0)])) ** 2))
        + np.sqrt(np.mean((_norm(responses[(0.0, 1.0)]) - _norm(responses[(1.0, 1.0)])) ** 2))
    ))

    # -- SK: sub presence (mean energy in 20-60 Hz across morph sweep) -----
    _sub_mask = (FREQS >= 20) & (FREQS <= 60)
    sub_energies = []
    for m in _OCC_STEPS:
        db_q0 = responses[(m, 0.0)]
        sub_energies.append(float(np.mean(db_q0[_sub_mask])))
    sub_presence = float(np.mean(sub_energies))
    sub_stability = -float(np.std(sub_energies))  # higher = more stable (negated std)

    # -- SK: spectral contrast (sub vs mid-high separation) ----------------
    _mid_mask = (FREQS >= 200) & (FREQS <= 2000)
    contrast_vals = []
    for m in _OCC_STEPS:
        db_q0 = responses[(m, 0.0)]
        sub_e = float(np.mean(db_q0[_sub_mask]))
        mid_e = float(np.mean(db_q0[_mid_mask]))
        contrast_vals.append(sub_e - mid_e)
    sub_contrast = float(np.mean(contrast_vals))

    # -- CDS: structural evolution (peak count change across morph) --------
    def _count_peaks(db, threshold_above_mean=3.0):
        mean_val = float(np.mean(db))
        count = 0
        for i in range(1, len(db) - 1):
            if db[i] > db[i-1] and db[i] > db[i+1] and db[i] > mean_val + threshold_above_mean:
                count += 1
        return count

    peaks_closed = _count_peaks(responses[(0.0, 0.0)])
    peaks_open = _count_peaks(responses[(1.0, 0.0)])
    # Reward bodies where the open state has more spectral peaks than closed
    peak_evolution = min(1.0, max(0.0, (peaks_open - peaks_closed) / 6.0))

    # -- CDS: spectral topology change (shape correlation between endpoints)
    db_closed_n = _norm(responses[(0.0, 0.0)])
    db_open_n = _norm(responses[(1.0, 0.0)])
    shape_corr = float(np.corrcoef(db_closed_n, db_open_n)[0, 1])
    # Low correlation = topologically different shapes (good for CDS)
    topology_shift = max(0.0, min(1.0, 1.0 - shape_corr))

    return {
        "gate_fail": gate_fail,
        "occupancy": occupancy,
        "trajectory": trajectory,
        "continuity": continuity,
        "flavor": flavor,
        "talkingness": talkingness,
        "ridge": ridge,
        "coherence": coherence,
        "dynamic_range": dynamic_range,
        "ruggedness": ruggedness,
        "morph_distance": morph_dist,
        "norm_morph_distance": norm_morph_dist,
        "spectral_tilt": tilt,
        "peak_db": max_peak,
        "sub_presence": sub_presence,
        "sub_stability": sub_stability,
        "sub_contrast": sub_contrast,
        "peak_evolution": peak_evolution,
        "topology_shift": topology_shift,
    }


# ---------------------------------------------------------------------------
# Stitch body builder (for vocal/drum/bass/pad/lead)
# ---------------------------------------------------------------------------

def _vowel_spec_from_params(name: str, vowel_idx: int, focus: float, weight: float) -> BodySpec:
    """Build a vowel BodySpec from optimizer parameters."""
    key = VOWEL_KEYS[int(vowel_idx) % len(VOWEL_KEYS)]
    f1, f2, f3, f4 = VOWELS[key]
    slots = [_empty_slot(a) for a in Actor.ALL]
    used = set()
    for freq in (f1, f2, f3, f4):
        for idx, actor in enumerate(Actor.ALL):
            if idx not in used and actor.freq_floor() <= freq <= actor.freq_ceiling():
                c, clr = ACTOR_CONTOUR.get(actor, (ContourKind.UNIT_CIRCLE, 0.5))
                state = SlotState(
                    place=freq_to_place(actor, freq), focus=focus, weight=weight,
                    contour=c, color=clr, hue=0.0, contour_override=None)
                slots[idx] = SlotSpec(actor=actor, state_a=state, state_b=state,
                                      pressure=PressureBehavior.TIGHTEN,
                                      compile_macro=CompileMacro.SINGLE, spread=0.0)
                used.add(idx)
                break
    return BodySpec(name=name, slots=slots, boost=4.0)


def _build_stitch(vowel_a_idx: int, vowel_b_idx: int,
                  focus: float, weight: float) -> Body:
    """Build a stitched body from two vowel indices + shared params."""
    sa = _vowel_spec_from_params("_a", vowel_a_idx, focus, weight)
    sb = _vowel_spec_from_params("_b", vowel_b_idx, focus, weight)
    ca = compile_body(sa)
    cb = compile_body(sb)
    corners = CornerArray(
        a=ca.corner(CornerName.A), b=ca.corner(CornerName.B),
        c=cb.corner(CornerName.C), d=cb.corner(CornerName.D),
    )
    return Body(name="opt", corners=corners, boost=4.0)


# ---------------------------------------------------------------------------
# 12-stage body builder (for 808/hihat/lofi/pluck/fx)
# ---------------------------------------------------------------------------

def _build_12stage(params: np.ndarray, category: str) -> Body:
    """Build a 12-stage body from a flat parameter vector.

    Per stage: [log_freq, radius, weight, notch_type(0-3)]
    Total: 12 * 4 = 48 params per morph endpoint, 96 total.
    params[0:48] = morph=0, params[48:96] = morph=1.
    """
    def _make_corner(p: np.ndarray, q_boost: float) -> CornerState:
        stages = []
        for i in range(12):
            base = i * 4
            freq = np.exp(p[base])  # log-freq space
            freq = max(20.0, min(SR * 0.48 - 50, freq))
            radius = min(0.998, max(0.30, p[base + 1]) + q_boost)
            weight = max(0.0, min(0.80, p[base + 2]))
            notch_type = int(p[base + 3]) % 4

            theta = TWO_PI * freq / SR
            r = radius
            a1 = -2.0 * r * math.cos(theta)
            val1 = -1.0 + weight * 0.85

            if notch_type == 0:  # pure
                val2, val3 = 0.0, 0.0
            elif notch_type == 1:  # unit_circle
                val2, val3 = ContourFamily.UNIT_CIRCLE.val2_val3(a1, r)
            elif notch_type == 2:  # interior
                zero_r = r * 0.6
                family = ContourFamily.interior_zero(zero_r, theta)
                val2, val3 = family.val2_val3(a1, r)
            else:  # offset notch
                offset = freq * 0.2
                zero_freq = max(20.0, min(SR * 0.48, freq + offset))
                zero_theta = TWO_PI * zero_freq / SR
                zero_a1 = -2.0 * math.cos(zero_theta)
                val2 = zero_a1 - a1
                val3 = r * r - 1.0

            stages.append(StageParams(a1=a1, r=r, val1=val1, val2=val2, val3=val3))

        return CornerState(stages=stages, boost=4.0)

    p_a = params[:48]
    p_b = params[48:96]

    # Interpolate for 4 corners
    corners = []
    for morph, q in [(0.0, 0.0), (0.0, 1.0), (1.0, 0.0), (1.0, 1.0)]:
        p = p_a * (1 - morph) + p_b * morph
        q_boost = q * 0.02
        corners.append(_make_corner(p, q_boost))

    return Body(name="opt", corners=CornerArray(
        a=corners[0], b=corners[1], c=corners[2], d=corners[3]), boost=4.0)


# ---------------------------------------------------------------------------
# SK builder: locked sub foundation + searched character layer
# ---------------------------------------------------------------------------

_SK_FOUNDATION_FREQS = [30.0, 45.0, 55.0, 70.0]  # 4 sub-band stages, constant across morph

_SK_CHAR_BOUNDS = [
    (100, 300),    # choke zone
    (200, 500),    # rip zone
    (400, 1200),   # rattle zone
    (800, 2000),   # cry zone
    (1500, 5000),  # stress zone
    (2500, 8000),  # fracture zone
    (4000, 12000), # air zone
    (6000, 16000), # sizzle zone
]


def _build_sk(params: np.ndarray) -> Body:
    """Build Speaker Knockerz body: locked sub + searched character.

    Layout:
      params[0:12]  = foundation: 4 stages × (log_freq, radius, weight)
      params[12:76] = character:  8 stages × 4 params × 2 morph endpoints
      params[76:84] = Q delta:    8 per-character-stage radius offsets

    Foundation stages are constant across morph and Q — the sub never moves.
    Character stages differ between M0 (dormant) and M100 (active).
    Q boosts character stage radii by the per-stage delta.
    """
    # -- Foundation: 4 stages, constant everywhere -------------------------
    foundation = []
    for i in range(4):
        base = i * 3
        freq = np.exp(params[base])
        freq = max(20.0, min(80.0, freq))
        radius = min(0.998, max(0.70, params[base + 1]))
        weight = max(0.0, min(0.50, params[base + 2]))
        theta = TWO_PI * freq / SR
        a1 = -2.0 * radius * math.cos(theta)
        foundation.append(StageParams(a1=a1, r=radius, val1=-1.0 + weight * 0.85,
                                       val2=0.0, val3=0.0))

    # -- Character: 8 stages × 2 endpoints --------------------------------
    def _make_char_stages(p: np.ndarray) -> list[StageParams]:
        stages = []
        for i in range(8):
            base = i * 4
            freq = np.exp(p[base])
            lo, hi = _SK_CHAR_BOUNDS[i]
            freq = max(float(lo), min(float(hi), freq))
            radius = min(0.998, max(0.30, p[base + 1]))
            weight = max(0.0, min(0.80, p[base + 2]))
            notch_type = int(p[base + 3]) % 4

            theta = TWO_PI * freq / SR
            r = radius
            a1 = -2.0 * r * math.cos(theta)
            val1 = -1.0 + weight * 0.85

            if notch_type == 0:
                val2, val3 = 0.0, 0.0
            elif notch_type == 1:
                val2, val3 = ContourFamily.UNIT_CIRCLE.val2_val3(a1, r)
            elif notch_type == 2:
                zero_r = r * 0.6
                family = ContourFamily.interior_zero(zero_r, theta)
                val2, val3 = family.val2_val3(a1, r)
            else:
                offset = freq * 0.2
                zero_freq = max(20.0, min(SR * 0.48, freq + offset))
                zero_theta = TWO_PI * zero_freq / SR
                zero_a1 = -2.0 * math.cos(zero_theta)
                val2 = zero_a1 - a1
                val3 = r * r - 1.0

            stages.append(StageParams(a1=a1, r=r, val1=val1, val2=val2, val3=val3))
        return stages

    char_m0 = _make_char_stages(params[12:44])   # 8 stages × 4 = 32 params
    char_m100 = _make_char_stages(params[44:76])  # 8 stages × 4 = 32 params
    q_deltas = params[76:84]                       # 8 per-stage radius offsets

    def _corner(morph: float, q: float) -> CornerState:
        stages = list(foundation)  # copy — sub is constant
        for i in range(8):
            s0 = char_m0[i]
            s1 = char_m100[i]
            # Interpolate character stages by morph
            a1 = s0.a1 * (1 - morph) + s1.a1 * morph
            r_base = s0.r * (1 - morph) + s1.r * morph
            r = min(0.998, r_base + q * max(0.0, min(0.15, q_deltas[i])))
            val1 = s0.val1 * (1 - morph) + s1.val1 * morph
            val2 = s0.val2 * (1 - morph) + s1.val2 * morph
            val3 = s0.val3 * (1 - morph) + s1.val3 * morph
            stages.append(StageParams(a1=a1, r=r, val1=val1, val2=val2, val3=val3))
        return CornerState(stages=stages, boost=4.0)

    return Body(name="opt", corners=CornerArray(
        a=_corner(0.0, 0.0), b=_corner(0.0, 1.0),
        c=_corner(1.0, 0.0), d=_corner(1.0, 1.0)), boost=4.0)


class SKProblem(Problem):
    """Speaker Knockerz: locked sub foundation + searched character layer."""

    def __init__(self):
        # 12 foundation + 32 char_m0 + 32 char_m100 + 8 q_delta = 84
        xl, xu = [], []
        # Foundation: 4 × (log_freq, radius, weight)
        for freq in _SK_FOUNDATION_FREQS:
            xl.extend([math.log(max(20, freq * 0.7)), 0.70, 0.0])
            xu.extend([math.log(min(80, freq * 1.3)), 0.92, 0.30])
        # Character M0 (dormant): 8 × (log_freq, radius, weight, notch_type)
        for lo, hi in _SK_CHAR_BOUNDS:
            xl.extend([math.log(lo), 0.05, 0.0, 0])
            xu.extend([math.log(hi), 0.30, 0.20, 3.99])
        # Character M100 (active): 8 × (log_freq, radius, weight, notch_type)
        for lo, hi in _SK_CHAR_BOUNDS:
            xl.extend([math.log(lo), 0.20, 0.03, 0])
            xu.extend([math.log(hi), 0.80, 0.60, 3.99])
        # Q deltas: 8 radius offsets
        xl.extend([0.0] * 8)
        xu.extend([0.10] * 8)

        self.objectives = CATEGORY_OBJECTIVES["speaker_knockerz"]
        super().__init__(n_var=84, n_obj=len(self.objectives),
                         n_ieq_constr=2,
                         xl=np.array(xl), xu=np.array(xu))

    def _evaluate(self, X, out, *args, **kwargs):
        F, G = [], []
        for x in X:
            try:
                body = _build_sk(x)
                metrics = _eval_body(body)
                row = [-metrics[obj] for obj in self.objectives]
                F.append(row)
                g1 = metrics["peak_db"] - 50.0
                gate = shipping_gate(body, "speaker_knockerz")
                if gate["passed"]:
                    g2 = -1.0
                else:
                    meas = gate["measurements"]
                    g2 = meas.get("worst_sub_drop_db", 6.0) - 6.0
                G.append([g1, g2])
            except Exception:
                F.append([0.0] * len(self.objectives))
                G.append([50.0, 50.0])
        out["F"] = np.array(F)
        out["G"] = np.array(G)


# ---------------------------------------------------------------------------
# Pymoo Problem definitions
# ---------------------------------------------------------------------------

CATEGORY_OBJECTIVES = {
    "vocal":  ["trajectory", "continuity", "flavor"],
    "drum":   ["trajectory", "continuity", "flavor"],
    "bass":   ["trajectory", "continuity", "flavor"],
    "pad":    ["trajectory", "continuity", "coherence"],
    "lead":   ["trajectory", "continuity", "flavor"],
    "808":    ["dynamic_range", "ruggedness", "morph_distance"],
    "hihat":  ["ridge", "dynamic_range", "morph_distance"],
    "lofi":   ["dynamic_range", "morph_distance", "spectral_tilt"],
    "pluck":  ["ridge", "dynamic_range", "morph_distance"],
    "fx":     ["ridge", "morph_distance", "ruggedness"],
    "speaker_knockerz": ["sub_contrast", "norm_morph_distance", "ridge"],
    "cul_de_sac":       ["topology_shift", "peak_evolution", "norm_morph_distance"],
}

STITCH_CATEGORIES = {"vocal", "drum", "bass", "pad", "lead"}


class StitchProblem(Problem):
    """Optimize stitched bodies with trajectory-first objectives.

    Constraints:
        g1: hard gate (peak > 50 dB or dead midpoint) → infeasible
        g2: occupancy ≥ 0.4 (at least 40 % of 5×5 bins alive)
    """

    def __init__(self, category: str):
        n_vowels = len(VOWEL_KEYS)
        xl = np.array([0, 0, 0.80, 0.40])
        xu = np.array([n_vowels - 0.01, n_vowels - 0.01, 0.99, 0.80])
        self.objectives = CATEGORY_OBJECTIVES[category]
        super().__init__(n_var=4, n_obj=len(self.objectives),
                         n_ieq_constr=2, xl=xl, xu=xu)

    def _evaluate(self, X, out, *args, **kwargs):
        F, G = [], []
        for x in X:
            try:
                body = _build_stitch(int(x[0]), int(x[1]), x[2], x[3])
                metrics = _eval_body(body)
                row = [-metrics[obj] for obj in self.objectives]
                F.append(row)
                g1 = 1.0 if metrics["gate_fail"] else -1.0
                g2 = 0.4 - metrics["occupancy"]
                G.append([g1, g2])
            except Exception:
                F.append([0.0] * len(self.objectives))
                G.append([1.0, 1.0])
        out["F"] = np.array(F)
        out["G"] = np.array(G)


_GATED_CATEGORIES = {"speaker_knockerz", "cul_de_sac"}


class TwelveStageProblem(Problem):
    """Optimize 12-stage bodies: search frequencies, radii, weights, notch types."""

    def __init__(self, category: str):
        # 96 params: 48 per morph endpoint (12 stages * 4 params each)
        # Per stage: [log_freq, radius, weight, notch_type]
        freq_bounds = self._freq_bounds(category)
        xl, xu = [], []
        for _ in range(2):  # two morph endpoints
            for i in range(12):
                xl.extend([math.log(freq_bounds[i][0]), 0.30, 0.03, 0])
                xu.extend([math.log(freq_bounds[i][1]), 0.95, 0.70, 3.99])
        self.objectives = CATEGORY_OBJECTIVES[category]
        self.category = category
        n_constr = 2 if category in _GATED_CATEGORIES else 0
        super().__init__(n_var=96, n_obj=len(self.objectives),
                         n_ieq_constr=n_constr,
                         xl=np.array(xl), xu=np.array(xu))

    def _freq_bounds(self, category: str) -> list[tuple[float, float]]:
        """Per-stage frequency bounds by category."""
        if category == "808":
            return [(20, 80)] * 4 + [(30, 500)] * 8
        elif category == "hihat":
            return [(400, 2000)] * 2 + [(1000, 6000)] * 4 + [(3000, 16000)] * 6
        elif category == "lofi":
            return [(60, 400)] * 3 + [(300, 3000)] * 3 + [(2000, 8000)] * 3 + [(5000, 18000)] * 3
        elif category == "pluck":
            return [(40, 400)] * 3 + [(200, 2000)] * 3 + [(1000, 8000)] * 3 + [(4000, 18000)] * 3
        elif category == "fx":
            return [(20, 200)] * 2 + [(100, 1000)] * 2 + [(500, 5000)] * 4 + [(3000, 18000)] * 4
        elif category == "speaker_knockerz":
            # Sub anchor (20-80), low-mid choke (100-400), mid rattle (400-1200),
            # presence cry (1k-5k), stress/fracture (4k-16k)
            return [(20, 80)] * 2 + [(100, 400)] * 2 + [(400, 1200)] * 2 + \
                   [(1000, 5000)] * 2 + [(2000, 8000)] * 2 + [(4000, 16000)] * 2
        elif category == "cul_de_sac":
            # Root tether (20-100), tube body (60-500), mid structure (300-3k),
            # comb region (1k-8k), shimmer (5k-18k)
            return [(20, 100)] * 2 + [(60, 500)] * 2 + [(300, 3000)] * 2 + \
                   [(1000, 8000)] * 2 + [(3000, 12000)] * 2 + [(5000, 18000)] * 2
        return [(20, 18000)] * 12

    def _evaluate(self, X, out, *args, **kwargs):
        gated = self.category in _GATED_CATEGORIES
        F, G = [], []
        for x in X:
            try:
                body = _build_12stage(x, self.category)
                metrics = _eval_body(body)
                row = [-metrics[obj] for obj in self.objectives]
                # Penalize extreme dynamic range for non-808 categories
                if self.category not in ("808", "speaker_knockerz") and metrics["dynamic_range"] > 100:
                    penalty = (metrics["dynamic_range"] - 100) * 0.05
                    row = [r + penalty for r in row]
                F.append(row)
                if gated:
                    # g1: peak ceiling — continuous overshoot above 50 dB (<=0 = feasible)
                    # P2K corpus p90 = 51.1 dB; 35 dB rejects 42% of real E-mu bodies
                    g1 = metrics["peak_db"] - 50.0
                    # g2: shipping gate — continuous violation magnitude
                    gate = shipping_gate(body, self.category)
                    if gate["passed"]:
                        g2 = -1.0
                    else:
                        meas = gate["measurements"]
                        g2 = meas.get("worst_sub_drop_db",
                             meas.get("worst_root_drop_db", 6.0)) - 6.0
                    G.append([g1, g2])
            except Exception:
                F.append([0.0] * len(self.objectives))
                if gated:
                    G.append([50.0, 50.0])
        out["F"] = np.array(F)
        if gated:
            out["G"] = np.array(G)


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def _plot_pareto(results, objectives, out_path):
    """Plot Pareto front in objective space."""
    F = -results.F  # un-negate for display
    fig, axes = plt.subplots(1, min(3, len(objectives) - 1), figsize=(5 * min(3, len(objectives) - 1), 4))
    if len(objectives) == 2:
        axes = [axes]
    for i, ax in enumerate(axes if hasattr(axes, '__iter__') else [axes]):
        if i + 1 < F.shape[1]:
            ax.scatter(F[:, 0], F[:, i + 1], c="steelblue", s=20, alpha=0.7)
            ax.set_xlabel(objectives[0])
            ax.set_ylabel(objectives[i + 1])
            ax.set_title(f"Pareto front")
            ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"  Pareto plot: {out_path}")


def _plot_responses(bodies: list[Body], names: list[str], out_path):
    """Overlay frequency responses of top bodies."""
    fig, ax = plt.subplots(figsize=(10, 5))
    for body, name in zip(bodies, names):
        enc = body.corners.interpolate(0.5, 0.5)
        db = cascade_response_db(enc, FREQS, SR)
        ax.semilogx(FREQS, db, label=name, alpha=0.7)
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Magnitude (dB)")
    ax.set_title("Top Candidates — Frequency Response @ M50 Q50")
    ax.legend(fontsize=7, ncol=2)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(20, SR / 2)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"  Response plot: {out_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def optimize_category(category: str, generations: int = 50, pop_size: int = 40) -> Path:
    """Run optimization for one category."""
    print(f"\n{'='*60}")
    print(f"  Optimizing: {category}")
    print(f"  Objectives: {CATEGORY_OBJECTIVES[category]}")
    print(f"  Generations: {generations}, Population: {pop_size}")
    print(f"{'='*60}")

    is_stitch = category in STITCH_CATEGORIES
    is_sk = category == "speaker_knockerz"
    if is_sk:
        problem = SKProblem()
    elif is_stitch:
        problem = StitchProblem(category)
    else:
        problem = TwelveStageProblem(category)

    algorithm = NSGA2(
        pop_size=pop_size,
        sampling=FloatRandomSampling(),
        crossover=SBX(prob=0.9, eta=15),
        mutation=PM(eta=20),
    )

    result = minimize(problem, algorithm, ("n_gen", generations),
                      seed=42, verbose=False)

    # Build Pareto-optimal bodies
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = VAULT_DIR / f"opt_{category}_{ts}"
    out_dir.mkdir(parents=True, exist_ok=True)

    bodies = []
    names = []
    metrics_list = []

    if result.X is None:
        print(f"\n  No feasible solutions found for {category}.")
        print(f"  Try more generations or larger population.")
        return out_dir

    for i, x in enumerate(result.X):
        if is_sk:
            body = _build_sk(x)
            name = f"Speaker_knockerz_{i:03d}"
        elif is_stitch:
            body = _build_stitch(int(x[0]), int(x[1]), x[2], x[3])
            va = VOWEL_KEYS[int(x[0]) % len(VOWEL_KEYS)]
            vb = VOWEL_KEYS[int(x[1]) % len(VOWEL_KEYS)]
            name = f"{category.capitalize()}_{va}_{vb}_{i:03d}"
        else:
            body = _build_12stage(x, category)
            name = f"{category.capitalize()}_{i:03d}"

        body_named = Body(name=name, corners=body.corners, boost=body.boost)
        metrics = _eval_body(body_named)

        path = out_dir / f"{name}.json"
        path.write_text(body_named.to_compiled_json())

        bodies.append(body_named)
        names.append(name)
        metrics_list.append(metrics)

    # --- Taste model scoring ---
    taste_scores = {}
    if _taste_scorer is not None:
        print(f"\n  Scoring {len(bodies)} survivors with taste model...")
        for body_obj, name in zip(bodies, names):
            try:
                tf = _extract_taste_features(body_obj)
                taste_scores[name] = _taste_scorer.score(tf)
            except Exception:
                taste_scores[name] = 0.0
    else:
        print("\n  (taste_model.json not found — skipping taste scoring)")

    # Sort by taste probability if available, else first objective
    obj0 = CATEGORY_OBJECTIVES[category][0]
    if taste_scores:
        ranked = sorted(zip(names, metrics_list),
                        key=lambda x: -taste_scores.get(x[0], 0.0))
    else:
        ranked = sorted(zip(names, metrics_list), key=lambda x: -x[1][obj0])

    print(f"\n  {len(bodies)} Pareto-optimal bodies saved to {out_dir}")
    if taste_scores:
        print(f"\n  Top 10 by taste probability:")
        for name, m in ranked[:10]:
            vals = "  ".join(f"{k}={m[k]:.2f}" for k in CATEGORY_OBJECTIVES[category])
            print(f"    {name}: taste={taste_scores[name]:.3f}  {vals}")
    else:
        print(f"\n  Top 10 by {obj0}:")
        for name, m in ranked[:10]:
            vals = "  ".join(f"{k}={m[k]:.2f}" for k in CATEGORY_OBJECTIVES[category])
            print(f"    {name}: {vals}")

    # Write taste ranking to file
    if taste_scores:
        ranking_path = out_dir / "taste_ranking.txt"
        with open(ranking_path, "w") as f:
            f.write(f"# Taste model ranking — {category}\n")
            f.write(f"# {len(bodies)} Pareto survivors, scored by taste_model.json\n\n")
            for i, (name, m) in enumerate(ranked):
                vals = "  ".join(f"{k}={m[k]:.2f}" for k in CATEGORY_OBJECTIVES[category])
                f.write(f"{i+1:3d}. {name}: taste={taste_scores[name]:.3f}  {vals}\n")
        print(f"  Taste ranking: {ranking_path}")

    # Plots
    _plot_pareto(result, CATEGORY_OBJECTIVES[category], out_dir / "pareto.png")
    top_bodies = [b for b, _ in sorted(zip(bodies, metrics_list),
                  key=lambda x: -taste_scores.get(x[0], -x[1][obj0]))[:10]]
    top_names = [n for n, _ in ranked[:10]]
    _plot_responses(top_bodies, top_names, out_dir / "responses.png")

    return out_dir


def main():
    parser = argparse.ArgumentParser(description="Forge optimizer")
    parser.add_argument("--category", default="vocal", help="Category to optimize")
    parser.add_argument("--generations", type=int, default=50)
    parser.add_argument("--pop", type=int, default=40)
    args = parser.parse_args()

    categories = list(CATEGORY_OBJECTIVES.keys()) if args.category == "all" else [args.category]

    for cat in categories:
        optimize_category(cat, args.generations, args.pop)


if __name__ == "__main__":
    main()
