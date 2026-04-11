"""Body analysis — motion, crossings, cancellations, profile.

Judges body quality like a Z-plane sound designer: morph trajectory
distance, spectral zero crossings, pole-zero proximity, composite profile.
"""
from __future__ import annotations

import math

import numpy as np

from pyruntime.body import Body
from pyruntime.constants import SR
from pyruntime.corner import CornerName, CornerState
from pyruntime.encode import EncodedCoeffs
from pyruntime.freq_response import cascade_response_db, freq_points
from pyruntime.stage_params import StageParams


def morph_trajectory_distance(body: Body) -> dict:
    """RMS spectral distance between M0 and M100 corners.

    Measures at both Q=0 and Q=1.  Higher = more morph motion.
    Zero = static filter (morph knob does nothing).

    Returns {"q0_rms_db": float, "q1_rms_db": float, "mean_rms_db": float}.
    """
    freqs = freq_points(sr=SR)
    sr = SR

    # M0_Q0 vs M100_Q0
    enc_m0_q0 = body.corners.interpolate(0.0, 0.0)
    enc_m100_q0 = body.corners.interpolate(1.0, 0.0)
    db_m0_q0 = cascade_response_db(enc_m0_q0, freqs, sr)
    db_m100_q0 = cascade_response_db(enc_m100_q0, freqs, sr)
    q0_rms = _rms_distance(db_m0_q0, db_m100_q0)

    # M0_Q1 vs M100_Q1
    enc_m0_q1 = body.corners.interpolate(0.0, 1.0)
    enc_m100_q1 = body.corners.interpolate(1.0, 1.0)
    db_m0_q1 = cascade_response_db(enc_m0_q1, freqs, sr)
    db_m100_q1 = cascade_response_db(enc_m100_q1, freqs, sr)
    q1_rms = _rms_distance(db_m0_q1, db_m100_q1)

    mean_rms = (q0_rms + q1_rms) / 2.0

    return {
        "q0_rms_db": float(q0_rms),
        "q1_rms_db": float(q1_rms),
        "mean_rms_db": float(mean_rms),
    }


def zero_crossing_count(encoded_stages: list, freqs: np.ndarray, sr: float) -> int:
    """Count where cascade magnitude response crosses 0 dB.

    More crossings = more complex spectral shape.
    """
    if not encoded_stages:
        return 0
    db = cascade_response_db(encoded_stages, freqs, sr)
    if len(db) < 2:
        return 0
    # Sign changes in the dB array (positive = above 0 dB, negative = below)
    signs = np.sign(db)
    # A crossing occurs where consecutive signs differ (and neither is zero)
    diffs = np.diff(signs)
    return int(np.count_nonzero(diffs))


def pole_zero_proximity(stage: StageParams) -> float:
    """Angular distance between pole and zero for a single stage.

    Close to 0 = cancellation.  Large = reinforcement.
    Passthrough stages (r < 0.01) return pi (not meaningful).
    """
    if stage.r < 0.01:
        return math.pi

    r = stage.r
    a1 = stage.a1

    # Pole angle: acos(-a1 / (2*r))
    pole_arg = -a1 / (2.0 * r)
    pole_arg = max(-1.0, min(1.0, pole_arg))
    pole_angle = math.acos(pole_arg)

    # Zero coefficients: b1 = a1 + val2, b2 = r^2 - val3
    b1 = a1 + stage.val2
    b2 = r * r - stage.val3

    if b2 <= 0.0:
        # Zero radius is zero or imaginary — no meaningful angle
        return math.pi

    zero_r = math.sqrt(b2)
    zero_arg = -b1 / (2.0 * zero_r)
    zero_arg = max(-1.0, min(1.0, zero_arg))
    zero_angle = math.acos(zero_arg)

    return abs(pole_angle - zero_angle)


def body_profile(body: Body) -> dict:
    """Composite analysis combining all metrics.

    Returns a dict with morph distance, per-corner crossing counts,
    stage proximity scores, active stage count, and spectral tilt.
    """
    freqs = freq_points(sr=SR)
    sr = SR

    # Morph trajectory
    morph_dist = morph_trajectory_distance(body)

    # Encoded corners for crossing counts
    enc_m0_q0 = body.corners.interpolate(0.0, 0.0)
    enc_m0_q100 = body.corners.interpolate(0.0, 1.0)
    enc_m100_q0 = body.corners.interpolate(1.0, 0.0)
    enc_m100_q100 = body.corners.interpolate(1.0, 1.0)

    crossings = {
        "m0_q0": zero_crossing_count(enc_m0_q0, freqs, sr),
        "m0_q100": zero_crossing_count(enc_m0_q100, freqs, sr),
        "m100_q0": zero_crossing_count(enc_m100_q0, freqs, sr),
        "m100_q100": zero_crossing_count(enc_m100_q100, freqs, sr),
    }

    # Stage proximity from raw StageParams at corner A (M0_Q0)
    corner_a = body.corners.corner(CornerName.A)
    raw_stages = corner_a.stages

    stage_prox = [pole_zero_proximity(s) for s in raw_stages]

    active = sum(1 for s in raw_stages if s.r > 0.01)

    # Spectral tilt at M0_Q0
    db_a = cascade_response_db(enc_m0_q0, freqs, sr)
    if len(db_a) >= 2:
        tilt = float(db_a[-1] - db_a[0])
    else:
        tilt = 0.0

    audit = midpoint_audit(body)

    return {
        "morph_distance": morph_dist,
        "crossings": crossings,
        "stage_proximity": stage_prox,
        "active_stages": active,
        "spectral_tilt_db": tilt,
        "midpoint_audit": audit,
    }


def midpoint_audit(body: Body, peak_limit_db: float = 40.0) -> dict:
    """Audit the morph trajectory for resonance spikes.

    The midpoint collapse: if stages are scrambled between corners,
    interpolating coefficients at 25%/50%/75% morph can cause poles
    to tear across the unit circle, producing catastrophic gain stacking.

    Tests morph at 0.0, 0.25, 0.5, 0.75, 1.0 for both Q=0 and Q=1.
    Returns pass/fail with peak dB at each test point.

    peak_limit_db: maximum allowed cascade peak (default 35 dB).
    """
    freqs = freq_points(sr=SR)
    sr = SR
    test_morphs = [0.0, 0.25, 0.5, 0.75, 1.0]
    test_qs = [0.0, 1.0]

    worst_peak = -999.0
    worst_point = (0.0, 0.0)
    points = []

    for m in test_morphs:
        for q in test_qs:
            enc = body.corners.interpolate(m, q)
            db = cascade_response_db(enc, freqs, sr)
            peak = float(np.max(db))
            points.append({"morph": m, "q": q, "peak_db": round(peak, 1)})
            if peak > worst_peak:
                worst_peak = peak
                worst_point = (m, q)

    passed = worst_peak <= peak_limit_db

    return {
        "passed": passed,
        "worst_peak_db": round(worst_peak, 1),
        "worst_point": {"morph": worst_point[0], "q": worst_point[1]},
        "peak_limit_db": peak_limit_db,
        "points": points,
    }


def _rms_distance(a: np.ndarray, b: np.ndarray) -> float:
    """RMS of element-wise difference between two dB arrays."""
    diff = a - b
    return float(np.sqrt(np.mean(diff * diff)))


def _normalize_db(db: np.ndarray) -> np.ndarray:
    """Shift a dB curve so its mean is 0. Removes gross level inflation."""
    return db - float(np.mean(db))


def morph_trajectory_distance_normalized(body: Body) -> dict:
    """Level-normalized morph distance.

    Same as morph_trajectory_distance but normalizes each endpoint's
    response curve to zero mean before computing RMS difference.
    This prevents candidates from winning by simply being louder
    at one endpoint.

    Returns raw and normalized distances for both.
    """
    freqs = freq_points(sr=SR)
    sr = SR

    enc_m0_q0 = body.corners.interpolate(0.0, 0.0)
    enc_m100_q0 = body.corners.interpolate(1.0, 0.0)
    db_m0_q0 = cascade_response_db(enc_m0_q0, freqs, sr)
    db_m100_q0 = cascade_response_db(enc_m100_q0, freqs, sr)

    raw_q0 = _rms_distance(db_m0_q0, db_m100_q0)
    norm_q0 = _rms_distance(_normalize_db(db_m0_q0), _normalize_db(db_m100_q0))

    return {
        "raw_rms_db": float(raw_q0),
        "normalized_rms_db": float(norm_q0),
    }


def dense_midpoint_audit(body: Body, peak_limit_db: float = 40.0,
                         n_morph: int = 21, n_q: int = 3) -> dict:
    """Dense sweep midpoint audit.

    Sweeps a fine morph grid (default 21 points: 0.0, 0.05, ..., 1.0)
    at multiple Q values. Reports worst-case peak anywhere in the sweep.
    """
    freqs = freq_points(sr=SR)
    sr = SR
    test_morphs = [i / (n_morph - 1) for i in range(n_morph)]
    test_qs = [i / max(1, n_q - 1) for i in range(n_q)]

    worst_peak = -999.0
    worst_point = (0.0, 0.0)

    for m in test_morphs:
        for q in test_qs:
            enc = body.corners.interpolate(m, q)
            db = cascade_response_db(enc, freqs, sr)
            peak = float(np.max(db))
            if peak > worst_peak:
                worst_peak = peak
                worst_point = (m, q)

    return {
        "passed": worst_peak <= peak_limit_db,
        "worst_peak_db": round(worst_peak, 1),
        "worst_point": {"morph": round(worst_point[0], 3), "q": round(worst_point[1], 3)},
        "peak_limit_db": peak_limit_db,
        "grid_size": n_morph * n_q,
    }


# ---------------------------------------------------------------------------
# Shipping invariant gates
#
# Hard pass/fail derived from canonical shipping body names only.
# Prose is non-normative; invariants are explicit, measurable assertions.
# ---------------------------------------------------------------------------

def _count_peaks_in_band(db: np.ndarray, freqs: np.ndarray,
                         lo: float, hi: float, threshold_above_mean: float = 3.0) -> int:
    """Count spectral peaks in a frequency band above local mean."""
    mask = (freqs >= lo) & (freqs <= hi)
    band = db[mask]
    if len(band) < 3:
        return 0
    mean_val = float(np.mean(band))
    count = 0
    for i in range(1, len(band) - 1):
        if band[i] > band[i-1] and band[i] > band[i+1] and band[i] > mean_val + threshold_above_mean:
            count += 1
    return count


def _band_energy(db: np.ndarray, freqs: np.ndarray, lo: float, hi: float) -> float:
    """Mean dB in a frequency band."""
    mask = (freqs >= lo) & (freqs <= hi)
    band = db[mask]
    return float(np.mean(band)) if len(band) > 0 else -120.0


def _band_peak(db: np.ndarray, freqs: np.ndarray, lo: float, hi: float) -> float:
    """Peak dB in a frequency band."""
    mask = (freqs >= lo) & (freqs <= hi)
    band = db[mask]
    return float(np.max(band)) if len(band) > 0 else -120.0


def _band_valley(db: np.ndarray, freqs: np.ndarray, lo: float, hi: float) -> float:
    """Valley dB in a frequency band."""
    mask = (freqs >= lo) & (freqs <= hi)
    band = db[mask]
    return float(np.min(band)) if len(band) > 0 else -120.0


# ---------------------------------------------------------------------------
# World Q contract — role bands, world mapping, hard gates
# ---------------------------------------------------------------------------

_ROLE_BANDS = [
    ("sub", 20.0, 200.0),
    ("low_mid", 200.0, 900.0),
    ("character", 900.0, 3000.0),
    ("air", 3000.0, 16000.0),
]

_BODY_WORLD = {
    "speaker_knockerz": "bass_pressure",
    "aluminum_siding": "brittle_air",
    "small_talk_ah_ee": "dual_formant",
    "cul_de_sac": "null_fracture",
}

_WORLD_Q_CONTRACTS = {
    "bass_pressure": {
        "min_dominant_peaks_qmax": 2,
        "min_occupied_bands_qmax": 2,
        "max_domination_ratio_qmax": 0.70,
        "hf_takeover_forbidden": True,
        "hf_takeover_hz": 8000.0,
    },
    "brittle_air": {
        "min_dominant_peaks_qmax": 1,
        "min_occupied_bands_qmax": 1,
        "max_domination_ratio_qmax": 0.80,
        "hf_takeover_forbidden": False,
    },
    "dual_formant": {
        "min_dominant_peaks_qmax": 2,
        "min_occupied_bands_qmax": 2,
        "max_domination_ratio_qmax": 0.65,
        "hf_takeover_forbidden": False,
    },
    "null_fracture": {
        "min_dominant_peaks_qmax": 1,
        "min_occupied_bands_qmax": 1,
        "max_domination_ratio_qmax": 0.85,
        "hf_takeover_forbidden": False,
    },
}

_MIN_Q_VITALITY_RMS_DB = 1.5


def _q_vitality(body: Body, test_morphs: list[float]) -> dict:
    """Measure how much the Q axis changes the response.

    RMS spectral distance between q=0 and q=1 at each morph position.
    A dead Q axis produces near-zero values.
    """
    freqs = freq_points(sr=SR)
    sr = SR
    deltas = []
    for m in test_morphs:
        db_q0 = cascade_response_db(body.corners.interpolate(m, 0.0), freqs, sr)
        db_q1 = cascade_response_db(body.corners.interpolate(m, 1.0), freqs, sr)
        delta = _rms_distance(db_q0, db_q1)
        deltas.append({"morph": m, "q_delta_rms_db": round(delta, 2)})
    max_delta = max(d["q_delta_rms_db"] for d in deltas)
    mean_delta = sum(d["q_delta_rms_db"] for d in deltas) / len(deltas)
    return {
        "max_q_delta_rms_db": round(max_delta, 2),
        "mean_q_delta_rms_db": round(mean_delta, 2),
        "per_morph": deltas,
    }


def _occupied_bands(db: np.ndarray, freqs: np.ndarray,
                    dynamic_range_db: float = 20.0) -> list[str]:
    """Return names of role bands with significant energy.

    A band is occupied if its mean energy is within dynamic_range_db
    of the loudest band. This catches both peaks and shelves.
    """
    band_energies = []
    for name, lo, hi in _ROLE_BANDS:
        energy = _band_energy(db, freqs, lo, hi)
        band_energies.append((name, energy))
    if not band_energies:
        return []
    max_energy = max(e for _, e in band_energies)
    return [name for name, e in band_energies if e >= max_energy - dynamic_range_db]


def _domination_ratio(db: np.ndarray, freqs: np.ndarray,
                      lo: float = 20.0, hi: float = 16000.0) -> float:
    """Approximate single-peak domination ratio.

    Fraction of total peak prominence concentrated in the strongest peak.
    High values (>0.7) indicate single-stage takeover.
    """
    mask = (freqs >= lo) & (freqs <= hi)
    band = db[mask]
    if len(band) < 5:
        return 1.0
    mean_val = float(np.mean(band))
    peak_heights = []
    for i in range(1, len(band) - 1):
        if band[i] > band[i - 1] and band[i] > band[i + 1] and band[i] > mean_val:
            peak_heights.append(float(band[i] - mean_val))
    if not peak_heights:
        return 1.0
    total = sum(peak_heights)
    if total <= 0:
        return 1.0
    return max(peak_heights) / total


def _canonical_shipping_name(body_name: str) -> str | None:
    token = body_name.strip().lower()
    mapping = {
        "speaker knockerz": "speaker_knockerz",
        "aluminum siding": "aluminum_siding",
        "small talk ah-ee": "small_talk_ah_ee",
        "cul-de-sac": "cul_de_sac",
    }
    return mapping.get(token)


def _top_peak_freqs(
    db: np.ndarray,
    freqs: np.ndarray,
    lo: float,
    hi: float,
    *,
    threshold_above_mean: float,
    max_peaks: int = 8,
) -> list[float]:
    mask = (freqs >= lo) & (freqs <= hi)
    f = freqs[mask]
    v = db[mask]
    if len(v) < 5:
        return []
    mean_v = float(np.mean(v))
    peaks: list[tuple[float, float]] = []
    for i in range(1, len(v) - 1):
        if v[i] > v[i - 1] and v[i] > v[i + 1] and v[i] > mean_v + threshold_above_mean:
            peaks.append((float(v[i]), float(f[i])))
    peaks.sort(key=lambda p: -p[0])
    return [freq for _, freq in peaks[:max_peaks]]


def shipping_gate(body: Body, body_name: str) -> dict:
    """Run the shipping gate using canonical body names.

    Returns {"passed": bool, "failures": [...], "measurements": {...}}.
    Unknown names fail closed.
    """
    freqs = freq_points(sr=SR)
    sr = SR
    test_morphs = [0.0, 0.25, 0.5, 0.75, 1.0]
    test_qs = [0.0, 0.5, 1.0]
    failures = []
    measurements: dict[str, object] = {}
    canonical = _canonical_shipping_name(body_name)
    measurements["canonical_name"] = canonical
    measurements["test_grid_morphs"] = test_morphs
    measurements["test_grid_qs"] = test_qs

    if canonical is None:
        failures.append(f"Unknown shipping body name '{body_name}'")
        return {
            "passed": False,
            "body": body_name,
            "canonical_body": None,
            "failures": failures,
            "measurements": measurements,
        }

    if canonical == "speaker_knockerz":
        # Name contract:
        # - sub anchor never disappears across morph/q
        # - danger zone adds upper-bass/low-mid bite
        ref_enc = body.corners.interpolate(0.0, 0.0)
        ref_db = cascade_response_db(ref_enc, freqs, sr)
        ref_sub = _band_peak(ref_db, freqs, 20.0, 60.0)
        measurements["ref_sub_peak_db"] = round(ref_sub, 1)

        worst_drop = 0.0
        worst_point = (0.0, 0.0)
        for m in test_morphs:
            for q in test_qs:
                enc = body.corners.interpolate(m, q)
                db = cascade_response_db(enc, freqs, sr)
                sub_peak = _band_peak(db, freqs, 20.0, 60.0)
                drop = ref_sub - sub_peak
                if drop > worst_drop:
                    worst_drop = drop
                    worst_point = (m, q)

        measurements["worst_sub_drop_db"] = round(worst_drop, 1)
        measurements["worst_sub_point"] = worst_point
        if worst_drop > 6.0:
            failures.append(
                f"Sub anchor failed: {worst_drop:.1f} dB drop at morph={worst_point[0]}, Q={worst_point[1]} (limit: 6 dB)")

        calm = cascade_response_db(body.corners.interpolate(0.0, 0.0), freqs, sr)
        danger = cascade_response_db(body.corners.interpolate(0.75, 0.0), freqs, sr)
        bite_gain = _band_energy(danger, freqs, 700.0, 2500.0) - _band_energy(calm, freqs, 700.0, 2500.0)
        measurements["danger_bite_gain_db"] = round(float(bite_gain), 2)
        if bite_gain < 1.5:
            failures.append(f"Danger bite failed: {bite_gain:.2f} dB (limit: 1.5 dB)")

    elif canonical == "aluminum_siding":
        # Name contract:
        # - 1k void remains scooped
        # - high-air zone can open
        worst_scoop = 999.0
        worst_point = (0.0, 0.0)
        for m in test_morphs:
            for q in test_qs:
                enc = body.corners.interpolate(m, q)
                db = cascade_response_db(enc, freqs, sr)
                void_energy = _band_energy(db, freqs, 800.0, 1200.0)
                hf_peak = _band_peak(db, freqs, 3000.0, 15000.0)
                scoop = hf_peak - void_energy
                if scoop < worst_scoop:
                    worst_scoop = scoop
                    worst_point = (m, q)

        measurements["worst_void_scoop_db"] = round(worst_scoop, 1)
        measurements["worst_void_point"] = worst_point
        if worst_scoop < 10.0:
            failures.append(
                f"1kHz void failed: only {worst_scoop:.1f} dB scoop at morph={worst_point[0]}, Q={worst_point[1]} (limit: 10 dB)")

        # Open-state air check.
        enc_open = body.corners.interpolate(1.0, 0.0)
        db_open = cascade_response_db(enc_open, freqs, sr)
        db_calm = cascade_response_db(body.corners.interpolate(0.0, 0.0), freqs, sr)
        air_gain = _band_energy(db_open, freqs, 7000.0, 16000.0) - _band_energy(db_calm, freqs, 7000.0, 16000.0)
        measurements["open_air_gain_db"] = round(float(air_gain), 2)
        if air_gain < 1.0:
            failures.append(f"Open-air gain failed: {air_gain:.2f} dB (limit: 1.0 dB)")

    elif canonical == "small_talk_ah_ee":
        # Name contract:
        # - two dominant formant peaks should be present and separated
        min_peak_count = 999
        max_peak_count = -1
        min_separation = 1e9
        worst_point = (0.0, 0.0)
        for m in [0.0, 0.5, 1.0]:
            for q in [0.0, 0.5, 1.0]:
                enc = body.corners.interpolate(m, q)
                db = cascade_response_db(enc, freqs, sr)
                peak_freqs = _top_peak_freqs(
                    db, freqs, 200.0, 5000.0,
                    threshold_above_mean=6.0,
                    max_peaks=5,
                )
                n_peaks = len(peak_freqs)
                if n_peaks < min_peak_count:
                    min_peak_count = n_peaks
                    worst_point = (m, q)
                max_peak_count = max(max_peak_count, n_peaks)
                if n_peaks >= 2:
                    sep = abs(peak_freqs[0] - peak_freqs[1])
                    min_separation = min(min_separation, sep)

        measurements["min_formant_peaks"] = int(min_peak_count)
        measurements["max_formant_peaks"] = int(max_peak_count)
        measurements["min_peaks_point"] = worst_point
        measurements["min_formant_separation_hz"] = round(float(min_separation if min_separation < 1e9 else 0.0), 1)
        if min_peak_count < 2:
            failures.append(
                f"Formant floor failed: only {min_peak_count} dominant peaks at morph={worst_point[0]}, Q={worst_point[1]} (need >= 2)")
        if max_peak_count > 4:
            failures.append(f"Formant ceiling failed: {max_peak_count} dominant peaks (limit: 4)")
        if min_separation < 250.0:
            failures.append(f"Formant separation failed: {min_separation:.1f} Hz (limit: 250 Hz)")

    elif canonical == "cul_de_sac":
        # Name contract:
        # - root tether persists
        # - a null event exists near midpoint
        # - fracture state exhibits comb-like multi-peak structure
        ref_enc = body.corners.interpolate(0.0, 0.0)
        ref_db = cascade_response_db(ref_enc, freqs, sr)
        ref_root = _band_peak(ref_db, freqs, 20.0, 200.0)
        measurements["ref_root_peak_db"] = round(ref_root, 1)

        worst_drop = 0.0
        worst_point = (0.0, 0.0)
        for m in test_morphs:
            for q in test_qs:
                enc = body.corners.interpolate(m, q)
                db = cascade_response_db(enc, freqs, sr)
                root_peak = _band_peak(db, freqs, 20.0, 200.0)
                drop = ref_root - root_peak
                if drop > worst_drop:
                    worst_drop = drop
                    worst_point = (m, q)

        measurements["worst_root_drop_db"] = round(worst_drop, 1)
        measurements["worst_root_point"] = worst_point
        if worst_drop > 6.0:
            failures.append(
                f"Root tether failed: {worst_drop:.1f} dB drop at morph={worst_point[0]}, Q={worst_point[1]} (limit: 6 dB)")

        calm = cascade_response_db(body.corners.interpolate(0.0, 0.0), freqs, sr)
        boundary = cascade_response_db(body.corners.interpolate(0.5, 0.0), freqs, sr)
        fracture = cascade_response_db(body.corners.interpolate(1.0, 0.0), freqs, sr)
        calm_energy = _band_energy(calm, freqs, 200.0, 8000.0)
        boundary_energy = _band_energy(boundary, freqs, 200.0, 8000.0)
        fracture_energy = _band_energy(fracture, freqs, 200.0, 8000.0)
        null_depth_calm = calm_energy - boundary_energy
        null_depth_fracture = fracture_energy - boundary_energy
        measurements["null_depth_vs_calm_db"] = round(float(null_depth_calm), 2)
        measurements["null_depth_vs_fracture_db"] = round(float(null_depth_fracture), 2)
        if null_depth_calm < 3.0 or null_depth_fracture < 3.0:
            failures.append(
                f"Null event failed: depth vs calm={null_depth_calm:.2f} dB, vs fracture={null_depth_fracture:.2f} dB (limit: 3 dB both)")

        comb_peaks = _count_peaks_in_band(fracture, freqs, 1200.0, 12000.0, threshold_above_mean=2.5)
        measurements["fracture_comb_peak_count"] = int(comb_peaks)
        if comb_peaks < 6:
            failures.append(f"Comb structure failed: only {comb_peaks} peaks at fracture state (limit: 6)")

    # ---- World Q contract gates (docs/world_q_behavior_gate_spec.md) ----
    world = _BODY_WORLD.get(canonical)
    if world is not None:
        contract = _WORLD_Q_CONTRACTS[world]

        # Q-axis vitality — is the Q knob doing anything?
        vitality = _q_vitality(body, test_morphs)
        measurements["q_vitality"] = vitality
        if vitality["max_q_delta_rms_db"] < _MIN_Q_VITALITY_RMS_DB:
            failures.append(
                f"Q AXIS IS DEAD: max delta {vitality['max_q_delta_rms_db']:.1f} dB RMS "
                f"(limit: {_MIN_Q_VITALITY_RMS_DB:.1f} dB)")

        # Per-world checks at q_max across morph positions
        worst_peaks_qmax = 999
        worst_peaks_morph = 0.0
        worst_bands_qmax = 999
        worst_bands_morph = 0.0
        worst_bands_names: list[str] = []
        worst_domination = 0.0
        worst_dom_morph = 0.0
        hf_takeover_detected = False
        hf_takeover_morph = 0.0

        for m in test_morphs:
            enc = body.corners.interpolate(m, 1.0)
            db_qmax = cascade_response_db(enc, freqs, sr)

            peaks = _top_peak_freqs(
                db_qmax, freqs, 20.0, 16000.0,
                threshold_above_mean=6.0, max_peaks=10)
            if len(peaks) < worst_peaks_qmax:
                worst_peaks_qmax = len(peaks)
                worst_peaks_morph = m

            bands = _occupied_bands(db_qmax, freqs)
            if len(bands) < worst_bands_qmax:
                worst_bands_qmax = len(bands)
                worst_bands_morph = m
                worst_bands_names = bands

            dom = _domination_ratio(db_qmax, freqs)
            if dom > worst_domination:
                worst_domination = dom
                worst_dom_morph = m

            if contract.get("hf_takeover_forbidden", False) and peaks:
                hf_limit = contract.get("hf_takeover_hz", 8000.0)
                if peaks[0] > hf_limit:
                    hf_takeover_detected = True
                    hf_takeover_morph = m

        measurements["q_contract"] = {
            "world": world,
            "worst_dominant_peaks_qmax": worst_peaks_qmax,
            "worst_peaks_morph": worst_peaks_morph,
            "worst_occupied_bands_qmax": worst_bands_qmax,
            "worst_bands_morph": worst_bands_morph,
            "worst_bands_names": worst_bands_names,
            "worst_domination_ratio_qmax": round(worst_domination, 3),
            "worst_domination_morph": worst_dom_morph,
            "hf_takeover": hf_takeover_detected,
        }

        min_peaks = contract["min_dominant_peaks_qmax"]
        if worst_peaks_qmax < min_peaks:
            failures.append(
                f"Q MAX COLLAPSED TO {worst_peaks_qmax} PEAK(S) "
                f"at morph={worst_peaks_morph:.2f} "
                f"({world} requires >= {min_peaks})")

        min_bands = contract["min_occupied_bands_qmax"]
        if worst_bands_qmax < min_bands:
            failures.append(
                f"Q LOST FREQUENCY BANDS: {worst_bands_qmax} occupied "
                f"[{', '.join(worst_bands_names)}] "
                f"at morph={worst_bands_morph:.2f} "
                f"({world} requires >= {min_bands})")

        max_dom = contract["max_domination_ratio_qmax"]
        if worst_domination > max_dom:
            failures.append(
                f"SINGLE PEAK DOMINATION: ratio {worst_domination:.2f} "
                f"at morph={worst_dom_morph:.2f} "
                f"({world} limit: {max_dom:.2f})")

        if hf_takeover_detected:
            failures.append(
                f"HF TAKEOVER at morph={hf_takeover_morph:.2f}: "
                f"dominant peak above {contract.get('hf_takeover_hz', 8000.0):.0f} Hz "
                f"(forbidden for {world})")

    return {
        "passed": len(failures) == 0,
        "body": body_name,
        "canonical_body": canonical,
        "failures": failures,
        "measurements": measurements,
    }
