"""P2k-trained stage-role vocabulary utilities.

Goal:
- Learn a structural role vocabulary from P2k bodies (not audio clip metrics).
- Score candidate bodies by distance to role targets.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np

from pyruntime.body import Body
from pyruntime.constants import SR
from pyruntime.corner import CornerName
from pyruntime.freq_response import freq_points, stage_response


CORNER_ORDER = (
    ("M0_Q0", CornerName.A),
    ("M0_Q100", CornerName.B),
    ("M100_Q0", CornerName.C),
    ("M100_Q100", CornerName.D),
)


@dataclass(frozen=True)
class StageFeature:
    peak_hz: float
    peak_db: float
    valley_hz: float
    valley_db: float
    range_db: float
    low_db: float
    lowmid_db: float
    mid_db: float
    high_db: float
    air_db: float
    tilt_db: float
    peak_count: int
    role: str
    passthrough: bool


def _band_mean(db: np.ndarray, freqs: np.ndarray, lo: float, hi: float) -> float:
    m = (freqs >= lo) & (freqs <= hi)
    if not np.any(m):
        return -120.0
    return float(np.mean(db[m]))


def _count_peaks(db: np.ndarray, freqs: np.ndarray, lo: float, hi: float, threshold_above_mean: float = 3.0) -> int:
    m = (freqs >= lo) & (freqs <= hi)
    v = db[m]
    if len(v) < 5:
        return 0
    mv = float(np.mean(v))
    count = 0
    for i in range(1, len(v) - 1):
        if v[i] > v[i - 1] and v[i] > v[i + 1] and v[i] > mv + threshold_above_mean:
            count += 1
    return int(count)


def _is_passthrough(enc) -> bool:
    return (
        abs(enc.c0 - 1.0) < 1e-6 and
        abs(enc.c1) < 1e-6 and
        abs(enc.c2) < 1e-6 and
        abs(enc.c3) < 1e-6 and
        abs(enc.c4) < 1e-6
    )


def _infer_role(
    *,
    peak_hz: float,
    range_db: float,
    valley_db: float,
    peak_count: int,
    passthrough: bool,
) -> str:
    if passthrough:
        return "passthrough"

    if peak_hz < 100.0:
        base = "sub_anchor"
    elif peak_hz < 300.0:
        base = "chest_mass"
    elif peak_hz < 900.0:
        base = "throat_choke"
    elif peak_hz < 2500.0:
        base = "formant_bite"
    elif peak_hz < 7000.0:
        base = "presence_shard"
    else:
        base = "air_edge"

    if range_db > 20.0 and valley_db < -10.0 and peak_count >= 2:
        return f"{base}_notch"
    return base


def stage_feature(enc, freqs: np.ndarray | None = None, sr: float = SR) -> StageFeature:
    if freqs is None:
        freqs = freq_points(n=512, sr=sr)
    h = stage_response(enc, freqs, sr)
    db = 20.0 * np.log10(np.maximum(np.abs(h), 1e-20))

    peak_idx = int(np.argmax(db))
    valley_idx = int(np.argmin(db))
    peak_hz = float(freqs[peak_idx])
    peak_db = float(db[peak_idx])
    valley_hz = float(freqs[valley_idx])
    valley_db = float(db[valley_idx])
    range_db = peak_db - valley_db
    low_db = _band_mean(db, freqs, 20.0, 120.0)
    lowmid_db = _band_mean(db, freqs, 120.0, 600.0)
    mid_db = _band_mean(db, freqs, 600.0, 2500.0)
    high_db = _band_mean(db, freqs, 2500.0, 9000.0)
    air_db = _band_mean(db, freqs, 9000.0, 18000.0)
    tilt_db = float(db[-1] - db[0]) if len(db) >= 2 else 0.0
    peak_count = _count_peaks(db, freqs, 100.0, 12000.0, threshold_above_mean=3.0)
    passthrough = _is_passthrough(enc)
    role = _infer_role(
        peak_hz=peak_hz,
        range_db=range_db,
        valley_db=valley_db,
        peak_count=peak_count,
        passthrough=passthrough,
    )

    return StageFeature(
        peak_hz=peak_hz,
        peak_db=peak_db,
        valley_hz=valley_hz,
        valley_db=valley_db,
        range_db=range_db,
        low_db=low_db,
        lowmid_db=lowmid_db,
        mid_db=mid_db,
        high_db=high_db,
        air_db=air_db,
        tilt_db=tilt_db,
        peak_count=peak_count,
        role=role,
        passthrough=passthrough,
    )


def body_signature(body: Body, sr: float = SR) -> dict:
    freqs = freq_points(n=512, sr=sr)
    out = {}
    for label, corner_name in CORNER_ORDER:
        enc = body.corners.corner(corner_name).encode()
        stages = []
        for i, s in enumerate(enc[:6], start=1):
            f = stage_feature(s, freqs=freqs, sr=sr)
            stages.append({
                "stage": i,
                "role": f.role,
                "passthrough": f.passthrough,
                "peak_hz": f.peak_hz,
                "peak_db": f.peak_db,
                "valley_hz": f.valley_hz,
                "valley_db": f.valley_db,
                "range_db": f.range_db,
                "low_db": f.low_db,
                "lowmid_db": f.lowmid_db,
                "mid_db": f.mid_db,
                "high_db": f.high_db,
                "air_db": f.air_db,
                "tilt_db": f.tilt_db,
                "peak_count": f.peak_count,
            })
        out[label] = stages
    return out


def _numeric_vector(stage: dict) -> np.ndarray:
    return np.array([
        float(stage["peak_hz"]),
        float(stage["peak_db"]),
        float(stage["valley_hz"]),
        float(stage["valley_db"]),
        float(stage["range_db"]),
        float(stage["low_db"]),
        float(stage["lowmid_db"]),
        float(stage["mid_db"]),
        float(stage["high_db"]),
        float(stage["air_db"]),
        float(stage["tilt_db"]),
        float(stage["peak_count"]),
    ], dtype=np.float64)


def role_distance(signature: dict, target_signature: dict) -> dict:
    """Distance score between two body signatures.

    Lower is better.
    """
    corner_scores = {}
    role_mismatch = 0
    l2_total = 0.0
    count = 0

    for label in signature.keys():
        stages_a = signature[label]
        stages_b = target_signature[label]
        corner_l2 = 0.0
        for sa, sb in zip(stages_a, stages_b):
            va = _numeric_vector(sa)
            vb = _numeric_vector(sb)
            # Lightly normalize major frequency axes.
            va[[0, 2]] /= 1000.0
            vb[[0, 2]] /= 1000.0
            d = float(np.linalg.norm(va - vb))
            corner_l2 += d
            l2_total += d
            count += 1
            if sa["role"] != sb["role"]:
                role_mismatch += 1
        corner_scores[label] = corner_l2 / max(1, len(stages_a))

    mean_l2 = l2_total / max(1, count)
    mismatch_rate = role_mismatch / max(1, count)
    # composite: role mismatch is primary, l2 refines tie-break.
    composite = mismatch_rate * 10.0 + mean_l2
    return {
        "composite_distance": composite,
        "role_mismatch_rate": mismatch_rate,
        "mean_l2": mean_l2,
        "corner_mean_l2": corner_scores,
        "stage_count": count,
    }


def summarize_role_vocabulary(signatures: Iterable[tuple[str, dict]]) -> dict:
    """Build aggregate vocabulary stats from `(name, signature)` rows."""
    by_role: dict[str, list[np.ndarray]] = {}
    by_stage: dict[str, dict[str, int]] = {}
    body_sequences: dict[str, dict[str, list[str]]] = {}

    for body_name, sig in signatures:
        body_sequences[body_name] = {}
        for label, stages in sig.items():
            seq = []
            by_stage.setdefault(label, {})
            for i, stage in enumerate(stages, start=1):
                role = stage["role"]
                seq.append(role)
                by_stage[label][role] = by_stage[label].get(role, 0) + 1
                by_role.setdefault(role, []).append(_numeric_vector(stage))
            body_sequences[body_name][label] = seq

    role_prototypes = {}
    for role, vectors in by_role.items():
        arr = np.vstack(vectors)
        med = np.median(arr, axis=0)
        role_prototypes[role] = {
            "peak_hz": float(med[0]),
            "peak_db": float(med[1]),
            "valley_hz": float(med[2]),
            "valley_db": float(med[3]),
            "range_db": float(med[4]),
            "low_db": float(med[5]),
            "lowmid_db": float(med[6]),
            "mid_db": float(med[7]),
            "high_db": float(med[8]),
            "air_db": float(med[9]),
            "tilt_db": float(med[10]),
            "peak_count": float(med[11]),
            "count": int(arr.shape[0]),
        }

    return {
        "role_prototypes": role_prototypes,
        "stage_role_distribution": by_stage,
        "body_role_sequences": body_sequences,
    }
