"""Utilities for working with compiled-v1 cartridge JSON without re-encoding.

Why this exists:
- `Body.from_dict()` can round-trip compiled-v1 coefficients through StageParams.
- That is useful for authoring, but it is not "truth-preserving" for extracted
  cartridges when you want exact c0..c4 behavior.

This module provides:
- A strict compiled-v1 detector
- Exact bilinear interpolation of c0..c4 directly from keyframes
"""

from __future__ import annotations

from typing import Any

from pyruntime.encode import EncodedCoeffs

_CORNERS = ("M0_Q0", "M0_Q100", "M100_Q0", "M100_Q100")


def is_compiled_v1(body_dict: dict[str, Any]) -> bool:
    if not isinstance(body_dict, dict):
        return False
    if body_dict.get("format") != "compiled-v1":
        return False
    keyframes = body_dict.get("keyframes")
    return isinstance(keyframes, list) and len(keyframes) >= 4


def _clamp01(x: float) -> float:
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    return float(x)


def _extract_keyframes(body_dict: dict[str, Any]) -> dict[str, list[dict[str, float]]]:
    keyframes = body_dict.get("keyframes")
    if not isinstance(keyframes, list):
        raise ValueError("compiled-v1 body missing keyframes list")

    by_label: dict[str, list[dict[str, float]]] = {}
    for kf in keyframes:
        if not isinstance(kf, dict):
            continue
        label = kf.get("label")
        if label not in _CORNERS:
            continue
        stages = kf.get("stages")
        if not isinstance(stages, list) or not stages:
            continue
        stage_rows: list[dict[str, float]] = []
        for s in stages:
            if not isinstance(s, dict):
                stage_rows.append({"c0": 1.0, "c1": 0.0, "c2": 0.0, "c3": 0.0, "c4": 0.0})
                continue
            stage_rows.append(
                {
                    "c0": float(s.get("c0", 1.0)),
                    "c1": float(s.get("c1", 0.0)),
                    "c2": float(s.get("c2", 0.0)),
                    "c3": float(s.get("c3", 0.0)),
                    "c4": float(s.get("c4", 0.0)),
                }
            )
        by_label[str(label)] = stage_rows

    missing = [c for c in _CORNERS if c not in by_label]
    if missing:
        raise ValueError(f"compiled-v1 body missing corners: {missing}")

    return by_label


def interpolate_compiled_v1_stages(
    body_dict: dict[str, Any],
    morph: float,
    q: float,
) -> list[EncodedCoeffs]:
    """Exact bilinear interpolation (Q then morph) over compiled-v1 c0..c4."""
    morph = _clamp01(morph)
    q = _clamp01(q)

    corners = _extract_keyframes(body_dict)

    stage_counts = [len(corners[c]) for c in _CORNERS]
    num_stages = min(stage_counts) if stage_counts else 0
    if num_stages <= 0:
        return []

    out: list[EncodedCoeffs] = []
    for i in range(num_stages):
        def interp(key: str) -> float:
            m0_q0 = corners["M0_Q0"][i][key]
            m0_q100 = corners["M0_Q100"][i][key]
            m100_q0 = corners["M100_Q0"][i][key]
            m100_q100 = corners["M100_Q100"][i][key]

            q_m0 = m0_q0 + (m0_q100 - m0_q0) * q
            q_m1 = m100_q0 + (m100_q100 - m100_q0) * q
            return q_m0 + (q_m1 - q_m0) * morph

        out.append(
            EncodedCoeffs(
                c0=interp("c0"),
                c1=interp("c1"),
                c2=interp("c2"),
                c3=interp("c3"),
                c4=interp("c4"),
            )
        )

    return out

