"""Weighted centroid shift bias for the pill sampler.

v1 contract (locked 2026-04-19):
  keep   =  1.00  strong pull toward
  maybe  =  0.35  weak pull toward
  reject = -0.20  mild push away

  shift formula: shift = Σ(w_i × (v_i − center)) / Σ(|w_i|)
  — not Σ(w_i × v_i) / Σ(|w_i|), which would treat negative weights as
    phantom anti-points rather than deviations from the current center.

  max shift per axis = 12.5% of current envelope width (per round)
  exploration fraction = 0.25  (25% of each batch drawn from unshifted envelope)
"""
from __future__ import annotations

import csv
import json
from copy import deepcopy
from pathlib import Path

VERDICT_WEIGHTS = {"keep": 1.00, "maybe": 0.35, "reject": -0.20}
MAX_SHIFT_FRAC = 0.125  # 12.5% of envelope width per axis
EXPLORE_FRAC = 0.25     # 25% of each batch uses the unshifted envelope

# Pill keyframe label → scope envelope corner key
KEYFRAME_TO_CORNER = {
    "M0_Q0":      "c000",
    "M100_Q0":    "c100",
    "M0_Q100":    "c010",
    "M100_Q100":  "c110",
}


def extract_bias_params(authored: dict) -> dict:
    """Per-keyframe (freq_mean, gain_mean) from an authored pill's sections.

    Stored in the compiled pill's provenance so future bias rounds can read
    interpretable values instead of biquad coefficients.
    """
    params = {}
    for label, payload in authored.get("keyframes", {}).items():
        sections = payload.get("template", {}).get("sections", [])
        if not sections:
            continue
        freq_mean = sum(s["low_freq"] for s in sections) / len(sections)
        gain_mean = sum(s["low_gain"] for s in sections) / len(sections)
        params[label] = {
            "freq_mean": round(freq_mean, 2),
            "gain_mean": round(gain_mean, 2),
        }
    return params


def pill_fingerprint(bias_params: dict) -> str:
    """Stable duplicate-suppression key from per-corner rounded freq/gain."""
    parts = [
        f"{label}:{round(bp['freq_mean'])}:{round(bp['gain_mean'])}"
        for label, bp in sorted(bias_params.items())
    ]
    return "|".join(parts)


def load_verdicts(csv_path: Path, scope: str) -> list[tuple[float, dict]]:
    """Load (weight, bias_params) pairs for *scope* from the verdicts CSV.

    Pills without provenance.bias_params are silently skipped — they predate
    the bias pipeline and carry no usable envelope information.
    """
    if not csv_path.exists():
        return []
    results = []
    with open(csv_path, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            if row.get("scope") != scope:
                continue
            weight = VERDICT_WEIGHTS.get(row.get("verdict"))
            if weight is None:
                continue
            fpath = Path(row.get("file", ""))
            if not fpath.exists():
                continue
            try:
                pill = json.loads(fpath.read_text(encoding="utf-8"))
            except Exception:
                continue
            bp = pill.get("provenance", {}).get("bias_params")
            if not bp:
                continue
            results.append((weight, bp))
    return results


def shift_envelope(scope_env: dict, verdicts: list[tuple[float, dict]]) -> tuple[dict, dict]:
    """Return (shifted_scope_env, shift_report).

    Shift formula (per corner, per axis):
        shift = Σ(w_i × (v_i − center)) / Σ(|w_i|)
        new_center = clamp(center + clamp(shift, ±max_shift), scope_bounds)

    shift_report keys: n_verdicts, corners → {freq/gain center_before/after/shift_pct}
    """
    if not verdicts:
        return deepcopy(scope_env), {"n_verdicts": 0, "corners": {}}

    shifted = deepcopy(scope_env)
    freq_width = scope_env["freq_max"] - scope_env["freq_min"]
    gain_width = scope_env["gain_max"] - scope_env["gain_min"]
    max_freq_shift = MAX_SHIFT_FRAC * freq_width
    max_gain_shift = MAX_SHIFT_FRAC * gain_width
    report_corners: dict[str, dict] = {}

    for kf_label, corner_key in KEYFRAME_TO_CORNER.items():
        if corner_key not in scope_env.get("corners", {}):
            continue
        corner = scope_env["corners"][corner_key]
        freq_c = corner["freq_center"]
        gain_c = corner["gain_center"]

        fdev = gdev = abs_w = 0.0
        for weight, bp in verdicts:
            kf = bp.get(kf_label)
            if kf is None:
                continue
            fdev += weight * (kf["freq_mean"] - freq_c)
            gdev += weight * (kf["gain_mean"] - gain_c)
            abs_w += abs(weight)

        if abs_w == 0.0:
            continue

        freq_shift = max(-max_freq_shift, min(max_freq_shift, fdev / abs_w))
        gain_shift = max(-max_gain_shift, min(max_gain_shift, gdev / abs_w))

        new_freq = max(scope_env["freq_min"],
                       min(scope_env["freq_max"], round(freq_c + freq_shift)))
        new_gain = max(scope_env["gain_min"],
                       min(scope_env["gain_max"], round(gain_c + gain_shift)))

        shifted["corners"][corner_key] = deepcopy(corner)
        shifted["corners"][corner_key]["freq_center"] = new_freq
        shifted["corners"][corner_key]["gain_center"] = new_gain

        report_corners[corner_key] = {
            "freq_center_before": freq_c,
            "freq_center_after": new_freq,
            "freq_shift_pct": round(abs(freq_shift) / freq_width * 100, 1) if freq_width else 0.0,
            "gain_center_before": gain_c,
            "gain_center_after": new_gain,
            "gain_shift_pct": round(abs(gain_shift) / gain_width * 100, 1) if gain_width else 0.0,
        }

    return shifted, {"n_verdicts": len(verdicts), "corners": report_corners}
