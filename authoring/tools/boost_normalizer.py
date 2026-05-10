"""forge/boost_normalizer.py — analytical per-corner boost normalization.

Reads a RawStage cartridge JSON, evaluates the 6-stage cascade peak |H(jw)|
over [20 Hz, 20000 Hz] at 39062.5 Hz authoring SR, and returns the boost
value that keeps cascade peak amplitude <= peak_amplitude_target.

Target headroom justification (from agc.rs):
  AGC_TABLE = [1.0001, 1.0001, 0.996, 0.990, 0.920, 0.500, ...]
  idx = (agc_gain * abs_sample) as u32 & 0xF
  The cliff at idx 4->5 drops gain from 0.920 -> 0.500 — halved in one step.
  idx 4 fires when agc_gain * abs_sample is in [4.0, 5.0).
  To stay clear of idx 5 (abs_sample >= 5.0 with agc_gain=1.0), we target
  abs_sample < 5.0 at steady state.  High-Q resonators ring up ~2x peak on
  attack transients, so target = 5.0 / 2.0 = 2.5 gives margin for both.
  20*log10(5.0) = +14.0 dB  (cliff)
  20*log10(2.5) = +7.96 dB  (target — call it +8 dB)
"""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np

AUTHORING_SR: float = 39062.5
EVAL_FREQS = np.geomspace(20.0, 20000.0, 2048)

# AGC cliff at idx=5 fires when agc_gain * abs_sample >= 5.0.
# With agc_gain=1.0 worst case and 2x transient overshoot budget:
PEAK_AMPLITUDE_TARGET: float = 2.5  # linear


def _compile_rawstage(rs: dict[str, Any]) -> dict[str, float] | None:
    """Compile one RawStage dict to DF2T coefficients.

    Uses the exact formula from emu_resonator.rs:
        c0 = 1.0 + val1
        c1 = a1  + val2
        c2 = r^2 - val3
        c3 = a1
        c4 = r^2

    Returns None for passthrough (r == 0 and a1 == 0).
    """
    a1 = rs["a1"]
    r = rs["r"]
    if r == 0.0 and a1 == 0.0:
        return None
    a2 = r * r
    return {
        "c0": 1.0 + rs["val1"],
        "c1": a1 + rs["val2"],
        "c2": a2 - rs["val3"],
        "c3": a1,
        "c4": a2,
    }


def _cascade_peak(coefs: list[dict[str, float]], freqs: np.ndarray) -> float:
    """Evaluate |H(jw)| cascade product over freqs, return peak linear amplitude."""
    w = 2.0 * np.pi * freqs / AUTHORING_SR
    e1 = np.exp(-1j * w)
    e2 = np.exp(-2j * w)
    mag = np.ones(len(freqs), dtype=np.float64)
    for c in coefs:
        num = c["c0"] + c["c1"] * e1 + c["c2"] * e2
        den = 1.0 + c["c3"] * e1 + c["c4"] * e2
        mag *= np.abs(num / den)
    return float(np.max(mag))


def compute_boosts_from_cartridge(
    json_path: Path,
    target: float = PEAK_AMPLITUDE_TARGET,
) -> dict[str, dict[str, float]]:
    """Parse cartridge JSON, evaluate each keyframe, return per-corner results.

    Returns a dict keyed by corner label with:
        peak   — raw cascade peak amplitude (no existing boost applied)
        boost  — proposed boost = target / peak
        check  — peak * boost (should be <= target * 1.001)
    """
    data = json.loads(json_path.read_text(encoding="utf-8"))
    results: dict[str, dict[str, float]] = {}

    for kf in data["keyframes"]:
        label: str = kf["label"]
        raw_stages = kf["stages"][:6]  # only first 6 active stages
        coefs = [c for rs in raw_stages if (c := _compile_rawstage(rs)) is not None]
        peak = _cascade_peak(coefs, EVAL_FREQS)
        boost = target / peak if peak > 0.0 else 1.0
        check = peak * boost
        assert check <= target * 1.001, (
            f"NO-FAKE-SUCCESS: {json_path.name} {label}: "
            f"peak={peak:.4f} boost={boost:.6f} check={check:.4f} > {target * 1.001:.4f}"
        )
        results[label] = {"peak": peak, "boost": boost, "check": check}

    return results


def compute_boosts_from_stages(
    coefs_per_corner: dict[str, list[dict[str, float]]],
    target: float = PEAK_AMPLITUDE_TARGET,
) -> dict[str, dict[str, float]]:
    """Evaluate pre-compiled DF2T coef lists per corner, return same result shape."""
    results: dict[str, dict[str, float]] = {}
    for label, coefs in coefs_per_corner.items():
        peak = _cascade_peak(coefs, EVAL_FREQS)
        boost = target / peak if peak > 0.0 else 1.0
        check = peak * boost
        assert check <= target * 1.001, (
            f"NO-FAKE-SUCCESS: {label}: "
            f"peak={peak:.4f} boost={boost:.6f} check={check:.4f} > {target * 1.001:.4f}"
        )
        results[label] = {"peak": peak, "boost": boost, "check": check}
    return results


def plot_boost_validation(
    body_name: str,
    corner_coefs: dict[str, list[dict[str, float]]],
    boosts: dict[str, dict[str, float]],
    out_path: Path,
) -> None:
    """Plot |H(jw)| in dB for all 4 corners, dashed=current (boost=1) solid=proposed.

    Cliff line at +14 dB (=20*log10(5.0)), target line at +7.96 dB (=20*log10(2.5)).
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    freqs = EVAL_FREQS
    cliff_db = 20.0 * math.log10(5.0)     # +14.0 dB
    target_db = 20.0 * math.log10(2.5)    # +7.96 dB

    corners = list(corner_coefs.keys())
    fig, axes = plt.subplots(2, 2, figsize=(16, 9), dpi=140)
    fig.patch.set_facecolor("#16191a")
    fig.suptitle(
        f"BOOST NORMALIZATION :: {body_name}\n"
        f"dashed=current (boost=1.0)  solid=proposed  "
        f"red=cliff(+14dB)  green=target(+8dB)",
        color="#ffba00", fontweight="bold", fontsize=11,
    )

    for ax, label in zip(axes.flat, corners):
        coefs = corner_coefs[label]
        info = boosts[label]
        boost = info["boost"]

        w = 2.0 * np.pi * freqs / AUTHORING_SR
        e1 = np.exp(-1j * w)
        e2 = np.exp(-2j * w)
        mag = np.ones(len(freqs), dtype=np.float64)
        for c in coefs:
            num = c["c0"] + c["c1"] * e1 + c["c2"] * e2
            den = 1.0 + c["c3"] * e1 + c["c4"] * e2
            mag *= np.abs(num / den)

        db_raw = 20.0 * np.log10(np.maximum(mag, 1e-12))
        db_boosted = db_raw + 20.0 * math.log10(max(boost, 1e-12))

        ax.set_facecolor("#0d0f10")
        ax.semilogx(freqs, db_raw, "--", color="#888888", lw=1.5,
                    label=f"boost=1.0 (raw)")
        ax.semilogx(freqs, db_boosted, "-", color="#22ddff", lw=2.2,
                    label=f"boost={boost:.4f}")
        ax.axhline(cliff_db, color="#ff4444", lw=1.4, linestyle="--",
                   label=f"cliff +{cliff_db:.1f} dB")
        ax.axhline(target_db, color="#44ee44", lw=1.4, linestyle="--",
                   label=f"target +{target_db:.1f} dB")
        ax.set_xscale("log")
        ax.set_xlim(20, 20000)
        ax.set_xlabel("Hz", color="white", fontsize=8)
        ax.set_ylabel("dB", color="white", fontsize=8)
        ax.set_title(
            f"{label}  peak={info['peak']:.3f}  boost={boost:.4f}  "
            f"peak*boost={info['check']:.3f}",
            color="#ffba00", fontsize=9, fontweight="bold",
        )
        ax.grid(True, which="both", alpha=0.18, color="white")
        ax.tick_params(colors="white", labelsize=7)
        for spine in ax.spines.values():
            spine.set_color("#3a3e42")
        ax.legend(facecolor="#22262a", edgecolor="#444", labelcolor="white",
                  fontsize=7, loc="upper right")

    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"  [boost_normalizer] saved plot -> {out_path}")
