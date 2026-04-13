"""Parity null: TRENCH DSP pipeline vs Talking Hedz hardware reference.

Pipeline (matches trenchwork_clean trench-core::engine::FilterEngine):
    raw ROM stage → (a1, r, val1, val2, val3, flag) → (b0,b1,b2,a1,a2)
      flag=0  → b0=1, b1=b2=0       (final stage lowpass reinterpretation)
      flag=1  → b0=1+val1, b1=a1+val2, b2=r²-val3
    → scipy.signal.sosfilt cascade (6 serial stages, no clamp)
    → per-sample AGC (16-entry table, verified against EmulatorX.dll)
    → best-fit lag + gain null vs hardware capture

The AGC is the missing runtime gain-staging piece. Its table and
algorithm are ported from trenchwork_clean/trench-core/src/agc.rs
(documented as verified against the EmulatorX.dll binary).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import soundfile as sf
from scipy.signal import sosfilt

REF_DIR = Path(r"C:/Users/hooki/trenchwork_clean/ref")
CART = Path(r"C:/Users/hooki/trenchwork_clean/cartridges/00_talking_hedz.json")

CORNERS = [
    ("M0_Q0",     "hedzm0q0.wav"),
    ("M0_Q100",   "hedzmorph0q100.wav"),
    ("M100_Q0",   "hedzmorph100q0.wav"),
    ("M100_Q100", "hedzmorph100q100.wav"),
]

# AGC_TABLE ported from trenchwork_clean/trench-core/src/agc.rs
AGC_TABLE = np.array(
    [1.0001, 1.0001, 0.996, 0.990, 0.920, 0.500, 0.200, 0.160,
     0.120, 0.120, 0.120, 0.120, 0.120, 0.120, 0.120, 0.120],
    dtype=np.float64,
)


def db(x: float) -> float:
    x = float(x)
    return -300.0 if x <= 0.0 else 20.0 * np.log10(x)


def mono(d: np.ndarray) -> np.ndarray:
    return d[:, 0] if d.ndim > 1 else d


def stage_coefficients(stage: dict) -> tuple[float, float, float, float, float]:
    """Raw ROM stage → (b0,b1,b2,a1,a2). Matches the RE stage-response model
    and the FilterEngine coefficient path."""
    a1 = float(stage["a1"])
    r = min(float(stage["r"]), 0.999999)  # stability clamp
    a2 = r * r
    if float(stage.get("flag", 1.0)) < 0.5:
        b0, b1, b2 = 1.0, 0.0, 0.0
    else:
        b0 = 1.0 + float(stage["val1"])
        b1 = a1 + float(stage["val2"])
        b2 = a2 - float(stage["val3"])
    return b0, b1, b2, a1, a2


def build_sos(stages: list[dict]) -> np.ndarray:
    rows = []
    for st in stages:
        b0, b1, b2, a1, a2 = stage_coefficients(st)
        rows.append([b0, b1, b2, 1.0, a1, a2])
    return np.asarray(rows, dtype=np.float64)


def apply_agc(samples: np.ndarray) -> np.ndarray:
    """Per-sample post-cascade soft limiter (EmulatorX.dll AGC)."""
    out = np.empty_like(samples)
    gain = 1.0
    for i, s in enumerate(samples):
        idx = int(gain * abs(s)) & 0xF
        new_gain = gain * AGC_TABLE[idx]
        gain = new_gain if new_gain < 1.0 else 1.0
        out[i] = s * gain
    return out


def apply_dc_block(samples: np.ndarray, sample_rate: float = 44100.0) -> np.ndarray:
    """One-pole DC blocker matching trench-core engine.rs DcBlocker."""
    r = 1.0 - (2.0 * np.pi * 20.0 / sample_rate)
    out = np.empty_like(samples)
    x_prev = 0.0
    y_prev = 0.0
    for i, x in enumerate(samples):
        y = x - x_prev + r * y_prev
        x_prev = x
        y_prev = y if np.isfinite(y) else 0.0
        out[i] = y_prev
    return out


def best_fit_null(pred: np.ndarray, ref: np.ndarray, lag_search: int = 512):
    n = min(len(pred), len(ref))
    pred = pred[:n].astype(np.float64)
    ref = ref[:n].astype(np.float64)
    best = (0, 0.0, float("inf"), 0.0, 0.0)  # lag, gain, rms, peak, ref_rms
    ref_rms_full = float(np.sqrt(np.mean(ref * ref)))
    for lag in range(-lag_search, lag_search + 1):
        if lag >= 0:
            p = pred[: n - lag]
            r = ref[lag:]
        else:
            p = pred[-lag:]
            r = ref[: n + lag]
        if len(p) < 1024:
            continue
        den = float((p * p).sum())
        if den <= 0:
            continue
        gain = float((p * r).sum() / den)
        res = r - gain * p
        rms = float(np.sqrt(np.mean(res * res)))
        if rms < best[2]:
            peak = float(np.max(np.abs(res)))
            ref_rms = float(np.sqrt(np.mean(r * r)))
            best = (lag, gain, rms, peak, ref_rms)
    return best


def corner_keyframe(body: dict, label: str) -> dict:
    kf = next((k for k in body["keyframes"] if k.get("label") == label), None)
    if kf is None:
        raise RuntimeError(f"missing corner {label}")
    return kf


def main() -> int:
    if not REF_DIR.exists() or not CART.exists():
        print(f"reference data missing; skip ({REF_DIR.exists()=} {CART.exists()=})")
        return 0

    body = json.loads(CART.read_text())
    dry = mono(sf.read(str(REF_DIR / "bypassed-pinknoise.wav"))[0])

    print(f"cartridge:  {CART.name}  ({body.get('name')})")
    print(f"dry input:  bypassed-pinknoise.wav  {len(dry)} samples")
    print(f"pipeline:   raw-ROM → SOS cascade → AGC → ×boost → lag+gain null")
    print()
    print(f"{'corner':10} {'boost':>6} {'lag':>5} {'gain':>8} "
          f"{'peak_res':>10} {'rms_res':>10} {'ref_rms':>10} {'rel_null':>10}")

    for label, wav in CORNERS:
        ref_path = REF_DIR / wav
        if not ref_path.exists():
            print(f"{label:10} missing {wav}")
            continue
        ref = mono(sf.read(str(ref_path))[0])
        kf = corner_keyframe(body, label)
        sos = build_sos(kf["stages"])
        boost = float(kf.get("boost", 1.0))
        cascaded = sosfilt(sos, dry.copy())
        pred = apply_agc(cascaded) * boost
        lag, gain, rms, peak, ref_rms = best_fit_null(pred, ref)
        rel = db(rms) - db(ref_rms)
        print(
            f"{label:10} {boost:6.2f} {lag:5d} {gain:8.4f} "
            f"{db(peak):+9.1f}dB {db(rms):+9.1f}dB "
            f"{db(ref_rms):+9.1f}dB {rel:+9.1f}dB"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
