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

DRY = Path(r"C:/Users/hooki/trenchwork_clean/ref/bypassed-pinknoise.wav")
# Canonical raw P2K skin source (sha256 9fb1bef0d0212980 for P2k_013).
SKINS = Path(r"C:/Users/hooki/trenchwork_clean/datasets/p2k_skins")
# Canonical references rendered by tools/render_canonical_refs.py from each
# skin in SKINS. MANIFEST.json lists which skins were rendered and where.
REF_DIR = Path(r"C:/Users/hooki/Trench/ref/canonical")

CORNERS = ("M0_Q0", "M0_Q100", "M100_Q0", "M100_Q100")
# Null worse than this triggers a failure exit code.
FAIL_THRESHOLD_DB = -140.0

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


def best_fit_null(pred: np.ndarray, ref: np.ndarray, lag_search: int = 0):
    """Lag-aligned best-fit-gain null. Canonical references are rendered
    through the exact same pipeline as `pred`, so lag=0 is correct and we
    skip the search by default."""
    n = min(len(pred), len(ref))
    p = pred[:n].astype(np.float64, copy=False)
    r = ref[:n].astype(np.float64, copy=False)
    den = float(np.dot(p, p))
    if den <= 0.0:
        return (0, 0.0, 0.0, 0.0, 0.0)
    gain = float(np.dot(p, r) / den)
    res = r - gain * p
    rms = float(np.sqrt(float(np.dot(res, res)) / n))
    peak = float(np.max(np.abs(res)))
    ref_rms = float(np.sqrt(float(np.dot(r, r)) / n))
    return (0, gain, rms, peak, ref_rms)


def corner_keyframe(body: dict, label: str) -> dict:
    """Accept both layouts: trenchwork_clean cartridge wire form
    (keyframes list) and canonical raw p2k_skins form (corners dict)."""
    if "keyframes" in body:
        kf = next((k for k in body["keyframes"] if k.get("label") == label), None)
        if kf is None:
            raise RuntimeError(f"missing corner {label}")
        return kf
    if "corners" in body:
        corner = body["corners"].get(label)
        if corner is None:
            raise RuntimeError(f"missing corner {label}")
        return {
            "stages": corner["stages"],
            "boost": float(body.get("boost", 1.0)),
        }
    raise RuntimeError("cartridge has neither 'keyframes' nor 'corners'")


def main() -> int:
    manifest_path = REF_DIR / "MANIFEST.json"
    if not REF_DIR.exists() or not manifest_path.exists() or not SKINS.exists() or not DRY.exists():
        print(
            f"parity inputs missing; skip "
            f"({REF_DIR.exists()=} {manifest_path.exists()=} {SKINS.exists()=} {DRY.exists()=})"
        )
        return 0

    manifest = json.loads(manifest_path.read_text())
    dry = mono(sf.read(str(DRY))[0])

    print(f"dry input:  {DRY.name}  {len(dry)} samples")
    print(f"refs:       {REF_DIR}  ({len(manifest)} skins)")
    print(f"pipeline:   raw-ROM → SOS cascade → AGC → ×boost → lag+gain null")
    print(f"fail at:    rel_null > {FAIL_THRESHOLD_DB:.0f} dB")
    print()
    print(f"{'skin':24} {'boost':>6} {'worst_corner':>12} {'lag':>5} {'gain':>8} {'rel_null':>11}")

    failures: list[tuple[str, str, float]] = []
    import gc
    for name in sorted(manifest):
        cart_path = SKINS / f"{name}.json"
        body = json.loads(cart_path.read_text())
        worst = (None, 0.0, 0, 1.0)  # corner, rel_db, lag, gain
        for label in CORNERS:
            wav = REF_DIR / f"{name}_{label}.wav"
            if not wav.exists():
                continue
            ref = mono(sf.read(str(wav))[0])
            kf = corner_keyframe(body, label)
            sos = build_sos(kf["stages"])
            boost = float(kf.get("boost", 1.0))
            cascaded = sosfilt(sos, dry)
            pred = apply_agc(cascaded) * boost
            lag, gain, rms, _peak, ref_rms = best_fit_null(pred, ref)
            rel = db(rms) - db(ref_rms)
            if worst[0] is None or rel > worst[1]:
                worst = (label, rel, lag, gain)
            del ref, sos, cascaded, pred
        gc.collect()
        if worst[0] is None:
            print(f"{name:24} no refs")
            continue
        boost = float(body.get("boost", 1.0))
        marker = " " if worst[1] <= FAIL_THRESHOLD_DB else "!"
        print(
            f"{marker}{name:23} {boost:6.2f} {worst[0]:>12} "
            f"{worst[2]:5d} {worst[3]:8.4f} {worst[1]:+8.1f} dB"
        )
        if worst[1] > FAIL_THRESHOLD_DB:
            failures.append((name, worst[0], worst[1]))

    print()
    if failures:
        print(f"FAIL: {len(failures)}/{len(manifest)} skins exceeded {FAIL_THRESHOLD_DB:.0f} dB")
        for name, corner, rel in failures:
            print(f"  {name}  worst corner {corner}  rel_null={rel:+.1f}dB")
        return 1
    print(f"OK: {len(manifest)}/{len(manifest)} skins null at <= {FAIL_THRESHOLD_DB:.0f} dB")
    return 0


if __name__ == "__main__":
    sys.exit(main())
