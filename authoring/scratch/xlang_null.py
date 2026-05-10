"""Null-test the Python port against Rust resampler dump."""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from rust_resampler import resample, NATIVE_SR

TARGET = Path(__file__).resolve().parent.parent.parent / "target"


def load_f32(path: Path) -> np.ndarray:
    return np.frombuffer(path.read_bytes(), dtype="<f4").astype(np.float32, copy=True)


def make_sine(sr: float, freq: float, n: int) -> np.ndarray:
    return np.sin(2.0 * math.pi * freq * np.arange(n, dtype=np.float64) / sr).astype(np.float32)


def null(pred: np.ndarray, ref: np.ndarray) -> tuple[float, float, int]:
    n = min(len(pred), len(ref))
    if n == 0:
        return 0.0, 0.0, 0
    p = pred[:n].astype(np.float64)
    r = ref[:n].astype(np.float64)
    den = float(np.dot(p, p))
    gain = float(np.dot(p, r) / den) if den > 0 else 0.0
    res = r - gain * p
    rel = float(np.sqrt(np.dot(res, res) / n))
    r_rms = float(np.sqrt(np.dot(r, r) / n))
    db = -300.0 if r_rms <= 0 or rel <= 0 else 20.0 * math.log10(rel / r_rms)
    return gain, db, n


def main() -> int:
    n = 4096

    # Host -> native
    inp = make_sine(44100.0, 1000.0, n)
    py_native = resample(inp, 44100.0, NATIVE_SR)
    rs_native = load_f32(TARGET / "resampler_dump_44k_to_native.f32")
    print(f"host->native: py={len(py_native)} rs={len(rs_native)}")
    gain, rel_db, n_used = null(py_native, rs_native)
    print(f"  null: gain={gain:.6f} rel_null={rel_db:+.1f} dB  over {n_used} samples")

    # Native -> host
    inp2 = make_sine(NATIVE_SR, 1000.0, n)
    py_host = resample(inp2, NATIVE_SR, 44100.0)
    rs_host = load_f32(TARGET / "resampler_dump_native_to_44k.f32")
    print(f"native->host: py={len(py_host)} rs={len(rs_host)}")
    gain, rel_db, n_used = null(py_host, rs_host)
    print(f"  null: gain={gain:.6f} rel_null={rel_db:+.1f} dB  over {n_used} samples")
    return 0


if __name__ == "__main__":
    sys.exit(main())
