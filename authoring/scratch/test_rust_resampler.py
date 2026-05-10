"""Self-test for tools/rust_resampler.py — mirrors the Rust unit tests
in resampler.rs and adds a roundtrip null threshold.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from rust_resampler import Resampler, resample, NATIVE_SR


def make_sine(sr: float, freq: float, n: int) -> np.ndarray:
    return np.sin(2.0 * math.pi * freq * np.arange(n, dtype=np.float64) / sr).astype(np.float32)


def rms(x: np.ndarray) -> float:
    x = np.asarray(x, dtype=np.float64)
    return float(np.sqrt(np.mean(x * x)))


def db(ratio: float) -> float:
    return -300.0 if ratio <= 0.0 else 20.0 * math.log10(ratio)


def test_silence() -> None:
    zeros = np.zeros(1024, dtype=np.float32)
    out = resample(zeros, 44100.0, NATIVE_SR)
    assert len(out) > 0
    assert float(np.max(np.abs(out))) < 1e-6, f"silence produced {np.max(np.abs(out))}"
    print(f"  silence_is_silence       OK  ({len(out)} samples, max={float(np.max(np.abs(out))):.2e})")


def test_roundtrip_44k() -> None:
    host_sr = 44100.0
    freq = 1000.0
    n = 4096
    inp = make_sine(host_sr, freq, n)
    rms_in = rms(inp[256:])
    native = resample(inp, host_sr, NATIVE_SR)
    out = resample(native, NATIVE_SR, host_sr)
    rms_out = rms(out[256:])
    ratio = rms_out / rms_in
    print(f"  roundtrip_44k_1kHz       rms_in={rms_in:.4f} rms_out={rms_out:.4f} ratio={ratio:.4f} ({db(ratio):+.3f} dB)")
    assert 0.891 <= ratio <= 1.122, f"out of ±1 dB: ratio={ratio}"


def test_roundtrip_48k() -> None:
    host_sr = 48000.0
    freq = 440.0
    n = 4096
    inp = make_sine(host_sr, freq, n)
    rms_in = rms(inp[256:])
    native = resample(inp, host_sr, NATIVE_SR)
    out = resample(native, NATIVE_SR, host_sr)
    rms_out = rms(out[256:])
    ratio = rms_out / rms_in
    print(f"  roundtrip_48k_440Hz      rms_in={rms_in:.4f} rms_out={rms_out:.4f} ratio={ratio:.4f} ({db(ratio):+.3f} dB)")
    assert 0.891 <= ratio <= 1.122, f"out of ±1 dB: ratio={ratio}"


def test_antialias() -> None:
    host_sr = 44100.0
    alias_freq = 22000.0
    n = 4096
    inp = make_sine(host_sr, alias_freq, n)
    rms_in = rms(inp[256:])
    native = resample(inp, host_sr, NATIVE_SR)
    rms_native = rms(native[256:])
    print(f"  antialias_22kHz          rms_in={rms_in:.4f} rms_native={rms_native:.4f} attenuation={db(rms_native/rms_in):+.1f} dB")
    assert rms_native < rms_in * 0.1


def test_reset() -> None:
    r = Resampler(44100.0, NATIVE_SR)
    loud = np.ones(128, dtype=np.float32)
    out1 = np.zeros(128, dtype=np.float32)
    r.process(loud, out1)
    r.reset()
    zeros = np.zeros(128, dtype=np.float32)
    out2 = np.zeros(128, dtype=np.float32)
    r.process(zeros, out2)
    mx = float(np.max(np.abs(out2)))
    assert mx < 1e-5, f"after reset: expected silence, got max={mx}"
    print(f"  reset_clears_state       OK  (max after reset={mx:.2e})")


def test_step_semantics() -> None:
    """output_count: 32 input @ 44100 -> NATIVE_SR should yield 27-30."""
    r = Resampler(44100.0, NATIVE_SR)
    est = r.output_count(32)
    assert 27 <= est <= 30, f"output_count out of range: {est}"
    print(f"  output_count             OK  (32 in -> {est} out)")


def main() -> int:
    print("rust_resampler self-test:")
    test_silence()
    test_roundtrip_44k()
    test_roundtrip_48k()
    test_antialias()
    test_reset()
    test_step_semantics()
    print("ALL OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
