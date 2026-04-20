"""Port of trench-core/src/resampler.rs to Python.

Bit-close translation of FilterEngine's polyphase FIR resampler so that
the canonical render chain in tools/parity_null.py uses the SAME FIR as
the runtime Rust engine. scipy.signal.resample_poly uses a Kaiser
window with ~23 taps per subfilter; the Rust engine uses 32 taps ×
512 phases Blackman-windowed sinc with per-phase DC normalisation.
Different magnitude responses near Nyquist produce a residual that
scales with filter resonance.

Invariants (must match resampler.rs):
    TAPS       = 32
    PHASES     = 512
    RING_LEN   = 64        (power of 2, > TAPS)
    wpos init  = TAPS / 2  (pre-filled zeros for causal window)
    rpos init  = 0.0
    step       = input_sr / output_sr
    cutoff     = min(output_sr / input_sr, 1.0)

All runtime math is f32. Table construction runs in f64, result cast
to f32 (matches `(proto[idx] * phase_scale) as f32` in Rust).

Public surface:
    class Resampler — stateful, mirrors resampler.rs `struct Resampler`
        .process(input, output) -> (consumed, produced)
        .reset()
        .output_count(n_input) -> int
    def resample(x, input_sr, output_sr) -> np.ndarray (f32)
"""
from __future__ import annotations

import math

import numpy as np

TAPS = 32
PHASES = 512
RING_LEN = 64
RING_MASK = RING_LEN - 1
NATIVE_SR = 39_062.5


def _build_sinc_table(cutoff: float) -> np.ndarray:
    """Build `(PHASES + 1) * TAPS` Blackman-windowed sinc table.

    Port of trench-core/src/resampler.rs::build_sinc_table. Prototype
    is evaluated in f64 (same as Rust) and the per-phase normalised
    polyphase rows are cast to f32 at the end (`as f32` in Rust).
    """
    n = TAPS * PHASES
    half_n = n / 2.0
    proto = np.zeros(n + 1, dtype=np.float64)
    for i in range(n + 1):
        t = (float(i) - half_n) / float(PHASES)
        if abs(t) < 1e-10:
            sinc = cutoff
        else:
            sinc = math.sin(math.pi * cutoff * t) / (math.pi * t)
        w = (
            0.42
            - 0.5 * math.cos(2.0 * math.pi * float(i) / float(n))
            + 0.08 * math.cos(4.0 * math.pi * float(i) / float(n))
        )
        proto[i] = sinc * w

    size = (PHASES + 1) * TAPS
    table = np.zeros(size, dtype=np.float32)
    for p in range(PHASES + 1):
        # Per-phase DC normalisation (f64 sum, same as Rust).
        phase_sum = 0.0
        for k in range(TAPS):
            idx = k * PHASES + p
            if idx <= n:
                phase_sum += proto[idx]
        if abs(phase_sum) > 1e-12:
            phase_scale = 1.0 / phase_sum
        else:
            phase_scale = 0.0
        for k in range(TAPS):
            idx = k * PHASES + p
            if idx <= n:
                # Cast f64 product to f32 — matches `as f32` in Rust.
                table[p * TAPS + k] = np.float32(proto[idx] * phase_scale)
    return table


class Resampler:
    """Polyphase FIR resampler. Straight port of resampler.rs."""

    def __init__(self, input_sr: float, output_sr: float) -> None:
        self.step = float(input_sr) / float(output_sr)
        cutoff = min(float(output_sr) / float(input_sr), 1.0)
        self.table = _build_sinc_table(cutoff)
        self.ring = np.zeros(RING_LEN, dtype=np.float32)
        # Pre-fill TAPS/2 zeros so rpos=0.0 has a full causal window.
        self.wpos = TAPS // 2
        self.rpos = 0.0

    def reset(self) -> None:
        self.ring.fill(0.0)
        self.wpos = TAPS // 2
        self.rpos = 0.0

    def output_count(self, n_input: int) -> int:
        return int(float(n_input) / self.step)

    def process(self, input_buf: np.ndarray, output_buf: np.ndarray) -> tuple[int, int]:
        """Feed `input_buf`, drain into `output_buf`. (consumed, produced)."""
        half = TAPS // 2
        in_idx = 0
        out_idx = 0
        in_len = len(input_buf)
        out_len = len(output_buf)
        ring = self.ring
        step = self.step

        while True:
            center = int(self.rpos)  # floor for non-negative; rpos never negative
            need = center + half

            while self.wpos <= need:
                if in_idx >= in_len:
                    return in_idx, out_idx
                ring[self.wpos & RING_MASK] = np.float32(input_buf[in_idx])
                self.wpos += 1
                in_idx += 1

            if out_idx >= out_len:
                return in_idx, out_idx

            output_buf[out_idx] = self._interpolate()
            out_idx += 1
            self.rpos += step

    def _interpolate(self) -> np.float32:
        # Match Rust's interpolate() exactly:
        #   let center = self.rpos.floor() as i64;
        #   let frac = (self.rpos - center as f64) as f32;
        #   let phase_f = frac * PHASES as f32;
        #   let phase_lo = phase_f as usize;
        #   let phase_hi = (phase_lo + 1).min(PHASES);
        #   let t = phase_f - phase_lo as f32;
        center = int(math.floor(self.rpos))
        frac = np.float32(self.rpos - float(center))
        phase_f = np.float32(frac * np.float32(PHASES))
        phase_lo = int(phase_f)  # truncation, matches `as usize`
        if phase_lo > PHASES:
            phase_lo = PHASES
        phase_hi = phase_lo + 1
        if phase_hi > PHASES:
            phase_hi = PHASES
        t = np.float32(phase_f - np.float32(phase_lo))

        base_lo = phase_lo * TAPS
        base_hi = phase_hi * TAPS
        half = TAPS // 2

        # Vectorised tap loop. All arithmetic stays in f32 — we construct
        # the coefficients in f32 and dot against the f32 ring slice.
        # Ring indices: (center - half + 1 + k) & MASK for k in 0..TAPS.
        idx = (np.arange(TAPS, dtype=np.int64) + (center - half + 1)) & RING_MASK
        taps = self.ring[idx]  # f32 slice
        c_lo = self.table[base_lo : base_lo + TAPS]
        c_hi = self.table[base_hi : base_hi + TAPS]
        # c = c_lo + (c_hi - c_lo) * t, all f32
        c = c_lo + (c_hi - c_lo) * t
        # Rust accumulates sample-by-sample in f32: `acc += x * c`. numpy's
        # dot on two f32 arrays reduces in f64 internally by default, which
        # would diverge from Rust. Use a manual f32 accumulation loop
        # written in numpy by casting the per-element products back to f32
        # and summing with np.add.reduce at f32 dtype.
        prod = (taps * c).astype(np.float32, copy=False)
        # Force left-to-right f32 accumulation to match Rust's sequential
        # `acc += x * c` (the f32 cascade of adds is not bit-identical to
        # a pairwise tree sum, but both fit in f32 so the difference is
        # <= 1 ULP per add; negligible at our threshold).
        acc = np.float32(0.0)
        for v in prod:
            acc = np.float32(acc + v)
        return acc


def resample(x: np.ndarray, input_sr: float, output_sr: float) -> np.ndarray:
    """Resample `x` from `input_sr` to `output_sr`. Returns f32 ndarray.

    Wraps Resampler.process over a whole buffer. Output length matches
    what Rust's `process` produces given unlimited output space —
    approximately `len(x) * output_sr / input_sr`, possibly ±1 sample.
    """
    r = Resampler(input_sr, output_sr)
    x_f32 = np.asarray(x, dtype=np.float32)
    # Conservative upper bound (matches resample_all helper in resampler.rs tests).
    max_out = int(float(len(x_f32)) / r.step + 2.0)
    out = np.zeros(max_out, dtype=np.float32)
    consumed, produced = r.process(x_f32, out)
    # Rust's process returns when output is full OR input is exhausted;
    # we sized output generously so we expect input to be exhausted.
    if consumed < len(x_f32):
        raise RuntimeError(
            f"resample: only consumed {consumed}/{len(x_f32)} samples "
            f"(output buffer undersized: max_out={max_out}, produced={produced})"
        )
    return out[:produced]
