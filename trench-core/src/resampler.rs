//! Zero-allocation polyphase FIR resampler for 39062.5 Hz ↔ host SR.
//!
//! The DF2T cascade is acoustically authoritative only at NATIVE_SR. This
//! module converts host-rate audio into and out of the cascade's native rate.
//! Every buffer is allocated in `Resampler::new`; `process` makes zero heap
//! calls, so it is safe to call from the audio thread.
//!
//! # Algorithm
//! Windowed-sinc polyphase FIR.  A `(PHASES + 1) × TAPS` table is built
//! once at construction.  At runtime, each output sample is produced by:
//!   1. Selecting the two polyphase phases that bracket the fractional read
//!      position (phase_lo, phase_lo + 1).
//!   2. Computing a TAPS-point dot-product against the ring buffer for each.
//!   3. Linearly blending the two results.
//!
//! # Latency
//! `TAPS / 2` input samples.  At 44100 Hz, ≈ 0.18 ms — negligible, and the
//! DAW can compensate if the plugin reports it.

/// E-mu native sample rate — the cascade's authoritative acoustic dimension.
pub const NATIVE_SR: f64 = 39_062.5;

/// Taps per polyphase phase. 32 gives ≥ 60 dB stopband at the 39062.5 ↔
/// 44100 and 39062.5 ↔ 48000 ratios with the Blackman window. 16 was too
/// short — the transition band covered 22 kHz, so anti-alias protection
/// collapsed into the stopband shoulder.
const TAPS: usize = 32;

/// Number of precomputed phases (sub-sample positions).  512 phases + linear
/// interpolation gives fractional-delay error < 1e-3 dB.
const PHASES: usize = 512;

/// Power-of-2 ring buffer.  Must be > TAPS.
const RING_LEN: usize = 64;
const RING_MASK: usize = RING_LEN - 1;

/// Windowed-sinc polyphase resampler.  Construct with `Resampler::new`; the
/// sinc table is heap-allocated there.  All audio-thread paths are
/// allocation-free.
pub struct Resampler {
    /// `(PHASES + 1) × TAPS` polyphase sinc table.
    /// Row `p` is the filter for fractional delay `p / PHASES`.
    table: Box<[f32]>,
    ring: [f32; RING_LEN],
    /// Absolute write position (monotonically increasing; taken mod RING_LEN).
    wpos: usize,
    /// Fractional read position (absolute units matching wpos).
    rpos: f64,
    /// Input samples consumed per output sample (= input_sr / output_sr).
    step: f64,
}

impl Resampler {
    /// Construct a resampler from `input_sr` to `output_sr`.
    ///
    /// Allocates the sinc table on the heap here; subsequent `process` calls
    /// allocate nothing.
    pub fn new(input_sr: f64, output_sr: f64) -> Self {
        let step = input_sr / output_sr;
        // Anti-aliasing cutoff: lower of the two Nyquists, normalised to
        // input Nyquist.  For upsamplers (step < 1) this is 1.0 (no
        // anti-aliasing needed); for downsamplers (step > 1) this band-limits
        // before decimation.
        let cutoff = (output_sr / input_sr).min(1.0);
        let table = build_sinc_table(cutoff);
        Self {
            table,
            ring: [0.0; RING_LEN],
            // Pre-fill TAPS/2 zero samples so the first read position (0.0)
            // has a full causal window available without waiting for input.
            wpos: TAPS / 2,
            rpos: 0.0,
            step,
        }
    }

    /// Feed `input` samples and drain resampled samples into `output`.
    ///
    /// Returns `(consumed, produced)`.  Production stops when either `output`
    /// is full or `input` is exhausted.  Callers that need a fixed number of
    /// output samples must loop until the output slice is filled.
    ///
    /// **No heap allocation.**
    pub fn process(&mut self, input: &[f32], output: &mut [f32]) -> (usize, usize) {
        let half = TAPS / 2;
        let mut in_idx = 0;
        let mut out_idx = 0;

        loop {
            // `center` = integer part of the current read position.
            // The symmetric window accesses ring positions
            // center - (half-1) .. center + half.
            // We need wpos > center + half, i.e. the lookahead half must be
            // in the ring.
            let center = self.rpos as usize;
            let need = center + half;

            while self.wpos <= need {
                if in_idx >= input.len() {
                    return (in_idx, out_idx);
                }
                self.ring[self.wpos & RING_MASK] = input[in_idx];
                self.wpos += 1;
                in_idx += 1;
            }

            if out_idx >= output.len() {
                return (in_idx, out_idx);
            }

            output[out_idx] = self.interpolate();
            out_idx += 1;
            self.rpos += self.step;
        }
    }

    /// Reset delay state.  Call on transport reset or when a gap in audio
    /// occurs.  Does not reallocate the sinc table.
    pub fn reset(&mut self) {
        self.ring = [0.0; RING_LEN];
        self.wpos = TAPS / 2;
        self.rpos = 0.0;
    }

    /// How many output samples `n_input` input samples will yield (floor).
    /// Use this to size the output slice when consuming a full input block.
    pub fn output_count(&self, n_input: usize) -> usize {
        // After consuming n_input, rpos advances by n_input * step.
        // The number of complete output samples equals the integer advances.
        ((n_input as f64) / self.step) as usize
    }

    #[inline]
    fn interpolate(&self) -> f32 {
        let center = self.rpos.floor() as i64;
        let frac = (self.rpos - center as f64) as f32;
        let half = (TAPS / 2) as i64;

        // Select the two polyphase phases bracketing `frac`.
        let phase_f = frac * PHASES as f32;
        let phase_lo = phase_f as usize;           // 0 .. PHASES-1
        let phase_hi = (phase_lo + 1).min(PHASES); // 1 .. PHASES
        let t = phase_f - phase_lo as f32;         // blend weight for phase_hi

        let base_lo = phase_lo * TAPS;
        let base_hi = phase_hi * TAPS;

        let mut acc = 0.0f32;
        for k in 0..TAPS {
            // Symmetric window: tap k spans position center-half+1+k.
            let pos = (center - half + 1 + k as i64) as usize & RING_MASK;
            let x = self.ring[pos];
            let c = self.table[base_lo + k] + (self.table[base_hi + k] - self.table[base_lo + k]) * t;
            acc += x * c;
        }
        acc
    }
}

// ---------------------------------------------------------------------------
// Sinc table construction (called once at Resampler::new, not on audio thread)
// ---------------------------------------------------------------------------

/// Build a `(PHASES + 1) × TAPS` Blackman-windowed sinc table.
///
/// Each row `p` is the polyphase FIR kernel for fractional delay `p / PHASES`.
/// The table is normalised so that row 0 (zero fractional delay) has unity
/// DC gain: any passband signal is reproduced at full amplitude.
fn build_sinc_table(cutoff: f64) -> Box<[f32]> {
    use std::f64::consts::PI;

    // Prototype filter length: TAPS × PHASES taps.
    let n = TAPS * PHASES;
    let half_n = (n as f64) / 2.0;

    // Compute windowed-sinc prototype. `t` is in *real* output-sample units
    // (one unit = one real sample), NOT prototype-sample units — the
    // prototype is oversampled by PHASES so each real sample spans PHASES
    // prototype samples. Expressing `t` in real units is what makes the
    // polyphase sub-filters correct fractional-delay FIRs.
    let mut proto = vec![0.0f64; n + 1];
    for i in 0..=n {
        let t = (i as f64 - half_n) / PHASES as f64;
        let sinc = if t.abs() < 1e-10 {
            cutoff
        } else {
            (PI * cutoff * t).sin() / (PI * t)
        };
        // Blackman window — better stopband than Hann for ≥80 dB rejection.
        let w = 0.42
            - 0.5 * (2.0 * PI * i as f64 / n as f64).cos()
            + 0.08 * (4.0 * PI * i as f64 / n as f64).cos();
        proto[i] = sinc * w;
    }

    // Polyphase decomposition: table[p][k] = proto[k * PHASES + p].
    //
    // Each polyphase sub-filter is a fractional-delay FIR in its own right,
    // so each must have unity DC gain. Normalising by phase 0 alone is
    // wrong because phase 0's sum varies wildly with `cutoff`: for
    // cutoff=1.0 the sinc zero-crossings at every integer t collapse phase
    // 0 to a one-tap delta with sum ≈ 1, while every other phase sums to
    // ≈ 1/PHASES; for cutoff<1 the ratio shifts again. Per-phase
    // normalisation is independent of cutoff and keeps every fractional
    // delay at unity passband gain.
    let size = (PHASES + 1) * TAPS;
    let mut table = vec![0.0f32; size];
    for p in 0..=PHASES {
        let phase_sum: f64 = (0..TAPS)
            .map(|k| {
                let idx = k * PHASES + p;
                if idx <= n { proto[idx] } else { 0.0 }
            })
            .sum();
        let phase_scale = if phase_sum.abs() > 1e-12 {
            1.0 / phase_sum
        } else {
            0.0
        };
        for k in 0..TAPS {
            let idx = k * PHASES + p;
            if idx <= n {
                table[p * TAPS + k] = (proto[idx] * phase_scale) as f32;
            }
        }
    }

    table.into_boxed_slice()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // Helpers ----------------------------------------------------------------

    fn make_sine(sr: f64, freq: f64, n: usize) -> Vec<f32> {
        (0..n)
            .map(|i| (2.0 * std::f64::consts::PI * freq * i as f64 / sr).sin() as f32)
            .collect()
    }

    fn rms(buf: &[f32]) -> f32 {
        let sum: f32 = buf.iter().map(|s| s * s).sum();
        (sum / buf.len() as f32).sqrt()
    }

    fn resample_all(r: &mut Resampler, input: &[f32]) -> Vec<f32> {
        // Conservative upper bound on output length.
        let max_out = (input.len() as f64 / r.step + 2.0) as usize;
        let mut out = vec![0.0f32; max_out];
        let mut in_off = 0;
        let mut out_off = 0;
        while in_off < input.len() {
            let (consumed, produced) = r.process(&input[in_off..], &mut out[out_off..]);
            in_off += consumed;
            out_off += produced;
            if consumed == 0 && produced == 0 {
                break;
            }
        }
        out.truncate(out_off);
        out
    }

    // Tests ------------------------------------------------------------------

    /// Silence in → silence out (no DC injection from filter).
    #[test]
    fn silence_is_silence() {
        let mut r = Resampler::new(44100.0, NATIVE_SR);
        let zeros = vec![0.0f32; 1024];
        let out = resample_all(&mut r, &zeros);
        assert!(!out.is_empty());
        for &s in &out {
            assert!(s.abs() < 1e-6, "silence produced {s}");
        }
    }

    /// Down then up should preserve a mid-band sine within ±1 dB.
    #[test]
    fn roundtrip_preserves_sine() {
        let host_sr = 44100.0;
        let freq = 1000.0;
        let n = 4096;

        let input = make_sine(host_sr, freq, n);
        let rms_in = rms(&input[256..]); // skip transient

        let mut down = Resampler::new(host_sr, NATIVE_SR);
        let native = resample_all(&mut down, &input);

        let mut up = Resampler::new(NATIVE_SR, host_sr);
        let output = resample_all(&mut up, &native);

        let rms_out = rms(&output[256..]);
        let ratio = rms_out / rms_in;
        assert!(
            (0.891..=1.122).contains(&ratio),
            "roundtrip gain {ratio:.4} outside ±1 dB (rms_in={rms_in:.4}, rms_out={rms_out:.4})"
        );
    }

    /// Frequencies above native Nyquist (~19.5 kHz) must be strongly
    /// attenuated in the downsampled output.
    #[test]
    fn antialias_blocks_above_native_nyquist() {
        let host_sr = 44100.0;
        let alias_freq = 22000.0; // above NATIVE_SR/2 ≈ 19531 Hz
        let n = 4096;

        let input = make_sine(host_sr, alias_freq, n);
        let rms_in = rms(&input[256..]);

        let mut down = Resampler::new(host_sr, NATIVE_SR);
        let native = resample_all(&mut down, &input);

        // Must be heavily attenuated in the native-rate stream.
        let rms_native = rms(&native[256..]);
        assert!(
            rms_native < rms_in * 0.1,
            "alias at {alias_freq} Hz not attenuated: rms_in={rms_in:.4}, rms_native={rms_native:.4}"
        );
    }

    /// 48 kHz ↔ NATIVE_SR roundtrip.
    #[test]
    fn roundtrip_48k() {
        let host_sr = 48000.0;
        let freq = 440.0;
        let n = 4096;

        let input = make_sine(host_sr, freq, n);
        let rms_in = rms(&input[256..]);

        let mut down = Resampler::new(host_sr, NATIVE_SR);
        let native = resample_all(&mut down, &input);

        let mut up = Resampler::new(NATIVE_SR, host_sr);
        let output = resample_all(&mut up, &native);

        let rms_out = rms(&output[256..]);
        let ratio = rms_out / rms_in;
        assert!(
            (0.891..=1.122).contains(&ratio),
            "48k roundtrip gain {ratio:.4} outside ±1 dB"
        );
    }

    /// Reset clears delay state.
    #[test]
    fn reset_clears_state() {
        let mut r = Resampler::new(44100.0, NATIVE_SR);
        // Feed a loud signal.
        let loud = vec![1.0f32; 128];
        resample_all(&mut r, &loud);
        // Reset, then silence should produce silence.
        r.reset();
        let zeros = vec![0.0f32; 128];
        let out = resample_all(&mut r, &zeros);
        for &s in &out {
            assert!(s.abs() < 1e-5, "after reset: expected silence, got {s}");
        }
    }

    /// `output_count` returns a conservative lower bound.
    #[test]
    fn output_count_is_conservative() {
        let r = Resampler::new(44100.0, NATIVE_SR);
        let estimated = r.output_count(32);
        // Actual will be estimated or estimated+1.
        assert!(estimated >= 27, "output_count too low: {estimated}");
        assert!(estimated <= 30, "output_count too high: {estimated}");
    }
}
