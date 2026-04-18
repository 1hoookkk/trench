//! Peak safety limiter (Task 10).
//!
//! A deliberately simple output-stage backstop for the TRENCH signal
//! chain. The point is not transparent mastering-grade limiting — it's
//! catching the rare runaway that escapes the Q AGC / cascade without
//! audibly pumping on normal material.
//!
//! ## Semantics
//!
//! - Ceiling: **-0.3 dBFS** (`CEILING_DB`). Linear ceiling is
//!   `10^(-0.3/20) ≈ 0.9661`.
//! - Attack: **instantaneous** (1 sample). Any peak above the ceiling
//!   immediately drops `current_gain` to whatever value places the
//!   instantaneous peak exactly at the ceiling.
//! - Release: **50 ms** (`RELEASE_MS`), one-pole exponential decay of
//!   `(1 - current_gain)` toward 0 — i.e. gain rising back toward 1.0.
//! - Stereo-linked: the same gain value multiplies both channels every
//!   sample, so the stereo image does not wander.
//! - **Transparent below ceiling**: when `|L|, |R| < ceiling_linear`
//!   *and* `current_gain == 1.0`, output is bit-identical to input.
//!   Specifically, `x * 1.0` followed by `x.clamp(-c, c)` where `|x|
//!   < c` returns the original bits, so the transparent path is exact.
//! - **Hard-clamp backstop**: after applying gain, the sample is
//!   clamped to `[-ceiling_linear, +ceiling_linear]`. This catches
//!   first-sample runaways that the release-state gain can't attenuate
//!   in a single step (e.g. the very first sample is +10.0 with
//!   `current_gain` still at 1.0).
//!
//! ## Contract
//!
//! - `process_stereo` never allocates.
//! - `reset()` sets `current_gain` back to 1.0 (unity).
//! - Sample-rate-correct release: `release_coeff = exp(-1 / (release_ms
//!   * sr / 1000))`. Verified by `sample_rate_independent_release_time`.

/// Peak limiter with a hardcoded -0.3 dBFS ceiling and 50 ms release.
pub struct SafetyLimiter {
    ceiling_lin: f32,
    release_coeff: f32,
    /// 1.0 at rest; drops to `ceiling / |peak|` when peak exceeds
    /// ceiling; releases back toward 1.0 via one-pole.
    current_gain: f32,
}

impl SafetyLimiter {
    /// Ceiling in dBFS. Public constant for test access; no runtime setter.
    pub const CEILING_DB: f32 = -0.3;
    /// Release time in milliseconds (one-pole time constant τ).
    pub const RELEASE_MS: f32 = 50.0;

    pub fn new(sample_rate: f32) -> Self {
        let sr = sample_rate.max(1.0);
        let release_samples = (Self::RELEASE_MS * 0.001 * sr).max(1.0);
        Self {
            ceiling_lin: 10f32.powf(Self::CEILING_DB / 20.0),
            release_coeff: (-1.0_f32 / release_samples).exp(),
            current_gain: 1.0,
        }
    }

    /// Process one stereo sample in-place. Allocation-free.
    #[inline]
    pub fn process_stereo(&mut self, l: &mut f32, r: &mut f32) {
        let peak = l.abs().max(r.abs());
        // Required gain to place instantaneous peak at the ceiling.
        let required_gain = if peak > self.ceiling_lin {
            self.ceiling_lin / peak
        } else {
            1.0
        };
        // Attack: instant drop. Release: one-pole toward 1.0.
        if required_gain < self.current_gain {
            self.current_gain = required_gain;
        } else {
            // current_gain <- 1 - (1 - current_gain) * release_coeff
            self.current_gain = 1.0 - (1.0 - self.current_gain) * self.release_coeff;
        }
        *l *= self.current_gain;
        *r *= self.current_gain;
        // Hard-clamp backstop: catches first-sample runaways where the
        // release-state gain had no prior opportunity to attenuate.
        *l = l.clamp(-self.ceiling_lin, self.ceiling_lin);
        *r = r.clamp(-self.ceiling_lin, self.ceiling_lin);
    }

    /// Reset internal gain state to unity. Leaves coefficients intact.
    pub fn reset(&mut self) {
        self.current_gain = 1.0;
    }
}
