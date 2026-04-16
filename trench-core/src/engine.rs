//! FilterEngine — thin runtime wrapper over the frozen `Cascade`.
//!
//! Order (matches `trenchwork_clean/trench-core/src/engine.rs` and the
//! `tools/parity_null.py` reference pipeline):
//! cascade -> AGC (with mix) -> output_gain (boost, ramped) -> DC blocker.
//!
//! The Trench `Cascade` already owns 12-stage biquad DSP, bilinear
//! interpolation, and per-sample coefficient ramping, so this module is
//! additive only: control-rate parameter dispatch, cascade-peak gain
//! ceiling, AGC, DC blocker.
//!
//! Doctrine: `Cascade` is frozen. This module does not touch its
//! internals. The cascade-peak estimator uses direct biquad form
//! `(b0, b1, b2, a1, a2)` — not the shifted minifloat-domain kernel from
//! the `trenchwork_clean` engine. See SESSION_STATE rev 4 for the split.

use crate::agc::agc_step;
use crate::cartridge::{Cartridge, CornerData};
#[cfg(test)]
use crate::cartridge::{NUM_COEFFS, NUM_STAGES};
use crate::cascade::{Cascade, BLOCK_SIZE};

/// First-order DC blocker (~20 Hz highpass).
///     y[n] = x[n] - x[n-1] + R * y[n-1]
#[derive(Debug, Clone, Copy)]
struct DcBlocker {
    x_prev: f32,
    y_prev: f32,
    r: f32,
}

impl DcBlocker {
    fn new(sample_rate: f64) -> Self {
        let fc = 20.0_f64;
        let r = 1.0 - (2.0 * std::f64::consts::PI * fc / sample_rate);
        Self {
            x_prev: 0.0,
            y_prev: 0.0,
            r: r as f32,
        }
    }

    #[inline(always)]
    fn process(&mut self, x: f32) -> f32 {
        let y = x - self.x_prev + self.r * self.y_prev;
        self.x_prev = x;
        self.y_prev = if y.is_finite() { y } else { 0.0 };
        self.y_prev
    }
}

/// Runtime debug toggles. All default on to match shipping behavior.
#[derive(Debug, Clone, Copy)]
pub struct DebugToggles {
    pub agc_enabled: bool,
    pub dc_block_enabled: bool,
    pub output_gain_clamp_enabled: bool,
}

impl Default for DebugToggles {
    fn default() -> Self {
        Self {
            agc_enabled: true,
            dc_block_enabled: true,
            output_gain_clamp_enabled: true,
        }
    }
}

pub struct FilterEngine {
    cascade: Cascade,
    cartridge: Option<Cartridge>,
    sample_rate: f64,
    output_gain: f32,
    target_output_gain: f32,
    delta_output_gain: f32,
    agc_gain: f32,
    agc_mix: f32,
    gain_ceiling: f32,
    dc_blocker: DcBlocker,
    pub debug: DebugToggles,
}

impl Default for FilterEngine {
    fn default() -> Self {
        Self::new()
    }
}

impl FilterEngine {
    pub fn new() -> Self {
        Self {
            cascade: Cascade::new(),
            cartridge: None,
            sample_rate: 44100.0,
            output_gain: 1.0,
            target_output_gain: 1.0,
            delta_output_gain: 0.0,
            agc_gain: 1.0,
            agc_mix: 1.0,
            gain_ceiling: 100.0,
            dc_blocker: DcBlocker::new(44100.0),
            debug: DebugToggles::default(),
        }
    }

    /// Reset runtime state and set the sample rate. Call before processing.
    pub fn prepare(&mut self, sample_rate: f64) {
        self.sample_rate = sample_rate;
        self.cascade.reset();
        self.output_gain = 1.0;
        self.target_output_gain = 1.0;
        self.delta_output_gain = 0.0;
        self.agc_gain = 1.0;
        self.dc_blocker = DcBlocker::new(sample_rate);
    }

    /// Install a cartridge. Does not reset biquad state (avoids clicks on
    /// live cartridge switch).
    pub fn load_cartridge(&mut self, cart: Cartridge) {
        self.cartridge = Some(cart);
    }

    /// AGC wet/dry mix. 0.0 = bypass, 1.0 = full AGC. The AGC state
    /// machine runs unconditionally so gain tracking stays warm.
    pub fn set_agc_mix(&mut self, mix: f32) {
        self.agc_mix = mix.clamp(0.0, 1.0);
    }

    /// Linear gain ceiling for the cascade-peak safety clamp. Default 100.0 (40 dB).
    pub fn set_gain_ceiling(&mut self, ceiling: f32) {
        self.gain_ceiling = ceiling.max(1.0);
    }

    /// Control-rate parameter dispatch for one sub-block of `chunk_size`
    /// samples. Interpolates the cartridge at (morph, q), hands the
    /// corner to the cascade (which ramps per-sample internally), and
    /// computes a ramped `output_gain` that folds in the per-corner boost
    /// and the peak-based gain ceiling clamp.
    fn set_parameters(&mut self, morph: f64, q: f64, chunk_size: usize) {
        let cart = match &self.cartridge {
            Some(c) => c,
            None => return,
        };

        let corner: CornerData = cart.interpolate(morph, q);
        let mut gain = cart.interpolate_boost(morph, q) as f32;

        if self.debug.output_gain_clamp_enabled {
            let peak = compute_cascade_peak(&corner, self.sample_rate);
            if peak * gain > self.gain_ceiling {
                gain *= self.gain_ceiling / (peak * gain);
            }
        }

        self.cascade.set_targets(&corner, chunk_size);
        // Boost is applied externally in `process_sample_inner` (after AGC,
        // matching the reference pipeline). Keep cascade internal boost
        // pinned to 1.0 so it doesn't double-apply.
        self.cascade.set_boost(1.0, chunk_size);

        self.target_output_gain = gain;
        self.delta_output_gain = (gain - self.output_gain) / chunk_size.max(1) as f32;
    }

    #[inline]
    fn process_sample_inner(&mut self, input: f32) -> f32 {
        // 1. Cascade (frozen) — includes per-sample coefficient ramp and
        //    internal boost (pinned to 1.0 here).
        let mut signal = self.cascade.tick(input);

        // 2. AGC with mix blend. State machine always runs.
        if self.debug.agc_enabled {
            let agc_out = agc_step(signal, &mut self.agc_gain);
            signal += (agc_out - signal) * self.agc_mix;
        }

        // 3. Output gain (per-corner boost × ceiling clamp, ramped).
        self.output_gain += self.delta_output_gain;
        signal *= self.output_gain;

        // 4. DC blocker — last in chain, matches C++ order.
        if self.debug.dc_block_enabled {
            signal = self.dc_blocker.process(signal);
        }

        signal
    }

    /// Snap ramped output gain to target at block boundary (prevent float
    /// drift). The cascade handles its own coefficient snap internally.
    fn snap_to_target(&mut self) {
        self.output_gain = self.target_output_gain;
        self.delta_output_gain = 0.0;
    }

    /// Process a mono f32 buffer in place. `morph` and `q` are control-rate
    /// intents held constant across the buffer; the engine ramps coefficients
    /// and output gain per-sample inside each 32-sample control block.
    pub fn process_block(&mut self, data: &mut [f32], morph: f64, q: f64) {
        if self.cartridge.is_none() {
            return;
        }
        let mut offset = 0;
        while offset < data.len() {
            let remaining = data.len() - offset;
            let chunk = remaining.min(BLOCK_SIZE);

            self.set_parameters(morph, q, chunk);

            for i in 0..chunk {
                data[offset + i] = self.process_sample_inner(data[offset + i]);
            }

            self.snap_to_target();
            offset += chunk;
        }
    }

    /// Returns whether the cascade encountered non-finite state since the
    /// last call. Forwarded from `Cascade::take_instability_flag`.
    pub fn take_instability_flag(&mut self) -> bool {
        self.cascade.take_instability_flag()
    }
}

/// Peak |H(e^jω)| across 1025 linearly-spaced bins from DC to Nyquist,
/// evaluated on the cascade target as direct-form biquads
/// `(b0, b1, b2, a1, a2)` = `(c0, c1, c2, c3, c4)`.
fn compute_cascade_peak(corner: &CornerData, sample_rate: f64) -> f32 {
    const NUM_BINS: usize = 1024;
    let nyquist = sample_rate / 2.0;
    let mut max_mag: f64 = 0.0;

    for bin in 0..=NUM_BINS {
        let freq = nyquist * (bin as f64 / NUM_BINS as f64);
        let omega = 2.0 * std::f64::consts::PI * freq / sample_rate;
        let cos_w = omega.cos();
        let cos_2w = (2.0 * omega).cos();
        let sin_w = omega.sin();
        let sin_2w = (2.0 * omega).sin();

        let mut cascade_mag_sq: f64 = 1.0;

        for stage in corner.iter() {
            // Direct biquad form — see trench-core doctrine and
            // SESSION_STATE rev 4 truth on `compiled-v1` kernel mapping.
            let b0 = stage[0];
            let b1 = stage[1];
            let b2 = stage[2];
            let a1 = stage[3];
            let a2 = stage[4];

            let num_real = b0 + b1 * cos_w + b2 * cos_2w;
            let num_imag = -(b1 * sin_w + b2 * sin_2w);
            let den_real = 1.0 + a1 * cos_w + a2 * cos_2w;
            let den_imag = -(a1 * sin_w + a2 * sin_2w);

            let den_mag_sq = den_real * den_real + den_imag * den_imag;
            if den_mag_sq > 1e-30 {
                let stage_mag_sq = (num_real * num_real + num_imag * num_imag) / den_mag_sq;
                cascade_mag_sq *= stage_mag_sq;
            }
        }

        let mag = cascade_mag_sq.sqrt();
        if mag > max_mag {
            max_mag = mag;
        }
    }

    max_mag as f32
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_passthrough_cartridge() -> Cartridge {
        let json = r#"{
            "format": "compiled-v1",
            "name": "passthrough",
            "sampleRate": 44100,
            "keyframes": [
                {"label": "M0_Q0",     "morph": 0.0, "q": 0.0,   "boost": 1.0,
                 "stages": [{"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},
                            {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},
                            {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},
                            {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},
                            {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},
                            {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0}]},
                {"label": "M0_Q100",   "morph": 0.0, "q": 1.0,   "boost": 1.0,
                 "stages": [{"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},
                            {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},
                            {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},
                            {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},
                            {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},
                            {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0}]},
                {"label": "M100_Q0",   "morph": 1.0, "q": 0.0,   "boost": 1.0,
                 "stages": [{"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},
                            {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},
                            {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},
                            {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},
                            {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},
                            {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0}]},
                {"label": "M100_Q100", "morph": 1.0, "q": 1.0,   "boost": 1.0,
                 "stages": [{"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},
                            {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},
                            {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},
                            {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},
                            {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},
                            {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0}]}
            ]
        }"#;
        Cartridge::from_json(json).expect("passthrough cartridge parses")
    }

    #[test]
    fn empty_engine_passes_buffer_through_unchanged() {
        let mut engine = FilterEngine::new();
        engine.prepare(44100.0);
        let mut buf = vec![0.1_f32, -0.2, 0.3, -0.4];
        let copy = buf.clone();
        engine.process_block(&mut buf, 0.0, 0.0);
        // No cartridge → process_block is a no-op.
        assert_eq!(buf, copy);
    }

    #[test]
    fn passthrough_cartridge_preserves_energy() {
        let mut engine = FilterEngine::new();
        engine.prepare(44100.0);
        engine.load_cartridge(make_passthrough_cartridge());

        // Impulse through passthrough stages + DC blocker. Energy should
        // survive the DC blocker to within a few percent for an impulse.
        let mut buf = vec![0.0_f32; 256];
        buf[0] = 1.0;
        let input_sum_sq: f32 = buf.iter().map(|&s| s * s).sum();

        engine.process_block(&mut buf, 0.5, 0.5);

        let output_sum_sq: f32 = buf.iter().map(|&s| s * s).sum();
        assert!(
            (output_sum_sq - input_sum_sq).abs() < 0.05,
            "impulse energy drifted: in={input_sum_sq} out={output_sum_sq}"
        );
        assert!(!engine.take_instability_flag());
    }

    #[test]
    fn gain_ceiling_clamps_runaway_boost() {
        // A passthrough-coefficient cartridge with a 1000x boost. With
        // gain_ceiling = 10.0 the engine should clamp the output gain
        // well under the raw 1000 target.
        let json = r#"{
            "format": "compiled-v1",
            "name": "loud",
            "sampleRate": 44100,
            "keyframes": [
                {"label": "M0_Q0",     "morph": 0.0, "q": 0.0, "boost": 1000.0,
                 "stages": [{"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},
                            {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},
                            {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},
                            {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},
                            {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},
                            {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0}]},
                {"label": "M0_Q100",   "morph": 0.0, "q": 1.0, "boost": 1000.0,
                 "stages": [{"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},
                            {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},
                            {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},
                            {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},
                            {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},
                            {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0}]},
                {"label": "M100_Q0",   "morph": 1.0, "q": 0.0, "boost": 1000.0,
                 "stages": [{"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},
                            {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},
                            {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},
                            {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},
                            {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},
                            {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0}]},
                {"label": "M100_Q100", "morph": 1.0, "q": 1.0, "boost": 1000.0,
                 "stages": [{"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},
                            {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},
                            {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},
                            {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},
                            {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},
                            {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0}]}
            ]
        }"#;
        let cart = Cartridge::from_json(json).unwrap();

        let mut engine = FilterEngine::new();
        engine.prepare(44100.0);
        engine.set_gain_ceiling(10.0);
        engine.debug.agc_enabled = false;
        engine.debug.dc_block_enabled = false;
        engine.load_cartridge(cart);

        let mut buf = vec![0.5_f32; 128];
        engine.process_block(&mut buf, 0.5, 0.5);

        let peak = buf.iter().map(|s| s.abs()).fold(0.0_f32, f32::max);
        assert!(
            peak <= 10.0 + 1e-3,
            "gain ceiling did not clamp: peak={peak}"
        );
        // And raw 0.5 * 1000 = 500 should definitely have been avoided.
        assert!(peak < 100.0, "peak={peak} too high, clamp not applied");
    }

    #[test]
    fn num_coeffs_matches_corner_shape() {
        // Compile-time sanity: corners are [[f64; 5]; NUM_STAGES].
        let _: [[f64; NUM_COEFFS]; NUM_STAGES] = [[0.0; NUM_COEFFS]; NUM_STAGES];
    }
}
