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
use crate::cartridge::{Cartridge, CornerData, SurfaceControlMode};
#[cfg(test)]
use crate::cartridge::{NUM_COEFFS, NUM_STAGES};
use crate::cascade::{Cascade, BLOCK_SIZE};
use crate::resampler::{Resampler, NATIVE_SR};

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
    /// Rate at which the cascade actually runs — NATIVE_SR when resampling, else sample_rate.
    cascade_sr: f64,
    target_x: f64,
    target_y: f64,
    target_z: f64,
    latched_y: Option<f64>,
    latched_z: Option<f64>,
    output_gain: f32,
    target_output_gain: f32,
    delta_output_gain: f32,
    agc_gain: f32,
    agc_mix: f32,
    gain_ceiling: f32,
    dc_blocker: DcBlocker,
    pub debug: DebugToggles,
    resampler_down: Option<Resampler>,
    resampler_up: Option<Resampler>,
    native_scratch: Vec<f32>,
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
            cascade_sr: 44100.0,
            target_x: 0.0,
            target_y: 0.0,
            target_z: 0.0,
            latched_y: None,
            latched_z: None,
            output_gain: 1.0,
            target_output_gain: 1.0,
            delta_output_gain: 0.0,
            agc_gain: 1.0,
            agc_mix: 1.0,
            gain_ceiling: 100.0,
            dc_blocker: DcBlocker::new(44100.0),
            debug: DebugToggles::default(),
            resampler_down: None,
            resampler_up: None,
            native_scratch: Vec::new(),
        }
    }

    /// Reset runtime state and set the sample rate. Call before processing (not on audio thread).
    pub fn prepare(&mut self, sample_rate: f64) {
        self.sample_rate = sample_rate;
        if (sample_rate - NATIVE_SR).abs() > 0.1 {
            self.cascade_sr = NATIVE_SR;
            self.resampler_down = Some(Resampler::new(sample_rate, NATIVE_SR));
            self.resampler_up = Some(Resampler::new(NATIVE_SR, sample_rate));
            self.native_scratch.resize(8192, 0.0);
        } else {
            self.cascade_sr = sample_rate;
            self.resampler_down = None;
            self.resampler_up = None;
        }
        self.cascade.reset();
        self.target_x = 0.0;
        self.target_y = 0.0;
        self.target_z = 0.0;
        self.latched_y = None;
        self.latched_z = None;
        self.output_gain = 1.0;
        self.target_output_gain = 1.0;
        self.delta_output_gain = 0.0;
        self.agc_gain = 1.0;
        self.dc_blocker = DcBlocker::new(self.cascade_sr);
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
    /// samples. Interpolates the cartridge at the active surface position, hands the
    /// corner to the cascade (which ramps per-sample internally), and
    /// computes a ramped `output_gain` that folds in the per-corner boost
    /// and the peak-based gain ceiling clamp.
    fn set_parameters_xyz(
        &mut self,
        x: f64,
        y: f64,
        z: f64,
        trigger_event: bool,
        chunk_size: usize,
    ) {
        let cart = match &self.cartridge {
            Some(c) => c,
            None => return,
        };

        self.target_x = x;
        self.target_y = y;
        self.target_z = z;

        let (resolved_y, resolved_z) = match cart.control_mode() {
            SurfaceControlMode::MorphQ => (y, 0.0),
            SurfaceControlMode::ModernLiveXyz => (y, z),
            SurfaceControlMode::LegacyLatchYz => {
                if trigger_event {
                    self.latched_y = Some(y);
                    self.latched_z = Some(z);
                }
                match (self.latched_y, self.latched_z) {
                    (Some(ly), Some(lz)) => (ly, lz),
                    _ => (y, z),
                }
            }
        };

        let corner: CornerData = cart.interpolate_xyz(x, resolved_y, resolved_z);
        let mut gain = cart.interpolate_boost_xyz(x, resolved_y, resolved_z) as f32;

        if self.debug.output_gain_clamp_enabled {
            let peak = compute_cascade_peak(&corner, self.cascade_sr);
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
        self.process_block_xyz_with_trigger(data, morph, q, 0.0, false);
    }

    pub fn process_block_xyz(&mut self, data: &mut [f32], x: f64, y: f64, z: f64) {
        self.process_block_xyz_with_trigger(data, x, y, z, false);
    }

    pub fn process_block_xyz_with_trigger(
        &mut self,
        data: &mut [f32],
        x: f64,
        y: f64,
        z: f64,
        trigger_event: bool,
    ) {
        if self.cartridge.is_none() {
            return;
        }

        if self.resampler_down.is_some() {
            // Resampled path: host SR → NATIVE_SR → cascade → host SR.
            // Use take/put to satisfy the borrow checker without allocating.
            let n_host = data.len();
            let m_max = ((n_host as f64 * NATIVE_SR / self.sample_rate).ceil() as usize + 4)
                .min(self.native_scratch.len());

            // Downsample host → native scratch.
            let produced_native = {
                let mut down = self.resampler_down.take().unwrap();
                let (_, n) = down.process(data, &mut self.native_scratch[..m_max]);
                self.resampler_down = Some(down);
                n
            };

            // Run cascade at NATIVE_SR on scratch buffer.
            let mut offset = 0;
            while offset < produced_native {
                let chunk = (produced_native - offset).min(BLOCK_SIZE);
                self.set_parameters_xyz(x, y, z, trigger_event && offset == 0, chunk);
                for i in 0..chunk {
                    let s = self.native_scratch[offset + i];
                    self.native_scratch[offset + i] = self.process_sample_inner(s);
                }
                self.snap_to_target();
                offset += chunk;
            }

            // Upsample native scratch → host output.
            let produced_host = {
                let mut up = self.resampler_up.take().unwrap();
                let (_, n) = up.process(&self.native_scratch[..produced_native], data);
                self.resampler_up = Some(up);
                n
            };

            // Fill any trailing samples the resampler couldn't produce yet.
            for s in &mut data[produced_host..] {
                *s = 0.0;
            }
            return;
        }

        // Direct path: host SR already equals cascade SR.
        let mut offset = 0;
        while offset < data.len() {
            let remaining = data.len() - offset;
            let chunk = remaining.min(BLOCK_SIZE);
            self.set_parameters_xyz(x, y, z, trigger_event && offset == 0, chunk);
            for i in 0..chunk {
                data[offset + i] = self.process_sample_inner(data[offset + i]);
            }
            self.snap_to_target();
            offset += chunk;
        }
    }

    /// Clear biquad state, resampler delay lines, and DC blocker without
    /// reallocating the sinc table. Call on transport reset.
    pub fn reset_dsp(&mut self) {
        self.cascade.reset();
        self.dc_blocker = DcBlocker::new(self.cascade_sr);
        self.agc_gain = 1.0;
        self.output_gain = 1.0;
        self.target_output_gain = 1.0;
        self.delta_output_gain = 0.0;
        if let Some(r) = &mut self.resampler_down {
            r.reset();
        }
        if let Some(r) = &mut self.resampler_up {
            r.reset();
        }
    }

    pub fn sample_rate(&self) -> f64 {
        self.sample_rate
    }

    /// Returns whether the cascade encountered non-finite state since the
    /// last call. Forwarded from `Cascade::take_instability_flag`.
    pub fn take_instability_flag(&mut self) -> bool {
        self.cascade.take_instability_flag()
    }

    pub fn latched_coordinates(&self) -> Option<(f64, f64)> {
        match (self.latched_y, self.latched_z) {
            (Some(y), Some(z)) => Some((y, z)),
            _ => None,
        }
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
        // Strip AGC and DC blocker so this tests only the cascade + resampler
        // path. AGC is a non-linear time-varying gain that ducks hard on an
        // impulse then recovers, which legitimately changes the output energy
        // by ~50%; that's audible-pipeline behaviour, not what this test is
        // meant to catch.
        engine.debug.agc_enabled = false;
        engine.debug.dc_block_enabled = false;
        engine.load_cartridge(make_passthrough_cartridge());

        // Band-limited input so the resampler's passband is tested rather
        // than the impulse energy redistributing across its stop-band. A
        // 1 kHz tone sits well inside passband for both 44100 and NATIVE_SR.
        let len = 1024;
        let mut buf: Vec<f32> = (0..len)
            .map(|i| {
                (2.0 * std::f64::consts::PI * 1000.0 * i as f64 / 44100.0).sin() as f32
            })
            .collect();
        let input_sum_sq: f32 = buf.iter().skip(128).map(|&s| s * s).sum();

        engine.process_block(&mut buf, 0.5, 0.5);

        let output_sum_sq: f32 = buf.iter().skip(128).map(|&s| s * s).sum();
        let ratio = output_sum_sq / input_sum_sq;
        assert!(
            (0.95..=1.05).contains(&ratio),
            "passband energy drifted: in={input_sum_sq} out={output_sum_sq} ratio={ratio}"
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
