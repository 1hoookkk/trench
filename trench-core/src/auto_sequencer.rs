//! AUTO sequencer — host-synced cycling through 4 snapshot slots.
//!
//! The shipping v1 AUTO mode (see memory `project_auto_mode`) cycles
//! through captured (morph, Q) snapshots in lockstep with host transport.
//! No generation, no reactivity: just a clock-driven index into up to 4
//! user-filled slots with an optional beat-relative glide between
//! consecutive slots.
//!
//! Transition vocabulary tracks the E-mu Function Generator lineage
//! (Lag Processor / DC Delay vs Linear): `glide_beats = 0.0` = **Hard
//! Snap** (jump instantly at the step boundary), `glide_beats > 0.0` =
//! **Lagged Glide** (slew linearly over that many beats, clamped to the
//! step length). Hardcoded conditional-jump logic and the heritage 63
//! transition shapes are out of scope.
//!
//! This module holds the reference math. The C++ plugin mirrors the same
//! evaluate_at_beats() logic in `trench-juce/plugin/source/AutoSequencer.h`
//! for runtime use — the two are kept small so drift stays visible.

pub const MAX_SLOTS: usize = 4;

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Slot {
    pub filled: bool,
    pub morph: f32,
    pub q: f32,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AutoOutput {
    pub morph: f32,
    pub q: f32,
}

/// Pure evaluator — stateless. All transport/host inputs are passed in.
///
/// `beats_position` is the host beat clock (ppq position) or a free-running
/// beat counter when the host isn't playing. `beats_per_step` comes from the
/// rate table. `glide_beats` is clamped to `[0, beats_per_step]` internally:
/// you cannot glide longer than the dwell. `glide_beats == 0.0` is the
/// Hard Snap case.
///
/// `base_morph` / `base_q` are the raw knob values used as fallback when
/// no slots are filled.
///
/// `just_enabled` = true on the first block after AUTO turns on; forces
/// the glide fraction to 1.0 so output snaps to the current slot with no
/// retroactive glide from a previous slot that wasn't actually playing.
pub fn evaluate_at_beats(
    beats_position: f64,
    beats_per_step: f64,
    glide_beats: f64,
    slots: &[Slot; MAX_SLOTS],
    base_morph: f32,
    base_q: f32,
    just_enabled: bool,
) -> AutoOutput {
    // Collect filled slot indices onto the stack.
    let mut filled: [usize; MAX_SLOTS] = [0; MAX_SLOTS];
    let mut filled_count: usize = 0;
    for (i, s) in slots.iter().enumerate() {
        if s.filled {
            filled[filled_count] = i;
            filled_count += 1;
        }
    }

    if filled_count == 0 {
        return AutoOutput { morph: base_morph, q: base_q };
    }
    if filled_count == 1 {
        let s = &slots[filled[0]];
        return AutoOutput { morph: s.morph, q: s.q };
    }

    // Step position: raw_step = beats / beats_per_step; cycle within filled_count.
    let bps = if beats_per_step > 0.0 { beats_per_step } else { 1.0 };
    let raw_step = beats_position / bps;
    let cycle = filled_count as f64;
    let step_mod = raw_step.rem_euclid(cycle);
    let step_index = step_mod.floor() as usize;
    let step_index = step_index.min(filled_count - 1);
    let intra_phase = step_mod - step_index as f64;

    let cur = &slots[filled[step_index]];
    let prev_idx = (step_index + filled_count - 1) % filled_count;
    let prev = &slots[filled[prev_idx]];

    // Glide fraction in [0, 1]. glide = 0 or just_enabled → Hard Snap.
    let glide_t = if just_enabled || glide_beats <= 0.0 {
        1.0_f64
    } else {
        let glide_phase = (glide_beats / bps).min(1.0);
        if glide_phase <= 0.0 {
            1.0
        } else if intra_phase < glide_phase {
            intra_phase / glide_phase
        } else {
            1.0
        }
    };

    let t = glide_t as f32;
    AutoOutput {
        morph: prev.morph + (cur.morph - prev.morph) * t,
        q: prev.q + (cur.q - prev.q) * t,
    }
}

/// Stateful wrapper: tracks free-run beat counter for when the host isn't
/// providing ppq. Single `f64` of mutable RT state — matches the plan.
pub struct AutoSequencer {
    free_run_beats: f64,
    was_enabled: bool,
}

impl Default for AutoSequencer {
    fn default() -> Self {
        Self::new()
    }
}

impl AutoSequencer {
    pub const fn new() -> Self {
        Self { free_run_beats: 0.0, was_enabled: false }
    }

    pub fn reset(&mut self) {
        self.free_run_beats = 0.0;
        self.was_enabled = false;
    }

    pub fn free_run_beats(&self) -> f64 {
        self.free_run_beats
    }

    /// Advance the free-run beat counter by one block's worth of samples.
    /// Returns the new beat position.
    pub fn advance_free_run(
        &mut self,
        num_samples: usize,
        sample_rate: f64,
        bpm: f64,
    ) -> f64 {
        if sample_rate > 0.0 && bpm > 0.0 && num_samples > 0 {
            let beats_delta = (num_samples as f64) * bpm / (60.0 * sample_rate);
            self.free_run_beats += beats_delta;
        }
        self.free_run_beats
    }

    /// Single-call evaluator for the audio thread.
    ///
    /// `ppq` = Some(host ppq) when the host provides it, else None and we
    /// use the internal free-run counter (which is advanced by this call).
    ///
    /// `enabled` is the AUTO toggle. On a false→true transition the
    /// first block returns the current slot with no glide from prev.
    pub fn evaluate(
        &mut self,
        enabled: bool,
        ppq: Option<f64>,
        num_samples: usize,
        sample_rate: f64,
        bpm: f64,
        beats_per_step: f64,
        glide_beats: f64,
        slots: &[Slot; MAX_SLOTS],
        base_morph: f32,
        base_q: f32,
    ) -> AutoOutput {
        if !enabled {
            self.was_enabled = false;
            return AutoOutput { morph: base_morph, q: base_q };
        }

        let just_enabled = !self.was_enabled;
        if just_enabled {
            self.free_run_beats = 0.0;
        }
        self.was_enabled = true;

        let beats = match ppq {
            Some(p) => p,
            None => {
                let b = self.free_run_beats;
                self.advance_free_run(num_samples, sample_rate, bpm);
                b
            }
        };

        evaluate_at_beats(
            beats,
            beats_per_step,
            glide_beats,
            slots,
            base_morph,
            base_q,
            just_enabled,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f32 = 1e-6;

    fn s(morph: f32, q: f32) -> Slot {
        Slot { filled: true, morph, q }
    }

    fn empty_slots() -> [Slot; MAX_SLOTS] {
        [Slot::default(); MAX_SLOTS]
    }

    #[test]
    fn step_boundary_math() {
        // 3 filled slots, beats_per_step = 1.0, ppq = 2.5.
        // raw_step = 2.5, step_mod = 2.5 % 3 = 2.5 → step=2, phase=0.5.
        let slots = [s(0.1, 0.9), s(0.2, 0.8), s(0.3, 0.7), Slot::default()];
        // Zero glide so we can infer step from output == slot[step] exactly.
        let out = evaluate_at_beats(2.5, 1.0, 0.0, &slots, 0.0, 0.0, false);
        assert!((out.morph - slots[2].morph).abs() < EPS);
        assert!((out.q - slots[2].q).abs() < EPS);
    }

    #[test]
    fn glide_phase() {
        // glide=0.5, bps=1.0 → glide_phase = 0.5.
        // intra=0.25 → glide_t = 0.5; intra=0.6 → glide_t = 1.0.
        let slots = [s(0.0, 0.0), s(1.0, 1.0), Slot::default(), Slot::default()];

        // intra = 0.25 at step 1 → beats = 1.25.
        // step 1 slot = (1.0,1.0), prev slot = (0.0,0.0), glide_t = 0.5 → out = 0.5.
        let out = evaluate_at_beats(1.25, 1.0, 0.5, &slots, 0.0, 0.0, false);
        assert!((out.morph - 0.5).abs() < EPS);
        assert!((out.q - 0.5).abs() < EPS);

        // intra = 0.6 → glide_t = 1.0 → out = cur = (1.0, 1.0).
        let out = evaluate_at_beats(1.6, 1.0, 0.5, &slots, 0.0, 0.0, false);
        assert!((out.morph - 1.0).abs() < EPS);
        assert!((out.q - 1.0).abs() < EPS);
    }

    #[test]
    fn slot_skip_when_empty() {
        // slots = {0=F, 1=E, 2=F, 3=E}, 2 filled.
        // step 0 → slot 0, step 1 → slot 2, step 2 → slot 0 (wrap).
        let mut slots = empty_slots();
        slots[0] = s(0.11, 0.19);
        slots[2] = s(0.22, 0.28);
        // Zero glide for clean step readout.

        // beats = 0.0 → step 0 → slot 0.
        let o = evaluate_at_beats(0.0, 1.0, 0.0, &slots, 0.0, 0.0, false);
        assert!((o.morph - 0.11).abs() < EPS);
        assert!((o.q - 0.19).abs() < EPS);

        // beats = 1.0 → step 1 → slot 2.
        let o = evaluate_at_beats(1.0, 1.0, 0.0, &slots, 0.0, 0.0, false);
        assert!((o.morph - 0.22).abs() < EPS);
        assert!((o.q - 0.28).abs() < EPS);

        // beats = 2.0 → step 2 → filled[2 % 2] = filled[0] → slot 0.
        let o = evaluate_at_beats(2.0, 1.0, 0.0, &slots, 0.0, 0.0, false);
        assert!((o.morph - 0.11).abs() < EPS);
        assert!((o.q - 0.19).abs() < EPS);
    }

    #[test]
    fn zero_filled_returns_manual() {
        let slots = empty_slots();
        let o = evaluate_at_beats(7.42, 0.5, 0.25, &slots, 0.4, 0.6, false);
        assert!((o.morph - 0.4).abs() < EPS);
        assert!((o.q - 0.6).abs() < EPS);
    }

    #[test]
    fn single_filled_slot_holds() {
        let mut slots = empty_slots();
        slots[1] = s(0.33, 0.77);
        // Any ppq, any glide → holds slot 1.
        for &b in &[0.0, 0.5, 1.5, 100.0] {
            let o = evaluate_at_beats(b, 1.0, 0.5, &slots, 0.0, 0.0, false);
            assert!((o.morph - 0.33).abs() < EPS, "beats={b}");
            assert!((o.q - 0.77).abs() < EPS, "beats={b}");
        }
    }

    #[test]
    fn glide_zero_is_hard_snap() {
        // glide_beats = 0 → glide_t = 1.0 at any intra_phase (Hard Snap).
        let slots = [s(0.0, 0.0), s(1.0, 1.0), Slot::default(), Slot::default()];
        for &b in &[1.0_f64, 1.1, 1.5, 1.99] {
            let o = evaluate_at_beats(b, 1.0, 0.0, &slots, 0.0, 0.0, false);
            assert!((o.morph - 1.0).abs() < EPS, "beats={b}");
            assert!((o.q - 1.0).abs() < EPS, "beats={b}");
        }
    }

    #[test]
    fn glide_at_full_beats_is_always_gliding() {
        // glide_beats == beats_per_step → glide_t == intra_phase.
        // Slot 0 = (0,0), Slot 1 = (1,1). Over step 1, output should equal intra_phase.
        let slots = [s(0.0, 0.0), s(1.0, 1.0), Slot::default(), Slot::default()];
        for &intra in &[0.0_f64, 0.25, 0.5, 0.75] {
            let beats = 1.0 + intra;
            let o = evaluate_at_beats(beats, 1.0, 1.0, &slots, 0.0, 0.0, false);
            let expected = intra as f32;
            assert!((o.morph - expected).abs() < 1e-5, "intra={intra}");
            assert!((o.q - expected).abs() < 1e-5, "intra={intra}");
        }
    }

    #[test]
    fn free_run_advance() {
        // bpm=120, sr=48000, num_samples=24000 → 0.5s → 1.0 beat.
        let mut seq = AutoSequencer::new();
        let b = seq.advance_free_run(24_000, 48_000.0, 120.0);
        assert!((b - 1.0).abs() < 1e-9, "b={b}");

        // Another 24000 samples → 2.0.
        let b = seq.advance_free_run(24_000, 48_000.0, 120.0);
        assert!((b - 2.0).abs() < 1e-9, "b={b}");

        // bpm=60, sr=100, num_samples=100 → 1.0s → 1.0 beat added → 3.0.
        let b = seq.advance_free_run(100, 100.0, 60.0);
        assert!((b - 3.0).abs() < 1e-9, "b={b}");
    }

    #[test]
    fn auto_just_enabled_snaps_to_first_filled_slot() {
        // 2 filled slots; first eval after enable should output slot 0 with
        // no fade-in from slot 1 (prev in the wraparound).
        let slots = [s(0.25, 0.75), s(0.9, 0.1), Slot::default(), Slot::default()];
        let mut seq = AutoSequencer::new();

        // Not yet enabled → returns base.
        let o = seq.evaluate(
            false, Some(0.0), 0, 48_000.0, 120.0, 1.0, 0.5, &slots, 0.05, 0.05,
        );
        assert!((o.morph - 0.05).abs() < EPS);
        assert!((o.q - 0.05).abs() < EPS);

        // First enabled block at ppq=0 → just_enabled=true → snap to slot 0,
        // no glide from slot 1 despite glide=0.5.
        let o = seq.evaluate(
            true, Some(0.0), 0, 48_000.0, 120.0, 1.0, 0.5, &slots, 0.0, 0.0,
        );
        assert!((o.morph - 0.25).abs() < EPS);
        assert!((o.q - 0.75).abs() < EPS);

        // Second block at ppq=0.1 (intra=0.1) still in step 0 → cur=slot0,
        // prev=slot1 (wrap). With just_enabled cleared and glide=0.5,
        // glide_phase=0.5, intra_phase=0.1 → glide_t = 0.2 →
        // out = prev + (cur - prev)*0.2 = slot1*0.8 + slot0*0.2.
        let o = seq.evaluate(
            true, Some(0.1), 0, 48_000.0, 120.0, 1.0, 0.5, &slots, 0.0, 0.0,
        );
        let expected_m = 0.9 * 0.8 + 0.25 * 0.2;
        let expected_q = 0.1 * 0.8 + 0.75 * 0.2;
        assert!((o.morph - expected_m).abs() < 1e-5, "m={}", o.morph);
        assert!((o.q - expected_q).abs() < 1e-5, "q={}", o.q);
    }

    #[test]
    fn disabled_returns_base_and_clears_first_cycle() {
        let slots = [s(0.25, 0.75), s(0.9, 0.1), Slot::default(), Slot::default()];
        let mut seq = AutoSequencer::new();

        // Enable, step forward, then disable.
        let _ = seq.evaluate(true, Some(0.0), 0, 48_000.0, 120.0, 1.0, 0.0, &slots, 0.0, 0.0);
        let _ = seq.evaluate(true, Some(2.0), 0, 48_000.0, 120.0, 1.0, 0.0, &slots, 0.0, 0.0);
        let o = seq.evaluate(false, Some(2.5), 0, 48_000.0, 120.0, 1.0, 0.0, &slots, 0.42, 0.58);
        assert!((o.morph - 0.42).abs() < EPS);
        assert!((o.q - 0.58).abs() < EPS);

        // Re-enable → should snap (just_enabled true again).
        let o = seq.evaluate(true, Some(0.0), 0, 48_000.0, 120.0, 1.0, 0.5, &slots, 0.0, 0.0);
        assert!((o.morph - 0.25).abs() < EPS);
        assert!((o.q - 0.75).abs() < EPS);
    }
}
