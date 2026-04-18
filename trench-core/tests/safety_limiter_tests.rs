//! Integration tests for `safety_limiter` (Task 10).
//!
//! Peak limiter with -0.3 dBFS ceiling, stereo-linked gain reduction,
//! instantaneous attack + 50 ms one-pole release. Transparent (bit-
//! identical) below the ceiling; hard-clamp backstop for extreme
//! first-sample runaways.

use trench_core::safety_limiter::SafetyLimiter;

const SR_48K: f32 = 48_000.0;
const SR_96K: f32 = 96_000.0;
const EPS: f32 = 1.0e-6;

fn ceiling_lin() -> f32 {
    10f32.powf(SafetyLimiter::CEILING_DB / 20.0)
}

/// Run the limiter in-place over paired L/R slices.
fn process_buffer(lim: &mut SafetyLimiter, l: &mut [f32], r: &mut [f32]) {
    assert_eq!(l.len(), r.len());
    for i in 0..l.len() {
        lim.process_stereo(&mut l[i], &mut r[i]);
    }
}

// ── 1. Transparent below ceiling ─────────────────────────────────────

#[test]
fn transparent_below_ceiling() {
    let mut lim = SafetyLimiter::new(SR_48K);

    // A mix of values all strictly below ceiling_linear ≈ 0.9661.
    let input_l: Vec<f32> = (0..1024)
        .map(|i| 0.5 * ((i as f32) * 0.01).sin())
        .collect();
    let input_r: Vec<f32> = (0..1024)
        .map(|i| 0.3 * ((i as f32) * 0.013).cos())
        .collect();
    let mut out_l = input_l.clone();
    let mut out_r = input_r.clone();

    process_buffer(&mut lim, &mut out_l, &mut out_r);

    for i in 0..input_l.len() {
        assert_eq!(
            out_l[i].to_bits(),
            input_l[i].to_bits(),
            "L[{i}] not bit-identical: in={} out={}",
            input_l[i], out_l[i]
        );
        assert_eq!(
            out_r[i].to_bits(),
            input_r[i].to_bits(),
            "R[{i}] not bit-identical: in={} out={}",
            input_r[i], out_r[i]
        );
    }
}

// ── 2. Peak above ceiling is clamped ─────────────────────────────────

#[test]
fn peak_above_ceiling_is_clamped() {
    let mut lim = SafetyLimiter::new(SR_48K);
    let c = ceiling_lin();

    // Baseline ~-20 dBFS with a single loud spike at sample 100.
    let n = 512;
    let mut l: Vec<f32> = (0..n).map(|i| 0.1 * ((i as f32) * 0.05).sin()).collect();
    let mut r: Vec<f32> = vec![0.0; n];
    l[100] = 0.99; // well above -0.3 dBFS ceiling

    process_buffer(&mut lim, &mut l, &mut r);

    let peak = l.iter().map(|x| x.abs()).fold(0.0_f32, f32::max);
    assert!(
        peak <= c + EPS,
        "output peak {peak} exceeds ceiling {c} (tol {EPS})"
    );
}

// ── 3. Stereo link preserves balance ─────────────────────────────────

#[test]
fn stereo_link_preserves_balance() {
    let mut lim = SafetyLimiter::new(SR_48K);

    // L has a loud transient; R is well below ceiling at the same
    // instant. Pre-seed a quiet signal on both so there's no dormant
    // gain state from earlier samples.
    let mut l = vec![0.05_f32; 64];
    let mut r = vec![0.02_f32; 64];
    let peak_idx = 32;
    let l_in = 0.99_f32;
    let r_in = 0.2_f32;
    l[peak_idx] = l_in;
    r[peak_idx] = r_in;

    process_buffer(&mut lim, &mut l, &mut r);

    let l_out = l[peak_idx];
    let r_out = r[peak_idx];
    let in_ratio = l_in / r_in;
    let out_ratio = l_out / r_out;
    assert!(
        (in_ratio - out_ratio).abs() < 1e-5,
        "stereo balance lost: in_ratio={in_ratio}, out_ratio={out_ratio}, \
         l_in={l_in}, r_in={r_in}, l_out={l_out}, r_out={r_out}"
    );
}

// ── 4. Release returns to unity ──────────────────────────────────────

#[test]
fn release_returns_to_unity() {
    let mut lim = SafetyLimiter::new(SR_48K);

    // One loud sample to engage gain reduction.
    let mut l0 = 0.99_f32;
    let mut r0 = 0.0_f32;
    lim.process_stereo(&mut l0, &mut r0);

    // Then ~200 ms (4× release) of silence.
    let n = (0.200 * SR_48K) as usize;
    let mut l = vec![0.0_f32; n];
    let mut r = vec![0.0_f32; n];
    process_buffer(&mut lim, &mut l, &mut r);

    // Probe gain via a small test sample: output/input should be ~1.0.
    let probe = 0.1_f32;
    let mut l_probe = probe;
    let mut r_probe = probe;
    lim.process_stereo(&mut l_probe, &mut r_probe);
    let gain_l = l_probe / probe;
    let gain_r = r_probe / probe;
    assert!(
        (gain_l - 1.0).abs() < 1e-3,
        "gain did not return to unity after release: gain_l={gain_l}"
    );
    assert!(
        (gain_r - 1.0).abs() < 1e-3,
        "gain did not return to unity after release: gain_r={gain_r}"
    );
}

// ── 5. Hard clamp on extreme input ───────────────────────────────────

#[test]
fn hard_clamp_on_extreme_input() {
    let mut lim = SafetyLimiter::new(SR_48K);
    let c = ceiling_lin();

    // An extreme first-sample runaway. Release-based gain state starts
    // at 1.0 so the output would naively be clamped to ceiling via the
    // hard-clamp backstop rather than by the gain-reduction path.
    let mut l = 10.0_f32;
    let mut r = -10.0_f32;
    lim.process_stereo(&mut l, &mut r);

    assert!(
        l.abs() <= c + EPS,
        "L not clamped: got {l}, ceiling {c}"
    );
    assert!(
        r.abs() <= c + EPS,
        "R not clamped: got {r}, ceiling {c}"
    );
}

// ── 6. Sample-rate-independent release time ──────────────────────────

#[test]
fn sample_rate_independent_release_time() {
    // Measure the gain envelope N samples into the release at both SRs,
    // where N corresponds to the same wall-clock time. The measured
    // gain values should be very close if the release coefficient
    // scales correctly with sample rate.

    fn gain_after_ms(sr: f32, ms: f32) -> f32 {
        let mut lim = SafetyLimiter::new(sr);
        // Trigger limiting with one loud sample.
        let mut l = 0.99_f32;
        let mut r = 0.0_f32;
        lim.process_stereo(&mut l, &mut r);
        // Run `ms` of silence to release.
        let n = (ms * 0.001 * sr) as usize;
        for _ in 0..n {
            let mut li = 0.0_f32;
            let mut ri = 0.0_f32;
            lim.process_stereo(&mut li, &mut ri);
        }
        // Probe.
        let probe = 0.1_f32;
        let mut lp = probe;
        let mut rp = probe;
        lim.process_stereo(&mut lp, &mut rp);
        lp / probe
    }

    let g48 = gain_after_ms(SR_48K, 25.0);
    let g96 = gain_after_ms(SR_96K, 25.0);

    // Gains at same wall time should be within ~10%.
    let rel_err = (g48 - g96).abs() / g48.max(1e-6);
    assert!(
        rel_err < 0.10,
        "release time not sample-rate independent: g48={g48}, g96={g96}, rel_err={rel_err}"
    );
}

// ── 7. Reset returns gain state to unity ─────────────────────────────

#[test]
fn reset_clears_gain_state() {
    let mut lim = SafetyLimiter::new(SR_48K);
    // Engage limiting.
    let mut l = 0.99_f32;
    let mut r = 0.0_f32;
    lim.process_stereo(&mut l, &mut r);
    // Reset.
    lim.reset();
    // Next sample should be transparent (below ceiling, unity gain).
    let probe = 0.1_f32;
    let mut lp = probe;
    let mut rp = probe;
    lim.process_stereo(&mut lp, &mut rp);
    assert_eq!(lp.to_bits(), probe.to_bits(), "reset did not restore unity gain on L");
    assert_eq!(rp.to_bits(), probe.to_bits(), "reset did not restore unity gain on R");
}
