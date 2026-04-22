//! Validates FilterEngine with the resampled path active (host SR ≠ NATIVE_SR).
//!
//! Every test that calls `engine.prepare(44100.0)` or `engine.prepare(48000.0)`
//! exercises the down-resampler (host→native) and up-resampler (native→host)
//! that were wired in the 2026-04-20 session.

use trench_core::{Cartridge, FilterEngine};

fn passthrough_cart() -> Cartridge {
    let json = r#"{
        "format": "compiled-v1",
        "name": "passthrough",
        "sampleRate": 44100,
        "keyframes": [
            {"label":"M0_Q0",     "morph":0.0,"q":0.0,  "boost":1.0,
             "stages":[{"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},
                       {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},
                       {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},
                       {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},
                       {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},
                       {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0}]},
            {"label":"M0_Q100",   "morph":0.0,"q":1.0,  "boost":1.0,
             "stages":[{"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},
                       {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},
                       {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},
                       {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},
                       {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},
                       {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0}]},
            {"label":"M100_Q0",   "morph":1.0,"q":0.0,  "boost":1.0,
             "stages":[{"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},
                       {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},
                       {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},
                       {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},
                       {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},
                       {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0}]},
            {"label":"M100_Q100", "morph":1.0,"q":1.0,  "boost":1.0,
             "stages":[{"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},
                       {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},
                       {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},
                       {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},
                       {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},
                       {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0}]}
        ]
    }"#;
    Cartridge::from_json(json).expect("passthrough cart")
}

fn make_sine(sr: f64, freq: f64, n: usize) -> Vec<f32> {
    (0..n)
        .map(|i| (2.0 * std::f64::consts::PI * freq * i as f64 / sr).sin() as f32)
        .collect()
}

fn rms(buf: &[f32]) -> f32 {
    if buf.is_empty() {
        return 0.0;
    }
    let sum: f32 = buf.iter().map(|s| s * s).sum();
    (sum / buf.len() as f32).sqrt()
}

/// Silence through a passthrough at 44100 must remain silence (no DC from
/// the resamplers).
#[test]
fn silence_remains_silence_44100() {
    let mut engine = FilterEngine::new();
    engine.prepare(44100.0);
    engine.load_cartridge(passthrough_cart());
    let mut buf = vec![0.0f32; 512];
    engine.process_block(&mut buf, 0.0, 0.0);
    for (i, &s) in buf.iter().enumerate() {
        assert!(s.abs() < 1e-5, "sample {i}: expected silence, got {s}");
    }
}

/// Same at 48 kHz.
#[test]
fn silence_remains_silence_48000() {
    let mut engine = FilterEngine::new();
    engine.prepare(48000.0);
    engine.load_cartridge(passthrough_cart());
    let mut buf = vec![0.0f32; 512];
    engine.process_block(&mut buf, 0.0, 0.0);
    for (i, &s) in buf.iter().enumerate() {
        assert!(s.abs() < 1e-5, "sample {i}: expected silence, got {s}");
    }
}

/// 1 kHz sine at 44100 through a passthrough. Energy must survive the
/// down-resample + cascade + up-resample within ±2 dB (steady-state).
#[test]
fn passthrough_preserves_1khz_energy_44100() {
    let n = 4096;
    let mut engine = FilterEngine::new();
    engine.prepare(44100.0);
    engine.debug.agc_enabled = false;
    engine.debug.dc_block_enabled = false;
    engine.load_cartridge(passthrough_cart());

    let mut buf = make_sine(44100.0, 1000.0, n);
    let rms_in = rms(&buf[256..]); // skip initial transient

    engine.process_block(&mut buf, 0.5, 0.5);

    let rms_out = rms(&buf[256..]);
    let ratio = rms_out / rms_in;
    assert!(
        (0.794..=1.259).contains(&ratio),
        "passthrough energy ratio {ratio:.4} outside ±2 dB (in={rms_in:.4}, out={rms_out:.4})"
    );
}

/// 440 Hz at 48 kHz, same gate.
#[test]
fn passthrough_preserves_440hz_energy_48000() {
    let n = 4096;
    let mut engine = FilterEngine::new();
    engine.prepare(48000.0);
    engine.debug.agc_enabled = false;
    engine.debug.dc_block_enabled = false;
    engine.load_cartridge(passthrough_cart());

    let mut buf = make_sine(48000.0, 440.0, n);
    let rms_in = rms(&buf[256..]);

    engine.process_block(&mut buf, 0.5, 0.5);

    let rms_out = rms(&buf[256..]);
    let ratio = rms_out / rms_in;
    assert!(
        (0.794..=1.259).contains(&ratio),
        "48k passthrough ratio {ratio:.4} outside ±2 dB"
    );
}

/// No instability flag on passthrough with resampling active.
#[test]
fn resampled_passthrough_no_instability() {
    let mut engine = FilterEngine::new();
    engine.prepare(44100.0);
    engine.load_cartridge(passthrough_cart());
    let mut buf = make_sine(44100.0, 440.0, 1024);
    engine.process_block(&mut buf, 0.5, 0.5);
    assert!(!engine.take_instability_flag());
}

/// reset_dsp() clears resampler delay lines — silence after reset.
#[test]
fn reset_dsp_clears_resamplers() {
    let mut engine = FilterEngine::new();
    engine.prepare(44100.0);
    engine.load_cartridge(passthrough_cart());

    // Feed a loud signal to load up the resampler ring buffers.
    let loud = vec![1.0f32; 256];
    let mut loud = loud;
    engine.process_block(&mut loud, 0.0, 0.0);

    // Reset, then silence should produce silence.
    engine.reset_dsp();
    let mut zeros = vec![0.0f32; 256];
    engine.process_block(&mut zeros, 0.0, 0.0);
    for (i, &s) in zeros.iter().enumerate() {
        assert!(
            s.abs() < 1e-4,
            "sample {i} after reset: expected silence, got {s}"
        );
    }
}

/// process_block still works with no cartridge loaded (no-op).
#[test]
fn no_cartridge_is_noop_at_44100() {
    let mut engine = FilterEngine::new();
    engine.prepare(44100.0);
    let original = vec![0.1f32, -0.2, 0.3, -0.4];
    let mut buf = original.clone();
    engine.process_block(&mut buf, 0.0, 0.0);
    assert_eq!(buf, original);
}

/// Hedz ROM cartridge loads and processes without instability at 44100 Hz.
#[test]
fn hedz_rom_no_instability_44100() {
    let mut engine = FilterEngine::new();
    engine.prepare(44100.0);
    engine.load_cartridge(Cartridge::hedz_rom());
    let mut buf = make_sine(44100.0, 440.0, 2048);
    engine.process_block(&mut buf, 0.5, 0.5);
    assert!(!engine.take_instability_flag());
}
