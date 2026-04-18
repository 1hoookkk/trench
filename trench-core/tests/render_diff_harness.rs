//! Render-diff harness — Hedz cascade + AGC baseline.
//!
//! Produces the canonical golden WAV from the Rust `FilterEngine` running
//! the baked Talking Hedz cartridge. The JUCE plugin's offline-render mode
//! (upcoming) must emit a bit-identical WAV for the same input; `tools/
//! render_diff.py` is the diff tool.
//!
//! **First run** (no committed fixture): generates
//! `tests/fixtures/hedz_reference_48k.wav`, panics with a message asking
//! you to commit the fixture and re-run.
//!
//! **Subsequent runs:** the test re-renders and asserts each output sample
//! matches the committed fixture bit-for-bit (`f32::to_bits` equality).
//! Drift in the Rust reference chain fails the test immediately.
//!
//! Input signal is a deterministic linear chirp (20 Hz → 20 kHz, 1 s,
//! amplitude 0.5, mono, 48 kHz) generated in-process; no binary input
//! fixture is required.
//!
//! Scope of this harness (step 2a of the render-diff doctrine):
//!   cascade → AGC (with mix) → output_gain → DC blocker
//! Post-cascade SPACE, pre-cascade Mackie drive, and `mod_fn` TRIG are
//! layered in once their Rust stages exist (plan Tasks 4, 6, 8) and are
//! folded into `FilterEngine`.

use hound::{SampleFormat, WavReader, WavSpec, WavWriter};
use std::f32::consts::TAU;
use std::fs;
use std::path::{Path, PathBuf};
use trench_core::{Cartridge, FilterEngine};

const SAMPLE_RATE: u32 = 48_000;
const DURATION_SAMPLES: usize = SAMPLE_RATE as usize; // 1.0 s
const INPUT_AMPLITUDE: f32 = 0.5;
const CHIRP_F0_HZ: f32 = 20.0;
const CHIRP_F1_HZ: f32 = 20_000.0;
const MORPH: f64 = 0.5;
const Q: f64 = 0.5;
const FIXTURE_RELPATH: &str = "tests/fixtures/hedz_reference_48k.wav";

/// Linear-frequency chirp, phase-accumulated sample-by-sample so the
/// phase integral matches exactly across runs and platforms. Output is
/// the same whether built in debug or release — no FMA, no SIMD.
fn generate_chirp() -> Vec<f32> {
    let sr = SAMPLE_RATE as f32;
    let n = DURATION_SAMPLES;
    let mut out = Vec::with_capacity(n);
    let mut phase: f32 = 0.0;
    for i in 0..n {
        let t = i as f32 / n as f32;
        let f = CHIRP_F0_HZ + (CHIRP_F1_HZ - CHIRP_F0_HZ) * t;
        phase += TAU * f / sr;
        if phase > TAU {
            phase -= TAU;
        }
        out.push(INPUT_AMPLITUDE * phase.sin());
    }
    out
}

fn render_rust_reference() -> Vec<f32> {
    let mut engine = FilterEngine::new();
    engine.prepare(SAMPLE_RATE as f64);
    engine.load_cartridge(Cartridge::hedz_rom());
    let mut buf = generate_chirp();
    engine.process_block(&mut buf, MORPH, Q);
    assert!(
        !engine.take_instability_flag(),
        "Hedz cascade raised instability flag during reference render"
    );
    buf
}

fn write_f32_wav(path: &Path, samples: &[f32]) {
    let spec = WavSpec {
        channels: 1,
        sample_rate: SAMPLE_RATE,
        bits_per_sample: 32,
        sample_format: SampleFormat::Float,
    };
    let mut w = WavWriter::create(path, spec).expect("create wav");
    for &s in samples {
        w.write_sample(s).expect("write sample");
    }
    w.finalize().expect("finalize wav");
}

fn read_f32_wav(path: &Path) -> (Vec<f32>, WavSpec) {
    let mut r = WavReader::open(path).expect("open wav");
    let spec = r.spec();
    let samples: Vec<f32> = r
        .samples::<f32>()
        .map(|s| s.expect("read sample"))
        .collect();
    (samples, spec)
}

#[test]
fn hedz_cascade_agc_render_matches_committed_fixture() {
    let samples = render_rust_reference();
    let fixture = PathBuf::from(FIXTURE_RELPATH);

    if !fixture.exists() {
        fs::create_dir_all(fixture.parent().unwrap()).expect("mkdir fixtures");
        write_f32_wav(&fixture, &samples);
        panic!(
            "render_diff_harness: generated new fixture at {}. Commit this WAV \
             (git add {}) and re-run the test. This failure is expected exactly \
             once — on the commit that introduces the harness.",
            fixture.display(),
            fixture.display()
        );
    }

    let (committed, spec) = read_f32_wav(&fixture);
    assert_eq!(spec.channels, 1, "fixture must be mono");
    assert_eq!(spec.sample_rate, SAMPLE_RATE, "fixture sample rate mismatch");
    assert_eq!(spec.bits_per_sample, 32);
    assert!(matches!(spec.sample_format, SampleFormat::Float));
    assert_eq!(
        samples.len(),
        committed.len(),
        "render length changed: re-render has {} samples, fixture has {}",
        samples.len(),
        committed.len()
    );

    let mut diffs = 0usize;
    let mut first_diff: Option<(usize, f32, f32)> = None;
    for (i, (&a, &b)) in samples.iter().zip(committed.iter()).enumerate() {
        if a.to_bits() != b.to_bits() {
            diffs += 1;
            if first_diff.is_none() {
                first_diff = Some((i, a, b));
            }
        }
    }
    if diffs != 0 {
        let (i, a, b) = first_diff.unwrap();
        panic!(
            "{diffs} sample(s) drifted vs committed fixture. First drift at \
             sample {i}: re-render={a} (bits {:#x}), fixture={b} (bits {:#x}). \
             The Rust reference chain changed — either revert the change or \
             regenerate the fixture (delete the WAV, re-run, commit).",
            a.to_bits(),
            b.to_bits()
        );
    }
}
