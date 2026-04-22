//! Parity gate: `compiled_talking_hedz` through FilterEngine at 44100 Hz
//! vs `canonical_wav_cal_talking_hedz` (see `FILTER_ARTIFACTS.json`).
//!
//! Authoritative Talking Hedz source is `cal_talking_hedz`
//! (`docs/calibration/Talking_Hedz.json`, RE'd 5-stage runtime model),
//! compiled to `compiled_talking_hedz`
//! (`trench-juce/cartridges/heritage/Talking_Hedz.json`).
//!
//! **Forbidden comparison**: `compiled_talking_hedz` must NOT be nulled
//! against `canonical_wav_p2k_013`, and `rom_hedz` must NOT be nulled
//! against `canonical_wav_cal_talking_hedz`. They are distinct filters
//! despite sharing the display name — see `forbidden_comparisons` in
//! the registry.
//!
//! Canonical render chain (see `ref/canonical/PROVENANCE.json`) —
//! AGC-OFF native-SR rebake:
//!   raw stage → resample host→native → SOS cascade @ 39062.5 Hz
//!             → × boost → resample native→host → float32 WAV
//!
//! FilterEngine here runs at 44100 Hz, resamples to native (39062.5 Hz),
//! runs the compiled-v1 DF2T cascade, and folds per-corner boost into the
//! ramped output gain. **AGC and DC blocker are both disabled** so the gate
//! validates cascade math + coefficient interpolation + resampler + boost.
//!
//! AGC parity is out of scope here on purpose: the AGC is a bit-masked
//! nonlinear state machine (`(gain * abs_sample) as u32 & 0xF` indexing a
//! 16-entry cliff-shaped table, see `agc.rs:12`). Rust runs it at f32, the
//! Python canonical chain at f64 — the index-flip trajectories diverge
//! after thousands of samples into totally different gain states. AGC
//! correctness is covered by the `agc.rs` unit tests (byte-identical port
//! of the EmulatorX.dll reference), not by this parity gate.
//!
//! Canonical refs are rendered through `tools/rust_resampler.py` — a
//! straight port of `trench-core/src/resampler.rs` — so FilterEngine and
//! the canonical chain share the SAME polyphase FIR at the host↔native
//! boundary. The ±256-sample lag search in `null_db_lag` is kept as
//! defense-in-depth against future resampler changes; in practice
//! observed lag is 0.
//!
//! Threshold −30 dB: compile-path (raw→compiled-v1 quantise) plus resampler
//! rounding. Audibility gate, not math gate.

use std::f32::consts::TAU;
use std::fs;
use std::path::PathBuf;

use trench_core::{Cartridge, FilterEngine};

// ─── paths ──────────────────────────────────────────────────────────────────

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("workspace parent")
        .to_path_buf()
}

fn hedz_json_path() -> PathBuf {
    repo_root()
        .join("trench-juce")
        .join("cartridges")
        .join("heritage")
        .join("Talking_Hedz.json")
}

// ─── input signal ────────────────────────────────────────────────────────────

const SR: f64 = 44100.0;
const SIGNAL_LEN: usize = 44100; // 1 s

fn make_chirp() -> Vec<f32> {
    let n = SIGNAL_LEN;
    let sr = SR as f32;
    let mut phase = 0.0f32;
    (0..n)
        .map(|i| {
            let t = i as f32 / n as f32;
            let f = 20.0 + (20_000.0 - 20.0) * t;
            phase += TAU * f / sr;
            if phase > TAU {
                phase -= TAU;
            }
            0.5 * phase.sin()
        })
        .collect()
}

// ─── WAV helpers ─────────────────────────────────────────────────────────────

fn read_wav_f32(path: &std::path::Path) -> Option<Vec<f32>> {
    let mut r = hound::WavReader::open(path).ok()?;
    let spec = r.spec();
    if spec.sample_format != hound::SampleFormat::Float || spec.bits_per_sample != 32 {
        return None;
    }
    let ch = spec.channels as usize;
    let samples: Vec<f32> = r
        .samples::<f32>()
        .filter_map(|s| s.ok())
        .step_by(ch) // take channel 0 only
        .collect();
    Some(samples)
}

// ─── engine render ────────────────────────────────────────────────────────────

fn render(cart: Cartridge, morph: f64, q: f64, input: &[f32]) -> Vec<f32> {
    let mut engine = FilterEngine::new();
    engine.prepare(SR);
    // Canonical refs are AGC-off native-SR rebakes — disable AGC + DC blocker
    // so the gate measures cascade + resampler + coefficient interpolation +
    // boost only. AGC correctness is covered by agc.rs unit tests.
    engine.debug.agc_enabled = false;
    engine.debug.dc_block_enabled = false;
    engine.load_cartridge(cart);
    // Chunk into DAW-sized host blocks. FilterEngine's internal
    // `native_scratch` is sized for ~8192 native samples (≈9252 host
    // samples at 44100/39062.5); passing the full 132 k-sample dry in
    // one shot clamps output to the scratch size and zero-fills the
    // remainder, defeating the parity gate. Real DAWs hand the engine
    // blocks of 64–2048 samples, so 4096 here still mirrors a realistic
    // call pattern while staying well under the scratch ceiling.
    const HOST_BLOCK: usize = 4096;
    let mut buf = input.to_vec();
    let mut offset = 0;
    while offset < buf.len() {
        let end = (offset + HOST_BLOCK).min(buf.len());
        engine.process_block(&mut buf[offset..end], morph, q);
        offset = end;
    }
    buf
}

// ─── metrics ──────────────────────────────────────────────────────────────────

/// Optimal-gain null with ±`max_lag` sample alignment search.
///
/// FilterEngine's 32-tap Rust polyphase resampler and the Python canonical
/// render chain's `scipy.signal.resample_poly` (default Kaiser polyphase)
/// both produce high-quality output but have different constant group
/// delays. A lag=0 null sees the two as uncorrelated even when the shapes
/// match; this search absorbs that constant offset.
///
/// Returns `(best_lag, gain, rel_db)` — the lag (in samples, `pred[i]`
/// aligned with `reference[i + lag]`), the optimal scalar gain at that
/// lag, and the residual RMS relative to the reference RMS in dB.
fn null_db_lag(pred: &[f32], reference: &[f32], max_lag: i64) -> (i64, f64, f64) {
    let mut best = (0i64, 0.0f64, f64::INFINITY);
    for lag in -max_lag..=max_lag {
        let (p_start, r_start) = if lag >= 0 {
            (0usize, lag as usize)
        } else {
            ((-lag) as usize, 0usize)
        };
        let n = pred
            .len()
            .saturating_sub(p_start)
            .min(reference.len().saturating_sub(r_start));
        if n < 1024 {
            continue;
        }
        let mut num = 0.0f64;
        let mut den = 0.0f64;
        for i in 0..n {
            let p = pred[p_start + i] as f64;
            let r = reference[r_start + i] as f64;
            num += p * r;
            den += p * p;
        }
        if den <= 0.0 {
            continue;
        }
        let gain = num / den;
        let mut res2 = 0.0f64;
        let mut r2 = 0.0f64;
        for i in 0..n {
            let p = pred[p_start + i] as f64;
            let r = reference[r_start + i] as f64;
            let residual = r - gain * p;
            res2 += residual * residual;
            r2 += r * r;
        }
        if r2 <= 0.0 {
            continue;
        }
        let rel = 10.0 * (res2 / r2).log10();
        if rel < best.2 {
            best = (lag, gain, rel);
        }
    }
    best
}

// ─── corner table ─────────────────────────────────────────────────────────────

const CORNERS: &[(&str, f64, f64)] = &[
    ("M0_Q0", 0.0, 0.0),
    ("M0_Q100", 0.0, 1.0),
    ("M100_Q0", 1.0, 0.0),
    ("M100_Q100", 1.0, 1.0),
];

// ─── tests ────────────────────────────────────────────────────────────────────

fn load_hedz() -> Option<Cartridge> {
    let path = hedz_json_path();
    if !path.exists() {
        eprintln!(
            "hedz_parity: skipping — Talking_Hedz.json not found at {}",
            path.display()
        );
        return None;
    }
    let json = fs::read_to_string(&path).expect("read Talking_Hedz.json");
    Some(Cartridge::from_json(&json).expect("parse Talking_Hedz.json"))
}

fn canonical_wav_path_named(prefix: &str, corner: &str) -> PathBuf {
    repo_root()
        .join("ref")
        .join("canonical")
        .join(format!("{prefix}_{corner}.wav"))
}

fn load_dry() -> Vec<f32> {
    let dry_path = PathBuf::from(r"C:/Users/hooki/trenchwork_clean/ref/bypassed-pinknoise.wav");
    if dry_path.exists() {
        read_wav_f32(&dry_path).unwrap_or_else(make_chirp)
    } else {
        make_chirp()
    }
}

/// Core null loop shared by the compile-path and runtime-identity gates.
fn run_parity_gate(gate_name: &str, cart: &Cartridge, canonical_prefix: &str, threshold_db: f64) {
    let input = load_dry();
    eprintln!(
        "{gate_name}: dry = {} samples, threshold = {threshold_db} dB",
        input.len()
    );

    let skip_transient = input.len() / 20;
    let mut any_fail = false;

    for &(label, morph, q) in CORNERS {
        let wav_path = canonical_wav_path_named(canonical_prefix, label);
        if !wav_path.exists() {
            eprintln!(
                "  {label:<12} — canonical WAV missing at {}",
                wav_path.display()
            );
            any_fail = true;
            continue;
        }
        let reference = match read_wav_f32(&wav_path) {
            Some(v) => v,
            None => {
                eprintln!("  {label:<12} — could not read WAV");
                any_fail = true;
                continue;
            }
        };

        let pred = render(cart.clone(), morph, q, &input);
        let n = pred.len().min(reference.len());
        let (lag, gain, rel) =
            null_db_lag(&pred[skip_transient..n], &reference[skip_transient..n], 256);
        let ok = rel <= threshold_db;
        eprintln!(
            "  {canonical_prefix:<25} {label:<12}  lag={lag:+4}  gain={gain:.4}  null={rel:+8.2} dB  {}",
            if ok { "OK" } else { "FAIL" }
        );
        if !ok {
            any_fail = true;
        }
    }

    assert!(
        !any_fail,
        "{gate_name}: one or more corners exceeded {threshold_db} dB threshold"
    );
}

/// Compile-path parity: compiled_talking_hedz (Rust DF2T kernel coeffs) nulled
/// against the Python canonical render of cal_talking_hedz (pole_freq_hz,
/// radius → on-the-fly a1 = -2r·cos(...) derivation). Measures how much the
/// compile step `calibration_re → compiled_v1` degrades the filter. Residual
/// lives in f32 quantisation of ~1e-4 coefficient differences through the
/// resonant cascade — amplified at M100_Q0 (narrowest bandwidth). Threshold
/// -30 dB: audibility gate, not math gate.
///
/// AGC + DC blocker disabled on the FilterEngine side (see `render()`); the
/// Python canonical also runs AGC-off. AGC correctness is covered by
/// `trench-core/src/agc.rs` unit tests (byte-identical EmulatorX.dll port).
/// The lag search in `null_db_lag` is kept as defense-in-depth against
/// future resampler changes; in practice observed lag is 0 because
/// `tools/rust_resampler.py` ports `trench-core/src/resampler.rs` verbatim.
#[test]
fn compile_path_parity_compiled_v1_vs_calibration() {
    let hedz = match load_hedz() {
        Some(c) => c,
        None => return,
    };
    run_parity_gate("compile_path_parity", &hedz, "cal_Talking_Hedz", -30.0);
}

/// Runtime identity: compiled_talking_hedz nulled against a Python canonical
/// rendered from the SAME compiled-v1 cartridge. Both sides use bit-identical
/// coefficients, the same ported polyphase FIR, AGC off, boost × 4.0.
///
/// Non-resonant corners (M0_Q0, M0_Q100) null at the f32 WAV noise floor
/// (≈ -180 dB) — proving the runtime path is bit-close to the reference
/// implementation when the filter state is well-behaved.
///
/// Resonant corners (M100_Q0, M100_Q100) have pole radii ≈ 0.999 (nearly
/// on the unit circle), where tiny numerical differences amplify through
/// the feedback loop. Rust now snaps the cascade coefficients and boost to
/// their exact targets at each 32-sample block boundary, so the static-corner
/// runtime converges to the Python reference's fixed-f64 behavior instead of
/// recirculating ulp-scale ramp residue forever.
///
/// Threshold -100 dB: significantly below audibility and comfortably above
/// the f32 WAV floor, so any future block-boundary drift shows up loudly
/// without touching the compile-path gate.
#[test]
fn runtime_identity_compiled_v1_both_sides() {
    let hedz = match load_hedz() {
        Some(c) => c,
        None => return,
    };
    run_parity_gate("runtime_identity", &hedz, "compiled_talking_hedz", -100.0);
}
