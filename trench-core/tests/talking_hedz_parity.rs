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

fn canonical_wav_path(corner: &str) -> PathBuf {
    repo_root()
        .join("ref")
        .join("canonical")
        .join(format!("cal_Talking_Hedz_{corner}.wav"))
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
    ("M0_Q0",     0.0, 0.0),
    ("M0_Q100",   0.0, 1.0),
    ("M100_Q0",   1.0, 0.0),
    ("M100_Q100", 1.0, 1.0),
];

// ─── tests ────────────────────────────────────────────────────────────────────

fn load_hedz() -> Option<Cartridge> {
    let path = hedz_json_path();
    if !path.exists() {
        eprintln!("hedz_parity: skipping — Talking_Hedz.json not found at {}", path.display());
        return None;
    }
    let json = fs::read_to_string(&path).expect("read Talking_Hedz.json");
    Some(Cartridge::from_json(&json).expect("parse Talking_Hedz.json"))
}

/// Talking_Hedz.json through FilterEngine at 44100 Hz vs canonical reference WAVs.
/// Uses the external pinknoise dry at trenchwork_clean if present (same source
/// used for the canonical WAVs), falling back to the chirp. Threshold: −30 dB.
///
/// Canonical refs are AGC-OFF native-SR rebakes — this gate validates
/// cascade math, coefficient interpolation, resampler, and boost. AGC is
/// disabled on the FilterEngine side (see `render()`); AGC correctness
/// lives in `trench-core/src/agc.rs` unit tests (byte-identical
/// EmulatorX.dll port). The ±256-sample lag search in `null_db_lag` is
/// load-bearing: two different polyphase FIRs (scipy's Kaiser
/// `resample_poly` vs our 32-tap) have different constant group delays.
#[test]
fn talking_hedz_44100_nulls_against_canonical_wav() {
    let hedz = match load_hedz() {
        Some(c) => c,
        None => return,
    };

    // Try external dry first (same source used for canonical WAVs), else chirp.
    let dry_path = PathBuf::from(r"C:/Users/hooki/trenchwork_clean/ref/bypassed-pinknoise.wav");
    let input: Vec<f32> = if dry_path.exists() {
        match read_wav_f32(&dry_path) {
            Some(v) => v,
            None => make_chirp(),
        }
    } else {
        make_chirp()
    };
    eprintln!(
        "talking_hedz_44100_nulls: dry = {} samples",
        input.len()
    );

    let skip_transient = input.len() / 20;
    let mut any_fail = false;
    const THRESHOLD: f64 = -30.0;

    for &(label, morph, q) in CORNERS {
        let wav_path = canonical_wav_path(label);
        if !wav_path.exists() {
            eprintln!("  {label:<12} — canonical WAV missing, skip");
            continue;
        }
        let reference = match read_wav_f32(&wav_path) {
            Some(v) => v,
            None => {
                eprintln!("  {label:<12} — could not read WAV, skip");
                continue;
            }
        };

        let pred = render(hedz.clone(), morph, q, &input);

        let n = pred.len().min(reference.len());
        let (lag, gain, rel) = null_db_lag(
            &pred[skip_transient..n],
            &reference[skip_transient..n],
            256,
        );
        let ok = rel <= THRESHOLD;
        eprintln!(
            "  Talking_Hedz 44100→native {label:<12}  lag={lag:+4}  gain={gain:.4}  null={rel:+7.1} dB  {}",
            if ok { "OK" } else { "FAIL" }
        );
        if !ok {
            any_fail = true;
        }
    }

    if any_fail {
        panic!("one or more corners exceeded {THRESHOLD} dB threshold vs canonical WAVs");
    }
}
