//! Rust-side parity gate: trench-core Cascade + AGC + boost vs python
//! reference WAVs in `ref/canonical/`.
//!
//! Loads raw P2K skin JSON from `datasets/p2k_skins/`, computes
//! direct-biquad coefficients in f64 (matching tools/parity_null.py's
//! `stage_coefficients_raw`), feeds them into a fresh `Cascade` for each
//! corner, processes the dry pink noise, applies the ported `agc_step`,
//! multiplies by the per-skin boost, and asserts the result nulls against
//! the python-rendered reference at better than -120 dB.
//!
//! Why -120 dB and not -140 dB (python's threshold): the rust pipeline
//! does the f32 cast at the cascade output and applies AGC + boost in
//! f32, while the python pipeline keeps f64 through AGC + boost and only
//! casts to f32 when writing the WAV. -120 dB is comfortable margin.
//!
//! Why raw skins instead of compiled cartridges: the cartridges in
//! `cartridges/p2k/` were quantized through an f32 round-trip in the
//! compile path (~6e-8 per coefficient), which catastrophically
//! amplifies in high-Q stages and pushes nulls to -60..-120 dB.
//! That is a separate finding about the compile path, not about cascade
//! math. The cascade itself should match python at near-f32 precision
//! when fed identical f64 coefficients.
//!
//! Coverage: 33 raw P2K skins (6-stage). Skips the 2 vocal alternates
//! (3-stage) for now — they need a padding strategy.

use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};

use serde::Deserialize;
use trench_core::{Cascade, CornerData, BLOCK_SIZE, NUM_STAGES};

/// f64 reimplementation of `trench_core::agc::agc_step` for the parity
/// test only. The shipping `agc_step` operates in f32 because the JUCE
/// runtime buffers are f32, but the python parity references are
/// rendered through an f64 AGC. AGC is path-dependent (the branch index
/// is `int(gain * |sample|) & 0xF`), so f32 quantization at any sample
/// where the cascade output exceeds ~1.0 produces a different gain
/// trajectory and tens of dB of residual against the f64 reference.
/// Comparing against an f64 inline mirror keeps the test focused on the
/// cascade math; the f32-vs-f64 AGC question lives in agc.rs.
const AGC_TABLE: [f64; 16] = [
    1.0001, 1.0001, 0.996, 0.990, 0.920, 0.500, 0.200, 0.160, 0.120, 0.120, 0.120, 0.120,
    0.120, 0.120, 0.120, 0.120,
];

#[inline(always)]
fn agc_step_f64(sample: f64, agc_gain: &mut f64) -> f64 {
    let abs_sample = sample.abs();
    let idx = ((*agc_gain * abs_sample) as u32 & 0xF) as usize;
    let new_gain = *agc_gain * AGC_TABLE[idx];
    *agc_gain = if new_gain < 1.0 { new_gain } else { 1.0 };
    sample * *agc_gain
}

const FAIL_THRESHOLD_DB: f64 = -120.0;
const CORNER_LABELS: [&str; 4] = ["M0_Q0", "M0_Q100", "M100_Q0", "M100_Q100"];

#[derive(Deserialize)]
struct RawSkin {
    #[serde(default)]
    boost: f64,
    #[serde(rename = "stageCount")]
    stage_count: Option<usize>,
    corners: BTreeMap<String, RawCorner>,
}

#[derive(Deserialize)]
struct RawCorner {
    stages: Vec<RawStage>,
}

#[derive(Deserialize)]
struct RawStage {
    a1: f64,
    r: f64,
    #[serde(default)]
    val1: f64,
    #[serde(default)]
    val2: f64,
    #[serde(default)]
    val3: f64,
    #[serde(default = "default_flag")]
    flag: f64,
}

fn default_flag() -> f64 {
    1.0
}

fn raw_stage_to_corner_row(s: &RawStage) -> [f64; 5] {
    let r = s.r.min(0.999999);
    let a2 = r * r;
    let (b0, b1, b2) = if s.flag < 0.5 {
        (1.0, 0.0, 0.0)
    } else {
        (1.0 + s.val1, s.a1 + s.val2, a2 - s.val3)
    };
    [b0, b1, b2, s.a1, a2]
}

fn build_corner(stages: &[RawStage]) -> Option<CornerData> {
    if stages.len() != NUM_STAGES {
        return None;
    }
    let mut corner = [[0.0f64; 5]; NUM_STAGES];
    for (i, s) in stages.iter().enumerate() {
        corner[i] = raw_stage_to_corner_row(s);
    }
    Some(corner)
}

fn trench_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("trench-core has a parent dir")
        .to_path_buf()
}

fn dry_path() -> PathBuf {
    PathBuf::from(r"C:/Users/hooki/trenchwork_clean/ref/bypassed-pinknoise.wav")
}

fn read_wav_mono(path: &Path) -> Vec<f32> {
    let mut reader = hound::WavReader::open(path)
        .unwrap_or_else(|e| panic!("opening {}: {e}", path.display()));
    let spec = reader.spec();
    let channels = spec.channels as usize;
    assert_eq!(
        spec.sample_format,
        hound::SampleFormat::Float,
        "expected float WAV at {}",
        path.display()
    );
    assert_eq!(spec.bits_per_sample, 32, "expected f32 WAV at {}", path.display());
    let samples: Vec<f32> = reader
        .samples::<f32>()
        .map(|s| s.expect("wav sample read"))
        .collect();
    if channels == 1 {
        samples
    } else {
        samples.into_iter().step_by(channels).collect()
    }
}

fn db(x: f64) -> f64 {
    if x <= 0.0 {
        -300.0
    } else {
        20.0 * x.log10()
    }
}

fn rms(v: &[f32]) -> f64 {
    if v.is_empty() {
        return 0.0;
    }
    let s: f64 = v.iter().map(|&x| (x as f64) * (x as f64)).sum();
    (s / v.len() as f64).sqrt()
}

/// Best-fit-gain null. Returns `(gain, rel_db)` where
/// `rel_db = 20*log10(rms(residual) / rms(ref))`.
fn null_db(pred: &[f32], reference: &[f32]) -> (f64, f64) {
    let n = pred.len().min(reference.len());
    let mut num = 0.0f64;
    let mut den = 0.0f64;
    for i in 0..n {
        num += (pred[i] as f64) * (reference[i] as f64);
        den += (pred[i] as f64) * (pred[i] as f64);
    }
    if den <= 0.0 {
        return (0.0, 0.0);
    }
    let gain = num / den;
    let mut residual = Vec::with_capacity(n);
    for i in 0..n {
        residual.push(((reference[i] as f64) - gain * (pred[i] as f64)) as f32);
    }
    let res_rms = rms(&residual);
    let ref_rms = rms(&reference[..n]);
    let rel = db(res_rms) - db(ref_rms);
    (gain, rel)
}

/// Render one corner of a raw skin through the rust pipeline.
///
/// Uses Cascade (f64 internally, f32 cast at output) followed by f64
/// AGC and f64 boost. The python reference is f64 throughout. The only
/// precision-loss point is the cascade's f32 output cast, which is then
/// promoted back to f64 before AGC. AGC is path-dependent, so even a
/// 1-ULP f32 difference at the cascade output can branch the gain
/// trajectory at corners where the cascade peaks above ~1.0.
fn render_corner(corner: &CornerData, boost: f64, dry: &[f32]) -> Vec<f32> {
    let mut cascade = Cascade::new();
    // Snap coefficients to target in 1 sample, then reset delay state
    // so the run starts with zero w1/w2 (matching python sosfilt's
    // zero-init).
    cascade.set_targets(corner, 1);
    let _ = cascade.tick(0.0);
    cascade.reset();
    // Re-apply set_targets so the rust ramp deltas are zero across the
    // run; coefs already equal target so deltas = 0.
    cascade.set_targets(corner, BLOCK_SIZE);
    cascade.set_boost(1.0, BLOCK_SIZE);

    let mut agc_gain = 1.0f64;
    let mut out = Vec::with_capacity(dry.len());
    for &x in dry {
        let cascaded = cascade.tick(x) as f64;
        let post_agc = agc_step_f64(cascaded, &mut agc_gain);
        out.push((post_agc * boost) as f32);
    }
    out
}

#[test]
fn rust_cascade_nulls_against_python_canonical_refs() {
    let root = trench_root();
    let skins_dir = root.join("datasets/p2k_skins");
    let ref_dir = root.join("ref/canonical");
    let dry = dry_path();

    if !skins_dir.exists() || !ref_dir.exists() || !dry.exists() {
        eprintln!(
            "canonical_parity: skipping (skins_dir={}, ref_dir={}, dry={})",
            skins_dir.exists(),
            ref_dir.exists(),
            dry.exists()
        );
        return;
    }

    let dry_samples = read_wav_mono(&dry);
    eprintln!("dry: {} samples", dry_samples.len());

    let mut entries: Vec<PathBuf> = fs::read_dir(&skins_dir)
        .expect("read skins dir")
        .filter_map(|e| e.ok().map(|e| e.path()))
        .filter(|p| p.extension().map(|s| s == "json").unwrap_or(false))
        .collect();
    entries.sort();
    eprintln!("raw skins: {}", entries.len());

    let mut failures: Vec<(String, String, f64)> = Vec::new();
    let mut worst_pass: f64 = -300.0;
    let mut covered = 0usize;
    let mut skipped_stage_count = 0usize;

    for path in &entries {
        let stem = path.file_stem().unwrap().to_string_lossy().to_string();
        let json = fs::read_to_string(path).expect("read skin");
        let skin: RawSkin = match serde_json::from_str(&json) {
            Ok(s) => s,
            Err(e) => {
                eprintln!("  {stem}: parse error: {e}");
                failures.push((stem, "parse".to_string(), 0.0));
                continue;
            }
        };

        let boost = if skin.boost == 0.0 { 1.0 } else { skin.boost };
        let stage_count = skin.stage_count.unwrap_or_else(|| {
            skin.corners
                .values()
                .next()
                .map(|c| c.stages.len())
                .unwrap_or(0)
        });
        if stage_count != NUM_STAGES {
            skipped_stage_count += 1;
            continue;
        }

        for label in CORNER_LABELS.iter() {
            let corner = match skin.corners.get(*label) {
                Some(c) => c,
                None => continue,
            };
            let corner_data = match build_corner(&corner.stages) {
                Some(cd) => cd,
                None => continue,
            };
            let ref_wav = ref_dir.join(format!("{stem}_{label}.wav"));
            if !ref_wav.exists() {
                continue;
            }
            let reference = read_wav_mono(&ref_wav);
            let pred = render_corner(&corner_data, boost, &dry_samples);
            let (gain, rel) = null_db(&pred, &reference);
            covered += 1;
            if rel > FAIL_THRESHOLD_DB {
                eprintln!(
                    " ! {stem:24} {label:>10}  gain={gain:.4}  rel_null={rel:+8.1} dB"
                );
                failures.push((stem.clone(), label.to_string(), rel));
            } else {
                if rel > worst_pass {
                    worst_pass = rel;
                }
                eprintln!(
                    "   {stem:24} {label:>10}  gain={gain:.4}  rel_null={rel:+8.1} dB"
                );
            }
        }
    }

    eprintln!();
    eprintln!(
        "rust parity: {covered} corners checked, {skipped_stage_count} skin(s) skipped (non-6-stage), \
         worst pass {worst_pass:+.1} dB, threshold {FAIL_THRESHOLD_DB:.0} dB"
    );
    if !failures.is_empty() {
        for (name, label, rel) in &failures {
            eprintln!("  FAIL: {name} {label} {rel:+.1} dB");
        }
        panic!(
            "{} corner(s) exceeded {} dB threshold",
            failures.len(),
            FAIL_THRESHOLD_DB
        );
    }
}
