//! Role distribution and interference-candidate diagnostic.
//!
//! Uses the newly-ported `Role::classify` and `StageFunction`/`InterferenceLink`
//! vocabulary to measure what the existing shipping candidates actually
//! contain at the coefficient level. Produces a JSON artifact under
//! `vault/_diagnostics/role_distribution.json` with, for every candidate:
//!
//! - extracted per-stage pole frequency and radius (from DF2T coefficients
//!   `c3 = -2 r cos θ`, `c4 = r²`)
//! - role classification per stage (Anchor / LowMid / Character / Air)
//! - interference-pair candidate count (stage pairs where both poles are
//!   high-Q and within 25% frequency spacing)
//!
//! The metric this is designed to expose: **how many of the shipping
//! candidates contain any cross-stage coordination structure at all.** A
//! candidate with zero interference-pair candidates is structurally flat
//! — the "basic shelf cascade" failure mode that the plots expose visually.

use std::fs;
use std::path::PathBuf;

use serde_json::json;
use trench_core::{Role, NUM_COEFFS, NUM_STAGES};

/// Extract (freq_hz, radius) for a stage from DF2T denominator coefficients.
///
/// DF2T denominator form is `1 + a1 z^-1 + a2 z^-2` with `a1 = -2 r cos(θ)`
/// and `a2 = r²`. The `c3 = a1`, `c4 = a2` mapping holds in trench-core's
/// kernel form (`y += -c3 y[-1] - c4 y[-2]`).
///
/// Returns `None` if the stage is a passthrough (all zero denom coefficients)
/// or if the coefficients don't yield a real pole.
fn extract_pole(stage: &[f64; NUM_COEFFS], sr: f64) -> Option<(f64, f64)> {
    let c3 = stage[3];
    let c4 = stage[4];
    if c4.abs() < 1e-12 && c3.abs() < 1e-12 {
        return None;
    }
    let r2 = c4;
    if r2 <= 0.0 || r2 >= 1.0 {
        return None;
    }
    let r = r2.sqrt();
    // a1 = -2 r cos(θ) → cos(θ) = -c3 / (2r)
    let cos_theta = -c3 / (2.0 * r);
    if !(-1.0..=1.0).contains(&cos_theta) {
        return None;
    }
    let theta = cos_theta.acos();
    let freq = theta * sr / (2.0 * std::f64::consts::PI);
    Some((freq, r))
}

/// Load the M100_Q100 corner's 6 active stages from a 12-stage forge body.
fn load_m100_q100_active_stages(path: &PathBuf) -> Vec<[f64; NUM_COEFFS]> {
    let json = fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("read {}: {e}", path.display()));
    let v: serde_json::Value = serde_json::from_str(&json)
        .unwrap_or_else(|e| panic!("parse {}: {e}", path.display()));
    let kfs = v["keyframes"].as_array().expect("keyframes");
    let kf = kfs
        .iter()
        .find(|k| k["label"] == "M100_Q100")
        .expect("missing M100_Q100 keyframe");
    let stages = kf["stages"].as_array().expect("stages");
    let mut out = Vec::with_capacity(NUM_STAGES);
    for s in stages.iter().take(NUM_STAGES) {
        out.push([
            s["c0"].as_f64().unwrap(),
            s["c1"].as_f64().unwrap(),
            s["c2"].as_f64().unwrap(),
            s["c3"].as_f64().unwrap(),
            s["c4"].as_f64().unwrap(),
        ]);
    }
    out
}

/// Candidate pairs: (i, j) where both poles are high-Q (r > 0.9) and
/// the pair's frequency spacing is within 25% of the lower freq.
fn count_interference_candidates(poles: &[(usize, f64, f64)]) -> (usize, Vec<(usize, usize, f64)>) {
    let mut pairs = Vec::new();
    for a in 0..poles.len() {
        for b in (a + 1)..poles.len() {
            let (_, fa, ra) = poles[a];
            let (_, fb, rb) = poles[b];
            if ra < 0.9 || rb < 0.9 {
                continue;
            }
            let lo = fa.min(fb);
            let hi = fa.max(fb);
            if lo <= 0.0 {
                continue;
            }
            let spacing_ratio = (hi - lo) / lo;
            if spacing_ratio < 0.25 {
                pairs.push((poles[a].0, poles[b].0, spacing_ratio));
            }
        }
    }
    (pairs.len(), pairs)
}

fn analyze_body(name: &str, path: &PathBuf, sr: f64) -> serde_json::Value {
    let stages = load_m100_q100_active_stages(path);
    let mut per_stage: Vec<serde_json::Value> = Vec::new();
    let mut poles: Vec<(usize, f64, f64)> = Vec::new();
    let mut role_counts = [0usize; 4];
    for (i, stage) in stages.iter().enumerate() {
        match extract_pole(stage, sr) {
            Some((freq, r)) => {
                let role = Role::classify(freq);
                match role {
                    Role::Anchor => role_counts[0] += 1,
                    Role::LowMid => role_counts[1] += 1,
                    Role::Character => role_counts[2] += 1,
                    Role::Air => role_counts[3] += 1,
                }
                per_stage.push(json!({
                    "stage": i,
                    "freq_hz": freq,
                    "pole_r": r,
                    "role": format!("{:?}", role),
                }));
                poles.push((i, freq, r));
            }
            None => {
                per_stage.push(json!({
                    "stage": i,
                    "freq_hz": null,
                    "pole_r": null,
                    "role": "None",
                }));
            }
        }
    }
    let (pair_count, pairs) = count_interference_candidates(&poles);
    let pair_details: Vec<serde_json::Value> = pairs
        .iter()
        .map(|(a, b, ratio)| json!({"a": a, "b": b, "spacing_ratio": ratio}))
        .collect();
    json!({
        "name": name,
        "source": path.display().to_string(),
        "stages": per_stage,
        "role_counts": {
            "Anchor": role_counts[0],
            "LowMid": role_counts[1],
            "Character": role_counts[2],
            "Air": role_counts[3],
        },
        "active_stages": poles.len(),
        "interference_candidates": pair_count,
        "interference_pairs": pair_details,
    })
}

#[test]
fn shipping_candidates_role_distribution() {
    let root = PathBuf::from("../vault/_shipping_finalists");
    let sr = 44100.0_f64;

    let mut bodies: Vec<(String, PathBuf)> = Vec::new();
    for ident in ["speaker_knockerz", "aluminum_siding", "small_talk", "cul_de_sac"] {
        let dir = root.join(ident);
        if !dir.exists() {
            continue;
        }
        let mut entries: Vec<_> = fs::read_dir(&dir)
            .unwrap()
            .filter_map(|e| e.ok())
            .filter(|e| {
                e.path()
                    .extension()
                    .map(|x| x == "json")
                    .unwrap_or(false)
            })
            .collect();
        entries.sort_by_key(|e| e.file_name());
        for e in entries {
            let p = e.path();
            let name = p.file_stem().unwrap().to_string_lossy().to_string();
            bodies.push((name, p));
        }
    }

    assert!(!bodies.is_empty(), "no shipping candidates found to diagnose");

    let mut report: Vec<serde_json::Value> = Vec::new();
    for (name, path) in &bodies {
        let body = analyze_body(name, path, sr);
        report.push(body);
    }

    let out = json!({
        "generated_by": "trench-core/tests/role_distribution_diagnostic.rs",
        "sr_hz": sr,
        "corner_probed": "M100_Q100",
        "interference_criteria": {
            "min_pole_r": 0.9,
            "max_spacing_ratio": 0.25,
        },
        "candidates": report,
    });

    let out_dir = PathBuf::from("../vault/_diagnostics");
    fs::create_dir_all(&out_dir).unwrap();
    let out_path = out_dir.join("role_distribution.json");
    fs::write(&out_path, serde_json::to_string_pretty(&out).unwrap()).unwrap();

    eprintln!("\n=== shipping candidate role distribution ===");
    for c in out["candidates"].as_array().unwrap() {
        eprintln!(
            "  {:<50}  active={}  roles={{A:{} L:{} C:{} Ai:{}}}  interference_pairs={}",
            c["name"].as_str().unwrap(),
            c["active_stages"].as_u64().unwrap(),
            c["role_counts"]["Anchor"].as_u64().unwrap(),
            c["role_counts"]["LowMid"].as_u64().unwrap(),
            c["role_counts"]["Character"].as_u64().unwrap(),
            c["role_counts"]["Air"].as_u64().unwrap(),
            c["interference_candidates"].as_u64().unwrap(),
        );
    }
    eprintln!("wrote {}", out_path.display());
}
