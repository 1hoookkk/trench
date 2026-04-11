//! Diagonal-motor trajectory test.
//!
//! Loads `cartridges/motors/snap_rip.json` and binds it to the live shipping
//! body `vault/_shipping_finalists/speaker_knockerz/speaker_knockerz__cand_04_fracture_hiss.json`.
//!
//! Walks two trajectories through the (morph, q) surface to the same end
//! point and asserts that the diagonal path samples coefficient territory
//! that the pure-morph path at q=0 cannot reach. This is the empirical floor
//! for the diagonal-gesture research bet identified in
//! `docs/modulation_exploration.md` ("Diagonal Gestures").

use std::fs;

use trench_core::{Cartridge, CornerData, ForceCurve, Motor, NUM_COEFFS, NUM_STAGES};

fn load_motor(path: &str) -> Motor {
    let json = fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("read motor {path}: {e}"));
    serde_json::from_str(&json)
        .unwrap_or_else(|e| panic!("parse motor {path}: {e}"))
}

/// Load a shipping body in the 12-stage forge keyframe format and extract
/// the 6 active stages (stages 0..6; stages 6..12 are passthrough).
///
/// The shipping forge emits 12-stage cartridges (6 active + 6 passthrough)
/// while `Cartridge::from_json` requires exactly 6 stages. This helper bridges
/// the two formats inline so the test binds to a real shipping body.
fn load_cartridge_active_only(path: &str) -> Cartridge {
    let json = fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("read cartridge {path}: {e}"));
    let v: serde_json::Value = serde_json::from_str(&json)
        .unwrap_or_else(|e| panic!("parse cartridge {path}: {e}"));

    let name = v["name"]
        .as_str()
        .unwrap_or_else(|| panic!("cartridge {path} missing 'name' string"))
        .to_string();
    let kfs = v["keyframes"]
        .as_array()
        .unwrap_or_else(|| panic!("cartridge {path} missing 'keyframes' array"));

    let extract = |label: &str| -> CornerData {
        let kf = kfs
            .iter()
            .find(|k| k["label"] == label)
            .unwrap_or_else(|| panic!("missing keyframe '{label}' in {path}"));
        let stages = kf["stages"]
            .as_array()
            .unwrap_or_else(|| panic!("keyframe '{label}' missing 'stages'"));
        assert!(
            stages.len() >= NUM_STAGES,
            "keyframe '{label}' has {} stages, need at least {NUM_STAGES} active",
            stages.len()
        );
        let mut corner = [[0.0_f64; NUM_COEFFS]; NUM_STAGES];
        for i in 0..NUM_STAGES {
            corner[i][0] = stages[i]["c0"].as_f64().expect("c0");
            corner[i][1] = stages[i]["c1"].as_f64().expect("c1");
            corner[i][2] = stages[i]["c2"].as_f64().expect("c2");
            corner[i][3] = stages[i]["c3"].as_f64().expect("c3");
            corner[i][4] = stages[i]["c4"].as_f64().expect("c4");
        }
        corner
    };

    Cartridge {
        name,
        corners: [
            extract("M0_Q0"),
            extract("M100_Q0"),
            extract("M0_Q100"),
            extract("M100_Q100"),
        ],
        boosts: [1.0; 4],
    }
}

/// Max absolute coefficient delta across all 6 stages and 5 coefficients.
fn max_abs_coeff_delta(
    a: &[[f64; NUM_COEFFS]; NUM_STAGES],
    b: &[[f64; NUM_COEFFS]; NUM_STAGES],
) -> f64 {
    let mut m = 0.0_f64;
    for stage in 0..NUM_STAGES {
        for c in 0..NUM_COEFFS {
            let d = (a[stage][c] - b[stage][c]).abs();
            if d > m {
                m = d;
            }
        }
    }
    m
}

#[test]
fn snap_rip_motor_loads_with_diagonal_vector() {
    let motor = load_motor("../cartridges/motors/snap_rip.json");
    assert_eq!(motor.name, "Snap-Rip");
    assert!(matches!(motor.curve, ForceCurve::Impulse));
    // Diagonal: both vector components must be non-zero, and Q must travel
    // opposite to morph (the "snap" — open the cascade while collapsing Q).
    assert!(motor.vector_morph.abs() > f32::EPSILON);
    assert!(motor.vector_q.abs() > f32::EPSILON);
    assert!(motor.vector_morph * motor.vector_q < 0.0,
        "Snap-Rip must travel diagonally with morph and q in opposition");
}

#[test]
fn snap_rip_position_at_traces_anti_diagonal() {
    let motor = load_motor("../cartridges/motors/snap_rip.json");
    // home (0, 1) → fully engaged (1, 0). Midpoint should be (0.5, 0.5).
    let (m0, q0) = motor.position_at(0.0);
    let (m_mid, q_mid) = motor.position_at(0.5);
    let (m1, q1) = motor.position_at(1.0);
    assert!((m0 - 0.0).abs() < 1e-6 && (q0 - 1.0).abs() < 1e-6);
    assert!((m_mid - 0.5).abs() < 1e-6 && (q_mid - 0.5).abs() < 1e-6);
    assert!((m1 - 1.0).abs() < 1e-6 && (q1 - 0.0).abs() < 1e-6);
}

#[test]
fn diagonal_trajectory_reaches_unreachable_coefficients() {
    let motor = load_motor("../cartridges/motors/snap_rip.json");
    let body = load_cartridge_active_only(
        "../vault/_shipping_finalists/speaker_knockerz/speaker_knockerz__cand_04_fracture_hiss.json",
    );

    // The diagonal path passes through the *interior* of the (m, q) surface.
    // The pure-morph path at q=0 walks the bottom edge only.
    // At any interior t ∈ (0, 1), the diagonal samples a point the bottom
    // edge cannot reach, so the coefficient sets must differ.
    //
    // We assert this at five interior samples and require the max absolute
    // coefficient delta to exceed a non-trivial floor at every sample.
    let samples = [0.1_f32, 0.25, 0.5, 0.75, 0.9];

    for &t in &samples {
        let (m_diag, q_diag) = motor.position_at(t);
        let diag_coeffs = body.interpolate(m_diag as f64, q_diag as f64);

        // Pure-morph path at q=0 sampled at the same morph value.
        let pure_morph_coeffs = body.interpolate(m_diag as f64, 0.0);

        let delta = max_abs_coeff_delta(&diag_coeffs, &pure_morph_coeffs);
        assert!(
            delta > 1e-4,
            "diagonal trajectory at t={t} (m={m_diag}, q={q_diag}) collapsed onto pure-morph path \
             at q=0; max abs coeff delta was {delta:.6e}. \
             Body's Q axis is degenerate or the motor vector is wrong."
        );
    }

    // Stronger: the diagonal midpoint (0.5, 0.5) must also differ from the
    // pure-Q midpoint at (0, 0.5), proving the diagonal samples a region
    // unreachable by either single-axis sweep.
    let (m_diag, q_diag) = motor.position_at(0.5);
    let diag_mid = body.interpolate(m_diag as f64, q_diag as f64);
    let pure_q_mid = body.interpolate(0.0, q_diag as f64);
    let delta_q = max_abs_coeff_delta(&diag_mid, &pure_q_mid);
    assert!(
        delta_q > 1e-4,
        "diagonal midpoint collapsed onto pure-Q path at m=0; max abs coeff delta was {delta_q:.6e}"
    );
}

#[test]
fn diagonal_trajectory_visits_distinct_coefficient_sets() {
    // The diagonal trajectory at successive t values should produce distinct
    // coefficient sets. If the path stalls (all samples identical), the body
    // is degenerate or the motor has zero vector. Reject either.
    let motor = load_motor("../cartridges/motors/snap_rip.json");
    let body = load_cartridge_active_only(
        "../vault/_shipping_finalists/speaker_knockerz/speaker_knockerz__cand_04_fracture_hiss.json",
    );

    let mut prev: Option<[[f64; NUM_COEFFS]; NUM_STAGES]> = None;
    let mut total_motion = 0.0_f64;
    for i in 0..=10 {
        let t = i as f32 / 10.0;
        let (m, q) = motor.position_at(t);
        let coeffs = body.interpolate(m as f64, q as f64);
        if let Some(p) = prev {
            total_motion += max_abs_coeff_delta(&coeffs, &p);
        }
        prev = Some(coeffs);
    }
    assert!(
        total_motion > 0.01,
        "diagonal trajectory through Speaker Knockerz C4 produced negligible motion \
         (total max-delta sum {total_motion:.6e}); body or motor is degenerate"
    );
}
