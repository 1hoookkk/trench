use serde::Deserialize;

/// Number of active stages in the cascade.
pub const NUM_STAGES: usize = 6;
/// Number of coefficients per stage (c0..c4).
pub const NUM_COEFFS: usize = 5;
/// Number of interpolation corners.
pub const NUM_CORNERS: usize = 4;

/// Passthrough stage coefficients: y = x (identity).
pub const PASSTHROUGH_COEFFS: [f64; NUM_COEFFS] = [1.0, 0.0, 0.0, 0.0, 0.0];

/// One corner of the interpolation surface: 6 stages × 5 coefficients.
pub type CornerData = [[f64; NUM_COEFFS]; NUM_STAGES];

// ── Array format (trench-core native): corners.M0_Q0 = [[c0,c1,c2,c3,c4], ...] ──

#[derive(Deserialize)]
struct ArrayCartridgeJson {
    version: String,
    name: String,
    corners: ArrayCornersJson,
}

#[derive(Deserialize)]
#[allow(non_snake_case)]
struct ArrayCornersJson {
    M0_Q0: CornerData,
    M100_Q0: CornerData,
    M0_Q100: CornerData,
    M100_Q100: CornerData,
}

// ── Keyframe format (forge output): keyframes[].stages = [{c0,c1,c2,c3,c4}, ...] ──

#[derive(Deserialize)]
struct KeyframeCartridgeJson {
    format: String,
    name: String,
    keyframes: Vec<KeyframeJson>,
}

#[derive(Deserialize)]
struct KeyframeJson {
    label: String,
    #[serde(default = "default_boost")]
    boost: f64,
    stages: Vec<StageCoeffsJson>,
}

fn default_boost() -> f64 {
    1.0
}

#[derive(Deserialize)]
struct StageCoeffsJson {
    c0: f64,
    c1: f64,
    c2: f64,
    c3: f64,
    c4: f64,
}

/// A loaded cartridge ready for runtime interpolation.
#[derive(Clone, Debug)]
pub struct Cartridge {
    pub name: String,
    /// Indexed as [corner_index][stage][coeff].
    /// Corner order: M0_Q0, M100_Q0, M0_Q100, M100_Q100.
    pub corners: [CornerData; NUM_CORNERS],
    /// Per-corner post-cascade gain (linear). Same corner order.
    pub boosts: [f64; NUM_CORNERS],
}

impl Cartridge {
    /// Load a cartridge from JSON. Accepts both formats:
    /// - Array format: `{"version":"compiled-v1", "corners":{"M0_Q0":[[...],...],...}}`
    /// - Keyframe format: `{"format":"compiled-v1", "keyframes":[{"label":"M0_Q0","stages":[{"c0":...},...]},...]}`
    pub fn from_json(json: &str) -> Result<Self, String> {
        let probe: serde_json::Value =
            serde_json::from_str(json).map_err(|e| format!("JSON parse error: {e}"))?;

        if probe.get("keyframes").is_some() {
            Self::from_keyframe_json(json)
        } else if probe.get("corners").is_some() {
            Self::from_array_json(json)
        } else {
            Err("JSON must have 'corners' or 'keyframes' key".to_string())
        }
    }

    fn from_array_json(json: &str) -> Result<Self, String> {
        let raw: ArrayCartridgeJson =
            serde_json::from_str(json).map_err(|e| format!("array format parse error: {e}"))?;

        if raw.version != "compiled-v1" {
            return Err(format!(
                "unsupported version: '{}', expected 'compiled-v1'",
                raw.version
            ));
        }

        Ok(Self {
            name: raw.name,
            corners: [
                raw.corners.M0_Q0,
                raw.corners.M100_Q0,
                raw.corners.M0_Q100,
                raw.corners.M100_Q100,
            ],
            boosts: [1.0; NUM_CORNERS],
        })
    }

    fn from_keyframe_json(json: &str) -> Result<Self, String> {
        let raw: KeyframeCartridgeJson =
            serde_json::from_str(json).map_err(|e| format!("keyframe format parse error: {e}"))?;

        if raw.format != "compiled-v1" {
            return Err(format!(
                "unsupported format: '{}', expected 'compiled-v1'",
                raw.format
            ));
        }

        let find_corner = |label: &str| -> Result<(CornerData, f64), String> {
            let kf = raw
                .keyframes
                .iter()
                .find(|kf| kf.label == label)
                .ok_or_else(|| format!("missing keyframe '{label}'"))?;

            let stage_count = kf.stages.len();
            let padded_stage_count = NUM_STAGES * 2;
            if stage_count != NUM_STAGES && stage_count != padded_stage_count {
                return Err(format!(
                    "keyframe '{}' has {} stages, expected {} or {} (6 active + 6 passthrough)",
                    label,
                    stage_count,
                    NUM_STAGES,
                    padded_stage_count,
                ));
            }

            let mut corner = [[0.0; NUM_COEFFS]; NUM_STAGES];
            for (i, s) in kf.stages.iter().take(NUM_STAGES).enumerate() {
                corner[i] = [s.c0, s.c1, s.c2, s.c3, s.c4];
            }

            // Accept padded forge output (12 stages) as long as the tail is passthrough.
            if stage_count == padded_stage_count {
                for (i, s) in kf.stages.iter().skip(NUM_STAGES).enumerate() {
                    let coeffs = [s.c0, s.c1, s.c2, s.c3, s.c4];
                    let mut ok = true;
                    for j in 0..NUM_COEFFS {
                        if (coeffs[j] - PASSTHROUGH_COEFFS[j]).abs() > 1.0e-10 {
                            ok = false;
                            break;
                        }
                    }
                    if !ok {
                        return Err(format!(
                            "keyframe '{}' stage {} is not passthrough (expected c0=1,c1-4=0)",
                            label,
                            NUM_STAGES + i + 1
                        ));
                    }
                }
            }
            Ok((corner, kf.boost))
        };

        let (c0, b0) = find_corner("M0_Q0")?;
        let (c1, b1) = find_corner("M100_Q0")?;
        let (c2, b2) = find_corner("M0_Q100")?;
        let (c3, b3) = find_corner("M100_Q100")?;

        Ok(Self {
            name: raw.name,
            corners: [c0, c1, c2, c3],
            boosts: [b0, b1, b2, b3],
        })
    }

    /// Build the hardcoded Talking Hedz cartridge from the baked
    /// consts in `emu_params`. One `String` allocation for the name at
    /// plugin init time — zero allocation, zero I/O, zero locks on
    /// the audio thread. See `tools/bake_hedz_const.py` and
    /// `emu_params.rs` for the provenance chain.
    pub fn hedz_rom() -> Self {
        Self {
            name: crate::emu_params::HEDZ_NAME.to_string(),
            corners: crate::emu_params::HEDZ_CORNERS,
            boosts: crate::emu_params::HEDZ_BOOSTS,
        }
    }

    /// Interpolate post-cascade boost (linear gain) at the given morph/q position.
    /// Q-axis first, then morph-axis (same order as coefficient interpolation).
    pub fn interpolate_boost(&self, morph: f64, q: f64) -> f64 {
        let q_m0 = self.boosts[0] + (self.boosts[2] - self.boosts[0]) * q;
        let q_m1 = self.boosts[1] + (self.boosts[3] - self.boosts[1]) * q;
        q_m0 + (q_m1 - q_m0) * morph
    }

    /// Interpolate coefficients for all 6 stages at the given morph/q position.
    /// Both morph and q are in [0.0, 1.0].
    /// Order: Q-axis first, then morph-axis (per spec).
    pub fn interpolate(&self, morph: f64, q: f64) -> CornerData {
        let mut result = [[0.0; NUM_COEFFS]; NUM_STAGES];

        let m0_q0 = &self.corners[0];
        let m100_q0 = &self.corners[1];
        let m0_q100 = &self.corners[2];
        let m100_q100 = &self.corners[3];

        for stage in 0..NUM_STAGES {
            for c in 0..NUM_COEFFS {
                // Q-axis first
                let q_m0 = m0_q0[stage][c] + (m0_q100[stage][c] - m0_q0[stage][c]) * q;
                let q_m1 = m100_q0[stage][c] + (m100_q100[stage][c] - m100_q0[stage][c]) * q;
                // Then morph-axis
                result[stage][c] = q_m0 + (q_m1 - q_m0) * morph;
            }
        }

        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn test_array_json() -> &'static str {
        r#"{
            "version": "compiled-v1",
            "name": "Test Body",
            "corners": {
                "M0_Q0":    [[1,0,0,0,0],[1,0,0,0,0],[1,0,0,0,0],[1,0,0,0,0],[1,0,0,0,0],[1,0,0,0,0]],
                "M100_Q0":  [[2,0,0,0,0],[2,0,0,0,0],[2,0,0,0,0],[2,0,0,0,0],[2,0,0,0,0],[2,0,0,0,0]],
                "M0_Q100":  [[3,0,0,0,0],[3,0,0,0,0],[3,0,0,0,0],[3,0,0,0,0],[3,0,0,0,0],[3,0,0,0,0]],
                "M100_Q100":[[4,0,0,0,0],[4,0,0,0,0],[4,0,0,0,0],[4,0,0,0,0],[4,0,0,0,0],[4,0,0,0,0]]
            }
        }"#
    }

    fn test_keyframe_json() -> &'static str {
        r#"{
            "format": "compiled-v1",
            "name": "Forge Body",
            "keyframes": [
                {"label": "M0_Q0", "boost": 4.0, "stages": [
                    {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},{"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},
                    {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},{"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},
                    {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},{"c0":1,"c1":0,"c2":0,"c3":0,"c4":0}
                ]},
                {"label": "M100_Q0", "boost": 4.0, "stages": [
                    {"c0":5,"c1":0,"c2":0,"c3":0,"c4":0},{"c0":5,"c1":0,"c2":0,"c3":0,"c4":0},
                    {"c0":5,"c1":0,"c2":0,"c3":0,"c4":0},{"c0":5,"c1":0,"c2":0,"c3":0,"c4":0},
                    {"c0":5,"c1":0,"c2":0,"c3":0,"c4":0},{"c0":5,"c1":0,"c2":0,"c3":0,"c4":0}
                ]},
                {"label": "M0_Q100", "boost": 4.0, "stages": [
                    {"c0":3,"c1":0,"c2":0,"c3":0,"c4":0},{"c0":3,"c1":0,"c2":0,"c3":0,"c4":0},
                    {"c0":3,"c1":0,"c2":0,"c3":0,"c4":0},{"c0":3,"c1":0,"c2":0,"c3":0,"c4":0},
                    {"c0":3,"c1":0,"c2":0,"c3":0,"c4":0},{"c0":3,"c1":0,"c2":0,"c3":0,"c4":0}
                ]},
                {"label": "M100_Q100", "boost": 4.0, "stages": [
                    {"c0":7,"c1":0,"c2":0,"c3":0,"c4":0},{"c0":7,"c1":0,"c2":0,"c3":0,"c4":0},
                    {"c0":7,"c1":0,"c2":0,"c3":0,"c4":0},{"c0":7,"c1":0,"c2":0,"c3":0,"c4":0},
                    {"c0":7,"c1":0,"c2":0,"c3":0,"c4":0},{"c0":7,"c1":0,"c2":0,"c3":0,"c4":0}
                ]}
            ]
        }"#
    }

    #[test]
    fn loads_array_format() {
        let cart = Cartridge::from_json(test_array_json()).unwrap();
        assert_eq!(cart.name, "Test Body");
        assert_eq!(cart.corners[0][0][0], 1.0);
    }

    #[test]
    fn loads_keyframe_format() {
        let cart = Cartridge::from_json(test_keyframe_json()).unwrap();
        assert_eq!(cart.name, "Forge Body");
        assert_eq!(cart.corners[0][0][0], 1.0);
        assert_eq!(cart.corners[1][0][0], 5.0);
        assert_eq!(cart.corners[0][5][0], 1.0);
    }

    #[test]
    fn keyframe_interpolation() {
        let cart = Cartridge::from_json(test_keyframe_json()).unwrap();
        // M0_Q0 c0=1, M100_Q0 c0=5, M0_Q100 c0=3, M100_Q100 c0=7
        // At midpoint: Q first → q_m0 = 1+(3-1)*0.5=2, q_m1 = 5+(7-5)*0.5=6
        // Then morph: 2+(6-2)*0.5=4
        let c = cart.interpolate(0.5, 0.5);
        assert!((c[0][0] - 4.0).abs() < 1e-10);
    }

    #[test]
    fn interpolate_corners() {
        let cart = Cartridge::from_json(test_array_json()).unwrap();

        let c = cart.interpolate(0.0, 0.0);
        assert!((c[0][0] - 1.0).abs() < 1e-10);

        let c = cart.interpolate(1.0, 1.0);
        assert!((c[0][0] - 4.0).abs() < 1e-10);

        let c = cart.interpolate(0.5, 0.5);
        assert!((c[0][0] - 2.5).abs() < 1e-10);
    }

    #[test]
    fn rejects_wrong_version() {
        let json = r#"{"version":"v2","name":"X","corners":{"M0_Q0":[[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0]],"M100_Q0":[[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0]],"M0_Q100":[[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0]],"M100_Q100":[[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0]]}}"#;
        assert!(Cartridge::from_json(json).is_err());
    }

    #[test]
    fn rejects_wrong_keyframe_format() {
        let json = r#"{"format":"v2","name":"X","keyframes":[]}"#;
        assert!(Cartridge::from_json(json).is_err());
    }

    #[test]
    fn rejects_wrong_stage_count() {
        let json = r#"{
            "format": "compiled-v1",
            "name": "Bad",
            "keyframes": [
                {"label": "M0_Q0", "stages": [
                    {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},
                    {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},
                    {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0}
                ]},
                {"label": "M100_Q0", "stages": []},
                {"label": "M0_Q100", "stages": []},
                {"label": "M100_Q100", "stages": []}
            ]
        }"#;
        assert!(Cartridge::from_json(json).is_err());
    }
}
