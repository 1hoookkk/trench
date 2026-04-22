use serde::Deserialize;

/// Optional drive-stage config (preceding cascade).
/// Public control is `slam` in [0, 1]. `input_gain_dB` remains as a
/// compatibility/authored field for older cartridges.
#[derive(Debug, Clone, Deserialize)]
pub struct DriveBlock {
    #[serde(default)]
    pub slam: Option<f32>,
    #[serde(rename = "input_gain_dB", default)]
    pub input_gain_db: f32,
    #[serde(default = "default_drive_model")]
    pub model: String,
}

fn default_drive_model() -> String {
    "desk_slam_v1".to_string()
}

impl Default for DriveBlock {
    fn default() -> Self {
        Self {
            slam: None,
            input_gain_db: 0.0,
            model: default_drive_model(),
        }
    }
}

impl DriveBlock {
    pub fn effective_slam(&self) -> f32 {
        self.slam
            .unwrap_or(self.input_gain_db / 60.0)
            .clamp(0.0, 1.0)
    }

    pub fn effective_input_gain_db(&self) -> f32 {
        self.effective_slam() * 60.0
    }
}

/// One segment of a TRENCH function-generator envelope.
#[derive(Debug, Clone, Deserialize)]
pub struct FnSegment {
    pub level: f32,
    pub time_ms: f32,
    pub shape: String,
    #[serde(default)]
    pub jump: Option<i32>,
}

/// Optional `mod_fn` block: a function-generator envelope + sync flags.
/// `key-sync` / `tempo-sync` are serialized as ints (0/1) in the E-mu style;
/// use [`ModFnBlock::key_sync`] / [`ModFnBlock::tempo_sync`] for booleans.
#[derive(Debug, Clone, Deserialize)]
pub struct ModFnBlock {
    pub segments: Vec<FnSegment>,
    #[serde(rename = "key-sync", default)]
    pub key_sync_int: i32,
    #[serde(rename = "tempo-sync", default)]
    pub tempo_sync_int: i32,
}

impl ModFnBlock {
    pub fn key_sync(&self) -> bool {
        self.key_sync_int != 0
    }
    pub fn tempo_sync(&self) -> bool {
        self.tempo_sync_int != 0
    }
}

/// Six-coefficient regression law (ITD or ILD). Feature order is the canonical
/// `[sin(az), sin(2*az), sin(3*az), sin(4*az), sin(5*az), sin(az)*abs(el)/30]`
/// from `docs/archive/qsound_spatial.md` §"Canonical Recon Model".
pub type LawCoeffs6 = [f32; 6];

/// Twelve-coefficient per-band regression law. Feature order is the canonical
/// `[1, log2(dist/0.25), log2(dist/0.25)^2, sin(az), cos(az), sin(2*az),
///   cos(2*az), sin(3*az), cos(3*az), el/30, sin(az)*el/30, cos(az)*el/30]`
/// from `docs/archive/qsound_spatial.md` §"Canonical Recon Model".
pub type BandLawCoeffs12 = [f32; 12];

/// Per-channel low/mid/high band coefficients.
#[derive(Debug, Clone, Deserialize)]
pub struct BandChannelCoeffs {
    pub low: BandLawCoeffs12,
    pub mid: BandLawCoeffs12,
    pub high: BandLawCoeffs12,
}

/// Left + right band coefficient matrices.
#[derive(Debug, Clone, Deserialize)]
pub struct BandCoeffs {
    pub l: BandChannelCoeffs,
    pub r: BandChannelCoeffs,
}

/// Typed QSound-style spatial profile. Units match the regression laws in
/// `docs/archive/qsound_spatial.md`: `azimuth` / `elevation` in radians,
/// `distance` in metres. All fields required when the block is present;
/// the block itself remains optional at the cartridge level.
#[derive(Debug, Clone, Deserialize)]
pub struct SpatialProfile {
    pub azimuth: f32,
    pub distance: f32,
    pub elevation: f32,
    pub itd_coeffs: LawCoeffs6,
    pub ild_coeffs: LawCoeffs6,
    pub band_coeffs: BandCoeffs,
}

/// Number of active stages in the cascade.
pub const NUM_STAGES: usize = 6;
/// Number of coefficients per stage (c0..c4).
pub const NUM_COEFFS: usize = 5;
/// Number of corners in the compiled-v1 bilinear surface.
pub const NUM_CORNERS: usize = 4;

/// Passthrough stage coefficients: y = x (identity).
pub const PASSTHROUGH_COEFFS: [f64; NUM_COEFFS] = [1.0, 0.0, 0.0, 0.0, 0.0];

/// One surface corner: 6 stages × 5 kernel-form coefficients.
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
///
/// One 4-corner bilinear filter surface: `M0_Q0`, `M100_Q0`, `M0_Q100`,
/// `M100_Q100`. Interpolation order is Q first, then morph.
#[derive(Clone, Debug)]
pub struct Cartridge {
    pub name: String,
    corners: [CornerData; NUM_CORNERS],
    boosts: [f64; NUM_CORNERS],
    /// Optional drive stage preceding the cascade. Absent block → inert default.
    pub drive: DriveBlock,
    /// Optional spatial profile (QSound-style). Absent → no spatial stage.
    pub spatial_profile: Option<SpatialProfile>,
    /// Optional function-generator envelope driving MORPH. Absent → no fn mod.
    pub mod_fn: Option<ModFnBlock>,
}

/// Lightweight probe for the three optional top-level blocks. Tolerant of
/// unknown fields and missing keys — any block absent stays `None`.
#[derive(Deserialize, Default)]
struct OptionalBlocksProbe {
    #[serde(default)]
    drive: Option<DriveBlock>,
    #[serde(default)]
    spatial_profile: Option<SpatialProfile>,
    #[serde(default)]
    mod_fn: Option<ModFnBlock>,
}

impl Cartridge {
    /// Load a cartridge from JSON. Accepts:
    /// - Array format: `{"version":"compiled-v1", "corners":{"M0_Q0":[[...],...],...}}`
    /// - Keyframe format: `{"format":"compiled-v1", "keyframes":[{"label":"M0_Q0","stages":[{"c0":...},...]},...]}`
    pub fn from_json(json: &str) -> Result<Self, String> {
        let probe: serde_json::Value =
            serde_json::from_str(json).map_err(|e| format!("JSON parse error: {e}"))?;

        let mut cart = if probe.get("keyframes").is_some() {
            Self::from_keyframe_json(json)?
        } else if probe.get("corners").is_some() {
            Self::from_array_json(json)?
        } else {
            return Err("JSON must have 'corners' or 'keyframes' key".to_string());
        };

        let blocks: OptionalBlocksProbe =
            serde_json::from_str(json).map_err(|e| format!("optional blocks parse error: {e}"))?;
        if let Some(d) = blocks.drive {
            cart.drive = d;
        }
        cart.spatial_profile = blocks.spatial_profile;
        cart.mod_fn = blocks.mod_fn;

        Ok(cart)
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
            drive: DriveBlock::default(),
            spatial_profile: None,
            mod_fn: None,
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
            parse_compiled_corner(label, &kf.stages, kf.boost)
        };

        let (c0, b0) = find_corner("M0_Q0")?;
        let (c1, b1) = find_corner("M100_Q0")?;
        let (c2, b2) = find_corner("M0_Q100")?;
        let (c3, b3) = find_corner("M100_Q100")?;

        Ok(Self {
            name: raw.name,
            corners: [c0, c1, c2, c3],
            boosts: [b0, b1, b2, b3],
            drive: DriveBlock::default(),
            spatial_profile: None,
            mod_fn: None,
        })
    }

    /// Build the hardcoded Talking Hedz cartridge from the baked
    /// consts in `hedz_rom`. One `String` allocation for the name at
    /// plugin init time — zero allocation, zero I/O, zero locks on
    /// the audio thread. See `tools/bake_hedz_const.py` and
    /// `hedz_rom.rs` for the provenance chain.
    pub fn hedz_rom() -> Self {
        Self {
            name: crate::hedz_rom::HEDZ_NAME.to_string(),
            corners: crate::hedz_rom::HEDZ_CORNERS,
            boosts: crate::hedz_rom::HEDZ_BOOSTS,
            drive: DriveBlock::default(),
            spatial_profile: None,
            mod_fn: None,
        }
    }

    pub fn corners(&self) -> &[CornerData; NUM_CORNERS] {
        &self.corners
    }

    pub fn boosts(&self) -> &[f64; NUM_CORNERS] {
        &self.boosts
    }

    /// Bilinear boost interpolation at the given morph/q position.
    pub fn interpolate_boost(&self, morph: f64, q: f64) -> f64 {
        let q_m0 = self.boosts[0] + (self.boosts[2] - self.boosts[0]) * q;
        let q_m1 = self.boosts[1] + (self.boosts[3] - self.boosts[1]) * q;
        q_m0 + (q_m1 - q_m0) * morph
    }

    /// Bilinear interpolation: produces coefficients for all 6 stages at the
    /// given morph/q position. Both morph and q are in [0.0, 1.0]. Order: Q
    /// first, then morph (SPEC.md §1, compiled-v1 path).
    pub fn interpolate(&self, morph: f64, q: f64) -> CornerData {
        let mut result = [[0.0; NUM_COEFFS]; NUM_STAGES];
        let m0_q0 = &self.corners[0];
        let m100_q0 = &self.corners[1];
        let m0_q100 = &self.corners[2];
        let m100_q100 = &self.corners[3];

        for stage in 0..NUM_STAGES {
            for c in 0..NUM_COEFFS {
                let q_m0 = m0_q0[stage][c] + (m0_q100[stage][c] - m0_q0[stage][c]) * q;
                let q_m1 = m100_q0[stage][c] + (m100_q100[stage][c] - m100_q0[stage][c]) * q;
                result[stage][c] = q_m0 + (q_m1 - q_m0) * morph;
            }
        }

        result
    }
}

fn parse_compiled_corner(
    label: &str,
    stages: &[StageCoeffsJson],
    boost: f64,
) -> Result<(CornerData, f64), String> {
    let stage_count = stages.len();
    let padded_stage_count = NUM_STAGES * 2;
    if stage_count != NUM_STAGES && stage_count != padded_stage_count {
        return Err(format!(
            "corner '{}' has {} stages, expected {} or {} (6 active + 6 passthrough)",
            label, stage_count, NUM_STAGES, padded_stage_count,
        ));
    }

    let mut corner = [[0.0; NUM_COEFFS]; NUM_STAGES];
    for (i, s) in stages.iter().take(NUM_STAGES).enumerate() {
        corner[i] = [s.c0, s.c1, s.c2, s.c3, s.c4];
    }

    if stage_count == padded_stage_count {
        for (i, s) in stages.iter().skip(NUM_STAGES).enumerate() {
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
                    "corner '{}' stage {} is not passthrough (expected c0=1,c1-4=0)",
                    label,
                    NUM_STAGES + i + 1
                ));
            }
        }
    }

    Ok((corner, boost))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn test_array_json() -> &'static str {
        r#"{
            "version": "compiled-v1",
            "name": "Test Filter",
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
            "name": "Forge Filter",
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
        assert_eq!(cart.name, "Test Filter");
        assert_eq!(cart.corners()[0][0][0], 1.0);
    }

    #[test]
    fn loads_keyframe_format() {
        let cart = Cartridge::from_json(test_keyframe_json()).unwrap();
        assert_eq!(cart.name, "Forge Filter");
        assert_eq!(cart.corners()[0][0][0], 1.0);
        assert_eq!(cart.corners()[1][0][0], 5.0);
        assert_eq!(cart.corners()[0][5][0], 1.0);
    }

    #[test]
    fn keyframe_interpolation() {
        let cart = Cartridge::from_json(test_keyframe_json()).unwrap();
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
