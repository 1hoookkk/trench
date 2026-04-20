use serde::Deserialize;
use serde_json::Value;

use crate::heritage_objects::HeritageObject;

pub const CUBE_CORNER_LABELS: [&str; 8] = [
    "c000", "c100", "c010", "c110", "c001", "c101", "c011", "c111",
];

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CubeAuthoringControlMode {
    ModernLiveXyz,
    LegacyLatchYz,
}

impl CubeAuthoringControlMode {
    pub fn parse(raw: &str) -> Result<Self, String> {
        match raw {
            "modern_live_xyz" => Ok(Self::ModernLiveXyz),
            "legacy_latch_yz" => Ok(Self::LegacyLatchYz),
            other => Err(format!("unsupported cube control mode: '{other}'")),
        }
    }

    pub fn as_str(&self) -> &'static str {
        match self {
            Self::ModernLiveXyz => "modern_live_xyz",
            Self::LegacyLatchYz => "legacy_latch_yz",
        }
    }
}

#[derive(Debug, Clone, Deserialize)]
struct CubeAuthoringJson {
    schema: String,
    id: String,
    name: String,
    provenance: Value,
    exactness: String,
    control_mode_default: String,
    axes: Value,
    legacy_behavior: Value,
    corners: std::collections::BTreeMap<String, HeritageObject>,
}

#[derive(Debug, Clone)]
pub struct CubeAuthoringDocument {
    pub id: String,
    pub name: String,
    pub provenance: Value,
    pub exactness: String,
    pub control_mode_default: CubeAuthoringControlMode,
    pub axes: Value,
    pub legacy_behavior: Value,
    pub corners: [HeritageObject; 8],
}

impl CubeAuthoringDocument {
    pub fn from_json(json: &str) -> Result<Self, String> {
        let raw: CubeAuthoringJson =
            serde_json::from_str(json).map_err(|e| format!("cube authoring parse error: {e}"))?;
        if raw.schema != "trench.authoring_path.cube.v1" {
            return Err(format!(
                "unsupported cube authoring schema: '{}'",
                raw.schema
            ));
        }
        let control_mode_default = CubeAuthoringControlMode::parse(&raw.control_mode_default)?;
        let corners = Self::resolve_corners(&raw.corners)?;
        Ok(Self {
            id: raw.id,
            name: raw.name,
            provenance: raw.provenance,
            exactness: raw.exactness,
            control_mode_default,
            axes: raw.axes,
            legacy_behavior: raw.legacy_behavior,
            corners,
        })
    }

    fn resolve_corners(
        corners: &std::collections::BTreeMap<String, HeritageObject>,
    ) -> Result<[HeritageObject; 8], String> {
        let mut resolved = Vec::with_capacity(8);
        for label in CUBE_CORNER_LABELS {
            let corner = corners
                .get(label)
                .ok_or_else(|| format!("cube authoring missing corner '{label}'"))?
                .clone();
            corner.validate()?;
            resolved.push(corner);
        }
        resolved
            .try_into()
            .map_err(|_| "cube authoring corner resolution failed".to_string())
    }

    pub fn normalized_control_mode(&self) -> &'static str {
        self.control_mode_default.as_str()
    }
}
