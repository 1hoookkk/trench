use serde::Deserialize;
use serde_json::Value;

fn default_sample_coordinate() -> f64 {
    0.0
}

#[derive(Debug, Clone, Copy, Deserialize, PartialEq)]
pub struct HeritageSamplePoint {
    #[serde(default = "default_sample_coordinate")]
    pub morph: f64,
    #[serde(default = "default_sample_coordinate")]
    pub q: f64,
}

impl Default for HeritageSamplePoint {
    fn default() -> Self {
        Self { morph: 0.0, q: 0.0 }
    }
}

#[derive(Debug, Clone, Deserialize, PartialEq)]
pub struct HeritageDesignerSection {
    pub index: usize,
    #[serde(rename = "type")]
    pub type_id: i32,
    pub low_freq: i32,
    pub low_gain: i32,
    pub high_freq: i32,
    pub high_gain: i32,
}

#[derive(Debug, Clone, Deserialize, PartialEq)]
pub struct HeritageDesignerTemplate {
    pub name: String,
    #[serde(default)]
    pub frequency: Option<f64>,
    #[serde(default)]
    pub gain: Option<f64>,
    pub sections: Vec<HeritageDesignerSection>,
}

fn default_native_poles_sr() -> f64 {
    44100.0
}

fn default_native_poles_boost() -> f64 {
    1.0
}

#[derive(Debug, Clone, Copy, Deserialize, PartialEq)]
pub struct NativePoleStage {
    pub freq_hz: f64,
    pub radius: f64,
}

#[derive(Debug, Clone, Deserialize, PartialEq)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum HeritageObject {
    MorphDesignerRef {
        template_name: String,
        #[serde(default)]
        sample: HeritageSamplePoint,
    },
    MorphDesignerInline {
        template: HeritageDesignerTemplate,
        #[serde(default)]
        sample: HeritageSamplePoint,
    },
    PeakShelfRef {
        path: String,
        #[serde(default)]
        sample: HeritageSamplePoint,
    },
    PeakShelfInline {
        document: Value,
        #[serde(default)]
        sample: HeritageSamplePoint,
    },
    NativePoles {
        #[serde(default = "default_native_poles_sr")]
        sr: f64,
        #[serde(default = "default_native_poles_boost")]
        boost: f64,
        stages: Vec<NativePoleStage>,
    },
}

impl HeritageObject {
    pub fn validate(&self) -> Result<(), String> {
        match self {
            Self::MorphDesignerRef { template_name, .. } => {
                if template_name.trim().is_empty() {
                    return Err("morph_designer_ref.template_name must be non-empty".to_string());
                }
                Ok(())
            }
            Self::MorphDesignerInline { template, .. } => {
                if template.sections.is_empty() {
                    return Err(
                        "morph_designer_inline.template.sections must not be empty".to_string()
                    );
                }
                Ok(())
            }
            Self::PeakShelfRef { path, .. } => {
                if path.trim().is_empty() {
                    return Err("peak_shelf_ref.path must be non-empty".to_string());
                }
                Ok(())
            }
            Self::PeakShelfInline { document, .. } => {
                reject_compiled_runtime_fields(document, "peak_shelf_inline.document")
            }
            Self::NativePoles { sr, stages, .. } => {
                if !sr.is_finite() || *sr <= 0.0 {
                    return Err(format!("native_poles.sr must be positive, got {sr}"));
                }
                if stages.is_empty() {
                    return Err("native_poles.stages must not be empty".to_string());
                }
                for (index, stage) in stages.iter().enumerate() {
                    if !stage.freq_hz.is_finite() || !stage.radius.is_finite() {
                        return Err(format!(
                            "native_poles.stages[{index}] freq_hz/radius must be finite"
                        ));
                    }
                }
                Ok(())
            }
        }
    }

    pub fn sample(&self) -> HeritageSamplePoint {
        match self {
            Self::MorphDesignerRef { sample, .. }
            | Self::MorphDesignerInline { sample, .. }
            | Self::PeakShelfRef { sample, .. }
            | Self::PeakShelfInline { sample, .. } => *sample,
            Self::NativePoles { .. } => HeritageSamplePoint::default(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn native_poles_parses_and_validates() {
        let json = r#"{
            "kind": "native_poles",
            "sr": 44100.0,
            "boost": 1.25,
            "stages": [
                {"freq_hz": 260.0, "radius": 0.995},
                {"freq_hz": 820.0, "radius": 0.993}
            ]
        }"#;
        let obj: HeritageObject = serde_json::from_str(json).unwrap();
        match &obj {
            HeritageObject::NativePoles { sr, boost, stages } => {
                assert_eq!(*sr, 44100.0);
                assert_eq!(*boost, 1.25);
                assert_eq!(stages.len(), 2);
                assert_eq!(stages[0].freq_hz, 260.0);
            }
            other => panic!("wrong variant: {other:?}"),
        }
        assert!(obj.validate().is_ok());
        assert_eq!(obj.sample(), HeritageSamplePoint::default());
    }

    #[test]
    fn native_poles_defaults_when_fields_omitted() {
        let json = r#"{
            "kind": "native_poles",
            "stages": [{"freq_hz": 440.0, "radius": 0.99}]
        }"#;
        let obj: HeritageObject = serde_json::from_str(json).unwrap();
        let HeritageObject::NativePoles { sr, boost, .. } = &obj else {
            panic!("wrong variant");
        };
        assert_eq!(*sr, 44100.0);
        assert_eq!(*boost, 1.0);
        assert!(obj.validate().is_ok());
    }

    #[test]
    fn native_poles_rejects_empty_stages() {
        let json = r#"{"kind": "native_poles", "stages": []}"#;
        let obj: HeritageObject = serde_json::from_str(json).unwrap();
        assert!(obj.validate().is_err());
    }

    #[test]
    fn native_poles_rejects_non_positive_sr() {
        let json = r#"{
            "kind": "native_poles",
            "sr": 0.0,
            "stages": [{"freq_hz": 440.0, "radius": 0.99}]
        }"#;
        let obj: HeritageObject = serde_json::from_str(json).unwrap();
        assert!(obj.validate().is_err());
    }
}

fn reject_compiled_runtime_fields(value: &Value, path: &str) -> Result<(), String> {
    match value {
        Value::Object(map) => {
            if let Some(format) = map.get("format").and_then(Value::as_str) {
                if format == "compiled-v1" {
                    return Err(format!("{path} must not embed compiled-v1 payloads"));
                }
            }
            if map.get("keyframes").is_some() {
                return Err(format!("{path} must not embed compiled keyframes"));
            }
            if map.get("representation").is_some() {
                return Err(format!(
                    "{path} must not embed runtime representation metadata"
                ));
            }
            for (key, nested) in map {
                if matches!(key.as_str(), "c0" | "c1" | "c2" | "c3" | "c4") {
                    return Err(format!("{path} must not embed runtime coefficients"));
                }
                reject_compiled_runtime_fields(nested, &format!("{path}.{key}"))?;
            }
            Ok(())
        }
        Value::Array(items) => {
            for (index, nested) in items.iter().enumerate() {
                reject_compiled_runtime_fields(nested, &format!("{path}[{index}]"))?;
            }
            Ok(())
        }
        _ => Ok(()),
    }
}
