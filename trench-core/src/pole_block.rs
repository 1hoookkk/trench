use std::error::Error;
use std::fmt;

use crate::cartridge::{CornerData, NUM_COEFFS, NUM_STAGES, PASSTHROUGH_COEFFS};

const CONJUGATE_ABS_TOLERANCE: f64 = 1.0e-9;
const CONJUGATE_REL_TOLERANCE: f64 = 1.0e-9;

/// One complex pole coordinate in the z-plane.
///
/// ROM-native pole data can be carried in cartesian form directly, avoiding any
/// lossy detour through authored Hz/Q abstractions before the runtime turns the
/// stage back into a biquad row for the cascade.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ComplexPole {
    pub re: f64,
    pub im: f64,
}

impl ComplexPole {
    pub const fn new(re: f64, im: f64) -> Self {
        Self { re, im }
    }

    pub fn from_polar(radius: f64, angle_radians: f64) -> Self {
        Self {
            re: radius * angle_radians.cos(),
            im: radius * angle_radians.sin(),
        }
    }

    pub const fn conjugate(self) -> Self {
        Self {
            re: self.re,
            im: -self.im,
        }
    }

    pub fn radius_squared(self) -> f64 {
        self.re.mul_add(self.re, self.im * self.im)
    }
}

/// One 2-pole resonant stage for the main cascade.
///
/// Numerator policy is the Morpheus-style pole-only audition path used by the
/// current cube tooling: `b = [gain, 0, 0]`, denominator driven by the complex
/// pole pair. If no explicit conjugate is supplied, the stage assumes the ROM
/// payload describes one half of a conjugate pair and mirrors it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PoleBlock {
    pub pole: ComplexPole,
    pub conjugate: Option<ComplexPole>,
    pub gain: f64,
}

impl PoleBlock {
    pub fn new(pole: ComplexPole) -> Self {
        Self {
            pole,
            conjugate: None,
            gain: 1.0,
        }
    }

    pub fn from_cartesian(re: f64, im: f64) -> Self {
        Self::new(ComplexPole::new(re, im))
    }

    pub fn from_pair(pole: ComplexPole, conjugate: ComplexPole) -> Self {
        Self {
            pole,
            conjugate: Some(conjugate),
            gain: 1.0,
        }
    }

    pub fn with_gain(mut self, gain: f64) -> Self {
        self.gain = gain;
        self
    }

    /// Compile this pole pair into direct-form kernel coefficients
    /// `[c0, c1, c2, c3, c4]` for the frozen DF2T cascade.
    pub fn try_coeffs(self) -> Result<[f64; NUM_COEFFS], PoleBlockError> {
        validate_finite("pole.re", self.pole.re)?;
        validate_finite("pole.im", self.pole.im)?;
        validate_finite("gain", self.gain)?;

        let conjugate = self.conjugate.unwrap_or_else(|| self.pole.conjugate());
        validate_finite("conjugate.re", conjugate.re)?;
        validate_finite("conjugate.im", conjugate.im)?;

        if !approx_eq(self.pole.re, conjugate.re) || !approx_eq(self.pole.im, -conjugate.im) {
            return Err(PoleBlockError::NotConjugate {
                pole: self.pole,
                conjugate,
            });
        }

        let radius_sq = self.pole.radius_squared();
        if radius_sq >= 1.0 {
            return Err(PoleBlockError::UnstablePole { radius_sq });
        }

        // For poles p1 and p2, denominator is:
        // 1 - (p1 + p2) z^-1 + (p1 * p2) z^-2
        let a1 = -(self.pole.re + conjugate.re);
        let a2 = self.pole.re * conjugate.re - self.pole.im * conjugate.im;

        validate_finite("a1", a1)?;
        validate_finite("a2", a2)?;

        Ok([self.gain, 0.0, 0.0, a1, a2])
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum PoleBlockError {
    NonFinite {
        field: &'static str,
        value: f64,
    },
    NotConjugate {
        pole: ComplexPole,
        conjugate: ComplexPole,
    },
    UnstablePole {
        radius_sq: f64,
    },
}

impl fmt::Display for PoleBlockError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NonFinite { field, value } => {
                write!(f, "{field} must be finite, got {value}")
            }
            Self::NotConjugate { pole, conjugate } => write!(
                f,
                "pole pair must be conjugates, got ({:.6}, {:.6}) and ({:.6}, {:.6})",
                pole.re, pole.im, conjugate.re, conjugate.im
            ),
            Self::UnstablePole { radius_sq } => {
                write!(f, "pole radius^2 must be < 1.0, got {radius_sq}")
            }
        }
    }
}

impl Error for PoleBlockError {}

#[derive(Debug, Clone, PartialEq)]
pub enum PoleBankError {
    TooManyStages {
        count: usize,
        max: usize,
    },
    Stage {
        index: usize,
        source: PoleBlockError,
    },
}

impl fmt::Display for PoleBankError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::TooManyStages { count, max } => {
                write!(
                    f,
                    "pole bank has {count} stages, max active stages is {max}"
                )
            }
            Self::Stage { index, source } => {
                write!(f, "pole bank stage {index} invalid: {source}")
            }
        }
    }
}

impl Error for PoleBankError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            Self::Stage { source, .. } => Some(source),
            Self::TooManyStages { .. } => None,
        }
    }
}

/// Compile a native pole-stage bank into one active `CornerData` payload for the
/// main cascade, padding unused stages with passthrough.
pub fn compile_pole_bank(stages: &[PoleBlock]) -> Result<CornerData, PoleBankError> {
    if stages.len() > NUM_STAGES {
        return Err(PoleBankError::TooManyStages {
            count: stages.len(),
            max: NUM_STAGES,
        });
    }

    let mut corner = [PASSTHROUGH_COEFFS; NUM_STAGES];
    for (index, stage) in stages.iter().copied().enumerate() {
        corner[index] = stage
            .try_coeffs()
            .map_err(|source| PoleBankError::Stage { index, source })?;
    }
    Ok(corner)
}

fn validate_finite(field: &'static str, value: f64) -> Result<(), PoleBlockError> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(PoleBlockError::NonFinite { field, value })
    }
}

fn approx_eq(a: f64, b: f64) -> bool {
    let delta = (a - b).abs();
    let scale = a.abs().max(b.abs()).max(1.0);
    delta <= CONJUGATE_ABS_TOLERANCE.max(CONJUGATE_REL_TOLERANCE * scale)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn single_cartesian_pole_compiles_to_resonator_kernel() {
        let radius = 0.92;
        let angle = 0.37;
        let pole = ComplexPole::from_polar(radius, angle);
        let coeffs = PoleBlock::new(pole).with_gain(0.75).try_coeffs().unwrap();

        let expected_a1 = -2.0 * radius * angle.cos();
        let expected_a2 = radius * radius;

        assert!((coeffs[0] - 0.75).abs() < 1.0e-12);
        assert_eq!(coeffs[1], 0.0);
        assert_eq!(coeffs[2], 0.0);
        assert!((coeffs[3] - expected_a1).abs() < 1.0e-12);
        assert!((coeffs[4] - expected_a2).abs() < 1.0e-12);
    }

    #[test]
    fn explicit_conjugate_pair_matches_single_pole_path() {
        let pole = ComplexPole::new(0.4, 0.2);
        let auto = PoleBlock::new(pole).try_coeffs().unwrap();
        let explicit = PoleBlock::from_pair(pole, pole.conjugate())
            .try_coeffs()
            .unwrap();

        assert_eq!(auto, explicit);
    }

    #[test]
    fn non_conjugate_pair_is_rejected() {
        let err = PoleBlock::from_pair(ComplexPole::new(0.4, 0.2), ComplexPole::new(0.41, -0.2))
            .try_coeffs()
            .unwrap_err();

        assert!(matches!(err, PoleBlockError::NotConjugate { .. }));
    }

    #[test]
    fn unstable_pole_is_rejected() {
        let err = PoleBlock::from_cartesian(1.0, 0.0)
            .try_coeffs()
            .unwrap_err();

        assert!(matches!(err, PoleBlockError::UnstablePole { .. }));
    }

    #[test]
    fn compile_pole_bank_pads_passthrough_stages() {
        let corner = compile_pole_bank(&[
            PoleBlock::from_cartesian(0.5, 0.1),
            PoleBlock::from_cartesian(0.2, 0.3).with_gain(0.8),
        ])
        .unwrap();

        assert_ne!(corner[0], PASSTHROUGH_COEFFS);
        assert_ne!(corner[1], PASSTHROUGH_COEFFS);
        for coeffs in corner.iter().skip(2) {
            assert_eq!(*coeffs, PASSTHROUGH_COEFFS);
        }
    }

    #[test]
    fn compile_pole_bank_rejects_stage_overflow() {
        let stages = vec![PoleBlock::from_cartesian(0.1, 0.1); NUM_STAGES + 1];
        let err = compile_pole_bank(&stages).unwrap_err();

        assert!(matches!(err, PoleBankError::TooManyStages { .. }));
    }
}
