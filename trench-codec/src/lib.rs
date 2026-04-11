//! 16-bit minifloat codec for TRENCH coefficient storage.
//!
//! Format: 4-bit exponent (bias 15) | 12-bit mantissa
//! Decode: ldexp((mantissa | 0x1000) * scale, exp - 15)
//! where scale = 2^-13 (normalizes the 13-bit significand to [0.5, 1.0))
//!
//! Sentinels:
//!   0x0000 = 0.0
//!   0xDFFF = passthrough gain (~0.25)
//!   0xFFFF = maximum constant (~1.0)

pub const SENTINEL_ZERO: u16 = 0x0000;
pub const SENTINEL_PASSTHROUGH: u16 = 0xDFFF;
pub const SENTINEL_MAX: u16 = 0xFFFF;

/// Passthrough gain value (decoded from 0xDFFF).
pub const PASSTHROUGH_GAIN: f64 = 8191.0 / 32768.0;

/// Decode a 16-bit minifloat to f64.
pub fn decode(raw: u16) -> f64 {
    match raw {
        SENTINEL_ZERO => 0.0,
        _ => {
            let exp = ((raw >> 12) & 0xF) as i32;
            let mantissa = (raw & 0xFFF) as u64;
            let significand = mantissa | 0x1000; // implicit leading 1
                                                 // significand * 2^(exp - 28)
                                                 // = significand * 2^(-13) * 2^(exp - 15)
            (significand as f64) * f64::from_bits(((1023_i64 + (exp as i64) - 28) as u64) << 52)
        }
    }
}

/// Encode an f64 value to 16-bit minifloat.
/// Clamps to representable range. Values <= 0.0 become SENTINEL_ZERO.
pub fn encode(value: f64) -> u16 {
    if value <= 0.0 || !value.is_finite() {
        return SENTINEL_ZERO;
    }

    // v = significand * 2^(exp_field - 28)
    // significand in [4096, 8191]
    // We need to find exp_field and mantissa.
    let bits = value.to_bits();
    let ieee_exp = ((bits >> 52) & 0x7FF) as i32 - 1023; // unbiased IEEE exponent
    let ieee_frac = (bits & 0x000F_FFFF_FFFF_FFFF) | 0x0010_0000_0000_0000; // 53-bit significand

    // value = ieee_frac * 2^(ieee_exp - 52)
    // We want: value = (mantissa | 0x1000) * 2^(exp_field - 28)
    // So: mantissa | 0x1000 = value / 2^(exp_field - 28)
    //
    // The significand (mantissa | 0x1000) is 13 bits in [4096, 8191].
    // value = frac * 2^(ieee_exp) where frac in [1.0, 2.0)
    // We want: significand in [4096, 8191] = [2^12, 2^13 - 1]
    // significand = frac * 4096.0
    // exp_field = ieee_exp + 28 - 12 = ieee_exp + 16

    // Extract top 13 bits of IEEE significand (which is 53 bits with implicit 1)
    let significand_13 = (ieee_frac >> (52 - 12)) as u16; // top 13 bits
    let mantissa = significand_13 & 0xFFF; // remove implicit 1
    let exp_field = ieee_exp + 16;

    if exp_field < 0 {
        return SENTINEL_ZERO; // underflow
    }
    if exp_field > 15 {
        return SENTINEL_MAX; // overflow (clamp)
    }

    ((exp_field as u16) << 12) | mantissa
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sentinel_zero_roundtrips() {
        assert_eq!(decode(SENTINEL_ZERO), 0.0);
        assert_eq!(encode(0.0), SENTINEL_ZERO);
    }

    #[test]
    fn sentinel_max_decodes_near_one() {
        let val = decode(SENTINEL_MAX);
        assert!(
            (val - 1.0).abs() < 0.001,
            "0xFFFF should decode near 1.0, got {val}"
        );
    }

    #[test]
    fn passthrough_decodes_correctly() {
        let val = decode(SENTINEL_PASSTHROUGH);
        assert!(
            (val - PASSTHROUGH_GAIN).abs() < 1e-10,
            "0xDFFF should decode to passthrough gain, got {val}"
        );
    }

    #[test]
    fn roundtrip_typical_values() {
        let test_values = [0.001, 0.01, 0.1, 0.25, 0.5, 0.75, 0.99];
        for &v in &test_values {
            let encoded = encode(v);
            let decoded = decode(encoded);
            let rel_err = ((decoded - v) / v).abs();
            assert!(
                rel_err < 0.001,
                "roundtrip failed for {v}: encoded={encoded:#06X}, decoded={decoded}, rel_err={rel_err}"
            );
        }
    }

    #[test]
    fn encode_negative_returns_zero() {
        assert_eq!(encode(-1.0), SENTINEL_ZERO);
    }
}
