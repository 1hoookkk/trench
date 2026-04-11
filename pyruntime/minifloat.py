"""Minifloat codec: 4-bit exponent (bias 15), 12-bit mantissa.

Transcribed from trench-core/src/minifloat.rs.
Decode: ldexp((mantissa | 0x1000) * scale, exp - 15).
"""
import math


def decode(raw: int) -> float:
    """Decode a u16 minifloat to f64."""
    if raw == 0:
        return 0.0
    exp = (raw >> 12) & 0xF
    mantissa = raw & 0xFFF
    return math.ldexp((mantissa | 0x1000), exp - 15)
