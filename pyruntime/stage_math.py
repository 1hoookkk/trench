"""Stage parameter construction from frequency/radius/zero specifications.

Transcribed from runtime/src/stage_math.rs.
"""
import math
from pyruntime.stage_params import StageParams
from pyruntime.constants import SR, TWO_PI


def resonator(freq_hz: float, radius: float, val1: float) -> StageParams:
    """Allpole resonator. val2=0, val3=0."""
    theta = TWO_PI * freq_hz / SR
    return StageParams(
        a1=-2.0 * radius * math.cos(theta),
        r=radius,
        val1=val1,
        val2=0.0,
        val3=0.0,
    )


def resonator_with_zero(
    freq_hz: float,
    radius: float,
    val1: float,
    zero_freq_hz: float,
    zero_radius: float,
) -> StageParams:
    """Resonator with explicit zero placement."""
    theta = TWO_PI * freq_hz / SR
    phi = TWO_PI * zero_freq_hz / SR

    a1 = -2.0 * radius * math.cos(theta)
    b1_target = -2.0 * zero_radius * math.cos(phi)
    b2_target = zero_radius * zero_radius

    val2 = b1_target - a1
    val3 = radius * radius - b2_target

    return StageParams(a1=a1, r=radius, val1=val1, val2=val2, val3=val3)


def zero_forced(freq_hz: float, radius: float, val1: float) -> StageParams:
    """Zero on unit circle at pole frequency. 30 dB deep null."""
    return resonator_with_zero(freq_hz, radius, val1, freq_hz, 1.0)


def zero_forced_offset(
    freq_hz: float, radius: float, val1: float, offset_semitones: float
) -> StageParams:
    """Zero on unit circle, offset from pole by semitones."""
    zero_freq = freq_hz * 2.0 ** (offset_semitones / 12.0)
    zero_freq = max(20.0, min(SR * 0.48, zero_freq))
    return resonator_with_zero(freq_hz, radius, val1, zero_freq, 1.0)
