"""Kernel-form encode: StageParams → EncodedCoeffs.

Transcribed from trench-core/src/encode.rs.
Critical: casts to float32 before computing (matches Rust as f32).

Direct format: c0=b0, c1=b1, c2=b2 (numerator), c3=a1, c4=a2 (denominator).
Passthrough: [1, 0, 0, 0, 0].
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from pyruntime.stage_params import StageParams

GUARD = 1e-12


@dataclass(frozen=True, slots=True)
class EncodedCoeffs:
    c0: float
    c1: float
    c2: float
    c3: float
    c4: float


def raw_to_encoded(stage: StageParams, flag: float = 1.0) -> EncodedCoeffs:
    """Convert StageParams to DF2T coefficients [c0..c4].

    Output maps directly to the TRENCH cascade:
        y  = c0*x + w1
        w1 = c1*x - c3*y + w2
        w2 = c2*x - c4*y

    So: c0=b0, c1=b1, c2=b2 (numerator), c3=a1, c4=a2 (denominator).

    Casts to float32 first to match Rust's f32 intermediate precision.
    """
    # Cast to f32 to match Rust pipeline
    a1 = float(np.float32(stage.a1))
    r = float(np.float32(stage.r))
    val1 = float(np.float32(stage.val1))
    val2 = float(np.float32(stage.val2))
    val3 = float(np.float32(stage.val3))

    if flag >= 0.5:
        # Resonator path
        b0 = 1.0 + val1
        b1 = a1 + val2
        b2 = r * r - val3

        c0 = b0
        c1 = b1
        c2 = b2
        c3 = a1
        c4 = r * r

        return EncodedCoeffs(c0=c0, c1=c1, c2=c2, c3=c3, c4=c4)
    else:
        # Lowpass path
        a1_bq = a1 - 2.0
        a2_bq = 1.0 - r
        b0_bq = (1.0 + a1_bq + a2_bq) / 4.0

        c0 = b0_bq
        c1 = 2.0 * b0_bq
        c2 = b0_bq
        c3 = a1_bq
        c4 = a2_bq

        return EncodedCoeffs(c0=c0, c1=c1, c2=c2, c3=c3, c4=c4)
