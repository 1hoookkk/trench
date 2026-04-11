"""Frequency response computation — numpy vectorized.

Evaluates the exact DF2T transfer function used by the Rust cascade:
    y  = c0*x + w1
    w1 = c1*x - c3*y + w2
    w2 = c2*x - c4*y

    H(z) = (c0 + c1*z^-1 + c2*z^-2) / (1 + c3*z^-1 + c4*z^-2)

Default sample rate is the authoring rate (39062.5 Hz).
All P2K bodies and forge output use this rate.
"""
import numpy as np
from pyruntime.encode import EncodedCoeffs
from pyruntime.constants import SR, FREQ_MIN, NUM_RESPONSE_POINTS, TWO_PI

SCOPE_FREQ_MAX = SR / 2.0 - 1.0


def freq_points(n: int = NUM_RESPONSE_POINTS, sr: float = SR) -> np.ndarray:
    """Log-spaced frequency points from 20 Hz to just under Nyquist."""
    f_max = sr / 2.0 - 1.0
    return np.logspace(np.log10(FREQ_MIN), np.log10(f_max), n)


def stage_response(enc: EncodedCoeffs, freqs: np.ndarray, sr: float = SR) -> np.ndarray:
    """Complex H(f) for one encoded stage.

    Direct DF2T transfer function — no coefficient translation.
    H(z) = (c0 + c1*z^-1 + c2*z^-2) / (1 + c3*z^-1 + c4*z^-2)
    """
    z_inv = np.exp(-1j * TWO_PI * freqs / sr)
    z_inv2 = z_inv * z_inv

    num = enc.c0 + enc.c1 * z_inv + enc.c2 * z_inv2
    den = 1.0 + enc.c3 * z_inv + enc.c4 * z_inv2

    return num / den


def cascade_response(
    stages: list, freqs: np.ndarray, sr: float = SR
) -> np.ndarray:
    """Complex cascade response: product of all stage responses."""
    h = np.ones(len(freqs), dtype=complex)
    for enc in stages:
        h *= stage_response(enc, freqs, sr)
    return h


def cascade_response_db(
    stages: list, freqs: np.ndarray, sr: float = SR
) -> np.ndarray:
    """Cascade magnitude in dB."""
    h = cascade_response(stages, freqs, sr)
    return 20.0 * np.log10(np.maximum(np.abs(h), 1e-20))
