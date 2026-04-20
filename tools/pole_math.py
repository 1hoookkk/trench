"""Pole → kernel-form biquad SOS row.

Canonical location for the Z-plane pole → coefficient conversion used by the
offline compile path. Callers: `tools/compile_cube.py` (native_poles corner
kind), `tools/cube_authoring/morpheus_to_compiled.py` (cube staging from
Morpheus vault), `tools/cube_authoring/edge_gate_native.py` (pole-space
edge gate).

Numerator is passthrough (b0=1, b1=0, b2=0); the pole pair yields the
denominator: a1 = -2r cos θ, a2 = r². θ = 2π · freq_hz / sr. Kernel-form
5-tuple matches the frozen DF2T cascade in `trench-core/src/cascade.rs`.
"""
from __future__ import annotations

import math


def pole_to_kernel_stage(
    freq_hz: float, radius: float, sr: float
) -> tuple[float, float, float, float, float]:
    """Complex conjugate pole pair → kernel-form (c0, c1, c2, c3, c4) SOS row."""
    theta = 2.0 * math.pi * (freq_hz / sr)
    a1 = -2.0 * radius * math.cos(theta)
    a2 = radius * radius
    return (1.0, 0.0, 0.0, a1, a2)
