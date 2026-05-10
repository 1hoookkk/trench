"""Canonical MorphLP/LPX lattice constants (proved RE contract)."""

from __future__ import annotations

MORPH_INDEX_SCALE = 16.0
MORPH_INDEX_MIN = 0
MORPH_INDEX_MAX = 15

# DAT_1806d73a0 / DAT_1806d73b0 (low 16-bit values consumed by compiler)
POLE_BASE_U16 = (0x1320, 0x1194, 0x0384, 0x00DC)   # 4896, 4500, 900, 220
POLE_SLOPE_U16 = (0x01BA, 0x01B8, 0x0195, 0x015E)  # 442, 440, 405, 350


def clamp_morph_index(idx: int) -> int:
    if idx < MORPH_INDEX_MIN:
        return MORPH_INDEX_MIN
    if idx > MORPH_INDEX_MAX:
        return MORPH_INDEX_MAX
    return idx


def morph_index(morph_norm: float, morph_offset: float = 0.0) -> int:
    """floor((morph_norm + morph_offset) * 16.0) clamped to [0,15]."""
    return clamp_morph_index(int((morph_norm + morph_offset) * MORPH_INDEX_SCALE))

