"""MorphLP zero-placement table — sanitized constant module.

The 16-entry CPhantomMorphLP zero table provides per-corner numerator
minifloat words (val1, val2, val3) for vocal-character Type 3 stages.
Unlike LPX, MorphLP gives **distinct zero coefficients at M0 and M100**:

- M0   (morph = 0): deep LP zeros (strong numerator character)
- M100 (morph = 1): zeros approach the passthrough sentinel (near all-pole)

That M0 → M100 dissolution is the mathematical source of the "talking"
character in vocal Type 3 bodies. Without per-corner zero variation the
morph reduces to pole-only sweep (LPX behavior).

## Layout

Each entry is (M0_v1, M0_v2, M0_v3, M100_v1, M100_v2, M100_v3) — six u16
minifloat values. Zero symmetry constraint: `v1 == v3` per corner
(conjugate-pair structure → b0 == b2 in the resulting biquad numerator).

Index 0 has the deepest M0 zeros and the most-dissolved M100 zeros (max
vocal divergence). Indices grow toward shallower M0 and more-present
M100 zeros — at index 14/15 (saturation) the divergence is smallest.

## Trust tier

Sanitized constant. No Ghidra addresses, decompiler symbols, or vendor
class names appear in this module — those are documented separately in
the RE vault. The numbers themselves are observable behavioural data
(filter-frequency-response signatures) and consequently do not carry
proprietary nomenclature.

## Usage

    from tools.morphlp_zeros import morphlp_corner_words

    m0_words  = morphlp_corner_words(idx=0, corner="M0")    # (v1, v2, v3)
    m100_words = morphlp_corner_words(idx=0, corner="M100")
"""
from __future__ import annotations

# 16 entries × (M0_v1, M0_v2, M0_v3, M100_v1, M100_v2, M100_v3)
_TABLE: tuple[tuple[int, int, int, int, int, int], ...] = (
    (0x11ED, 0x94FC, 0x11ED, 0xDFFC, 0xFFFD, 0xDFFC),  #  0 — max divergence
    (0x11EE, 0x90FC, 0x11ED, 0xE2FC, 0xFDFD, 0xE2FC),  #  1
    (0x11ED, 0x8DFC, 0x11ED, 0xE5FC, 0xFC7D, 0xE5FC),  #  2
    (0x11EE, 0x8AFC, 0x11ED, 0xE8FC, 0xFAFD, 0xE8FC),  #  3
    (0x11ED, 0x87FC, 0x11ED, 0xEBFC, 0xF97D, 0xEBFC),  #  4
    (0x11ED, 0x83FC, 0x11ED, 0xEDFC, 0xF87D, 0xEDFC),  #  5
    (0x11EE, 0x80FC, 0x11ED, 0xF07D, 0xF77D, 0xF07D),  #  6
    (0x0FEE, 0x7DFC, 0x0FEF, 0xF17D, 0xF67D, 0xF17D),  #  7
    (0x0FEF, 0x7AFC, 0x0FEF, 0xF27D, 0xF47D, 0xF27D),  #  8
    (0x0FEE, 0x77FC, 0x0FEF, 0xF3FD, 0xF27D, 0xF3FD),  #  9
    (0x0FEE, 0x73FC, 0x0FEF, 0xF47D, 0xF17D, 0xF47D),  # 10
    (0x0DEE, 0x70FC, 0x0DEF, 0xF4FD, 0xF07D, 0xF4FD),  # 11
    (0x0DED, 0x6DFC, 0x0DEF, 0xF5FD, 0xEFFC, 0xF5FD),  # 12
    (0x0BEF, 0x6AFC, 0x0BEF, 0xF67D, 0xEEFC, 0xF67D),  # 13
    (0x0BEF, 0x62FC, 0x0BEF, 0xF6FD, 0xECFC, 0xF6FD),  # 14
    (0x0BEF, 0x62FC, 0x0BEF, 0xF6FD, 0xECFC, 0xF6FD),  # 15 — saturation
)

NUM_ENTRIES = len(_TABLE)
PASSTHROUGH_SENTINEL = 0xDFFF


def morphlp_corner_words(idx: int, corner: str) -> tuple[int, int, int]:
    """Return (val1, val2, val3) numerator minifloat words for one corner
    of MorphLP entry `idx`. `corner` is "M0" or "M100"."""
    if not 0 <= idx < NUM_ENTRIES:
        raise ValueError(f"MorphLP index {idx} out of range [0, {NUM_ENTRIES})")
    row = _TABLE[idx]
    if corner == "M0":
        return (row[0], row[1], row[2])
    if corner == "M100":
        return (row[3], row[4], row[5])
    raise ValueError(f"corner must be 'M0' or 'M100', got {corner!r}")


def assert_symmetry(tol: int = 2) -> None:
    """Sanity check: val1 ≈ val3 within each corner (conjugate-pair structure).
    The RE doc claims strict equality but the firmware ROM has small LSB
    perturbations: most entries 0 or ±1, entry 12 has ±2. Treated as
    quantization artifact: tolerance defaults to ±2."""
    for idx, row in enumerate(_TABLE):
        m0_v1, _, m0_v3, m100_v1, _, m100_v3 = row
        if abs(m0_v1 - m0_v3) > tol:
            raise AssertionError(f"entry {idx}: M0 v1=0x{m0_v1:04X} v3=0x{m0_v3:04X}")
        if abs(m100_v1 - m100_v3) > tol:
            raise AssertionError(f"entry {idx}: M100 v1=0x{m100_v1:04X} v3=0x{m100_v3:04X}")


if __name__ == "__main__":
    assert_symmetry()
    print(f"MorphLP table: {NUM_ENTRIES} entries, symmetry OK")
    for idx in range(NUM_ENTRIES):
        m0 = morphlp_corner_words(idx, "M0")
        m100 = morphlp_corner_words(idx, "M100")
        print(
            f"  [{idx:2d}] M0=(0x{m0[0]:04X}, 0x{m0[1]:04X}, 0x{m0[2]:04X})  "
            f"M100=(0x{m100[0]:04X}, 0x{m100[1]:04X}, 0x{m100[2]:04X})"
        )
