"""P2K corner splice: combine corners from two bodies into novel morph trajectories.

Transcribed from runtime/src/splice.rs.

Corner layout: A=M0_Q0, B=M0_Q100, C=M100_Q0, D=M100_Q100.
Morph axis (A/B → C/D) traverses between two filter voices.
Q axis (A/C → B/D) applies pressure.
"""
from __future__ import annotations

from enum import Enum

from pyruntime.corner import CornerArray, CornerName


class SpliceMode(Enum):
    REST_TO_REST = "RestToRest"
    REST_TO_MORPHED = "RestToMorphed"
    MORPHED_TO_MORPHED = "MorphedToMorphed"
    CROSS = "Cross"


def splice_corners(
    body_a: CornerArray, body_b: CornerArray, mode: SpliceMode
) -> CornerArray:
    """Splice two bodies into a novel morph trajectory.

    Takes corners from body_a and body_b according to the splice mode
    and assembles them into a new 4-corner surface.
    """
    if mode == SpliceMode.REST_TO_REST:
        # A's resting → B's resting
        a = body_a.corner(CornerName.A)
        b = body_a.corner(CornerName.B)
        c = body_b.corner(CornerName.A)
        d = body_b.corner(CornerName.B)
    elif mode == SpliceMode.REST_TO_MORPHED:
        # A's resting → B's morphed
        a = body_a.corner(CornerName.A)
        b = body_a.corner(CornerName.B)
        c = body_b.corner(CornerName.C)
        d = body_b.corner(CornerName.D)
    elif mode == SpliceMode.MORPHED_TO_MORPHED:
        # A's morphed → B's morphed
        a = body_a.corner(CornerName.C)
        b = body_a.corner(CornerName.D)
        c = body_b.corner(CornerName.C)
        d = body_b.corner(CornerName.D)
    elif mode == SpliceMode.CROSS:
        # Diagonal cross through combined space
        a = body_a.corner(CornerName.A)
        b = body_b.corner(CornerName.B)
        c = body_b.corner(CornerName.C)
        d = body_a.corner(CornerName.D)
    else:
        raise ValueError(f"Unknown splice mode: {mode}")

    return CornerArray(a=a, b=b, c=c, d=d)
