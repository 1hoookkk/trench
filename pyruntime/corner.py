"""Named corners, CornerArray, and bilinear interpolation.

Transcribed from runtime/src/corner.rs.
"""
from __future__ import annotations
from enum import IntEnum
from dataclasses import dataclass
from pyruntime.stage_params import StageParams
from pyruntime.encode import EncodedCoeffs, raw_to_encoded
from pyruntime.constants import NUM_BODY_STAGES


class CornerName(IntEnum):
    A = 0  # M0_Q0
    B = 1  # M0_Q100
    C = 2  # M100_Q0
    D = 3  # M100_Q100

    def json_key(self) -> str:
        return ["M0_Q0", "M0_Q100", "M100_Q0", "M100_Q100"][self.value]

    def morph(self) -> float:
        return 0.0 if self in (CornerName.A, CornerName.B) else 1.0

    def q(self) -> float:
        return 0.0 if self in (CornerName.A, CornerName.C) else 1.0


@dataclass
class CornerState:
    stages: list  # list[StageParams], length NUM_BODY_STAGES
    boost: float = 4.0
    _pre_encoded: list | None = None  # Heritage bypass — authoritative EncodedCoeffs

    def encode(self) -> list:
        """Encode all stages to kernel-form, padded to NUM_BODY_STAGES."""
        if self._pre_encoded is not None:
            enc = list(self._pre_encoded)
        else:
            enc = [raw_to_encoded(s) for s in self.stages]
        # Pad to NUM_BODY_STAGES with passthrough
        passthrough = EncodedCoeffs(c0=1.0, c1=0.0, c2=0.0, c3=0.0, c4=0.0)
        while len(enc) < NUM_BODY_STAGES:
            enc.append(passthrough)
        return enc


class CornerArray:
    """The complete 4-corner surface for a body."""

    def __init__(self, a: CornerState, b: CornerState, c: CornerState, d: CornerState):
        self._corners = [a, b, c, d]

    def corner(self, name: CornerName) -> CornerState:
        return self._corners[name.value]

    def names(self) -> list[CornerName]:
        return [CornerName.A, CornerName.B, CornerName.C, CornerName.D]

    def interpolate(self, morph: float, q: float) -> list:
        """Bilinear interpolation at (morph, q). Q first, then morph.

        Returns list of EncodedCoeffs for all stages.
        """
        a_enc = self._corners[0].encode()
        b_enc = self._corners[1].encode()
        c_enc = self._corners[2].encode()
        d_enc = self._corners[3].encode()

        out = []
        for i in range(NUM_BODY_STAGES):
            # Q first, then morph (frozen invariant)
            q_m0 = _lerp_enc(a_enc[i], b_enc[i], q)
            q_m100 = _lerp_enc(c_enc[i], d_enc[i], q)
            out.append(_lerp_enc(q_m0, q_m100, morph))
        return out


def _lerp_enc(a: EncodedCoeffs, b: EncodedCoeffs, t: float) -> EncodedCoeffs:
    return EncodedCoeffs(
        c0=a.c0 + (b.c0 - a.c0) * t,
        c1=a.c1 + (b.c1 - a.c1) * t,
        c2=a.c2 + (b.c2 - a.c2) * t,
        c3=a.c3 + (b.c3 - a.c3) * t,
        c4=a.c4 + (b.c4 - a.c4) * t,
    )
