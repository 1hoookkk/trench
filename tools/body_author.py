"""Body author — 6 typed roles → compiled body.

The simplest authoring model: pick 6 roles with frequencies,
the compiler handles type/gain/radius/zeros.

Roles:
  anchor(freq)           — locked, spectral mass, type 1, high radius
  character(f_lo, f_hi)  — sweeps with morph, type 2 formant
  suppressor(freq)       — cancels a band, type 1, gain below unity
  diffuser(freq)         — broadband scatter, type 1, low radius
  boundary(freq)         — shelf/tilt edge, type 3

Usage:
    python tools/body_author.py
"""
from __future__ import annotations

import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np

from pyruntime.body import Body
from pyruntime.constants import SR, TWO_PI
from pyruntime.designer_compile import (
    DesignerCell, make_four_corner_template, compile_four_corner_to_body,
)
from pyruntime.freq_response import cascade_response_db, freq_points
from pyruntime.forge_generator import solve_c4_surface, apply_gain_budget

FREQS = freq_points()
VAULT_DIR = Path(__file__).parent.parent / "vault" / "_authored"

# Freq → packed integer (inverse of fw_freq_value for SR=44100)
# fw: freq_value = (220 * freq_param) >> 7 + 18
# So: freq_param = ((freq_value - 18) << 7) / 220
# And freq_value → Hz is via the semitone table, but simplified:
#   Hz ≈ 72 * 2^(freq_param / 18.3)  (empirical fit to the table)


def _hz_to_param(hz: float) -> int:
    """Convert Hz to 0-127 freq parameter."""
    if hz < 72:
        return 0
    param = 18.3 * math.log2(hz / 72.0)
    return max(0, min(127, int(round(param))))


def _hz_to_param_f(hz: float) -> float:
    """Float version for interpolation."""
    if hz < 72:
        return 0.0
    return max(0.0, min(127.0, 18.3 * math.log2(hz / 72.0)))


# =============================================================================
# Roles
# =============================================================================

@dataclass
class Stage:
    """Compiled stage ready for 4-corner template."""
    a: DesignerCell  # M0_Q0
    b: DesignerCell  # M0_Q100
    c: DesignerCell  # M100_Q0
    d: DesignerCell  # M100_Q100


def anchor(freq_hz: float) -> Stage:
    """Locked spectral mass. Type 1 gain=127 (max bandwidth resonator)."""
    p = _hz_to_param(freq_hz)
    return Stage(
        a=DesignerCell(type=1, freq=p, gain=127),
        b=DesignerCell(type=1, freq=p, gain=127),
        c=DesignerCell(type=1, freq=p, gain=127),
        d=DesignerCell(type=1, freq=p, gain=127),
    )


def character(f_lo: float, f_hi: float) -> Stage:
    """Sweeps with morph. Type 1 gain=127 (resonator). Heritage-proven gain level."""
    p_lo = _hz_to_param(f_lo)
    p_hi = _hz_to_param(f_hi)
    p_mid = (p_lo + p_hi) // 2
    p_q_lo = p_lo + (p_mid - p_lo) // 3
    p_q_hi = p_hi - (p_hi - p_mid) // 3
    return Stage(
        a=DesignerCell(type=1, freq=p_lo, gain=127),
        b=DesignerCell(type=1, freq=p_q_lo, gain=127),
        c=DesignerCell(type=1, freq=p_hi, gain=127),
        d=DesignerCell(type=1, freq=p_q_hi, gain=127),
    )


def suppressor(freq_hz: float) -> Stage:
    """Cancels a band. Type 1 gain=0 (minimum bandwidth — the anti-resonator)."""
    p = _hz_to_param(freq_hz)
    return Stage(
        a=DesignerCell(type=1, freq=p, gain=0),
        b=DesignerCell(type=1, freq=p, gain=0),
        c=DesignerCell(type=1, freq=p, gain=0),
        d=DesignerCell(type=1, freq=p, gain=0),
    )


def diffuser(freq_hz: float) -> Stage:
    """Broadband scatter. Type 2 (fixed radius formant)."""
    p = _hz_to_param(freq_hz)
    return Stage(
        a=DesignerCell(type=2, freq=p, gain=0),
        b=DesignerCell(type=2, freq=p, gain=0),
        c=DesignerCell(type=2, freq=p, gain=0),
        d=DesignerCell(type=2, freq=p, gain=0),
    )


def boundary(freq_hz: float) -> Stage:
    """Shelf/tilt at spectral edge. Type 3 shaped, gain=127."""
    p = _hz_to_param(freq_hz)
    return Stage(
        a=DesignerCell(type=3, freq=p, gain=127),
        b=DesignerCell(type=3, freq=p, gain=127),
        c=DesignerCell(type=3, freq=p, gain=127),
        d=DesignerCell(type=3, freq=p, gain=127),
    )


# =============================================================================
# Compiler: roles → body
# =============================================================================

def compile_roles(name: str, stages: list[Stage], target_db: float = 36.0) -> Body:
    """Compile 6 roles into a body through the heritage firmware path + c4 solver."""
    assert len(stages) <= 6, f"Max 6 stages, got {len(stages)}"

    corners = {"M0_Q0": [], "M0_Q100": [], "M100_Q0": [], "M100_Q100": []}
    for s in stages:
        corners["M0_Q0"].append(s.a)
        corners["M0_Q100"].append(s.b)
        corners["M100_Q0"].append(s.c)
        corners["M100_Q100"].append(s.d)

    # Pad with bypass
    while len(corners["M0_Q0"]) < 6:
        for k in corners:
            corners[k].append(DesignerCell(type=0, freq=0, gain=0))

    template = make_four_corner_template(name, corners)
    body = compile_four_corner_to_body(template)

    # Normalize cascade gain — without this, output is ~100 dB too quiet
    body = solve_c4_surface(body, target_db=target_db)
    body = apply_gain_budget(body)
    return body


def describe(body: Body) -> str:
    """Quick MCP-style characterization."""
    _VM = (FREQS >= 200) & (FREQS <= 5000)
    pts = [0.0, 0.5, 1.0]
    lines = []
    for m in pts:
        db = cascade_response_db(body.corners.interpolate(m, 0.5), FREQS, SR)
        pk = float(np.max(db))
        vl = float(np.min(db))
        dyn = pk - vl

        # talkingness
        vdb = db[_VM]
        mean = float(np.mean(db))
        talk_peaks = 0
        if len(vdb) > 2:
            for i in range(1, len(vdb) - 1):
                if vdb[i] > vdb[i-1] and vdb[i] > vdb[i+1] and vdb[i] > mean + 3:
                    talk_peaks += 1
        talk = min(1.0, talk_peaks / 4.0)

        # centroid
        lin = 10.0 ** (db / 20.0)
        d = float(np.sum(lin))
        cent = float(np.sum(FREQS * lin) / d) if d > 0 else 1000

        lines.append(f"  m={m:.1f}: peak={pk:>6.1f}  dyn={dyn:>5.1f}  talk={talk:.2f}  cent={cent:>6.0f} Hz")

    return "\n".join(lines)


# =============================================================================
# The 4 flagship bodies — authored as roles
# =============================================================================

def make_speaker_knockerz():
    return compile_roles("Speaker Knockerz", [
        anchor(45),                    # sub — the floor
        anchor(75),                    # sub harmonic
        character(200, 800),           # choke → rip
        character(400, 2000),          # rattle → cry
        character(800, 4500),          # stress → fracture
        suppressor(6000),              # HF ceiling
    ])


def make_aluminum_siding():
    return compile_roles("Aluminum Siding", [
        suppressor(1000),              # permanent midrange void
        boundary(250),                 # warmth floor
        character(5000, 8000),         # shimmer → glass
        character(7000, 11000),        # glass → tear
        character(9000, 14000),        # tear → shatter
        character(12000, 17000),       # near-Nyquist sizzle
    ])


def make_small_talk():
    return compile_roles("Small Talk Ah-Ee", [
        character(700, 350),           # F1: Ah → Ee (drops)
        character(1200, 2300),         # F2: Ah → Ee (rises)
        character(2500, 3000),         # F3: brightness
        character(3500, 3700),         # F4: nasality
        anchor(150),                   # warmth floor
        character(500, 1500),          # throat color
    ])


def make_cul_de_sac():
    return compile_roles("Cul-De-Sac", [
        anchor(100),                   # hum tether — never moves
        character(350, 1500),          # pipe → shard 1
        character(400, 3500),          # pipe → shard 2
        character(450, 6000),          # pipe → shard 3
        character(500, 10000),         # pipe → shard 4
        character(550, 15000),         # pipe → shard 5
    ])


ALL = {
    "speaker_knockerz": make_speaker_knockerz,
    "aluminum_siding": make_aluminum_siding,
    "small_talk": make_small_talk,
    "cul_de_sac": make_cul_de_sac,
}


def main():
    VAULT_DIR.mkdir(parents=True, exist_ok=True)

    for key, factory in ALL.items():
        body = factory()
        print(f"\n{'='*50}")
        print(f"  {body.name}")
        print(f"{'='*50}")
        print(describe(body))

        out = VAULT_DIR / f"{key}.json"
        with open(out, "w") as f:
            f.write(body.to_compiled_json(provenance="body-author"))
        print(f"  → {out}")


if __name__ == "__main__":
    main()
