#!/usr/bin/env python3
"""forge_cli.py — Smart Forge end-to-end Q-law CLI + Wizard.

Default invocation (no args) opens THE OPERATOR — an interactive,
ASCII-driven cartridge press. Pick a recipe, pick an arc, pick a policy,
and the forge solves, gates, and stamps a cartridge.

Subcommands:
    (no args)            interactive wizard (default)
    throw [arc]          one-shot: forge an arc with sane defaults
    list                 show available recipes and arcs
    audit <cartridge>    re-run the gate on an existing cartridge

Pipeline (invariant): plan → solve → decode → gate → export.
No external scripts called; this is the consolidated forge.

Doctrinal targets baked in (Rosetta Stone):
    r        = 0.80              (UNIFORM_RADIUS)
    val1     = -0.11             (stage_gain = 0.89)
    val3     = -0.36 @ M0        (zero on unit circle)
    val3     = +0.64 @ M100      (zero dissolved to origin)
    c4       = r²  shared        (gate I2)
    c0       < 0.90              (gate I3)
"""
from __future__ import annotations

import argparse
import json
import math
import random
import struct
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from rich.align import Align
from rich.columns import Columns
from rich.console import Console, Group
from rich.live import Live
from rich.panel import Panel
from rich.prompt import Confirm, IntPrompt, Prompt
from rich.table import Table
from rich.text import Text


# ─── doctrinal constants ──────────────────────────────────────────────────
AUTHORING_SR = 39062.5
NYQUIST_LIMIT_HZ = AUTHORING_SR * 0.5
ACTIVE_STAGES = 6
TOTAL_STAGES = 12

UNIFORM_RADIUS = 0.9999    # right against the f32 edge of self-oscillation.
                            # r=0.9999 → c4≈0.9998. Per-stage pole gain ≈
                            # 1/(1-r²) ≈ 5000× ≈ +74 dB. Pole magnitude
                            # 0.9999 < 1.0 → STAB still passes; any closer
                            # would round to 1.0 in f32 and explode.
                            # Q axis is NOT authored — only M0 and M100 frames
                            # are forged; the cartridge's Q corners duplicate
                            # their morph-paired frame so M0_Q0 == M0_Q100 and
                            # M100_Q0 == M100_Q100. Q stays a runtime knob.
STAGE_GAIN = 1.0           # RULE 5 — stage_gain MUST be ~1.0 (val1 ≈ 0).
                            # Pumping c0 in a 6-stage serial DF2T cascade
                            # multiplies broadband gain into a volume bomb;
                            # it does NOT build vocal peaks. Peak power comes
                            # exclusively from radius (c4 ≥ 0.90) and pole
                            # stacking (c3). Zeros (c1, c2) carve the notches.
EQ_ZERO_OFFSET_SEMITONES = 5.0   # fallback for unknown roles

# Per-role zero placement, calibrated against reference/calibration/Talking_Hedz.json.
# Calibration plants zeros in a "treble basket" (4500–8500 Hz) regardless of
# pole position — every low/mid stage has its zero parked up there, so the
# zeros cooperate to carve a fixed upper-mid notch zone instead of each
# stage trying to carve its own neighbourhood. HF_NYQ is the exception: pole
# is already near Nyquist, so its zero drops an octave BELOW the pole.
#
# mode "fixed_hz"          → zero at this absolute frequency
# mode "semitone_offset"   → zero at pole_hz × 2^(value/12)
ROLE_ZERO_PLACEMENT: dict[str, dict[str, Any]] = {
    # morph_lp lattice (3 stages)
    "LP_POLE":  {"mode": "fixed_hz",        "value": 4600.0},   # treble-basket scaffolding zero
    "LATTICE":  {"mode": "lattice_law",     "value": None},     # zeros copied from prev stage
    "BYPASS_X": {"mode": "passthrough",     "value": None},     # structural dummy
    # talking_formant calibration roles from Talking_Hedz.json.
    # Interior-zero peak resonators: poles move anatomically, zeros stay in
    # the treble basket instead of following each pole locally.
    "LOWMID":  {"mode": "fixed_hz",        "value": 4608.0},
    "BITE_LO": {"mode": "fixed_hz",        "value": 4852.0},
    "BITE_HI": {"mode": "fixed_hz",        "value": 5237.0},
    "AIR":     {"mode": "fixed_hz",        "value": 8392.0},
    # peak_shelf (synthetic — keep semitone offsets so the sweep tracks)
    "PEAK":    {"mode": "semitone_offset", "value": +5.0},
    "RES":     {"mode": "semitone_offset", "value": +7.0},
    "WIDE":    {"mode": "semitone_offset", "value": +9.0},
    "ANCHOR":  {"mode": "fixed_hz",        "value": 7300.0},
    "HF_NYQ":  {"mode": "semitone_offset", "value": -12.0},
}
LP_ZERO_HZ = NYQUIST_LIMIT_HZ * 0.985
HP_ZERO_HZ = max(1.0, AUTHORING_SR * 5e-4)
# ── Type 3 "split-code" frequency compression (FUN_1802c6590) ──────────────
# Compile-time hack for vocal/special stages (stage type byte = 3). At the
# fully-negative morph extreme (M0), any authored frequency byte above 0xDB
# gets force-clamped down to a hardcoded anchor — semitone 57 ≈ 220 Hz —
# locking the vowel anchor point so the cascade can't drift into chaos at
# sweep extremes. Runtime never sees this; coefficients are baked.
TYPE3_FREQ_CEILING_HZ = 4500.0   # Hz proxy for "byte > 0xDB"
TYPE3_ANCHOR_HZ       = 220.0    # Hz at semitone 57 (heritage table)


def _apply_type3_clamp(pole_hz: float, frame_key: str, *, vocal: bool) -> float:
    """At M0 extreme, vocal stages with f > ceiling get clamped to anchor."""
    if not vocal or frame_key != "M0":
        return pole_hz
    if pole_hz > TYPE3_FREQ_CEILING_HZ:
        return TYPE3_ANCHOR_HZ
    return pole_hz


M0_ZERO_RADIUS = 1.0       # zeros sit ON the unit circle — violent carves
M100_ZERO_RADIUS = 0.85    # zeros stay strong at M100 too, so the EE shape
                            # is recognisable (not just a dissolved shelf).
                            # Zeros still pull inward from 1.0 → 0.85 across
                            # the morph, so M0 is harsher than M100, but M100
                            # keeps real formant carving. Calibration's
                            # interior zeros sit at zero_r ≈ 0.7–0.85 — this
                            # is the doctrinally faithful soft-carve zone.

VAL1_TOLERANCE = 0.05      # I3 reframed: val1 = c0 - 1 must satisfy |val1| ≤ 0.05.
                            # Old rule "c0 < 0.90" was inverted — it punished
                            # the doctrinal target (c0 ≈ 1.0) and forced peaks
                            # to come from radius alone. The new check enforces
                            # stage_gain ≈ 1.0; deviation is a "volume bomb".

CORNERS = ("M0_Q0", "M0_Q100", "M100_Q0", "M100_Q100")
CORNER_TO_FRAME = {
    "M0_Q0": "M0", "M0_Q100": "M0",
    "M100_Q0": "M100", "M100_Q100": "M100",
}
CORNER_MORPH_Q = {
    "M0_Q0": (0.0, 0.0), "M0_Q100": (0.0, 1.0),
    "M100_Q0": (1.0, 0.0), "M100_Q100": (1.0, 1.0),
}
PASSTHROUGH_STAGE = {"c0": 1.0, "c1": 0.0, "c2": 0.0, "c3": 0.0, "c4": 0.0}

# Default cartridge output dir
DEFAULT_OUT_DIR = Path("cartridges/factory/generated/qlaw")


# ─── recipes & arc library ────────────────────────────────────────────────
@dataclass(frozen=True)
class StageRole:
    role: str
    landmark: str
    is_boundary: bool
    is_notch: bool


RECIPES: dict[str, dict[str, Any]] = {
    # ── Family A: MorphDesigner — vocal multi-formant (Talking_Hedz lineage) ──
    # 6-stage cascade. Roles match reference/calibration/Talking_Hedz.json
    # s0..s5: HF_Nyquist / LowMid / FormantBite / FormantBite / Air / Anchor.
    # All interior-zero peak resonators. ANCHOR is the foundation.
    "talking_formant": {
        "blurb": "vocal multi-formant (Talking_Hedz lineage, 6-stage MorphDesigner)",
        "arc_kind": "vowel",
        "topology": "cascade_6",
        "stages": (
            StageRole("HF_NYQ",  "air",        is_boundary=False, is_notch=False),
            StageRole("LOWMID",  "throat_f1",  is_boundary=False, is_notch=False),
            StageRole("BITE_LO", "mouth_f2",   is_boundary=False, is_notch=False),
            StageRole("BITE_HI", "bite",       is_boundary=False, is_notch=False),
            StageRole("AIR",     "presence",   is_boundary=False, is_notch=False),
            StageRole("ANCHOR",  "floor",      is_boundary=False, is_notch=False),  # foundation
        ),
    },
    # ── Family B: MorphLP/LPX — resonant LP sweep (3-stage lattice) ─────────
    # FUN_1802c59b0 inner compiler. Stage 0: LP biquad. Stage 1: lattice law
    # (c1=S0.c3, c2=S0.c4, c0=1.0). Stage 2: bypass with cross-reference.
    # Slots 3-11 forced passthrough. active_stages = 3.
    "morph_lp": {
        "blurb": "resonant LP sweep (MorphLP/LPX 3-stage cross-coupled lattice)",
        "arc_kind": "sweep",
        "topology": "lattice_3",
        "stages": (
            StageRole("LP_POLE",  "sweep_freq", is_boundary=False, is_notch=False),
            StageRole("LATTICE",  "sweep_freq", is_boundary=False, is_notch=False),
            StageRole("BYPASS_X", "floor",      is_boundary=False, is_notch=False),
        ),
    },
    # ── peak_shelf — heritage 6-stage cascade for filter sweeps ─────────────
    "peak_shelf": {
        "blurb": "heritage peak/shelf morph (6-stage cascade)",
        "arc_kind": "sweep",
        "topology": "cascade_6",
        "stages": (
            StageRole("HF_NYQ",  "ceiling",     is_boundary=False, is_notch=False),
            StageRole("PEAK",    "sweep_freq",  is_boundary=False, is_notch=False),
            StageRole("RES",     "sweep_freq",  is_boundary=False, is_notch=False),
            StageRole("SHELF",   "sweep_freq",  is_boundary=False, is_notch=True),
            StageRole("WIDE",    "sweep_freq",  is_boundary=False, is_notch=False),
            StageRole("ANCHOR",  "floor",       is_boundary=False, is_notch=False),
        ),
    },
}
PHALANX_ANCHOR = {
    "talking_formant": "LOWMID",   # F1 zone is the Path-B anchor for vocal cascade
    "peak_shelf":      "PEAK",
    # morph_lp uses lattice topology — phalanx anchor doesn't apply.
}

SYNTHETIC: dict[str, dict[str, Any]] = {
    "nasal_anti": {"sources": ("throat_f1", "mouth_f2"), "method": "log_mid"},
}

# Vowel formants (Hz) — F1/F2/F3/F4. Cribbed from docs/sonic_tables/tables.json
# but inlined so the CLI stays self-contained.
VOWELS: dict[str, dict[str, int]] = {
    "ah":    {"f1": 700, "f2": 1220, "f3": 2600, "f4": 3400},
    "ee":    {"f1": 342, "f2": 2322, "f3": 3000, "f4": 3657},
    "oo":    {"f1": 350, "f2": 1250, "f3": 2200, "f4": 3400},
    "er":    {"f1": 490, "f2": 1350, "f3": 1690, "f4": 3400},
    "aw":    {"f1": 570, "f2": 840,  "f3": 2410, "f4": 3400},
    "ih":    {"f1": 427, "f2": 2034, "f3": 2684, "f4": 3400},
    "eh":    {"f1": 580, "f2": 1840, "f3": 2480, "f4": 3400},
    "ae":    {"f1": 660, "f2": 1720, "f3": 2410, "f4": 3400},
    "uh":    {"f1": 640, "f2": 1190, "f3": 2390, "f4": 3400},
    "schwa": {"f1": 500, "f2": 1500, "f3": 2500, "f4": 3400},
}


@dataclass(frozen=True)
class Arc:
    """One arc = one frame_map + the recipe that consumes it.

    There is exactly ONE authoring schema (trench.ui.frame_map.v1). Vowel
    arcs and sweep arcs both produce that schema; they only differ in which
    landmarks live on each frame. The recipe binding chooses which six-stage
    role table to apply to those landmarks.
    """
    name: str
    blurb: str
    recipe: str
    frames: dict[str, dict[str, Any]]


def _vowel_arc(name: str, frm: str, to: str, blurb: str,
               floor_a: int = 80, floor_b: int = 110, air: int = 7000) -> Arc:
    a, b = VOWELS[frm], VOWELS[to]
    return Arc(
        name=name, blurb=blurb, recipe="talking_formant",
        frames={
            "A": {
                "floor":     {"hz": floor_a},
                "throat_f1": {"hz": a["f1"]},
                "mouth_f2":  {"hz": a["f2"]},
                "bite":      {"hz": a["f3"]},
                "presence":  {"hz": a["f4"]},
                "air":       {"hz": air},
            },
            "B": {
                "floor":     {"hz": floor_b},
                "throat_f1": {"hz": b["f1"]},
                "mouth_f2":  {"hz": b["f2"]},
                "bite":      {"hz": b["f3"]},
                "presence":  {"hz": b["f4"]},
                "air":       {"hz": air},
            },
        },
    )


def _sweep_arc(name: str, freq_a: int, freq_b: int, blurb: str,
               floor: int = 50, ceiling: int = 18000,
               recipe: str = "morph_lp") -> Arc:
    """Sweep arcs default-bind to morph_lp (3-stage lattice). Override via
    recipe='peak_shelf' for the heritage 6-stage cascade variant."""
    return Arc(
        name=name, blurb=blurb, recipe=recipe,
        frames={
            "A": {
                "ceiling":    {"hz": ceiling},
                "sweep_freq": {"hz": freq_a},
                "floor":      {"hz": floor},
            },
            "B": {
                "ceiling":    {"hz": ceiling},
                "sweep_freq": {"hz": freq_b},
                "floor":      {"hz": floor},
            },
        },
    )


ARCS: tuple[Arc, ...] = (
    # vocal arcs (recipe = talking_formant)
    _vowel_arc("ah_to_ee", "ah", "ee", "open mouth → closed front"),
    _vowel_arc("ee_to_ah", "ee", "ah", "closed front → open mouth"),
    _vowel_arc("oo_to_ee", "oo", "ee", "rounded back → closed front (belch)"),
    _vowel_arc("ee_to_oo", "ee", "oo", "closed front → rounded back"),
    _vowel_arc("er_to_aw", "er", "aw", "mid-rhotic → open back"),
    _vowel_arc("aw_to_er", "aw", "er", "open back → mid-rhotic"),
    _vowel_arc("ah_to_oo", "ah", "oo", "open → rounded"),
    _vowel_arc("oo_to_ah", "oo", "ah", "rounded → open"),
    _vowel_arc("ih_to_ae", "ih", "ae", "high lax → low front"),
    _vowel_arc("ae_to_ih", "ae", "ih", "low front → high lax"),
    # sweep arcs (recipe = peak_shelf)
    _sweep_arc("220_to_8k",         220,  8000, "classic low → high mid-band sweep"),
    _sweep_arc("bass_to_air",        80, 12000, "sub bass → airy top"),
    _sweep_arc("chest_to_presence", 250,  5000, "chest body → vocal presence"),
    _sweep_arc("rumble_to_sting",    60,  6000, "subbass → sting"),
    _sweep_arc("low_to_nyquist",    100, 18000, "wide-open low → high"),
    _sweep_arc("high_to_low",      8000,   220, "reverse sweep — high → low"),
)

POLICIES = {
    "A": ("anatomical_scatter", "distributed poles → vocal body"),
    "B": ("resonance_phalanx",  "stacked poles → spectral laser"),
}


# ─── ASCII art (Fallout: New Vegas Trench Trading Co. — cheerful end times) ───
TITLE_ART = r"""
   ████████╗██████╗ ███████╗███╗   ██╗ ██████╗██╗  ██╗
   ╚══██╔══╝██╔══██╗██╔════╝████╗  ██║██╔════╝██║  ██║
      ██║   ██████╔╝█████╗  ██╔██╗ ██║██║     ███████║
      ██║   ██╔══██╗██╔══╝  ██║╚██╗██║██║     ██╔══██║
      ██║   ██║  ██║███████╗██║ ╚████║╚██████╗██║  ██║
      ╚═╝   ╚═╝  ╚═╝╚══════╝╚═╝  ╚═══╝ ╚═════╝╚═╝  ╚═╝
        ╔═══════════════════════════════════════════╗
        ║   T R A D I N G   C O.  —  est. 20XX     ║
        ║   "fresh cartridges, hot off the anvil"  ║
        ╚═══════════════════════════════════════════╝
"""

OPERATOR_IDLE = r"""
       ╔════════════════════════╗
       ║   ░░░░░░░░░░░░░░░░░░   ║
       ║    ▄▀▀▀▀▀▀▀▀▀▀▀▀▀▄     ║
       ║   █   ◕        ◕  █    ║
       ║   █       ▼       █    ║
       ║   █    ╲_____╱    █    ║
       ║    ▀▄▄▄▄▄▄▄▄▄▄▄▄▄▀     ║
       ║      ╔═══════╗         ║
       ║   ═══╣ APRON ╠═══      ║
       ║      ╚═══════╝         ║
       ║   ═══════════════════  ║
       ║    [ GENERAL GOODS ]   ║
       ╚════════════════════════╝
            ─ OL' CHET ─
       Trench Trading Co.
"""

OPERATOR_THINKING = r"""
       ╔════════════════════════╗
       ║   ░░░░░░░░░░░░░░░░░░   ║
       ║    ▄▀▀▀▀▀▀▀▀▀▀▀▀▀▄     ║
       ║   █   •        •  █    ║
       ║   █       ▽       █    ║
       ║   █    ╲ ─── ╱    █    ║
       ║    ▀▄▄▄▄▄▄▄▄▄▄▄▄▄▀     ║
       ║      ╔═══════╗         ║
       ║   ═══╣ APRON ╠═══      ║
       ║      ╚═══════╝         ║
       ║   ═══════════════════  ║
       ║   ".......hmm......."  ║
       ╚════════════════════════╝
       *scratchin' the chin*
"""

OPERATOR_FORGING = r"""
       ╔════════════════════════╗
       ║   ✦  ° . ° · ° .  ✦    ║
       ║    ▄▀▀▀▀▀▀▀▀▀▀▀▀▀▄     ║
       ║   █   ●        ●  █    ║
       ║   █       ▼       █    ║
       ║   █    ▔▔▔▔▔▔▔    █    ║
       ║    ▀▄▄▄▄▄▄▄▄▄▄▄▄▄▀     ║
       ║      ╔═══════╗         ║
       ║   ═══╣ APRON ╠═══      ║
       ║      ╚═══════╝         ║
       ║       ▌▌  ▐▐           ║
       ║   ╔═════════════╗      ║
       ║   ║ ░  ANVIL  ░ ║      ║
       ║   ╚═════════════╝      ║
       ╚════════════════════════╝
        *clang clang clang!*
"""

OPERATOR_PASSED = r"""
       ╔════════════════════════╗
       ║  ★ ★ ★  PIPING HOT  ★ ★║
       ║    ▄▀▀▀▀▀▀▀▀▀▀▀▀▀▄     ║
       ║   █   ◕        ◕  █    ║
       ║   █       ▽       █    ║
       ║   █    ╲ ‿‿‿ ╱    █    ║
       ║    ▀▄▄▄▄▄▄▄▄▄▄▄▄▄▀     ║
       ║      ╔═══════╗         ║
       ║   ═══╣ APRON ╠═══      ║
       ║      ╚═══════╝         ║
       ║   ┌───────────────┐    ║
       ║   │  CARTRIDGE !  │    ║
       ║   └───────────────┘    ║
       ╚════════════════════════╝
       "She's a beaut, partner!"
"""

OPERATOR_REJECTED = r"""
       ╔════════════════════════╗
       ║   ✗ ✗  REJECT  ✗ ✗     ║
       ║    ▄▀▀▀▀▀▀▀▀▀▀▀▀▀▄     ║
       ║   █   ●        ●  █    ║
       ║   █       ▽       █    ║
       ║   █    ╲ ︵︵︵ ╱    █    ║
       ║    ▀▄▄▄▄▄▄▄▄▄▄▄▄▄▀     ║
       ║      ╔═══════╗         ║
       ║   ═══╣ APRON ╠═══      ║
       ║      ╚═══════╝         ║
       ║      ╳ ╳ ╳ ╳ ╳         ║
       ║   *back to the bench*  ║
       ╚════════════════════════╝
       "Whoo-ee, math went sideways."
"""

FORGE_STAGES = (
    ("SOLVING THE Q-LAW",   "░", "yellow"),
    ("STACKING POLES",       "▒", "bright_yellow"),
    ("CARVING ZEROS",        "▓", "bright_red"),
    ("RUNNING GATE",         "█", "bright_green"),
)

PASS_BANNER = r"""
 ██████╗  █████╗ ████████╗███████╗    ██████╗  █████╗ ███████╗███████╗
██╔════╝ ██╔══██╗╚══██╔══╝██╔════╝    ██╔══██╗██╔══██╗██╔════╝██╔════╝
██║  ███╗███████║   ██║   █████╗      ██████╔╝███████║███████╗███████╗
██║   ██║██╔══██║   ██║   ██╔══╝      ██╔═══╝ ██╔══██║╚════██║╚════██║
╚██████╔╝██║  ██║   ██║   ███████╗    ██║     ██║  ██║███████║███████║
 ╚═════╝ ╚═╝  ╚═╝   ╚═╝   ╚══════╝    ╚═╝     ╚═╝  ╚═╝╚══════╝╚══════╝
"""

FAIL_BANNER = r"""
 ██████╗  █████╗ ████████╗███████╗    ███████╗ █████╗ ██╗██╗
██╔════╝ ██╔══██╗╚══██╔══╝██╔════╝    ██╔════╝██╔══██╗██║██║
██║  ███╗███████║   ██║   █████╗      █████╗  ███████║██║██║
██║   ██║██╔══██║   ██║   ██╔══╝      ██╔══╝  ██╔══██║██║██║
╚██████╔╝██║  ██║   ██║   ███████╗    ██║     ██║  ██║██║███████╗
 ╚═════╝ ╚═╝  ╚═╝   ╚═╝   ╚══════╝    ╚═╝     ╚═╝  ╚═╝╚═╝╚══════╝
"""


# ─── helpers ──────────────────────────────────────────────────────────────
class ForgeError(Exception):
    pass


def f32(value: float) -> float:
    return struct.unpack("<f", struct.pack("<f", float(value)))[0]


def bark_to_hz(bark: float) -> float:
    return 1960.0 * (bark + 0.53) / (26.28 - bark)


def hz_to_bark(hz: float) -> float:
    return 26.81 * hz / (1960.0 + hz) - 0.53


def _frame_landmark_hz(frame: dict[str, Any], landmark: str, frame_name: str) -> float:
    if landmark in SYNTHETIC:
        spec = SYNTHETIC[landmark]
        sources = spec["sources"]
        method = spec["method"]
        for src in sources:
            if src not in frame:
                raise ForgeError(
                    f"frame {frame_name!r}: synthetic '{landmark}' needs source '{src}'"
                )
        freqs = [_frame_landmark_hz(frame, src, frame_name) for src in sources]
        if method == "log_mid":
            product = 1.0
            for f in freqs:
                product *= f
            return product ** (1.0 / len(freqs))
        raise ForgeError(f"unknown synthetic method {method!r}")

    if landmark not in frame:
        raise ForgeError(f"frame {frame_name!r} missing landmark '{landmark}'")
    value = frame[landmark]
    if not isinstance(value, dict):
        raise ForgeError(f"frame {frame_name!r}.{landmark} must be an object")
    if "hz" in value:
        hz = float(value["hz"])
    elif "bark" in value:
        hz = bark_to_hz(float(value["bark"]))
    else:
        raise ForgeError(f"frame {frame_name!r}.{landmark} needs hz or bark")
    if not math.isfinite(hz) or hz <= 0.0 or hz >= NYQUIST_LIMIT_HZ:
        raise ForgeError(f"frame {frame_name!r}.{landmark} freq out of range: {hz}")
    return hz


def _semitone_offset_hz(base_hz: float, semitones: float) -> float:
    out = base_hz * (2.0 ** (semitones / 12.0))
    return min(out, NYQUIST_LIMIT_HZ * 0.985)


def arc_to_frame_map(arc: Arc) -> dict[str, Any]:
    """Wrap arc frames in the canonical schema. ONE authoring schema for everything."""
    return {
        "schema": "trench.ui.frame_map.v1",
        "name": arc.name,
        "category": "VOW" if arc.recipe == "talking_formant" else "PSM",
        "boost": 1.0,
        "source_affinity": ["lead synth stabs", "one-shots"],
        "frames": arc.frames,
    }


def find_arc(name: str) -> Arc:
    for a in ARCS:
        if a.name == name:
            return a
    raise ForgeError(
        f"no arc {name!r}; known: {', '.join(a.name for a in ARCS)}"
    )


# ─── plan / solve / decode / gate ─────────────────────────────────────────
@dataclass
class PlannedStage:
    role: str
    landmark: str
    is_boundary: bool
    is_notch: bool
    hz_m0: float
    hz_m100: float


def plan(frame_map: dict[str, Any], recipe_name: str) -> tuple[list[PlannedStage], dict[str, Any]]:
    if recipe_name not in RECIPES:
        raise ForgeError(f"unknown recipe {recipe_name!r}; known: {sorted(RECIPES)}")
    if frame_map.get("schema") != "trench.ui.frame_map.v1":
        raise ForgeError(
            f"expected schema trench.ui.frame_map.v1, got {frame_map.get('schema')!r}"
        )
    frames = frame_map.get("frames", {})
    if "A" not in frames or "B" not in frames:
        raise ForgeError("frame_map needs frames.A and frames.B")

    stages: list[PlannedStage] = []
    for role in RECIPES[recipe_name]["stages"]:
        hz_a = _frame_landmark_hz(frames["A"], role.landmark, "A")
        hz_b = _frame_landmark_hz(frames["B"], role.landmark, "B")
        stages.append(PlannedStage(
            role=role.role, landmark=role.landmark,
            is_boundary=role.is_boundary, is_notch=role.is_notch,
            hz_m0=hz_a, hz_m100=hz_b,
        ))

    meta = {
        "recipe": recipe_name,
        "name": frame_map.get("name", "unnamed"),
        "category": frame_map.get("category", "PROG"),
        "source_affinity": frame_map.get("source_affinity", []),
    }
    return stages, meta


def _zero_hz_for_role(stage: PlannedStage, frame_key: str) -> float:
    pole_hz = stage.hz_m0 if frame_key == "M0" else stage.hz_m100
    if stage.is_notch:
        return pole_hz  # zero on unit circle AT pole = pure notch
    spec = ROLE_ZERO_PLACEMENT.get(stage.role)
    if spec is None:
        out = _semitone_offset_hz(pole_hz, EQ_ZERO_OFFSET_SEMITONES)
    elif spec["mode"] == "fixed_hz":
        out = float(spec["value"])
    else:  # semitone_offset
        out = pole_hz * (2.0 ** (float(spec["value"]) / 12.0))
    return max(20.0, min(out, NYQUIST_LIMIT_HZ * 0.985))


def _solve_stage(
    stage: PlannedStage,
    frame_key: str,
    *,
    pole_hz: float,
    radius: float,
    stage_gain: float,
    zero_radius: float,
    is_anchor_phalanx: bool,
) -> dict[str, float]:
    theta_pole = 2.0 * math.pi * pole_hz / AUTHORING_SR
    a1 = f32(-2.0 * radius * math.cos(theta_pole))
    radius_q = f32(radius)
    r2 = radius_q * radius_q

    if is_anchor_phalanx or zero_radius == 0.0:
        b1_target = 0.0
        b2_target = 0.0
    else:
        zero_hz = _zero_hz_for_role(stage, frame_key)
        theta_z = 2.0 * math.pi * zero_hz / AUTHORING_SR
        b1_target = -2.0 * zero_radius * math.cos(theta_z)
        b2_target = zero_radius * zero_radius

    val1 = f32(stage_gain - 1.0)
    val2 = f32(b1_target - a1)
    val3 = f32(r2 - b2_target)

    return {
        "c0": f32(1.0 + val1),
        "c1": f32(a1 + val2),
        "c2": f32(r2 - val3),
        "c3": f32(a1),
        "c4": f32(r2),
    }


def _recipe_for_stages(stages: list[PlannedStage]) -> str:
    seq = tuple(s.role for s in stages)
    for name, recipe in RECIPES.items():
        if tuple(r.role for r in recipe["stages"]) == seq:
            return name
    raise ForgeError("unrecognised stage sequence")


def _build_stage_biquad(pole_hz: float, zero_hz: float, zero_radius: float,
                        stage_gain: float, radius: float) -> dict[str, float]:
    """Compute c0..c4 for a single explicit_zero biquad."""
    theta_pole = 2.0 * math.pi * pole_hz / AUTHORING_SR
    a1 = f32(-2.0 * radius * math.cos(theta_pole))
    radius_q = f32(radius)
    r2 = radius_q * radius_q

    if zero_radius == 0.0:
        b1_target = 0.0
        b2_target = 0.0
    else:
        theta_z = 2.0 * math.pi * zero_hz / AUTHORING_SR
        b1_target = -2.0 * zero_radius * math.cos(theta_z)
        b2_target = zero_radius * zero_radius

    val1 = f32(stage_gain - 1.0)
    val2 = f32(b1_target - a1)
    val3 = f32(r2 - b2_target)
    return {
        "c0": f32(1.0 + val1),
        "c1": f32(a1 + val2),
        "c2": f32(r2 - val3),
        "c3": f32(a1),
        "c4": f32(r2),
    }


def _solve_lattice_3(stages: list[PlannedStage], boost: float) -> dict[str, Any]:
    """MorphLP/LPX 3-stage cross-coupled lattice (FUN_1802c59b0).

    Stage 0: independent LP biquad with treble-basket scaffolding zeros.
    Stage 1: LATTICE LAW — c1=S0.c3, c2=S0.c4, c0=1.0; own pole at S1 freq.
    Stage 2: bypass with cross-reference (passthrough).
    Slots 3-11: passthrough.

    Type 3 clamp applies at M0 to vocal stages above the 0xDB freq ceiling.
    """
    if len(stages) < 3:
        raise ForgeError("lattice_3 requires 3 stage roles (LP_POLE, LATTICE, BYPASS_X)")

    keyframes: list[dict[str, Any]] = []
    for corner in CORNERS:
        frame_key = CORNER_TO_FRAME[corner]
        zero_radius = M0_ZERO_RADIUS if frame_key == "M0" else M100_ZERO_RADIUS

        # ── Stage 0 — LP_POLE (independent biquad) ──────────────────────────
        s0_role = stages[0]
        s0_pole_hz = s0_role.hz_m0 if frame_key == "M0" else s0_role.hz_m100
        s0_pole_hz = _apply_type3_clamp(s0_pole_hz, frame_key, vocal=True)
        s0_spec = ROLE_ZERO_PLACEMENT.get(s0_role.role, {"mode": "semitone_offset", "value": 0.0})
        if s0_spec["mode"] == "fixed_hz":
            s0_zero_hz = float(s0_spec["value"])
        else:  # semitone_offset
            s0_zero_hz = s0_pole_hz * (2.0 ** (float(s0_spec["value"] or 0.0) / 12.0))
        s0_zero_hz = max(20.0, min(s0_zero_hz, NYQUIST_LIMIT_HZ * 0.985))
        s0 = _build_stage_biquad(
            s0_pole_hz, s0_zero_hz, zero_radius, STAGE_GAIN, UNIFORM_RADIUS,
        )

        # ── Stage 1 — LATTICE LAW ───────────────────────────────────────────
        # Stage 1's NUMERATOR equals Stage 0's DENOMINATOR: c1=S0.c3, c2=S0.c4
        # Own pole/denominator independent.
        s1_role = stages[1]
        s1_pole_hz = s1_role.hz_m0 if frame_key == "M0" else s1_role.hz_m100
        s1_pole_hz = _apply_type3_clamp(s1_pole_hz, frame_key, vocal=True)
        theta_s1 = 2.0 * math.pi * s1_pole_hz / AUTHORING_SR
        a1_s1 = f32(-2.0 * UNIFORM_RADIUS * math.cos(theta_s1))
        radius_q = f32(UNIFORM_RADIUS)
        r2_s1 = radius_q * radius_q
        s1 = {
            "c0": f32(1.0),         # lattice law: stage_gain locked to unity
            "c1": s0["c3"],         # lattice law: c1 == S0.c3 (numerator zeros == prev poles)
            "c2": s0["c4"],         # lattice law: c2 == S0.c4
            "c3": a1_s1,            # own denominator pole
            "c4": r2_s1,            # own denominator radius²
        }

        # ── Stage 2 — BYPASS_X (structural dummy) ──────────────────────────
        s2 = dict(PASSTHROUGH_STAGE)

        compiled = [s0, s1, s2]
        while len(compiled) < TOTAL_STAGES:
            compiled.append(dict(PASSTHROUGH_STAGE))

        morph, q = CORNER_MORPH_Q[corner]
        keyframes.append({
            "label": corner, "morph": morph, "q": q,
            "boost": boost, "stages": compiled,
        })
    return {"keyframes": keyframes}


def solve(stages: list[PlannedStage], policy: str, boost: float = 1.0) -> dict[str, Any]:
    """Dispatch by recipe topology."""
    recipe = _recipe_for_stages(stages)
    topology = RECIPES[recipe].get("topology", "cascade_6")

    if topology == "lattice_3":
        # Lattice doesn't take policy — it's a fixed cross-coupled topology
        return _solve_lattice_3(stages, boost)

    # ── cascade_6 (existing path) ────────────────────────────────────────────
    if policy not in {"A", "B"}:
        raise ForgeError("policy must be 'A' or 'B'")
    if policy == "B":
        anchor_role = PHALANX_ANCHOR.get(recipe)
        if anchor_role is None:
            raise ForgeError("phalanx anchor not configured for this recipe")
    else:
        anchor_role = None

    keyframes: list[dict[str, Any]] = []
    for corner in CORNERS:
        frame_key = CORNER_TO_FRAME[corner]
        zero_radius = M0_ZERO_RADIUS if frame_key == "M0" else M100_ZERO_RADIUS

        if policy == "B":
            anchor_stage = next(s for s in stages if s.role == anchor_role)
            shared_pole_hz = anchor_stage.hz_m0 if frame_key == "M0" else anchor_stage.hz_m100

        compiled: list[dict[str, float]] = []
        for stage in stages:
            if policy == "A":
                pole_hz = stage.hz_m0 if frame_key == "M0" else stage.hz_m100
                is_anchor_phalanx = False
            else:
                pole_hz = shared_pole_hz
                is_anchor_phalanx = (stage.role == anchor_role)
            compiled.append(_solve_stage(
                stage, frame_key,
                pole_hz=pole_hz, radius=UNIFORM_RADIUS, stage_gain=STAGE_GAIN,
                zero_radius=zero_radius, is_anchor_phalanx=is_anchor_phalanx,
            ))
        while len(compiled) < TOTAL_STAGES:
            compiled.append(dict(PASSTHROUGH_STAGE))

        morph, q = CORNER_MORPH_Q[corner]
        keyframes.append({
            "label": corner, "morph": morph, "q": q, "boost": boost, "stages": compiled,
        })
    return {"keyframes": keyframes}


@dataclass
class HeritageView:
    a1: float; r: float; val1: float; val2: float; val3: float


def decode_heritage(stage: dict[str, float]) -> HeritageView:
    c0, c1, c2, c3, c4 = (float(stage[f"c{i}"]) for i in range(5))
    if c3 == 0.0 and c4 == 0.0:
        return HeritageView(0.0, 0.0, 0.0, 0.0, 0.0)
    return HeritageView(
        a1=c3, r=math.sqrt(max(c4, 0.0)),
        val1=c0 - 1.0, val2=c1 - c3, val3=c4 - c2,
    )


@dataclass
class GateResult:
    passed: bool
    failures: list[str] = field(default_factory=list)


def _stab_check(s: dict, label: str, stage_idx: int) -> str | None:
    """Return failure msg if pole magnitude ≥ 1.0, else None."""
    if not (s["c3"] or s["c4"]):
        return None
    disc = s["c3"] ** 2 - 4.0 * s["c4"]
    if disc < 0:
        pole_mag = math.sqrt(s["c4"])
    else:
        r1 = (-s["c3"] + math.sqrt(disc)) / 2.0
        r2 = (-s["c3"] - math.sqrt(disc)) / 2.0
        pole_mag = max(abs(r1), abs(r2))
    if pole_mag >= 1.0:
        return f"{label} STAB: stage {stage_idx} pole_mag={pole_mag:.4f} ≥ 1.0"
    return None


def _gate_cascade_6(cartridge: dict[str, Any]) -> GateResult:
    """Original 6-stage cascade gate: 6 active, shared c4, val1≈0, stable."""
    failures: list[str] = []
    for kf in cartridge["keyframes"]:
        label = kf["label"]
        stages = kf["stages"]
        if len(stages) != TOTAL_STAGES:
            failures.append(f"{label}: expected {TOTAL_STAGES} stages, got {len(stages)}")
            continue
        actives = [s for s in stages if s["c3"] != 0.0 or s["c4"] != 0.0]
        if len(actives) != ACTIVE_STAGES:
            failures.append(f"{label} I1: {len(actives)} active stages (need {ACTIVE_STAGES})")
            continue
        c4_ref = actives[0]["c4"]
        for i, s in enumerate(stages):
            if (s["c3"] or s["c4"]) and s["c4"] != c4_ref:
                failures.append(f"{label} I2: stage {i} c4 drifts ({s['c4']!r} vs {c4_ref!r})")
        for i, s in enumerate(stages):
            if not (s["c3"] or s["c4"]):
                continue
            val1 = s["c0"] - 1.0
            if abs(val1) > VAL1_TOLERANCE:
                failures.append(
                    f"{label} I3 (Rule 5): stage {i} val1={val1:+.4f} "
                    f"(c0={s['c0']:.4f}) — pumping c0 stacks broadband gain"
                )
        for i, s in enumerate(stages):
            msg = _stab_check(s, label, i)
            if msg:
                failures.append(msg)
    return GateResult(passed=not failures, failures=failures)


def _gate_lattice_3(cartridge: dict[str, Any]) -> GateResult:
    """Lattice gate: enforces FUN_1802c59b0 cross-coupled topology.

    L1: Slots 0 and 1 must be active resonators.
    L2: Slot 2 must be passthrough (structural dummy / "bypass with cross-ref").
    L3: Slots 3-11 must be passthrough (active_stages = 3 max).
    L4: LATTICE LAW — Stage 1's c1 == Stage 0's c3, c2 == c4. Bit-identical f32.
    L5: Both Stage 0 and Stage 1 have c0 ≈ 1.0 (Rule 5 / Lattice convention).
    L6: Both stages stable (pole_mag < 1).
    """
    failures: list[str] = []
    for kf in cartridge["keyframes"]:
        label = kf["label"]
        stages = kf["stages"]
        if len(stages) != TOTAL_STAGES:
            failures.append(f"{label}: expected {TOTAL_STAGES} stages, got {len(stages)}")
            continue
        s0, s1, s2 = stages[0], stages[1], stages[2]

        # L1: slots 0 and 1 active
        if not (s0["c3"] or s0["c4"]):
            failures.append(f"{label} L1: slot 0 (LP_POLE) is inactive")
        if not (s1["c3"] or s1["c4"]):
            failures.append(f"{label} L1: slot 1 (LATTICE) is inactive")

        # L2: slot 2 passthrough
        if s2["c3"] != 0.0 or s2["c4"] != 0.0:
            failures.append(f"{label} L2: slot 2 (BYPASS_X) must be passthrough — "
                            f"got c3={s2['c3']!r} c4={s2['c4']!r}")

        # L3: slots 3-11 passthrough
        for i in range(3, TOTAL_STAGES):
            s = stages[i]
            if s["c3"] != 0.0 or s["c4"] != 0.0:
                failures.append(f"{label} L3: slot {i} must be passthrough (active_stages = 3)")

        # L4: LATTICE LAW — f32 bit-identical
        if s1["c1"] != s0["c3"]:
            failures.append(
                f"{label} L4 LATTICE: c1[S1]={s1['c1']!r} != c3[S0]={s0['c3']!r}"
            )
        if s1["c2"] != s0["c4"]:
            failures.append(
                f"{label} L4 LATTICE: c2[S1]={s1['c2']!r} != c4[S0]={s0['c4']!r}"
            )

        # L5: c0 ≈ 1.0 on active stages
        for idx, s in ((0, s0), (1, s1)):
            val1 = s["c0"] - 1.0
            if abs(val1) > VAL1_TOLERANCE:
                failures.append(
                    f"{label} L5 (Rule 5): slot {idx} val1={val1:+.4f} (c0={s['c0']:.4f})"
                )

        # L6: stability
        for i, s in enumerate((s0, s1)):
            msg = _stab_check(s, label, i)
            if msg:
                failures.append(msg)
    return GateResult(passed=not failures, failures=failures)


def gate(cartridge: dict[str, Any]) -> GateResult:
    expected = set(CORNERS)
    found = {kf["label"] for kf in cartridge["keyframes"]}
    if expected - found:
        return GateResult(False, [f"missing corners: {sorted(expected - found)}"])

    topology = cartridge.get("qlaw", {}).get("topology", "cascade_6")
    if topology == "lattice_3":
        return _gate_lattice_3(cartridge)
    return _gate_cascade_6(cartridge)


def build_cartridge(meta: dict[str, Any], solved: dict[str, Any], policy: str) -> dict[str, Any]:
    recipe = meta["recipe"]
    topology = RECIPES[recipe].get("topology", "cascade_6")
    return {
        "format": "compiled-v1",
        "name": meta["name"],
        "provenance": "forge_cli",
        "sampleRate": AUTHORING_SR,
        "qlaw": {
            "topology": topology,
            "policy": policy,
            "policy_name": POLICIES[policy][0],
            "uniform_radius": UNIFORM_RADIUS,
            "stage_gain": STAGE_GAIN,
            "m0_zero_radius": M0_ZERO_RADIUS,
            "m100_zero_radius": M100_ZERO_RADIUS,
            "type3_freq_ceiling_hz": TYPE3_FREQ_CEILING_HZ,
            "type3_anchor_hz": TYPE3_ANCHOR_HZ,
        },
        "recipe": recipe,
        "category": meta.get("category", "PROG"),
        "source_affinity": meta.get("source_affinity", []),
        "keyframes": solved["keyframes"],
    }


# ─── presentation ─────────────────────────────────────────────────────────
def show_title(console: Console, cinematic: bool = False) -> None:
    if not cinematic:
        console.print(Text(TITLE_ART, style="bold red"))
        return
    # cinematic intro — reveal lines top-down with a brief delay
    lines = TITLE_ART.splitlines()
    rendered: list[Text] = []
    with Live(Group(*rendered), console=console, refresh_per_second=24,
              transient=False) as live:
        for line in lines:
            rendered.append(Text(line, style="bold red"))
            live.update(Group(*rendered))
            time.sleep(0.05)


def show_operator(console: Console, art: str, line: str, style: str = "bold yellow1") -> None:
    console.print(Align.left(Text(art, style="bold gold3")))
    console.print(Align.left(Panel(Text(line, style=style), border_style="gold3",
                                   padding=(0, 2))))


def render_recipe_panel() -> Panel:
    body = Text()
    for i, (key, spec) in enumerate(RECIPES.items(), 1):
        body.append(f"  [{i}] ", style="bold yellow3")
        body.append(f"{key:18s}", style="bold yellow1")
        body.append(f"{spec['blurb']}\n", style="grey70")
    return Panel(body, title="[bold gold3]TODAY'S SPECIALS[/bold gold3] [grey50](recipes)[/grey50]",
                 border_style="gold3")


def render_arc_panel() -> Panel:
    body = Text()
    last_recipe = None
    for i, arc in enumerate(ARCS, 1):
        if arc.recipe != last_recipe:
            if last_recipe is not None:
                body.append("\n")
            kind_label = "VOCAL ARCS  (talking_formant)" if arc.recipe == "talking_formant" \
                else "FILTER SWEEPS  (peak_shelf)"
            body.append(f"  ── {kind_label} ──\n", style="bold gold3")
            last_recipe = arc.recipe
        body.append(f"  [{i:>2}] ", style="bold yellow3")
        body.append(f"{arc.name:22s}", style="bold yellow1")
        body.append(f"{arc.blurb}\n", style="grey70")
    return Panel(body, title="[bold gold3]TODAY'S CARTRIDGES[/bold gold3] [grey50](pick a route)[/grey50]",
                 border_style="gold3")


def render_policy_panel() -> Panel:
    body = Text()
    for key, (name, blurb) in POLICIES.items():
        body.append(f"  [{key}] ", style="bold yellow3")
        body.append(f"{name:24s}", style="bold yellow1")
        body.append(f"{blurb}\n", style="grey70")
    return Panel(body, title="[bold gold3]WORKBENCH MODES[/bold gold3] [grey50](policy)[/grey50]",
                 border_style="gold3")


def render_anvil_animation(console: Console, total_ticks: int = 16) -> None:
    """Live-rendered forge animation: cycling sparks + step labels."""
    spark_pool = "✦✧·.°*"
    rng = random.Random(0xC0FFEE)

    def make_frame(i: int) -> Panel:
        stage_idx = (i * len(FORGE_STAGES)) // total_ticks
        label, fill, color = FORGE_STAGES[min(stage_idx, len(FORGE_STAGES) - 1)]
        spark_top = "".join(rng.choice(spark_pool + "       ") for _ in range(40))
        spark_mid = "".join(rng.choice(spark_pool + "         ") for _ in range(40))
        bar_filled = (fill * ((i % 20) + 1)).ljust(20, " ")
        body = (
            f"  {spark_top}\n"
            f"      ╱│╲  ╱│╲     ╱│╲\n"
            f"   {spark_mid}\n"
            f"   ╔════════════════════╗\n"
            f"   ║{bar_filled}║\n"
            f"   ╚════════════════════╝\n"
            f"        ▀▀▀▀▀▀▀▀▀▀▀\n"
            f"         {label}"
        )
        return Panel(Text(body, style=f"bold {color}"),
                     border_style=color, padding=(0, 2))

    with Live(make_frame(0), console=console, refresh_per_second=12,
              transient=True) as live:
        for i in range(1, total_ticks + 1):
            time.sleep(0.08)
            live.update(make_frame(i))


ROLE_GLYPH = {
    # talking_formant — calibration role names from Talking_Hedz.json
    "HF_NYQ":  ("⌜", "blue"),           # spectral ceiling
    "LOWMID":  ("◉", "bright_red"),     # F1 zone body
    "BITE_LO": ("◉", "bright_yellow"),  # F2 zone formant
    "BITE_HI": ("◉", "bright_green"),   # F3 / bite formant
    "AIR":     ("◉", "cyan"),           # presence
    "ANCHOR":  ("⊙", "magenta"),        # foundation — sub-low body weight
    # peak_shelf
    "PEAK":    ("◉", "bright_yellow"),
    "RES":     ("◉", "bright_red"),
    "SHELF":   ("◇", "magenta"),
    "WIDE":    ("◉", "cyan"),
}


def _stage_pole_zero(st: dict[str, float]) -> tuple[complex, complex] | None:
    """Return (pole, zero) for an active stage, both as complex (real, imag).
    Returns None for passthrough."""
    c0, c1, c2, c3, c4 = (float(st[f"c{i}"]) for i in range(5))
    if c3 == 0.0 and c4 == 0.0:
        return None
    # Pole: roots of 1 + c3·z⁻¹ + c4·z⁻² = 0  →  z² + c3·z + c4 = 0
    disc_p = c3 * c3 - 4.0 * c4
    if disc_p < 0:
        # complex conjugate pair on circle of radius √c4
        re = -c3 / 2.0
        im = math.sqrt(-disc_p) / 2.0
        pole = complex(re, im)
    else:
        r1 = (-c3 + math.sqrt(disc_p)) / 2.0
        r2 = (-c3 - math.sqrt(disc_p)) / 2.0
        pole = complex(max(abs(r1), abs(r2)), 0.0)
    # Zero: roots of c0 + c1·z⁻¹ + c2·z⁻² = 0  →  c0·z² + c1·z + c2 = 0
    if abs(c0) < 1e-12:
        zero = complex(0, 0)
    elif c2 == 0.0 and c1 == 0.0:
        zero = complex(0, 0)
    else:
        disc_z = c1 * c1 - 4.0 * c0 * c2
        if disc_z < 0:
            re = -c1 / (2.0 * c0)
            im = math.sqrt(-disc_z) / (2.0 * c0)
            zero = complex(re, im)
        else:
            r1 = (-c1 + math.sqrt(disc_z)) / (2.0 * c0)
            r2 = (-c1 - math.sqrt(disc_z)) / (2.0 * c0)
            zero = complex(max(abs(r1), abs(r2)) if abs(r1) != abs(r2) else r1, 0.0)
    return pole, zero


def render_zplane(stages: list[PlannedStage], cartridge: dict[str, Any], corner: str,
                  width: int = 39, height: int = 19) -> Panel:
    """ASCII Z-plane: poles (●) and zeros (○) on the complex plane with unit circle."""
    cx, cy = width // 2, height // 2
    rx, ry = cx - 1, cy - 1
    grid = [[(" ", "dim") for _ in range(width)] for _ in range(height)]

    # background: real axis, imag axis, unit circle
    for x in range(width):
        grid[cy][x] = ("─", "grey37")
    for y in range(height):
        grid[y][cx] = ("│", "grey37")
    grid[cy][cx] = ("┼", "grey50")
    # unit circle dots
    for deg in range(0, 360, 4):
        t = math.radians(deg)
        x = int(round(cx + math.cos(t) * rx))
        y = int(round(cy - math.sin(t) * ry))
        if 0 <= x < width and 0 <= y < height and grid[y][x][0] in (" ", "─", "│", "┼"):
            grid[y][x] = ("·", "grey50")

    # poles + zeros per stage
    kf = next(k for k in cartridge["keyframes"] if k["label"] == corner)
    for idx, st in enumerate(kf["stages"][: len(stages)]):
        pz = _stage_pole_zero(st)
        if pz is None:
            continue
        pole, zero = pz
        glyph, color = ROLE_GLYPH.get(stages[idx].role, ("●", "white"))
        # plot pole conjugate pair
        for sign in (1, -1):
            px = int(round(cx + pole.real * rx))
            py = int(round(cy - sign * pole.imag * ry))
            if 0 <= px < width and 0 <= py < height:
                grid[py][px] = ("●", f"bold {color}")
        # plot zero conjugate pair (only if not at origin)
        if abs(zero.imag) > 0.02 or abs(zero.real) > 0.02:
            for sign in (1, -1):
                zx = int(round(cx + zero.real * rx))
                zy = int(round(cy - sign * zero.imag * ry))
                if 0 <= zx < width and 0 <= zy < height:
                    if grid[zy][zx][0] not in ("●",):  # don't stomp poles
                        grid[zy][zx] = ("○", color)

    text = Text()
    for row in grid:
        for ch, sty in row:
            text.append(ch, style=sty)
        text.append("\n")

    legend = Text()
    legend.append("● pole  ", style="bold red")
    legend.append("○ zero  ", style="cyan")
    legend.append("· unit circle", style="grey50")

    return Panel(
        Group(text, Align.center(legend)),
        title=f"[bold cyan]Z-PLANE :: {corner}[/bold cyan]",
        border_style="cyan",
        padding=(0, 1),
    )


def _composite_magnitude_db(cartridge: dict[str, Any], corner: str,
                            freqs_hz: list[float]) -> list[float]:
    kf = next(k for k in cartridge["keyframes"] if k["label"] == corner)
    boost = float(kf["boost"])
    out: list[float] = []
    for f in freqs_hz:
        w = 2.0 * math.pi * f / AUTHORING_SR
        e1r, e1i = math.cos(-w), math.sin(-w)
        e2r, e2i = math.cos(-2 * w), math.sin(-2 * w)
        mag_total = boost
        for st in kf["stages"]:
            c0, c1, c2, c3, c4 = (float(st[f"c{i}"]) for i in range(5))
            if c3 == 0.0 and c4 == 0.0:
                continue
            nr = c0 + c1 * e1r + c2 * e2r
            ni = c1 * e1i + c2 * e2i
            dr = 1.0 + c3 * e1r + c4 * e2r
            di = c3 * e1i + c4 * e2i
            n_mag = math.hypot(nr, ni)
            d_mag = math.hypot(dr, di)
            mag_total *= n_mag / max(d_mag, 1e-9)
        out.append(20.0 * math.log10(max(mag_total, 1e-9)))
    return out


def render_magnitude_bars(stages: list[PlannedStage], cartridge: dict[str, Any],
                          corner: str, n_bins: int = 22) -> Panel:
    """Horizontal log-frequency magnitude render with role markers."""
    freqs = [40.0 * (18000.0 / 40.0) ** (i / (n_bins - 1)) for i in range(n_bins)]
    db = _composite_magnitude_db(cartridge, corner, freqs)
    max_db = max(db)
    min_db = min(db)
    span = max(max_db - min_db, 1.0)

    # decide which freq bin each role lives in (for marker glyphs)
    role_at_bin: dict[int, str] = {}
    role_color_at_bin: dict[int, str] = {}
    frame_key = "M0" if corner.startswith("M0") else "M100"
    for stage in stages:
        if stage.role == "ANTI":
            target_hz = stage.hz_m0 if frame_key == "M0" else stage.hz_m100
        elif stage.role in ("F1", "F2", "BITE"):
            target_hz = stage.hz_m0 if frame_key == "M0" else stage.hz_m100
        elif stage.role == "LP":
            continue
        elif stage.role == "HP":
            continue
        else:
            continue
        # find nearest bin
        best_i = 0
        best_d = float("inf")
        for i, fb in enumerate(freqs):
            d = abs(math.log(fb) - math.log(target_hz))
            if d < best_d:
                best_d, best_i = d, i
        glyph, color = ROLE_GLYPH[stage.role]
        role_at_bin[best_i] = stage.role
        role_color_at_bin[best_i] = color

    text = Text()
    bar_width = 40
    for i, (f, d) in enumerate(zip(freqs, db)):
        normalised = (d - min_db) / span
        filled = int(round(normalised * bar_width))
        # marker
        marker = "│"
        marker_style = "grey50"
        if i in role_at_bin:
            role = role_at_bin[i]
            marker = ROLE_GLYPH[role][0]
            marker_style = f"bold {role_color_at_bin[i]}"

        text.append(f"{f:6.0f} Hz  ", style="grey62")
        text.append(marker, style=marker_style)
        text.append(" ")
        # gradient color across the bar
        if filled > 0:
            # cool-to-hot color depending on dB level
            level = (d - min_db) / span
            if level < 0.33:
                bar_style = "blue"
            elif level < 0.66:
                bar_style = "yellow"
            else:
                bar_style = "bright_red"
            text.append("█" * filled, style=bar_style)
        text.append("·" * (bar_width - filled), style="grey23")
        text.append(f"  {d:+6.1f} dB", style="grey70")
        if i in role_at_bin:
            text.append(f"  ◀ {role_at_bin[i]}", style=f"bold {role_color_at_bin[i]}")
        text.append("\n")

    return Panel(
        text,
        title=f"[bold cyan]MAGNITUDE :: {corner}[/bold cyan]   "
              f"[grey50](peaks/notches labelled by role)[/grey50]",
        border_style="cyan",
        padding=(0, 1),
    )


_PLOT_ROLE_COLORS = {
    # talking_formant calibration roles
    "ANCHOR":  "#cc55ff",   # foundation — magenta
    "LOWMID":  "#ff4444",   # F1-zone body — red
    "BITE_LO": "#ffaa22",   # F2 formant — orange
    "BITE_HI": "#55ff55",   # F3 / bite — green
    "AIR":     "#22ddff",   # presence — cyan
    "HF_NYQ":  "#5599ff",   # spectral ceiling — blue
    # peak_shelf
    "PEAK":    "#ffaa22",
    "RES":     "#ff4444",
    "SHELF":   "#cc55ff",
    "WIDE":    "#22ddff",
}


def plot_magnitude_curve(arc: Arc, stages: list[PlannedStage],
                         cartridge: dict[str, Any], out_path: Path) -> Path:
    """Write a matplotlib magnitude plot for the cartridge.

    Two panels (M0_Q0, M100_Q0) with composite |H(f)| in dB on log-Hz axis,
    role markers vertically lined and labelled by role colour.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    freqs = np.geomspace(40.0, 18000.0, 1024)
    fig, axes = plt.subplots(2, 1, figsize=(11.5, 7.5), dpi=150, sharex=True)
    fig.suptitle(
        f"TRENCH TRADING CO.   ::   {cartridge['name']}   "
        f"(Path {cartridge['qlaw']['policy']} — {cartridge['qlaw']['policy_name']})",
        color="#ffba00", fontweight="bold",
    )
    fig.patch.set_facecolor("#16191a")

    for ax, corner in zip(axes, ("M0_Q0", "M100_Q0")):
        ax.set_facecolor("#0d0f10")
        kf = next(k for k in cartridge["keyframes"] if k["label"] == corner)
        boost = float(kf["boost"])
        w = 2.0 * np.pi * freqs / AUTHORING_SR
        e1 = np.exp(-1j * w)
        e2 = np.exp(-2j * w)
        composite = np.ones_like(freqs)
        for st in kf["stages"]:
            c0, c1, c2, c3, c4 = (float(st[f"c{i}"]) for i in range(5))
            if c3 == 0.0 and c4 == 0.0:
                continue
            num = c0 + c1 * e1 + c2 * e2
            den = 1.0 + c3 * e1 + c4 * e2
            composite *= np.abs(num / den)
        composite_db = 20.0 * np.log10(np.maximum(composite * boost, 1e-9))
        ax.semilogx(freqs, composite_db, color="#ffba00", lw=2.0,
                    label="composite |H(f)|")

        # mark role frequencies
        frame_key = "M0" if corner.startswith("M0") else "M100"
        for stage in stages:
            target = stage.hz_m0 if frame_key == "M0" else stage.hz_m100
            color = _PLOT_ROLE_COLORS.get(stage.role, "#cccccc")
            ax.axvline(target, color=color, alpha=0.45, linestyle="--", lw=1.0)
            ax.annotate(
                stage.role,
                xy=(target, 0), xycoords=("data", "axes fraction"),
                xytext=(0, 6), textcoords="offset points",
                color=color, fontsize=9, ha="center", fontweight="bold",
            )

        ax.set_title(corner, color="#ffba00", fontweight="bold")
        ax.set_ylabel("magnitude (dB)", color="white")
        ax.grid(True, which="both", alpha=0.18, color="white")
        ax.tick_params(colors="white")
        for spine in ax.spines.values():
            spine.set_color("#444444")
        ax.legend(loc="lower left", framealpha=0.4, facecolor="#22262a",
                  edgecolor="#444444", labelcolor="white")

    axes[-1].set_xlabel("frequency (Hz)", color="white")
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, facecolor=fig.get_facecolor())
    plt.close(fig)
    return out_path


def render_morph_story(arc_label: str, stages: list[PlannedStage],
                       cartridge: dict[str, Any]) -> Panel:
    """Plain-English narrative of the cartridge stages, calibration-grounded."""
    body = Text()
    body.append("• ", style="grey50"); body.append("ARC: ", style="grey62")
    body.append(f"{arc_label}\n", style="bold white")
    body.append("• ", style="grey50"); body.append("STAGES (M0 / M100):\n", style="grey62")
    for s in stages:
        glyph, color = ROLE_GLYPH.get(s.role, ("•", "white"))
        body.append(f"    {glyph} {s.role:8s}", style=f"bold {color}")
        body.append(f" {s.hz_m0:6.0f} Hz", style="grey70")
        body.append(" → ", style="grey50")
        body.append(f"{s.hz_m100:6.0f} Hz", style="grey70")
        body.append(f"   {s.landmark}\n", style="grey50")
    body.append("• ", style="grey50"); body.append("FOUNDATION:   ", style="grey62")
    anchor = next((s for s in stages if s.role == "ANCHOR"), None)
    if anchor:
        body.append(f"ANCHOR at {anchor.hz_m0:.0f}→{anchor.hz_m100:.0f} Hz "
                    f"(sub-low body weight)\n", style="bold magenta")
    else:
        body.append("none\n", style="grey50")
    body.append("• ", style="grey50"); body.append("MORPH:        ", style="grey62")
    body.append("zeros DEEP at M0 → SOFT at M100\n", style="cyan")
    body.append("• ", style="grey50"); body.append("POLICY:       ", style="grey62")
    body.append(f"{cartridge['qlaw']['policy_name']}\n", style="bold green")
    return Panel(body, title="[bold yellow]WHAT YOU'RE HOLDING[/bold yellow]",
                 border_style="yellow", padding=(1, 2))


def render_cartridge_table(stages: list[PlannedStage], cartridge: dict[str, Any]) -> Group:
    chunks: list[Any] = []
    for kf in cartridge["keyframes"]:
        actives = [s for s in kf["stages"] if s["c3"] or s["c4"]]
        c4 = actives[0]["c4"] if actives else 0.0
        tbl = Table(
            title=f"{kf['label']}   c4={c4:.7f}   boost={kf['boost']}",
            title_style="bold cyan",
            show_lines=False, expand=False,
        )
        tbl.add_column("S", style="magenta", no_wrap=True)
        tbl.add_column("role", style="white", no_wrap=True)
        tbl.add_column("c0", justify="right")
        tbl.add_column("c1", justify="right")
        tbl.add_column("c2", justify="right")
        tbl.add_column("c3", justify="right")
        tbl.add_column("c4", justify="right")
        tbl.add_column("a1", justify="right", style="dim")
        tbl.add_column("r", justify="right", style="dim")
        tbl.add_column("val1", justify="right", style="dim")
        tbl.add_column("val2", justify="right", style="dim")
        tbl.add_column("val3", justify="right", style="dim")
        for i, st in enumerate(kf["stages"]):
            if not (st["c3"] or st["c4"]):
                continue
            role = stages[i].role if i < len(stages) else "OFF"
            h = decode_heritage(st)
            tbl.add_row(
                f"S{i}", role,
                f"{st['c0']:.4f}", f"{st['c1']:+.4f}", f"{st['c2']:.4f}",
                f"{st['c3']:+.4f}", f"{st['c4']:.4f}",
                f"{h.a1:+.4f}", f"{h.r:.4f}",
                f"{h.val1:+.4f}", f"{h.val2:+.4f}", f"{h.val3:+.4f}",
            )
        chunks.append(tbl)
    return Group(*chunks)


# ─── wizard ───────────────────────────────────────────────────────────────
def wizard(console: Console) -> int:
    console.clear()
    show_title(console, cinematic=True)
    show_operator(
        console, OPERATOR_IDLE,
        "Well howdy, partner! What kinda cartridge can I rustle up for ya today?",
    )

    # 1) arc — the arc carries its recipe; one pick, no schema gymnastics.
    console.print(render_arc_panel())
    arc_idx = IntPrompt.ask(
        Text("Pick yer route", style="bold yellow3"),
        choices=[str(i) for i in range(1, len(ARCS) + 1)],
        default="1",
    )
    arc = ARCS[arc_idx - 1]
    recipe = arc.recipe

    # 3) policy
    show_operator(console, OPERATOR_THINKING,
                  "Now then — body or laser? Both fine choices, mind ya.")
    console.print(render_policy_panel())
    policy = Prompt.ask(
        Text("Workbench mode", style="bold yellow3"),
        choices=["A", "B"], default="A",
    )

    # 4) confirm
    summary = Text()
    summary.append("recipe : ", style="grey50"); summary.append(f"{recipe}\n", style="bold yellow1")
    summary.append("route  : ", style="grey50"); summary.append(f"{arc.name}  ", style="bold yellow1")
    summary.append(f"({arc.frm} → {arc.to})\n", style="gold3")
    summary.append("mode   : ", style="grey50"); summary.append(f"{policy} ", style="bold yellow1")
    summary.append(f"— {POLICIES[policy][0]}\n", style="gold3")
    console.print(Panel(summary, title="[bold gold3]RING 'ER UP?[/bold gold3]",
                        border_style="gold3", padding=(1, 2)))
    if not Confirm.ask(Text("Forge it", style="bold yellow3"), default=True):
        console.print(Text("Suit yerself, partner. Come back any time.", style="yellow3"))
        return 1

    # 5) forge animation + actual work
    show_operator(console, OPERATOR_FORGING,
                  "Hammerin' on the anvil... this'll just take a sec...")
    render_anvil_animation(console, total_ticks=20)

    frame_map = arc_to_frame_map(arc)
    stages, meta = plan(frame_map, recipe)
    solved = solve(stages, policy)
    cartridge = build_cartridge(meta, solved, policy)
    gate_result = gate(cartridge)

    # 6) reveal
    if gate_result.passed:
        console.print(Text(PASS_BANNER, style="bold green3"))
        show_operator(console, OPERATOR_PASSED,
                      "She's a beaut, partner! Wrap her up for ya right now.",
                      style="bold green3")
        # Visual dashboard: morph story + Z-plane (M0 vs M100) + magnitude bars
        console.print(render_morph_story(
            f"{arc.name}  ({arc.blurb})", stages, cartridge,
        ))
        console.print(Columns([
            render_zplane(stages, cartridge, "M0_Q0"),
            render_zplane(stages, cartridge, "M100_Q0"),
        ], padding=(0, 1), equal=True, expand=False))
        console.print(render_magnitude_bars(stages, cartridge, "M0_Q0"))
        console.print(render_magnitude_bars(stages, cartridge, "M100_Q0"))

        out_path = DEFAULT_OUT_DIR / f"{arc.name}_path{policy}.json"
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(json.dumps(cartridge, indent=2) + "\n", encoding="utf-8")

        plot_path = out_path.with_suffix(".magnitude.png")
        try:
            plot_magnitude_curve(arc, stages, cartridge, plot_path)
            plot_line = f"   plot:     {plot_path}"
        except Exception as exc:  # matplotlib import or write failure
            plot_path = None
            plot_line = f"   plot:     [skipped — {exc}]"

        wrap = Text()
        wrap.append(f"   cart:     {out_path}\n", style="bold yellow1")
        wrap.append(plot_line + "\n", style="bold yellow1" if plot_path else "yellow3")
        wrap.append("\n   wrapped up nice. take a look at the curve, partner —",
                    style="green3")
        wrap.append("\n   if she sings right to ya, she's yours.", style="green3")
        console.print(Panel(wrap, title="[bold green3]WRAPPED & READY[/bold green3]",
                            border_style="green3", padding=(1, 2)))

        # Pop the plot open with the OS default viewer so the user can audition.
        if plot_path is not None and sys.platform == "win32":
            try:
                import os
                os.startfile(str(plot_path))  # type: ignore[attr-defined]
            except Exception:
                pass
        return 0

    console.print(Text(FAIL_BANNER, style="bold red"))
    show_operator(console, OPERATOR_REJECTED,
                  "Whoo-ee, math went sideways on me. Back to the bench.",
                  style="bold red")
    body = Text("\n".join(f"  ✗  {f}" for f in gate_result.failures), style="red")
    console.print(Panel(body, title="[bold red]GATE FAIL — NO SALE[/bold red]",
                        border_style="red"))
    return 3


# ─── one-shot subcommands ─────────────────────────────────────────────────
def cmd_throw(console: Console, arc_name: str | None, recipe: str | None, policy: str,
              boost: float, out: Path | None) -> int:
    arc = find_arc(arc_name) if arc_name else ARCS[0]
    # Recipe is bound to the arc by default; --recipe only overrides if given
    recipe = recipe or arc.recipe
    frame_map = arc_to_frame_map(arc)
    stages, meta = plan(frame_map, recipe)
    solved = solve(stages, policy, boost=boost)
    cartridge = build_cartridge(meta, solved, policy)
    gate_result = gate(cartridge)

    show_title(console)
    console.print(render_cartridge_table(stages, cartridge))
    if not gate_result.passed:
        body = Text("\n".join(f"  ✗  {f}" for f in gate_result.failures), style="red")
        console.print(Panel(body, title="[bold red]GATE FAIL[/bold red]", border_style="red"))
        return 3

    out = out or (DEFAULT_OUT_DIR / f"{arc.name}_path{policy}.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(cartridge, indent=2) + "\n", encoding="utf-8")

    plot_path = out.with_suffix(".magnitude.png")
    try:
        plot_magnitude_curve(arc, stages, cartridge, plot_path)
    except Exception:
        plot_path = None

    body = Text()
    body.append(f"PASS — {out}\n", style="bold green3")
    if plot_path:
        body.append(f"plot — {plot_path}", style="bold green3")
    console.print(Panel(body, border_style="green3"))
    return 0


def cmd_list(console: Console) -> int:
    show_title(console)
    console.print(render_recipe_panel())
    console.print(render_arc_panel())
    console.print(render_policy_panel())
    return 0


def cmd_audit(console: Console, path: Path) -> int:
    try:
        cart = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        console.print(Panel(Text(f"read error: {exc}", style="red"),
                            title="[bold red]forge_cli audit[/bold red]"))
        return 1
    result = gate(cart)
    if result.passed:
        console.print(Panel(Text(f"GATE PASS  —  {path}", style="bold green"),
                            border_style="green"))
        return 0
    body = Text("\n".join(f"  ✗  {f}" for f in result.failures), style="red")
    console.print(Panel(body, title=f"[bold red]GATE FAIL — {path}[/bold red]",
                        border_style="red"))
    return 3


# ─── entry ────────────────────────────────────────────────────────────────
def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(
        prog="forge_cli",
        description="TRENCH spectral forge — interactive cartridge press.",
    )
    sub = parser.add_subparsers(dest="cmd")

    sub.add_parser("wizard", help="interactive wizard (default if no args)")

    p_throw = sub.add_parser("throw", help="forge an arc with sane defaults")
    p_throw.add_argument("arc", nargs="?", default=None,
                         help=f"arc name (default {ARCS[0].name})")
    p_throw.add_argument("--recipe", default=None, choices=sorted(RECIPES),
                         help="override the arc's bound recipe (rare)")
    p_throw.add_argument("--policy", default="A", choices=["A", "B"])
    p_throw.add_argument("--boost", type=float, default=1.0)
    p_throw.add_argument("--out", type=Path, default=None)

    sub.add_parser("list", help="show recipes / arcs / policies")

    p_audit = sub.add_parser("audit", help="re-run the gate on an existing cartridge")
    p_audit.add_argument("path", type=Path)

    args = parser.parse_args(argv)

    # Windows legacy console can't render unicode box chars in cp1252;
    # force UTF-8 on stdout/stderr before rich writes anything.
    for stream_name in ("stdout", "stderr"):
        stream = getattr(sys, stream_name, None)
        if stream is not None and hasattr(stream, "reconfigure"):
            try:
                stream.reconfigure(encoding="utf-8")
            except Exception:
                pass

    console = Console(soft_wrap=False, force_terminal=True, legacy_windows=False)

    try:
        if args.cmd in (None, "wizard"):
            return wizard(console)
        if args.cmd == "throw":
            return cmd_throw(console, args.arc, args.recipe, args.policy,
                             args.boost, args.out)
        if args.cmd == "list":
            return cmd_list(console)
        if args.cmd == "audit":
            return cmd_audit(console, args.path)
    except ForgeError as exc:
        console.print(Panel(Text(f"forge error: {exc}", style="red"),
                            border_style="red"))
        return 2
    except KeyboardInterrupt:
        console.print(Text("\ninterrupted.", style="yellow"))
        return 130
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
