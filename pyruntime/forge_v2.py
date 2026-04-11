"""
forge_v2.py — Semantic Blueprint Compiler for TRENCH

Translates physical authoring intents (clusters, crossings, dominance)
directly into 4-corner bodies via the Hz/radius path (bypasses heritage
compiler's radius ceiling). The c4 solver handles shared b0 per corner.

Authoring: zones, clusters, crossings, Q compression law.
Output: compiled-v1 JSON, plugin-loadable.
"""

import json
import math
import os
import sys
from dataclasses import dataclass, field
from typing import List, Dict, Optional

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np

from pyruntime.corner import CornerName, CornerState, CornerArray
from pyruntime.encode import EncodedCoeffs, raw_to_encoded
from pyruntime.stage_params import StageParams
from pyruntime.body import Body
from pyruntime.constants import SR, NUM_BODY_STAGES

# Import c4 solver from forge_generator
from pyruntime.forge_generator import solve_c4_surface, apply_gain_budget, cascade_peak_db


# =============================================================================
# ZONE → Hz MAPPING
# =============================================================================

ZONES_HZ = {
    "sub": 60,
    "mud": 200,
    "boxy": 500,
    "bite": 2500,
    "air": 6000,
    "ceiling": 10000,
}

# Radius by role — P2K-grade values
ROLE_RADIUS = {
    "monster": 0.999,
    "suppressor": 0.990,
    "anchor": 0.997,
    "cluster": 0.997,
    "dormant": 0.0,
}

# Q gain modulation: how role affects radius at Q=100
ROLE_Q_RADIUS_MOD = {
    "monster": 0.0,         # stays at max
    "suppressor": -0.03,    # loosens
    "anchor": 0.0,
    "cluster": 0.002,       # tightens slightly
    "dormant": 0.0,
}


# =============================================================================
# SCHEMA
# =============================================================================

@dataclass
class StageDef:
    zone: str
    hz: float
    role: str
    type: int = 1  # unused in Hz/radius path, kept for documentation


@dataclass
class ClusterDef:
    name: str
    stages: List[int]
    morph_direction: str       # 'up', 'down', 'stationary'
    morph_target_zone: str
    internal_spread: str       # 'converge', 'diverge', 'cross'


@dataclass
class CrossingDef:
    pair: List[int]
    morph_region: str          # 'early', 'mid', 'late'


@dataclass
class QDominance:
    stage: int
    role: str


@dataclass
class QLawDef:
    compression: str           # 'extreme', 'moderate', 'none'
    center_hz: float
    dominance: List[QDominance]


@dataclass
class SemanticBody:
    name: str
    sentence: str
    stages: List[StageDef]
    clusters: List[ClusterDef]
    crossings: List[CrossingDef]
    q_law: QLawDef


# =============================================================================
# COMPILER: SEMANTIC → 4-CORNER (Hz, radius)
# =============================================================================

def _clamp_hz(hz: float) -> float:
    return max(20.0, min(SR / 2.0 - 1.0, hz))


def _clamp_r(r: float) -> float:
    return max(0.0, min(0.999, r))


def _stage_to_params(hz: float, radius: float) -> StageParams:
    """Hz/radius → StageParams with zeros at origin (all-pole resonator)."""
    theta = 2.0 * math.pi * _clamp_hz(hz) / SR
    r = _clamp_r(radius)
    a1 = -2.0 * r * math.cos(theta)
    # Zeros at origin: val2 = -a1, val3 = r²
    return StageParams(a1=a1, r=r, val1=0.0, val2=-a1, val3=r * r)


class SemanticCompiler:
    def __init__(self, body: SemanticBody):
        self.body = body
        self.n = len(body.stages)
        # 4 corners: lists of (hz, radius)
        self.m0_q0 = [(0.0, 0.0)] * self.n
        self.m100_q0 = [(0.0, 0.0)] * self.n
        self.m0_q100 = [(0.0, 0.0)] * self.n
        self.m100_q100 = [(0.0, 0.0)] * self.n

    def compile(self) -> Body:
        # 1. M0_Q0 — the Anchor
        for i, stg in enumerate(self.body.stages):
            r = ROLE_RADIUS.get(stg.role, 0.990)
            self.m0_q0[i] = (stg.hz, r)

        # 2. M100_Q0 — apply clusters and crossings
        self.m100_q0 = list(self.m0_q0)  # copy

        for cluster in self.body.clusters:
            target_hz = ZONES_HZ.get(cluster.morph_target_zone, 1000)
            num = len(cluster.stages)
            for idx, stg_idx in enumerate(cluster.stages):
                hz_base, r_base = self.m0_q0[stg_idx]

                if cluster.internal_spread == "converge":
                    hz_new = target_hz
                elif cluster.internal_spread == "diverge":
                    spread = (idx - num / 2) * 400  # ±400 Hz spread
                    hz_new = target_hz + spread
                elif cluster.internal_spread == "cross":
                    spread = (num / 2 - idx) * 500  # invert order, ±500 Hz
                    hz_new = target_hz + spread
                else:
                    hz_new = target_hz

                self.m100_q0[stg_idx] = (_clamp_hz(hz_new), r_base)

        # Apply forced crossings
        for cross in self.body.crossings:
            a, b = cross.pair
            _, ra = self.m100_q0[a]
            _, rb = self.m100_q0[b]
            hz_a = self.m100_q0[a][0]
            hz_b = self.m100_q0[b][0]
            ha0 = self.m0_q0[a][0]
            hb0 = self.m0_q0[b][0]
            # If they haven't crossed, swap destinations
            if (ha0 < hb0 and hz_a <= hz_b) or (ha0 > hb0 and hz_a >= hz_b):
                self.m100_q0[a] = (hz_b, ra)
                self.m100_q0[b] = (hz_a, rb)

        # 3. Q corners — compression toward center + dominance mods
        comp = {"extreme": 0.8, "moderate": 0.4, "none": 0.0}
        comp_factor = comp.get(self.body.q_law.compression, 0.0)
        center = self.body.q_law.center_hz
        dom_map = {d.stage: d.role for d in self.body.q_law.dominance}

        for i in range(self.n):
            role = dom_map.get(i, self.body.stages[i].role)
            r_mod = ROLE_Q_RADIUS_MOD.get(role, 0.0)

            # M0_Q100
            hz0, r0 = self.m0_q0[i]
            hz_q = hz0 + (center - hz0) * comp_factor
            self.m0_q100[i] = (_clamp_hz(hz_q), _clamp_r(r0 + r_mod))

            # M100_Q100
            hz100, r100 = self.m100_q0[i]
            hz100_q = hz100 + (center - hz100) * comp_factor
            self.m100_q100[i] = (_clamp_hz(hz100_q), _clamp_r(r100 + r_mod))

        # 4. Build Body from 4 corners
        corner_data = {
            CornerName.A: self.m0_q0,
            CornerName.B: self.m0_q100,
            CornerName.C: self.m100_q0,
            CornerName.D: self.m100_q100,
        }

        corner_map = {}
        for cn, stage_list in corner_data.items():
            stages = []
            encoded = []
            for hz, r in stage_list:
                if r < 0.01:  # dormant
                    stages.append(StageParams.passthrough())
                    encoded.append(EncodedCoeffs(c0=1.0, c1=0.0, c2=0.0, c3=0.0, c4=0.0))
                else:
                    sp = _stage_to_params(hz, r)
                    stages.append(sp)
                    encoded.append(raw_to_encoded(sp))

            while len(stages) < NUM_BODY_STAGES:
                stages.append(StageParams.passthrough())
                encoded.append(EncodedCoeffs(c0=1.0, c1=0.0, c2=0.0, c3=0.0, c4=0.0))

            corner_map[cn] = CornerState(stages=stages, boost=4.0, _pre_encoded=encoded)

        ca = CornerArray(
            a=corner_map[CornerName.A],
            b=corner_map[CornerName.B],
            c=corner_map[CornerName.C],
            d=corner_map[CornerName.D],
        )
        return Body(name=self.body.name, corners=ca, boost=4.0)


# =============================================================================
# THE 3 TEST BODIES
# =============================================================================

BODIES = [
    SemanticBody(
        name="Magma Rising",
        sentence="Sub cluster erupts to bite range. Stages cross. Q compresses to mud.",
        stages=[
            StageDef("sub", 40, "cluster"),
            StageDef("sub", 60, "cluster"),
            StageDef("sub", 80, "cluster"),
            StageDef("mud", 150, "cluster"),
            StageDef("air", 8000, "suppressor"),
            StageDef("ceiling", 12000, "dormant"),
        ],
        clusters=[
            ClusterDef("magma", [0, 1, 2, 3], "up", "bite", "cross"),
        ],
        crossings=[
            CrossingDef([1, 3], "mid"),
        ],
        q_law=QLawDef(
            compression="extreme",
            center_hz=200,
            dominance=[
                QDominance(0, "monster"),
                QDominance(4, "suppressor"),
            ],
        ),
    ),
    SemanticBody(
        name="Glass Ceiling",
        sentence="HF wall drops to midrange and shatters. Q pinches to a laser.",
        stages=[
            StageDef("boxy", 600, "anchor"),
            StageDef("air", 7000, "cluster"),
            StageDef("air", 7500, "cluster"),
            StageDef("air", 8000, "cluster"),
            StageDef("air", 8500, "cluster"),
            StageDef("ceiling", 11000, "suppressor"),
        ],
        clusters=[
            ClusterDef("glass", [1, 2, 3, 4], "down", "bite", "diverge"),
        ],
        crossings=[],
        q_law=QLawDef(
            compression="extreme",
            center_hz=3000,
            dominance=[
                QDominance(1, "monster"),
                QDominance(0, "suppressor"),
            ],
        ),
    ),
    SemanticBody(
        name="Two Mouths",
        sentence="Dual formants crossing over. Q chokes them to nasal center.",
        stages=[
            StageDef("boxy", 400, "cluster"),
            StageDef("boxy", 500, "cluster"),
            StageDef("bite", 2000, "cluster"),
            StageDef("bite", 2500, "cluster"),
            StageDef("mud", 250, "anchor"),
            StageDef("air", 6000, "suppressor"),
        ],
        clusters=[
            ClusterDef("mouth_a", [0, 1], "up", "bite", "converge"),
            ClusterDef("mouth_b", [2, 3], "down", "boxy", "converge"),
        ],
        crossings=[
            CrossingDef([1, 2], "mid"),
        ],
        q_law=QLawDef(
            compression="moderate",
            center_hz=1000,
            dominance=[
                QDominance(0, "monster"),
                QDominance(2, "monster"),
                QDominance(5, "suppressor"),
            ],
        ),
    ),
]


# =============================================================================
# CLI
# =============================================================================

def main():
    out_dir = "vault/sift_queue"
    os.makedirs(out_dir, exist_ok=True)

    for i, semantic_body in enumerate(BODIES):
        compiler = SemanticCompiler(semantic_body)
        body = compiler.compile()

        # Apply c4 solver + gain budget
        body = solve_c4_surface(body, target_db=36.0)
        body = apply_gain_budget(body)

        # Report
        peaks = cascade_peak_db(body)
        from pyruntime.corner import CornerName
        c4s = {}
        for cn in CornerName:
            enc = body.corners.corner(cn).encode()
            for e in enc:
                if not (abs(e.c0 - 2.0) < 0.01 and abs(e.c4 - 1.0) < 0.01):
                    c4s[cn.json_key()] = e.c4
                    break

        safe_name = semantic_body.name.replace(" ", "_")
        path = os.path.join(out_dir, f"V2_{safe_name}.json")
        with open(path, "w") as f:
            f.write(body.to_compiled_json(provenance="forge-v2-semantic"))

        c4_str = "  ".join(f"{k}={v:.4f}" for k, v in c4s.items())
        peak_str = max(peaks.values())
        print(f"  {semantic_body.name} → {path}")
        print(f"    c4: {c4_str}")
        print(f"    peak: {peak_str:+.1f} dB")
        print(f"    \"{semantic_body.sentence}\"")


if __name__ == "__main__":
    main()
