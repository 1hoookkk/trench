"""Body authoring session model.

A session tracks the authoring of one body through checkpoints,
takes, and gate results. Checkpoints are morph landmarks (not stages).
The compiler maps checkpoints to stage parameters.

Session file: vault/sessions/{name}.session.json
"""
from __future__ import annotations

import json
import math
import os
import time
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Optional

from pyruntime.body import Body
from pyruntime.constants import SR, TWO_PI
from pyruntime.corner import CornerArray, CornerState
from pyruntime.stage_math import resonator, resonator_with_zero
from pyruntime.stage_params import StageParams
from pyruntime.encode import raw_to_encoded
from pyruntime.freq_response import cascade_response_db, freq_points
from pyruntime.analysis import (
    midpoint_audit, dense_midpoint_audit,
    morph_trajectory_distance_normalized, shipping_gate,
)

VAULT_DIR = Path(__file__).parent.parent / "vault"
SESSION_DIR = VAULT_DIR / "sessions"
LIVE_PATH = Path(__file__).parent.parent / "trench_live.json"
NUM_STAGES = 6


# ── Checkpoint ──────────────────────────────────────────────────────────

@dataclass
class Checkpoint:
    """One morph landmark in the body arc.

    Not a stage. A checkpoint is a named beat in the morph timeline
    with perceptual directives that compile down to stage states.
    """
    name: str               # e.g. "Vault", "Choke", "Cry"
    morph: float             # 0.0 - 1.0 position in the sweep
    # Perceptual controls (each -1.0 to +1.0, 0 = neutral)
    brightness: float = 0.0  # negative=darker, positive=brighter
    intensity: float = 0.0   # negative=weaker, positive=stronger
    width: float = 0.0       # negative=tighter, positive=wider
    violence: float = 0.0    # negative=safer, positive=more violent
    # Anchor frequencies — the spectral landmarks at this checkpoint
    anchors: list[float] = field(default_factory=list)  # Hz values


@dataclass
class Take:
    """One saved version of the body."""
    version: int
    timestamp: float
    note: str = ""
    checkpoints: list[Checkpoint] = field(default_factory=list)
    gates: dict = field(default_factory=dict)
    body_path: str = ""


@dataclass
class Session:
    """Full authoring session for one body."""
    name: str
    world: str = ""
    brief: str = ""
    current_checkpoints: list[Checkpoint] = field(default_factory=list)
    takes: list[Take] = field(default_factory=list)
    active_take: int = 0

    def save(self):
        SESSION_DIR.mkdir(parents=True, exist_ok=True)
        path = SESSION_DIR / f"{self.name}.session.json"
        with open(path, "w") as f:
            json.dump(asdict(self), f, indent=2)
        return str(path)

    @staticmethod
    def load(name: str) -> Session:
        path = SESSION_DIR / f"{name}.session.json"
        with open(path) as f:
            d = json.load(f)
        s = Session(name=d["name"], world=d.get("world", ""),
                    brief=d.get("brief", ""))
        s.current_checkpoints = [Checkpoint(**c) for c in d.get("current_checkpoints", [])]
        s.takes = []
        for t in d.get("takes", []):
            cps = [Checkpoint(**c) for c in t.get("checkpoints", [])]
            s.takes.append(Take(version=t["version"], timestamp=t["timestamp"],
                                note=t.get("note", ""), checkpoints=cps,
                                gates=t.get("gates", {}), body_path=t.get("body_path", "")))
        s.active_take = d.get("active_take", 0)
        return s

    @staticmethod
    def list_sessions() -> list[str]:
        SESSION_DIR.mkdir(parents=True, exist_ok=True)
        return [f.stem.replace(".session", "")
                for f in SESSION_DIR.glob("*.session.json")]


# ── Checkpoint compiler ─────────────────────────────────────────────────

def _default_anchors(n: int = 6) -> list[float]:
    """Log-spaced default anchor frequencies."""
    return [80 * (2 ** (i * 1.2)) for i in range(n)]


def _compile_checkpoint_state(cp: Checkpoint, n_stages: int = NUM_STAGES) -> list[StageParams]:
    """Compile one checkpoint's directives into stage parameters.

    Maps perceptual controls to pole/zero placement:
    - brightness shifts anchor frequencies
    - intensity controls radius
    - width controls gain (val1)
    - violence adds zeros and pushes radius higher
    """
    anchors = cp.anchors if cp.anchors else _default_anchors(n_stages)
    # Pad or trim to n_stages
    while len(anchors) < n_stages:
        anchors.append(anchors[-1] * 1.5 if anchors else 1000.0)
    anchors = anchors[:n_stages]

    stages = []
    for i, base_hz in enumerate(anchors):
        # Brightness: shift frequency
        semitone_shift = cp.brightness * 6.0  # ±6 semitones at full
        hz = base_hz * (2.0 ** (semitone_shift / 12.0))
        hz = max(20.0, min(SR * 0.48, hz))

        # Intensity: control radius
        base_r = 0.88
        r = base_r + cp.intensity * 0.08  # 0.80 to 0.96
        r = max(0.30, min(0.996, r))

        # Width: control gain (val1)
        base_gain = -0.45
        val1 = base_gain + cp.width * 0.25  # -0.70 to -0.20
        val1 = max(-0.95, min(-0.05, val1))

        # Violence: add zero if > 0.3, push radius if > 0.5
        if cp.violence > 0.3:
            zero_hz = hz * (1.3 + cp.violence * 0.5)
            zero_hz = max(20.0, min(SR * 0.48, zero_hz))
            zero_r = 0.5 + cp.violence * 0.3
            if cp.violence > 0.5:
                r = min(0.996, r + (cp.violence - 0.5) * 0.06)
            stages.append(resonator_with_zero(hz, r, val1, zero_hz, zero_r))
        else:
            stages.append(resonator(hz, r, val1))

    return stages


def compile_session_body(checkpoints: list[Checkpoint]) -> Body:
    """Compile checkpoints into a body.

    Uses the first checkpoint as Frame A (morph=0)
    and the last checkpoint as Frame B (morph=1).
    Intermediate checkpoints shape the trajectory but
    the 4-corner system only stores endpoints.
    """
    if not checkpoints:
        raise ValueError("No checkpoints")

    sorted_cps = sorted(checkpoints, key=lambda c: c.morph)
    cp_a = sorted_cps[0]   # closest to morph=0
    cp_b = sorted_cps[-1]  # closest to morph=1

    stages_a = _compile_checkpoint_state(cp_a)
    stages_b = _compile_checkpoint_state(cp_b)

    # Pad to 6
    while len(stages_a) < NUM_STAGES:
        stages_a.append(StageParams.passthrough())
    while len(stages_b) < NUM_STAGES:
        stages_b.append(StageParams.passthrough())

    corner_a = CornerState(stages=stages_a, boost=4.0)
    corner_c = CornerState(stages=stages_b, boost=4.0)
    # Q collapsed for now (A=B, C=D)
    corners = CornerArray(a=corner_a, b=corner_a, c=corner_c, d=corner_c)
    return Body(name="session", corners=corners, boost=4.0)


def run_gates(body: Body, body_name: str = "") -> dict:
    """Run all safety gates on a body."""
    mid = midpoint_audit(body)
    dense = dense_midpoint_audit(body, peak_limit_db=50.0, n_morph=21, n_q=3)
    norm = morph_trajectory_distance_normalized(body)

    gates = {
        "midpoint": {"passed": mid["passed"], "worst_peak_db": mid["worst_peak_db"]},
        "dense": {"passed": dense["passed"], "worst_peak_db": dense["worst_peak_db"]},
        "morph_distance": norm,
    }

    if body_name:
        sg = shipping_gate(body, body_name)
        gates["shipping"] = {"passed": sg["passed"], "failures": sg["failures"],
                             "measurements": sg["measurements"]}

    return gates


def save_take(session: Session, note: str = "") -> Take:
    """Compile current checkpoints, run gates, save as a new take."""
    body = compile_session_body(session.current_checkpoints)
    gates = run_gates(body, session.name)

    version = len(session.takes) + 1
    body_path = str(VAULT_DIR / f"{session.name}_v{version}.json")
    Path(body_path).parent.mkdir(parents=True, exist_ok=True)
    with open(body_path, "w") as f:
        f.write(body.to_compiled_json(provenance=f"session-v{version}"))

    # Write to live for instant audition
    with open(LIVE_PATH, "w") as f:
        f.write(body.to_compiled_json(provenance="session-live"))

    take = Take(
        version=version,
        timestamp=time.time(),
        note=note,
        checkpoints=list(session.current_checkpoints),
        gates=gates,
        body_path=body_path,
    )
    session.takes.append(take)
    session.active_take = version
    session.save()
    return take


def push_live(session: Session):
    """Compile current checkpoints and write to trench_live.json."""
    body = compile_session_body(session.current_checkpoints)
    with open(LIVE_PATH, "w") as f:
        f.write(body.to_compiled_json(provenance="session-live"))
    return body
