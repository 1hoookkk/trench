"""FastAPI authoring runtime.

Active surfaces:
- `/`: unified morph designer + sift + analysis
- `/designer`: legacy four-corner designer
- `/sift`: legacy sift surface
- `/terrain`: 3D waypoint mesh capture surface
- `/author`: guided editor (direct stage authoring)
- `/phonemes`: corpus-grounded phoneme atlas browser
"""
from __future__ import annotations

import hashlib
import glob
import json
import math
import os
import random
import re
import statistics
import time
from datetime import datetime, timezone
from pathlib import Path

from fastapi import FastAPI, HTTPException, WebSocket, WebSocketDisconnect
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse
from pydantic import BaseModel, Field
from scipy.spatial import Delaunay, Voronoi, QhullError

from pyruntime.body import Body, list_vault_bodies, load_vault_body
from pyruntime.compiled_v1 import is_compiled_v1, interpolate_compiled_v1_stages
from pyruntime.constants import NUM_BODY_STAGES, SR
from pyruntime.corner import CornerArray, CornerState
from pyruntime.designer_compile import (
    compile_four_corner_to_body,
    legacy_sections_to_four_corner,
    make_four_corner_template,
)
from pyruntime.encode import EncodedCoeffs, raw_to_encoded
from pyruntime.freq_response import cascade_response_db, freq_points
from pyruntime.render import render_body, render_from_body
from pyruntime.splice import SpliceMode, splice_corners
from pyruntime.stage_params import StageParams
from pyruntime.target import (
    build_landmark_target,
    build_nasal_target,
    build_vowel_target,
    build_morph_target,
    build_composite_target,
    build_bell_target,
    build_electronic_target,
    build_target,
    build_weighted_target,
    get_bell_names,
    get_consonant_keys,
    get_electronic_keys,
    get_instrument_keys,
    get_landmark_names,
    get_moog_keys,
    get_nasal_keys,
    get_singer_keys,
    get_vowel_keys,
    load_sonic_tables,
)
from pyruntime.macro_compile import SlotState, compile_body, compile_slot_at
from pyruntime.analysis import body_profile, midpoint_audit, morph_trajectory_distance
from pyruntime.validator import validate as validate_corners


app = FastAPI(title="TRENCH Authoring Runtime")
app.add_middleware(CORSMiddleware, allow_origins=["*"], allow_methods=["*"], allow_headers=["*"])

STATIC_DIR = os.path.join(os.path.dirname(__file__), "static")
LIVE_PATH = os.path.join(os.path.dirname(__file__), "..", "trench_live.json")
VAULT_DIR = os.path.join(os.path.dirname(__file__), "..", "vault")
PHONEMES_DIR = os.path.join(VAULT_DIR, "_phonemes")
TRIAGE_DIR = os.path.join(VAULT_DIR, "_triage")
SAFE_NAME = re.compile(r"^[A-Za-z0-9_]+$")
_SAFE_P2K_ID = re.compile(r"^P2k_\\d{3}$")
_SAFE_PHONEME = re.compile(r"^PH_\\d{3}$")
VAULT_ROOT = Path(VAULT_DIR).resolve()
TERRAIN_CACHE_DIR = Path(__file__).parent / ".terrain_cache"
TERRAIN_STABILITY_CACHE = TERRAIN_CACHE_DIR / "stability.json"

PASSTHROUGH_ENC = EncodedCoeffs(c0=1.0, c1=0.0, c2=0.0, c3=0.0, c4=0.0)
CORNER_KEYS = ("M0_Q0", "M0_Q100", "M100_Q0", "M100_Q100")
_terrain_seed_cache: dict | None = None
_terrain_ws_clients: set[WebSocket] = set()


class ResponseRequest(BaseModel):
    body: dict
    morph: float = 0.5
    q: float = 0.5


class ExportRequest(BaseModel):
    body: dict
    name: str


class BakeRequest(BaseModel):
    body: dict
    provenance: str = "pyruntime"


class RenderRequest(BaseModel):
    body: dict | None = None
    morph: float = 0.5
    q: float = 0.5
    mackie_amount: float = 0.5
    erode_amount: float = 0.0
    corrode_amount: float = 0.0
    duration: float = 2.0


class DesignerSectionInput(BaseModel):
    type: int = 0
    low_freq: int = 0
    low_gain: int = 0
    high_freq: int = 0
    high_gain: int = 0


class DesignerCellInput(BaseModel):
    type: int = 0
    freq: int = 0
    gain: int = 0


class DesignerRequest(BaseModel):
    name: str = "Untitled"
    boost: float = 4.0
    corners: dict[str, list[DesignerCellInput]] = Field(default_factory=dict)
    sections: list[DesignerSectionInput] = Field(default_factory=list)


class DesignerRenderRequest(DesignerRequest):
    morph: float = 0.5
    q: float = 0.5
    mackie_amount: float = 0.5
    erode_amount: float = 0.0
    corrode_amount: float = 0.0
    duration: float = 2.0


class BatchRequest(BaseModel):
    base_sections: list = Field(default_factory=list)
    count: int = 10
    name_prefix: str = "candidate"
    boost: float = 4.0


# ── Direct pole-zero authoring ───────────────────────────────────────────

class DirectStageInput(BaseModel):
    """One stage: pole freq/radius, optional zero freq/radius."""
    pole_hz: float = 400.0
    radius: float = 0.95
    gain: float = 0.0       # val1: -1.0 (silent) to 0.0 (unity)
    zero_hz: float | None = None
    zero_radius: float | None = None

class DirectBodyRequest(BaseModel):
    """6 stages × 2 frames (A=closed, B=open). Direct Hz authoring."""
    name: str = "Untitled"
    frame_a: list[DirectStageInput] = Field(default_factory=list)
    frame_b: list[DirectStageInput] = Field(default_factory=list)
    morph: float = 0.0
    q: float = 0.0

class DirectRenderRequest(DirectBodyRequest):
    mackie_amount: float = 0.5
    duration: float = 3.0


ROLE_SPECS: list[dict] = [
    {
        "key": "anchor",
        "label": "Anchor",
        "band": (35.0, 90.0),
        "default_freq": 52.0,
        "radius_range": (0.58, 0.97),
        "gain_bias": 0.05,
        "q_push": 0.01,
        "description": "Sub/body anchor. Holds the floor.",
    },
    {
        "key": "chest",
        "label": "Chest",
        "band": (90.0, 320.0),
        "default_freq": 170.0,
        "radius_range": (0.70, 0.98),
        "gain_bias": 0.02,
        "q_push": 0.02,
        "description": "Low resonance mass and pressure.",
    },
    {
        "key": "mouth",
        "label": "Mouth",
        "band": (280.0, 1800.0),
        "default_freq": 920.0,
        "radius_range": (0.78, 0.992),
        "gain_bias": 0.0,
        "q_push": 0.03,
        "description": "Main identity and vowel center.",
    },
    {
        "key": "bite",
        "label": "Bite",
        "band": (1200.0, 4200.0),
        "default_freq": 2400.0,
        "radius_range": (0.70, 0.992),
        "gain_bias": 0.04,
        "q_push": 0.07,
        "description": "Upper-mid edge and brass bite.",
    },
    {
        "key": "air",
        "label": "Air",
        "band": (3500.0, 12000.0),
        "default_freq": 6800.0,
        "radius_range": (0.58, 0.96),
        "gain_bias": 0.06,
        "q_push": 0.08,
        "description": "Top sheen, sibilance, glass.",
    },
    {
        "key": "scar",
        "label": "Scar",
        "band": (250.0, 9000.0),
        "default_freq": 1900.0,
        "radius_range": (0.46, 0.94),
        "gain_bias": 0.0,
        "q_push": 0.12,
        "description": "Notch and betrayal lane.",
    },
]
ROLE_SPEC_BY_KEY = {spec["key"]: spec for spec in ROLE_SPECS}
_role_seed_cache: dict | None = None


class RoleLaneInput(BaseModel):
    enabled: bool = True
    freq_a: float | None = None
    freq_b: float | None = None
    amount: float = 0.6
    width: float = 0.5


class RoleBodyRequest(BaseModel):
    name: str = "Untitled"
    lanes: dict[str, RoleLaneInput] = Field(default_factory=dict)
    aggression: float = 0.45
    morph: float = 0.0
    q: float = 0.0


class RoleRenderRequest(RoleBodyRequest):
    mackie_amount: float = 0.5
    duration: float = 3.0


class TargetRequest(BaseModel):
    name: str
    source: str
    key: str


class SpliceRequest(BaseModel):
    name: str
    body_a: str
    body_b: str
    mode: str = "RestToMorphed"


class MorphTargetRequest(BaseModel):
    name: str
    source_a: str
    key_a: str
    source_b: str
    key_b: str


class CompositeTargetRequest(BaseModel):
    name: str
    slots: list[dict]


class AnalyzeRequest(BaseModel):
    body: dict


class TerrainPointRequest(BaseModel):
    waypoint_weights: dict[str, float] = Field(default_factory=dict)
    aggression: float = 0.0


class TerrainCureRequest(BaseModel):
    start: TerrainPointRequest
    end: TerrainPointRequest
    q_multiplier: float = 1.5
    name: str = "terrain_capture"


def _load_body(d: dict) -> Body:
    try:
        return Body.from_dict(d)
    except ValueError as e:
        raise HTTPException(400, f"Body parse error: {e}") from e


def _designer_template_from_request(req: DesignerRequest):
    if req.corners:
        missing = [key for key in CORNER_KEYS if key not in req.corners]
        if missing:
            raise HTTPException(400, f"Designer request missing corners: {missing}")
        corner_payload = {
            key: [
                {"type": cell.type, "freq": cell.freq, "gain": cell.gain}
                for cell in req.corners[key][:6]
            ]
            for key in CORNER_KEYS
        }
        return make_four_corner_template(req.name, corner_payload)

    if req.sections:
        tuples = [
            (section.type, section.low_freq, section.low_gain, section.high_freq, section.high_gain)
            for section in req.sections[:6]
        ]
        return legacy_sections_to_four_corner(req.name, tuples)

    return make_four_corner_template(
        req.name,
        {key: [{"type": 0, "freq": 0, "gain": 0} for _ in range(6)] for key in CORNER_KEYS},
    )


def _compile_designer_body(req: DesignerRequest) -> Body:
    template = _designer_template_from_request(req)
    return compile_four_corner_to_body(template, boost=req.boost)


def _terrain_safe_name(name: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_]+", "_", name.strip())
    cleaned = re.sub(r"_+", "_", cleaned).strip("_")
    return cleaned or "terrain_capture"


def _default_role_lane(spec: dict) -> RoleLaneInput:
    return RoleLaneInput(
        enabled=True,
        freq_a=spec["default_freq"],
        freq_b=spec["default_freq"],
        amount=0.58 if spec["key"] != "scar" else 0.30,
        width=0.42 if spec["key"] != "air" else 0.55,
    )


def _radius_to_role_width(radius: float, spec: dict) -> float:
    lo, hi = spec["radius_range"]
    if hi <= lo:
        return 0.5
    return max(0.0, min(1.0, 1.0 - ((radius - lo) / (hi - lo))))


def _raw_stage_to_role_metrics(stage: dict, spec: dict) -> dict:
    stage_params = StageParams(
        a1=float(stage["a1"]),
        r=float(stage["r"]),
        val1=float(stage["val1"]),
        val2=float(stage["val2"]),
        val3=float(stage["val3"]),
    )
    return {
        "freq_hz": max(20.0, min(18000.0, stage_params.freq_hz())),
        "amount": max(0.0, min(1.0, 1.0 + float(stage["val1"]))),
        "width": _radius_to_role_width(float(stage["r"]), spec),
        "zero_energy": math.sqrt(float(stage["val2"]) ** 2 + float(stage["val3"]) ** 2),
    }


def _pick_role_stage_metrics(stages: list[dict], spec: dict) -> dict:
    band_lo, band_hi = spec["band"]
    default_freq = spec["default_freq"]
    best: dict | None = None
    best_score = float("inf")
    for stage in stages[:6]:
        metrics = _raw_stage_to_role_metrics(stage, spec)
        freq = metrics["freq_hz"]
        dist = abs(math.log(max(freq, 20.0)) - math.log(default_freq))
        if freq < band_lo:
            dist += (band_lo - freq) / band_lo
        elif freq > band_hi:
            dist += (freq - band_hi) / band_hi
        if spec["key"] == "scar":
            dist -= metrics["zero_energy"] * 0.15
        if dist < best_score:
            best = metrics
            best_score = dist
    if best is not None:
        return best
    return {
        "freq_hz": default_freq,
        "amount": 0.55,
        "width": 0.45,
        "zero_energy": 0.0,
    }


def _build_role_seed_payload() -> dict:
    global _role_seed_cache
    if _role_seed_cache is not None:
        return _role_seed_cache

    templates: list[dict] = [
        {
            "key": "neutral_vowel",
            "label": "Neutral vowel",
            "lanes": {
                "anchor": {"enabled": True, "freq_a": 52.0, "freq_b": 52.0, "amount": 0.62, "width": 0.40},
                "chest": {"enabled": True, "freq_a": 165.0, "freq_b": 150.0, "amount": 0.54, "width": 0.42},
                "mouth": {"enabled": True, "freq_a": 700.0, "freq_b": 340.0, "amount": 0.70, "width": 0.32},
                "bite": {"enabled": True, "freq_a": 1220.0, "freq_b": 2322.0, "amount": 0.64, "width": 0.30},
                "air": {"enabled": True, "freq_a": 3400.0, "freq_b": 3657.0, "amount": 0.50, "width": 0.54},
                "scar": {"enabled": False, "freq_a": 1800.0, "freq_b": 2100.0, "amount": 0.20, "width": 0.58},
            },
        },
        {
            "key": "hostile_brass",
            "label": "Hostile brass",
            "lanes": {
                "anchor": {"enabled": True, "freq_a": 58.0, "freq_b": 60.0, "amount": 0.56, "width": 0.46},
                "chest": {"enabled": True, "freq_a": 220.0, "freq_b": 250.0, "amount": 0.50, "width": 0.44},
                "mouth": {"enabled": True, "freq_a": 740.0, "freq_b": 860.0, "amount": 0.60, "width": 0.35},
                "bite": {"enabled": True, "freq_a": 1300.0, "freq_b": 2500.0, "amount": 0.82, "width": 0.24},
                "air": {"enabled": True, "freq_a": 4200.0, "freq_b": 7600.0, "amount": 0.62, "width": 0.48},
                "scar": {"enabled": True, "freq_a": 1800.0, "freq_b": 2500.0, "amount": 0.42, "width": 0.44},
            },
        },
    ]

    p2k_templates: list[dict] = []
    for body in _get_p2k_bodies():
        corners = body.get("corners", {})
        corner_a = corners.get("M0_Q0", {}).get("stages", [])
        corner_c = corners.get("M100_Q0", {}).get("stages", [])
        if not corner_a or not corner_c:
            continue
        lanes: dict[str, dict] = {}
        for spec in ROLE_SPECS:
            metrics_a = _pick_role_stage_metrics(corner_a, spec)
            metrics_c = _pick_role_stage_metrics(corner_c, spec)
            lanes[spec["key"]] = {
                "enabled": True if spec["key"] != "scar" else (metrics_a["zero_energy"] > 0.12 or metrics_c["zero_energy"] > 0.12),
                "freq_a": round(max(spec["band"][0], min(spec["band"][1], metrics_a["freq_hz"])), 2),
                "freq_b": round(max(spec["band"][0], min(spec["band"][1], metrics_c["freq_hz"])), 2),
                "amount": round(max(0.08, min(1.0, (metrics_a["amount"] + metrics_c["amount"]) * 0.5)), 3),
                "width": round(max(0.05, min(1.0, (metrics_a["width"] + metrics_c["width"]) * 0.5)), 3),
            }
        p2k_templates.append(
            {
                "key": body.get("name", "unknown"),
                "label": body.get("name", "unknown"),
                "lanes": lanes,
            }
        )

    _role_seed_cache = {
        "roles": [
            {
                "key": spec["key"],
                "label": spec["label"],
                "band": [spec["band"][0], spec["band"][1]],
                "default_freq": spec["default_freq"],
                "description": spec["description"],
            }
            for spec in ROLE_SPECS
        ],
        "templates": templates,
        "p2k": p2k_templates,
    }
    return _role_seed_cache


def _resolve_role_lanes(raw_lanes: dict[str, RoleLaneInput]) -> dict[str, RoleLaneInput]:
    resolved: dict[str, RoleLaneInput] = {}
    for spec in ROLE_SPECS:
        lane = raw_lanes.get(spec["key"])
        if lane is None:
            resolved[spec["key"]] = _default_role_lane(spec)
            continue
        freq_a = lane.freq_a if lane.freq_a is not None else spec["default_freq"]
        freq_b = lane.freq_b if lane.freq_b is not None else freq_a
        resolved[spec["key"]] = RoleLaneInput(
            enabled=bool(lane.enabled),
            freq_a=max(spec["band"][0], min(spec["band"][1], float(freq_a))),
            freq_b=max(spec["band"][0], min(spec["band"][1], float(freq_b))),
            amount=max(0.0, min(1.0, float(lane.amount))),
            width=max(0.0, min(1.0, float(lane.width))),
        )
    return resolved


def _role_stage_input(spec: dict, lane: RoleLaneInput, frame: str, aggression: float, q_hot: float, variant: dict) -> DirectStageInput:
    freq = float(lane.freq_a if frame == "a" else lane.freq_b)
    freq = max(spec["band"][0], min(spec["band"][1], freq))
    min_r, max_r = spec["radius_range"]
    radius = max_r - lane.width * (max_r - min_r)
    radius += q_hot * aggression * spec["q_push"] * variant.get("q_scale", 1.0)
    radius = max(min_r, min(0.997, radius))

    gain = -0.92 + lane.amount * (0.80 + spec["gain_bias"])
    gain += q_hot * aggression * 0.16 * variant.get("gain_scale", 1.0)
    gain = max(-0.95, min(-0.04, gain))

    zero_hz: float | None = None
    zero_radius: float | None = None

    if spec["key"] == "scar" and lane.enabled:
        if lane.amount > 0.18 or (aggression * q_hot) > 0.10:
            zero_hz = max(spec["band"][0], min(spec["band"][1], freq * variant["scar_ratio"]))
            zero_radius = max(0.12, min(1.0, 0.15 + aggression * 0.55 + lane.amount * 0.22 + q_hot * 0.12))
    elif spec["key"] == "air" and aggression > 0.58 and q_hot > 0.5 and lane.amount > 0.34 and variant.get("air_zero", False):
        zero_hz = max(1200.0, min(spec["band"][1], freq * 0.78))
        zero_radius = max(0.12, min(0.85, 0.18 + aggression * 0.35))
    elif spec["key"] == "anchor" and aggression > 0.74 and q_hot > 0.72 and variant.get("anchor_zero", False):
        zero_hz = max(20.0, min(220.0, freq * 0.68))
        zero_radius = max(0.10, min(0.52, 0.10 + aggression * 0.28))

    return DirectStageInput(
        pole_hz=freq,
        radius=radius,
        gain=gain,
        zero_hz=zero_hz,
        zero_radius=zero_radius,
    )


def _build_role_corner(lanes: dict[str, RoleLaneInput], frame: str, aggression: float, q_hot: float, variant: dict) -> CornerState:
    compiled: list[tuple[float, StageParams]] = []
    for spec in ROLE_SPECS:
        lane = lanes[spec["key"]]
        if not lane.enabled:
            continue
        stage_input = _role_stage_input(spec, lane, frame, aggression, q_hot, variant)
        compiled.append((stage_input.pole_hz, _direct_stage(stage_input)))
    compiled.sort(key=lambda item: item[0])
    stages = [item[1] for item in compiled[:NUM_BODY_STAGES]]
    while len(stages) < NUM_BODY_STAGES:
        stages.append(StageParams.passthrough())
    return CornerState(stages=stages, boost=4.0)


def _encoded_rms_diff(enc_a, enc_b) -> float:
    freqs = freq_points(sr=SR)
    db_a = cascade_response_db(enc_a, freqs, SR).tolist()
    db_b = cascade_response_db(enc_b, freqs, SR).tolist()
    if not db_a:
        return 0.0
    total = 0.0
    for value_a, value_b in zip(db_a, db_b):
        delta = float(value_a) - float(value_b)
        total += delta * delta
    return math.sqrt(total / len(db_a))


def _role_body_metrics(body: Body) -> dict:
    audit = midpoint_audit(body)
    morph = morph_trajectory_distance(body)
    q_diff = _encoded_rms_diff(body.corners.interpolate(0.5, 0.0), body.corners.interpolate(0.5, 1.0))
    floor_freqs = freq_points(sr=SR)
    floor_db = cascade_response_db(body.corners.interpolate(0.0, 0.0), floor_freqs, SR).tolist()
    sub_peak = max((float(db) for freq, db in zip(floor_freqs.tolist(), floor_db) if 25.0 <= freq <= 90.0), default=-120.0)
    score = morph["mean_rms_db"] * 11.0 + q_diff * 9.0 + sub_peak * 0.15
    if audit["passed"]:
        score += 60.0
    else:
        score -= max(0.0, float(audit["worst_peak_db"]) - 55.0) * 2.0
    return {
        "audit": audit,
        "morph_distance": morph,
        "q_diff_rms_db": q_diff,
        "sub_peak_db": sub_peak,
        "score": score,
    }


def _compile_role_body(req: RoleBodyRequest) -> tuple[Body, dict]:
    lanes = _resolve_role_lanes(req.lanes)
    aggression = max(0.0, min(1.0, float(req.aggression)))
    variants = [
        {"scar_ratio": 0.78, "q_scale": 1.0, "gain_scale": 1.0, "air_zero": False, "anchor_zero": False},
        {"scar_ratio": 0.68, "q_scale": 1.08, "gain_scale": 1.03, "air_zero": False, "anchor_zero": False},
        {"scar_ratio": 0.84, "q_scale": 1.12, "gain_scale": 1.06, "air_zero": True, "anchor_zero": False},
        {"scar_ratio": 0.73, "q_scale": 1.18, "gain_scale": 1.10, "air_zero": False, "anchor_zero": True},
    ]
    best_body: Body | None = None
    best_metrics: dict | None = None
    for variant in variants:
        body = Body(
            name=req.name,
            corners=CornerArray(
                a=_build_role_corner(lanes, "a", aggression, 0.0, variant),
                b=_build_role_corner(lanes, "a", aggression, 1.0, variant),
                c=_build_role_corner(lanes, "b", aggression, 0.0, variant),
                d=_build_role_corner(lanes, "b", aggression, 1.0, variant),
            ),
            boost=4.0,
        )
        metrics = _role_body_metrics(body)
        if best_metrics is None or metrics["score"] > best_metrics["score"]:
            best_body = body
            best_metrics = metrics

    assert best_body is not None and best_metrics is not None
    return best_body, best_metrics


def _terrain_label(entry: dict, fallback: str) -> str:
    if "label" in entry:
        return str(entry["label"])
    if "name" in entry:
        return str(entry["name"])
    return fallback


def _terrain_formants(source: str, key: str, entry: dict) -> list[float]:
    if source == "vowel":
        return [float(entry[fk]) for fk in ("f1", "f2", "f3", "f4") if entry.get(fk) is not None]
    if source == "nasal":
        return [float(entry[fk]) for fk in ("f1", "f2", "anti1", "anti2") if entry.get(fk) is not None]
    if source == "consonant":
        return [float(formant["freq"]) for formant in entry.get("formants", [])]
    if source == "instrument":
        return [float(peak["freq"]) for peak in entry.get("peaks", [])]
    if source == "moog":
        return [float(entry["fc_default"])]
    if source == "singer":
        return [float(entry["f3"]), float(entry["f4"]), float(entry["f5"])]
    return [float(key)]


def _terrain_waypoint_records() -> list[dict]:
    tables = load_sonic_tables()
    keys = tables.get("_meta", {}).get("waypoint_keys", [])
    source_order = [
        ("vowel", "vowels"),
        ("nasal", "nasals"),
        ("consonant", "consonants"),
        ("instrument", "instruments"),
        ("moog", "moog"),
        ("singer", "singer"),
    ]
    records: list[dict] = []

    for key in keys:
        found = False
        for source, table_key in source_order:
            table = tables.get(table_key, {})
            if key not in table:
                continue
            entry = table[key]
            formants = _terrain_formants(source, key, entry)
            records.append(
                {
                    "key": key,
                    "source": source,
                    "category": source,
                    "label": _terrain_label(entry, key.upper()),
                    "formants": formants,
                    "dominant_freq": float(formants[0]) if formants else 1000.0,
                    "secondary_freq": float(formants[1]) if len(formants) > 1 else float(formants[0]) if formants else 1000.0,
                }
            )
            found = True
            break
        if not found:
            raise HTTPException(500, f"Terrain waypoint '{key}' is not present in sonic tables")

    return records


def _clamp01(value: float) -> float:
    return max(0.0, min(1.0, value))


def _log_norm(value: float, lo: float, hi: float) -> float:
    if lo <= 0.0 or hi <= 0.0 or hi <= lo:
        return 0.5
    value = max(lo, min(hi, value))
    return (math.log(value) - math.log(lo)) / (math.log(hi) - math.log(lo))


def _position_waypoints(records: list[dict]) -> None:
    vowels = [record for record in records if record["category"] == "vowel"]
    all_dom = [record["dominant_freq"] for record in records]
    all_second = [record["secondary_freq"] for record in records]
    dom_lo = min(all_dom) if all_dom else 80.0
    dom_hi = max(all_dom) if all_dom else 8000.0
    sec_lo = min(all_second) if all_second else 120.0
    sec_hi = max(all_second) if all_second else 12000.0

    vowel_f1_lo = min(record["formants"][0] for record in vowels) if vowels else 250.0
    vowel_f1_hi = max(record["formants"][0] for record in vowels) if vowels else 900.0
    vowel_f2_lo = min(record["formants"][1] for record in vowels if len(record["formants"]) > 1) if vowels else 800.0
    vowel_f2_hi = max(record["formants"][1] for record in vowels if len(record["formants"]) > 1) if vowels else 2500.0

    for record in records:
        category = record["category"]
        dominant = record["dominant_freq"]
        secondary = record["secondary_freq"]

        if category == "vowel" and len(record["formants"]) >= 2:
            f1 = float(record["formants"][0])
            f2 = float(record["formants"][1])
            x = 0.12 + 0.76 * _log_norm(f2, vowel_f2_lo, vowel_f2_hi)
            z = 0.18 + 0.60 * (1.0 - _log_norm(f1, vowel_f1_lo, vowel_f1_hi))
        elif category == "nasal":
            x = 0.10 + 0.55 * _log_norm(secondary, sec_lo, sec_hi)
            z = 0.10 + 0.16 * _log_norm(dominant, dom_lo, dom_hi)
        elif category == "consonant":
            centroid = sum(record["formants"]) / max(1, len(record["formants"]))
            x = 0.36 + 0.52 * _log_norm(centroid, sec_lo, sec_hi)
            z = 0.68 + 0.18 * _log_norm(dominant, dom_lo, sec_hi)
        elif category == "instrument":
            x = 0.22 + 0.58 * _log_norm(dominant, dom_lo, sec_hi)
            z = 0.40 + 0.24 * _log_norm(secondary, sec_lo, sec_hi)
        elif category == "moog":
            x = 0.58 + 0.22 * _log_norm(dominant, dom_lo, sec_hi)
            z = 0.10 + 0.05 * _log_norm(dominant, dom_lo, sec_hi)
        elif category == "singer":
            cluster = sum(record["formants"]) / max(1, len(record["formants"]))
            x = 0.44 + 0.18 * _log_norm(cluster, sec_lo, sec_hi)
            z = 0.82
        else:
            x = 0.15 + 0.70 * _log_norm(dominant, dom_lo, sec_hi)
            z = 0.15 + 0.70 * _log_norm(secondary, sec_lo, sec_hi)

        record["x"] = _clamp01(x)
        record["z"] = _clamp01(z)

    # Relax overlap without destroying the broad biome layout.
    for _ in range(28):
        for i in range(len(records)):
            for j in range(i + 1, len(records)):
                dx = records[j]["x"] - records[i]["x"]
                dz = records[j]["z"] - records[i]["z"]
                dist = math.hypot(dx, dz)
                if dist <= 1e-6 or dist >= 0.07:
                    continue
                push = (0.07 - dist) * 0.30
                nx = dx / dist
                nz = dz / dist
                records[i]["x"] = _clamp01(records[i]["x"] - nx * push)
                records[i]["z"] = _clamp01(records[i]["z"] - nz * push)
                records[j]["x"] = _clamp01(records[j]["x"] + nx * push)
                records[j]["z"] = _clamp01(records[j]["z"] + nz * push)
        for record in records:
            record["x"] = max(0.06, min(0.94, record["x"]))
            record["z"] = max(0.06, min(0.94, record["z"]))


def _voronoi_edges(points: list[list[float]]) -> list[list[list[float]]]:
    if len(points) < 3:
        return []
    try:
        vor = Voronoi(points)
    except QhullError:
        return []

    edges: list[list[list[float]]] = []
    for ridge in vor.ridge_vertices:
        if len(ridge) != 2 or -1 in ridge:
            continue
        a = vor.vertices[ridge[0]]
        b = vor.vertices[ridge[1]]
        edges.append(
            [
                [_clamp01(float(a[0])), _clamp01(float(a[1]))],
                [_clamp01(float(b[0])), _clamp01(float(b[1]))],
            ]
        )
    return edges


def _weights_for_point(seed: dict, x: float, z: float) -> dict[str, float]:
    points = seed["points"]
    delaunay: Delaunay | None = seed["delaunay"]
    coords = [x, z]

    if delaunay is not None:
        simplex = int(delaunay.find_simplex([coords])[0])
        if simplex >= 0:
            simplex_vertices = delaunay.simplices[simplex]
            transform = delaunay.transform[simplex]
            delta = [coords[0] - transform[2][0], coords[1] - transform[2][1]]
            bary = [
                transform[0][0] * delta[0] + transform[0][1] * delta[1],
                transform[1][0] * delta[0] + transform[1][1] * delta[1],
            ]
            bary.append(1.0 - bary[0] - bary[1])
            weights: dict[str, float] = {}
            for local_idx, vertex_idx in enumerate(simplex_vertices):
                weights[seed["waypoints"][int(vertex_idx)]["key"]] = max(0.0, float(bary[local_idx]))
            total = sum(weights.values())
            if total > 0.0:
                return {key: value / total for key, value in weights.items()}

    # Fallback outside the convex hull: inverse-distance blend of 3 nearest points.
    dists = []
    for index, point in enumerate(points):
        dx = x - point[0]
        dz = z - point[1]
        dist = math.hypot(dx, dz)
        dists.append((dist, index))
    dists.sort(key=lambda item: item[0])
    nearest = dists[:3]
    if nearest and nearest[0][0] <= 1e-8:
        waypoint = seed["waypoints"][nearest[0][1]]
        return {waypoint["key"]: 1.0}
    inv_total = sum(1.0 / max(item[0], 1e-5) for item in nearest)
    return {
        seed["waypoints"][index]["key"]: (1.0 / max(dist, 1e-5)) / inv_total
        for dist, index in nearest
    }


def _terrain_sources_from_weights(weights: dict[str, float], seed: dict) -> list[dict]:
    source_map = {waypoint["key"]: waypoint for waypoint in seed["waypoints"]}
    sources: list[dict] = []
    for key, weight in weights.items():
        if weight <= 0.0 or key not in source_map:
            continue
        waypoint = source_map[key]
        sources.append({"source": waypoint["source"], "key": key, "weight": float(weight)})
    sources.sort(key=lambda entry: entry["weight"], reverse=True)
    return sources


def _slot_state_at_morph(slot, morph: float) -> SlotState:
    contour = slot.state_a.contour if morph < 0.5 else slot.state_b.contour
    contour_override = slot.state_a.contour_override if morph < 0.5 else slot.state_b.contour_override
    return SlotState(
        place=slot.state_a.place + morph * (slot.state_b.place - slot.state_a.place),
        focus=slot.state_a.focus + morph * (slot.state_b.focus - slot.state_a.focus),
        weight=slot.state_a.weight + morph * (slot.state_b.weight - slot.state_a.weight),
        contour=contour,
        color=slot.state_a.color + morph * (slot.state_b.color - slot.state_a.color),
        hue=slot.state_a.hue + morph * (slot.state_b.hue - slot.state_a.hue),
        contour_override=contour_override,
    )


def _compile_point_from_spec(spec, pressure: float, morph: float = 0.5) -> CornerState:
    stages = [StageParams.passthrough()] * NUM_BODY_STAGES
    stage_idx = 0
    for slot in spec.slots:
        state = _slot_state_at_morph(slot, morph)
        compiled = compile_slot_at(slot, state, _clamp01(pressure))
        if stage_idx < NUM_BODY_STAGES:
            stages[stage_idx] = compiled.primary
            stage_idx += 1
        if compiled.secondary is not None and stage_idx < NUM_BODY_STAGES:
            stages[stage_idx] = compiled.secondary
            stage_idx += 1
    return CornerState(stages=stages, boost=spec.boost)


def _terrain_body_from_points(
    name: str,
    start_sources: list[dict],
    start_aggression: float,
    end_sources: list[dict],
    end_aggression: float,
    q_multiplier: float,
) -> Body:
    start_floor_spec = build_weighted_target(f"{name}_start_floor", start_sources, aggression=start_aggression)
    end_floor_spec = build_weighted_target(f"{name}_end_floor", end_sources, aggression=end_aggression)
    start_roof_aggr = _clamp01(start_aggression * q_multiplier)
    end_roof_aggr = _clamp01(end_aggression * q_multiplier)
    start_roof_spec = build_weighted_target(f"{name}_start_roof", start_sources, aggression=start_roof_aggr)
    end_roof_spec = build_weighted_target(f"{name}_end_roof", end_sources, aggression=end_roof_aggr)

    corners = CornerArray(
        a=_compile_point_from_spec(start_floor_spec, pressure=start_aggression),
        b=_compile_point_from_spec(start_roof_spec, pressure=start_roof_aggr),
        c=_compile_point_from_spec(end_floor_spec, pressure=end_aggression),
        d=_compile_point_from_spec(end_roof_spec, pressure=end_roof_aggr),
    )
    boost = (start_floor_spec.boost + end_floor_spec.boost + start_roof_spec.boost + end_roof_spec.boost) / 4.0
    return Body(name=name, corners=corners, boost=boost)


def _stability_score_from_audit(audit: dict) -> float:
    peak = float(audit.get("worst_peak_db", 120.0))
    if peak <= 35.0:
        return 1.0
    if peak >= 70.0:
        return 0.0
    return _clamp01(1.0 - ((peak - 35.0) / 35.0))


def _vault_positions(seed: dict) -> list[dict]:
    freqs = freq_points(sr=SR)
    dom_lo = min(waypoint["dominant_freq"] for waypoint in seed["waypoints"])
    dom_hi = max(waypoint["dominant_freq"] for waypoint in seed["waypoints"])
    bodies: list[dict] = []
    for name in list_vault_bodies()[:64]:
        try:
            body = load_vault_body(name)
            enc = body.corners.interpolate(0.0, 0.0)
            db = cascade_response_db(enc, freqs, SR)
            peak_index = int(db.argmax())
            peak_freq = float(freqs[peak_index])
            profile = body_profile(body)
            x = 0.10 + 0.80 * _log_norm(peak_freq, dom_lo, dom_hi)
            z = 0.18 + 0.60 * _clamp01(profile["morph_distance"]["mean_rms_db"] / 18.0)
            bodies.append({"name": name, "x": round(x, 4), "z": round(z, 4)})
        except Exception:
            continue
    return bodies


def _build_stability_map(seed: dict, grid_w: int = 24, grid_h: int = 14) -> dict:
    TERRAIN_CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cache_key = hashlib.sha1(
        json.dumps(
            {
                "waypoints": [(waypoint["key"], waypoint["x"], waypoint["z"]) for waypoint in seed["waypoints"]],
                "grid": [grid_w, grid_h],
                "version": 2,
            },
            sort_keys=True,
        ).encode("utf-8")
    ).hexdigest()

    if TERRAIN_STABILITY_CACHE.exists():
        try:
            with open(TERRAIN_STABILITY_CACHE) as handle:
                cached = json.load(handle)
            if cached.get("cache_key") == cache_key:
                return cached["map"]
        except Exception:
            pass

    data: list[float] = []
    for gy in range(grid_h):
        z = gy / max(1, grid_h - 1)
        for gx in range(grid_w):
            x = gx / max(1, grid_w - 1)
            weights = _weights_for_point(seed, x, z)
            sources = _terrain_sources_from_weights(weights, seed)
            if not sources:
                data.append(0.0)
                continue
            scores = []
            for aggression in (0.2, 0.55, 0.9):
                try:
                    probe = _terrain_body_from_points(
                        name="terrain_probe",
                        start_sources=sources,
                        start_aggression=aggression,
                        end_sources=sources,
                        end_aggression=aggression,
                        q_multiplier=1.35,
                    )
                    scores.append(_stability_score_from_audit(midpoint_audit(probe, peak_limit_db=55.0)))
                except Exception:
                    scores.append(0.0)
            data.append(round(min(scores) if scores else 0.0, 4))

    stability_map = {"grid_w": grid_w, "grid_h": grid_h, "data": data}
    with open(TERRAIN_STABILITY_CACHE, "w") as handle:
        json.dump({"cache_key": cache_key, "map": stability_map}, handle, indent=2)
    return stability_map


def _terrain_seed() -> dict:
    global _terrain_seed_cache
    if _terrain_seed_cache is not None:
        return _terrain_seed_cache

    waypoints = _terrain_waypoint_records()
    _position_waypoints(waypoints)
    points = [[waypoint["x"], waypoint["z"]] for waypoint in waypoints]
    try:
        delaunay = Delaunay(points) if len(points) >= 3 else None
    except QhullError:
        delaunay = None
    triangles = delaunay.simplices.tolist() if delaunay is not None else []

    seed = {
        "waypoints": [
            {
                **waypoint,
                "x": round(float(waypoint["x"]), 4),
                "z": round(float(waypoint["z"]), 4),
            }
            for waypoint in waypoints
        ],
        "points": points,
        "triangles": triangles,
        "voronoi_edges": _voronoi_edges(points),
        "delaunay": delaunay,
    }
    seed["stability_map"] = _build_stability_map(seed)
    seed["vault_bodies"] = _vault_positions(seed)
    _terrain_seed_cache = seed
    return seed


async def _terrain_broadcast(event: dict) -> None:
    dead: list[WebSocket] = []
    for client in list(_terrain_ws_clients):
        try:
            await client.send_json(event)
        except Exception:
            dead.append(client)
    for client in dead:
        _terrain_ws_clients.discard(client)


@app.get("/")
def serve_root():
    return FileResponse(os.path.join(STATIC_DIR, "index.html"))


@app.get("/designer")
def serve_designer():
    return FileResponse(os.path.join(STATIC_DIR, "index.html"))


@app.get("/workbench")
def serve_workbench():
    return FileResponse(os.path.join(STATIC_DIR, "workbench.html"))


@app.get("/sift")
def serve_sift():
    return FileResponse(os.path.join(STATIC_DIR, "sift.html"))


@app.get("/forge")
def serve_forge():
    return FileResponse(os.path.join(STATIC_DIR, "forge.html"))


@app.get("/terrain")
def serve_terrain():
    return FileResponse(os.path.join(STATIC_DIR, "terrain.html"))


@app.get("/phonemes")
def serve_phonemes():
    return FileResponse(os.path.join(STATIC_DIR, "phonemes.html"))


@app.get("/vaultlab")
def serve_vaultlab():
    return FileResponse(os.path.join(STATIC_DIR, "vaultlab.html"))


@app.get("/phonemes/inventory")
def phonemes_inventory():
    inv_path = Path(PHONEMES_DIR) / "p2k_phoneme_inventory.json"
    if not inv_path.is_file():
        raise HTTPException(404, f"Phoneme inventory not found: {inv_path}")
    with open(inv_path, "r", encoding="utf-8") as f:
        return json.load(f)


@app.get("/phonemes/atlas/{phoneme}")
def phonemes_atlas(phoneme: str):
    if not _SAFE_PHONEME.match(phoneme):
        raise HTTPException(400, "Invalid phoneme key (expected PH_###)")
    atlas_path = Path(PHONEMES_DIR) / "p2k_corner_atlas.jsonl"
    if not atlas_path.is_file():
        raise HTTPException(404, f"Phoneme atlas not found: {atlas_path}")

    records: list[dict] = []
    with open(atlas_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                rec = json.loads(line)
            except Exception:
                continue
            if rec.get("phoneme") == phoneme:
                records.append(rec)

    return {"phoneme": phoneme, "count": len(records), "records": records}


@app.get("/phonemes/body/{body_id}")
def phonemes_body(body_id: str):
    if not _SAFE_P2K_ID.match(body_id):
        raise HTTPException(400, "Invalid body id (expected P2k_###)")
    atlas_path = Path(PHONEMES_DIR) / "p2k_corner_atlas.jsonl"
    if not atlas_path.is_file():
        raise HTTPException(404, f"Phoneme atlas not found: {atlas_path}")

    records: list[dict] = []
    with open(atlas_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                rec = json.loads(line)
            except Exception:
                continue
            if rec.get("id") == body_id:
                records.append(rec)

    return {"id": body_id, "count": len(records), "records": records}


@app.get("/role-editor.js")
def serve_role_editor_js():
    return FileResponse(os.path.join(STATIC_DIR, "role_editor.js"))


@app.get("/role-editor.css")
def serve_role_editor_css():
    return FileResponse(os.path.join(STATIC_DIR, "role_editor.css"))


@app.get("/health")
def health():
    return {"status": "ok"}


@app.get("/live-response")
def live_response():
    if not os.path.exists(LIVE_PATH):
        return {"freqs": [], "db": [], "label": ""}
    try:
        with open(LIVE_PATH) as f:
            live_data = json.load(f)
        freqs = freq_points(sr=SR)
        if is_compiled_v1(live_data):
            enc = interpolate_compiled_v1_stages(live_data, morph=0.5, q=0.5)
        else:
            body = Body.from_dict(live_data)
            enc = body.corners.interpolate(0.5, 0.5)
        db = cascade_response_db(enc, freqs, SR)
        return {"freqs": freqs.tolist(), "db": db.tolist(), "label": live_data.get("name", "")}
    except Exception:
        return {"freqs": [], "db": [], "label": ""}


@app.post("/response")
def response_endpoint(req: ResponseRequest):
    freqs = freq_points(sr=SR)
    if is_compiled_v1(req.body):
        enc = interpolate_compiled_v1_stages(req.body, morph=req.morph, q=req.q)
    else:
        body = _load_body(req.body)
        enc = body.corners.interpolate(req.morph, req.q)
    db = cascade_response_db(enc, freqs, SR)
    return {"freqs": freqs.tolist(), "db": db.tolist()}


@app.post("/export")
def export_endpoint(req: ExportRequest):
    if not SAFE_NAME.match(req.name):
        raise HTTPException(400, "Name must be alphanumeric + underscore only")
    os.makedirs(VAULT_DIR, exist_ok=True)
    path = os.path.join(VAULT_DIR, f"{req.name}.json")
    body = _load_body(req.body)
    with open(path, "w") as f:
        f.write(body.to_compiled_json())
    return {"path": path}


@app.post("/bake")
def bake_endpoint(req: BakeRequest):
    body = _load_body(req.body)
    return json.loads(body.to_compiled_json(provenance=req.provenance))


@app.post("/render")
def render_endpoint(req: RenderRequest):
    return render_body(
        body_dict=req.body,
        morph=req.morph,
        q=req.q,
        mackie_amount=req.mackie_amount,
        erode_amount=req.erode_amount,
        corrode_amount=req.corrode_amount,
        duration=req.duration,
    )


@app.post("/designer")
def designer_endpoint(req: DesignerRequest):
    return json.loads(_compile_designer_body(req).to_json())


@app.post("/designer/render")
def designer_render_endpoint(req: DesignerRenderRequest):
    body = _compile_designer_body(req)
    return render_from_body(
        body,
        req.morph,
        req.q,
        req.mackie_amount,
        req.erode_amount,
        req.corrode_amount,
        req.duration,
    )


@app.post("/designer/response")
def designer_response_endpoint(req: DesignerRenderRequest):
    body = _compile_designer_body(req)
    enc = body.corners.interpolate(req.morph, req.q)
    freqs = freq_points(sr=SR)
    return {"freqs": freqs.tolist(), "db": cascade_response_db(enc, freqs, SR).tolist()}


@app.post("/designer/live")
def designer_live_endpoint(req: DesignerRequest):
    body = _compile_designer_body(req)
    compiled = body.to_compiled_json(provenance="designer-live")
    with open(LIVE_PATH, "w") as f:
        f.write(compiled)
    return {"status": "ok", "name": body.name}


_candidate_queue: list[dict] = []
_candidate_index: int = 0
_p2k_bodies: list[dict] | None = None


def _get_p2k_bodies() -> list[dict]:
    global _p2k_bodies
    if _p2k_bodies is not None:
        return _p2k_bodies

    bodies: list[dict] = []
    skin_dir = os.path.join(os.path.dirname(__file__), "..", "datasets", "p2k_skins")
    for path in sorted(glob.glob(os.path.join(skin_dir, "P2k_*.json"))):
        with open(path) as f:
            bodies.append(json.load(f))
    _p2k_bodies = bodies
    return bodies


@app.post("/sift/generate")
def sift_generate(req: BatchRequest):
    """Generate a candidate queue for rapid ear triage."""
    global _candidate_queue, _candidate_index
    _candidate_queue = []
    _candidate_index = 0

    bodies = _get_p2k_bodies()
    if len(bodies) < 2:
        raise HTTPException(500, "Need at least 2 P2K bodies in datasets/p2k_skins/")

    corner_labels = ["M0_Q0", "M0_Q100", "M100_Q0", "M100_Q100"]

    def body_to_corners(bdata: dict, boost: float) -> list[CornerState]:
        out: list[CornerState] = []
        for label in corner_labels:
            raw_stages = bdata["corners"][label]["stages"]
            stages: list[StageParams] = []
            pre_encoded: list[EncodedCoeffs] = []
            for stage in raw_stages[:6]:
                sp = StageParams(
                    a1=stage["a1"],
                    r=stage["r"],
                    val1=stage["val1"],
                    val2=stage["val2"],
                    val3=stage["val3"],
                )
                stages.append(sp)
                pre_encoded.append(raw_to_encoded(sp, flag=stage.get("flag", 1.0)))
            while len(stages) < NUM_BODY_STAGES:
                stages.append(StageParams.passthrough())
                pre_encoded.append(PASSTHROUGH_ENC)
            out.append(CornerState(stages=stages, boost=boost, _pre_encoded=pre_encoded))
        return out

    def perturb_corner(corner: CornerState, semitones: float) -> CornerState:
        ratio = 2.0 ** (semitones / 12.0)
        stages: list[StageParams] = []
        pre_encoded: list[EncodedCoeffs] = []
        for sp in corner.stages[:6]:
            if sp.r > 0.01:
                freq_hz = math.acos(max(-1.0, min(1.0, -sp.a1 / (2 * sp.r)))) * 39062.5 / (2 * math.pi)
                new_freq = max(20.0, min(18000.0, freq_hz * ratio))
                theta = 2.0 * math.pi * new_freq / 39062.5
                new_sp = StageParams(
                    a1=-2.0 * sp.r * math.cos(theta),
                    r=sp.r,
                    val1=sp.val1,
                    val2=sp.val2,
                    val3=sp.val3,
                )
            else:
                new_sp = sp
            stages.append(new_sp)
            pre_encoded.append(raw_to_encoded(new_sp, flag=1.0))
        while len(stages) < NUM_BODY_STAGES:
            stages.append(StageParams.passthrough())
            pre_encoded.append(PASSTHROUGH_ENC)
        return CornerState(stages=stages, boost=corner.boost, _pre_encoded=pre_encoded)

    source_pairs = [
        ("vowel", "vowel"),
        ("vowel", "nasal"),
        ("vowel", "bell"),
        ("nasal", "landmark"),
    ]
    source_key_getters = {
        "vowel": get_vowel_keys,
        "nasal": get_nasal_keys,
        "landmark": get_landmark_names,
        "bell": get_bell_names,
    }

    for i in range(req.count):
        strategy = random.choice(["splice", "perturb", "cross", "target"])
        a_body = random.choice(bodies)
        b_body = random.choice(bodies)
        name = f"{req.name_prefix}_{i:03d}"

        if strategy == "target":
            try:
                src_a, src_b = random.choice(source_pairs)
                key_a = random.choice(source_key_getters[src_a]())
                key_b = random.choice(source_key_getters[src_b]())
                spec = build_morph_target(name, src_a, key_a, src_b, key_b)
                target_corners = compile_body(spec)
                body = Body(name=name, corners=target_corners, boost=spec.boost)
                _candidate_queue.append(json.loads(body.to_compiled_json(provenance="sift-target")))
            except (ValueError, IndexError):
                pass
            continue

        if strategy == "splice":
            a_corners = body_to_corners(a_body, req.boost)
            b_corners = body_to_corners(b_body, req.boost)
            corners = [a_corners[0], a_corners[1], b_corners[2], b_corners[3]]
        elif strategy == "perturb":
            shift = random.uniform(-5.0, 5.0)
            base = body_to_corners(a_body, req.boost)
            corners = [perturb_corner(corner, shift) for corner in base]
        else:
            a_corners = body_to_corners(a_body, req.boost)
            b_corners = body_to_corners(b_body, req.boost)
            corners = []
            for idx in range(4):
                ac = a_corners[idx]
                bc = b_corners[idx]
                stages = list(ac.stages[:3]) + list(bc.stages[3:6]) + list(ac.stages[6:])
                pre = list(ac._pre_encoded[:3]) + list(bc._pre_encoded[3:6]) + list(ac._pre_encoded[6:])
                corners.append(CornerState(stages=stages, boost=req.boost, _pre_encoded=pre))

        body = Body(
            name=name,
            corners=CornerArray(a=corners[0], b=corners[1], c=corners[2], d=corners[3]),
            boost=req.boost,
        )
        _candidate_queue.append(json.loads(body.to_compiled_json(provenance=f"sift-{strategy}")))

    if _candidate_queue:
        with open(LIVE_PATH, "w") as f:
            json.dump(_candidate_queue[0], f)

    return {"count": len(_candidate_queue), "current": 0, "name": _candidate_queue[0].get("name", "") if _candidate_queue else ""}


@app.post("/sift/next")
def sift_next():
    global _candidate_index
    if not _candidate_queue:
        raise HTTPException(400, "No candidates. Call /sift/generate first.")
    _candidate_index = (_candidate_index + 1) % len(_candidate_queue)
    with open(LIVE_PATH, "w") as f:
        json.dump(_candidate_queue[_candidate_index], f)
    return {"current": _candidate_index, "count": len(_candidate_queue), "name": _candidate_queue[_candidate_index].get("name", "")}


@app.post("/sift/prev")
def sift_prev():
    global _candidate_index
    if not _candidate_queue:
        raise HTTPException(400, "No candidates. Call /sift/generate first.")
    _candidate_index = (_candidate_index - 1) % len(_candidate_queue)
    with open(LIVE_PATH, "w") as f:
        json.dump(_candidate_queue[_candidate_index], f)
    return {"current": _candidate_index, "count": len(_candidate_queue), "name": _candidate_queue[_candidate_index].get("name", "")}


@app.post("/sift/save")
def sift_save():
    if not _candidate_queue or _candidate_index >= len(_candidate_queue):
        raise HTTPException(400, "No current candidate.")
    candidate = _candidate_queue[_candidate_index]
    name = candidate.get("name", f"sift_{_candidate_index:03d}")
    os.makedirs(VAULT_DIR, exist_ok=True)
    path = os.path.join(VAULT_DIR, f"{name}.json")
    with open(path, "w") as f:
        json.dump(candidate, f, indent=2)
    return {"saved": name, "path": path}


@app.post("/sift/trash")
def sift_trash():
    global _candidate_queue, _candidate_index
    if not _candidate_queue:
        raise HTTPException(400, "No candidates.")
    _candidate_queue.pop(_candidate_index)
    if not _candidate_queue:
        return {"count": 0, "current": 0, "name": ""}
    _candidate_index = _candidate_index % len(_candidate_queue)
    with open(LIVE_PATH, "w") as f:
        json.dump(_candidate_queue[_candidate_index], f)
    return {"current": _candidate_index, "count": len(_candidate_queue), "name": _candidate_queue[_candidate_index].get("name", "")}


@app.get("/sift/status")
def sift_status():
    name = ""
    if _candidate_queue and _candidate_index < len(_candidate_queue):
        name = _candidate_queue[_candidate_index].get("name", "")
    return {"current": _candidate_index, "count": len(_candidate_queue), "name": name}


@app.get("/sift/current")
def sift_current():
    """Return the current sift candidate as a full body dict for audition."""
    if not _candidate_queue or _candidate_index >= len(_candidate_queue):
        raise HTTPException(400, "No current candidate.")
    return _candidate_queue[_candidate_index]


@app.get("/sonic-tables")
def sonic_tables():
    return {
        "vowels": get_vowel_keys(),
        "nasals": get_nasal_keys(),
        "landmarks": get_landmark_names(),
        "bells": get_bell_names(),
        "electronic": get_electronic_keys(),
        "consonants": get_consonant_keys(),
        "instruments": get_instrument_keys(),
        "moog": get_moog_keys(),
        "singer": get_singer_keys(),
    }


@app.get("/sonic-guides")
def sonic_guides():
    tables = load_sonic_tables()
    guides: list[dict] = []

    for key, entry in tables.get("vowels", {}).items():
        for formant_name in ("f1", "f2", "f3", "f4"):
            freq = entry.get(formant_name)
            if freq is None:
                continue
            guides.append({
                "category": "vowels",
                "group": key,
                "label": f"{key.upper()} {formant_name.upper()}",
                "freq_hz": float(freq),
            })

    for key, entry in tables.get("nasals", {}).items():
        for name in ("f1", "f2", "anti1", "anti2"):
            freq = entry.get(name)
            if freq is None:
                continue
            guides.append({
                "category": "nasals",
                "group": key,
                "label": f"{key.upper()} {name.upper()}",
                "freq_hz": float(freq),
            })

    for entry in tables.get("landmarks", {}).get("entries", []):
        guides.append({
            "category": "landmarks",
            "group": "landmarks",
            "label": entry["name"],
            "freq_hz": float(entry["freq_hz"]),
        })

    for key, entry in tables.get("instruments", {}).items():
        for index, peak in enumerate(entry.get("peaks", []), start=1):
            guides.append({
                "category": "instruments",
                "group": key,
                "label": f"{key} P{index}",
                "freq_hz": float(peak["freq"]),
            })

    for key, entry in tables.get("consonants", {}).items():
        for index, formant in enumerate(entry.get("formants", []), start=1):
            guides.append({
                "category": "consonants",
                "group": key,
                "label": f"{key.upper()} F{index}",
                "freq_hz": float(formant["freq"]),
            })

    for key, entry in tables.get("moog", {}).items():
        guides.append({
            "category": "moog",
            "group": key,
            "label": key,
            "freq_hz": float(entry["fc_default"]),
        })

    for key, entry in tables.get("singer", {}).items():
        for name in ("f3", "f4", "f5"):
            guides.append({
                "category": "singer",
                "group": key,
                "label": f"{key} {name.upper()}",
                "freq_hz": float(entry[name]),
            })

    guides.sort(key=lambda guide: guide["freq_hz"])
    return {"guides": guides}


@app.post("/target")
def target_endpoint(req: TargetRequest):
    try:
        spec = build_target(req.name, req.source, req.key)
    except ValueError as e:
        raise HTTPException(
            400,
            f"{e}. Use vowel, nasal, landmark, bell, electronic, consonant, instrument, moog, or singer.",
        ) from e

    corners = compile_body(spec)
    body = Body(name=req.name, corners=corners, boost=spec.boost)
    return json.loads(body.to_json())


@app.post("/morph-target")
def morph_target_endpoint(req: MorphTargetRequest):
    if not SAFE_NAME.match(req.name):
        raise HTTPException(400, "Name must be alphanumeric + underscore only")
    try:
        spec = build_morph_target(req.name, req.source_a, req.key_a, req.source_b, req.key_b)
        corners = compile_body(spec)
        body = Body(name=req.name, corners=corners, boost=spec.boost)
        return json.loads(body.to_json())
    except ValueError as e:
        raise HTTPException(400, str(e)) from e


@app.post("/composite-target")
def composite_target_endpoint(req: CompositeTargetRequest):
    if not SAFE_NAME.match(req.name):
        raise HTTPException(400, "Name must be alphanumeric + underscore only")
    try:
        spec = build_composite_target(req.name, req.slots)
        corners = compile_body(spec)
        body = Body(name=req.name, corners=corners, boost=spec.boost)
        return json.loads(body.to_json())
    except ValueError as e:
        raise HTTPException(400, str(e)) from e


@app.post("/analyze")
def analyze_endpoint(req: AnalyzeRequest):
    body = _load_body(req.body)
    return body_profile(body)


@app.get("/terrain/seed")
def terrain_seed():
    seed = _terrain_seed()
    return {
        "waypoints": seed["waypoints"],
        "triangles": seed["triangles"],
        "voronoi_edges": seed["voronoi_edges"],
        "stability_map": seed["stability_map"],
        "vault_bodies": seed["vault_bodies"],
    }


@app.post("/terrain/preview")
def terrain_preview(req: TerrainPointRequest):
    seed = _terrain_seed()
    if not req.waypoint_weights:
        raise HTTPException(400, "terrain preview requires waypoint_weights")

    sources = _terrain_sources_from_weights(req.waypoint_weights, seed)
    if not sources:
        raise HTTPException(400, "terrain preview resolved to zero valid waypoint weights")

    spec = build_weighted_target("terrain_preview", sources, aggression=req.aggression)
    point = _compile_point_from_spec(spec, pressure=req.aggression)
    enc = point.encode()
    freqs = freq_points(sr=SR)
    db = cascade_response_db(enc, freqs, SR)
    peak_db = float(max(db)) if len(db) else -120.0
    dominant = max(sources, key=lambda entry: entry["weight"])
    dominant_waypoint = next(
        waypoint for waypoint in seed["waypoints"] if waypoint["key"] == dominant["key"]
    )

    return {
        "freqs": freqs.tolist(),
        "db": db.tolist(),
        "label": dominant_waypoint["label"],
        "source": dominant_waypoint["source"],
        "aggression": _clamp01(req.aggression),
        "peak_db": round(peak_db, 2),
        "stable": peak_db <= 55.0,
        "weights": sources,
    }


@app.post("/terrain/cure")
async def terrain_cure(req: TerrainCureRequest):
    seed = _terrain_seed()
    start_sources = _terrain_sources_from_weights(req.start.waypoint_weights, seed)
    end_sources = _terrain_sources_from_weights(req.end.waypoint_weights, seed)
    if not start_sources:
        raise HTTPException(400, "terrain cure start resolved to zero waypoint weights")
    if not end_sources:
        raise HTTPException(400, "terrain cure end resolved to zero waypoint weights")

    safe_name = _terrain_safe_name(req.name)
    body = _terrain_body_from_points(
        name=safe_name,
        start_sources=start_sources,
        start_aggression=req.start.aggression,
        end_sources=end_sources,
        end_aggression=req.end.aggression,
        q_multiplier=req.q_multiplier,
    )
    audit = midpoint_audit(body, peak_limit_db=55.0)
    trajectory = morph_trajectory_distance(body)
    compiled = json.loads(body.to_compiled_json(provenance="terrain-cure"))

    written = False
    if audit["passed"]:
        with open(LIVE_PATH, "w") as handle:
            json.dump(compiled, handle, indent=2)
        written = True
        await _terrain_broadcast({"type": "body_updated", "name": safe_name})

    return {
        "name": safe_name,
        "written": written,
        "body": compiled,
        "audit": audit,
        "trajectory": trajectory,
    }


@app.post("/terrain/load/{name}")
async def terrain_load(name: str):
    if not SAFE_NAME.match(name):
        raise HTTPException(400, "Name must be alphanumeric + underscore only")
    try:
        body = load_vault_body(name)
    except FileNotFoundError as e:
        raise HTTPException(404, str(e)) from e

    compiled = json.loads(body.to_compiled_json(provenance="terrain-vault"))
    with open(LIVE_PATH, "w") as handle:
        json.dump(compiled, handle, indent=2)
    await _terrain_broadcast({"type": "body_updated", "name": name})
    return {"status": "ok", "name": name, "body": compiled}


@app.websocket("/ws/live")
async def terrain_live_socket(websocket: WebSocket):
    await websocket.accept()
    _terrain_ws_clients.add(websocket)
    try:
        while True:
            await websocket.receive_text()
    except WebSocketDisconnect:
        _terrain_ws_clients.discard(websocket)
    except Exception:
        _terrain_ws_clients.discard(websocket)


_SUGGEST_PAIRS = [
    # Cross-source morphs — maximal timbral contrast
    ("vowel", "ee", "nasal", "nasal_m", "bright vowel to nasal — F2 collision with anti-formant"),
    ("vowel", "ah", "bell", "Freiburg Hosanna", "open vowel to bell partials — formant-to-harmonic transition"),
    ("vowel", "oo", "electronic", "acid", "dark vowel to acid sweep — low formants meet resonant climb"),
    ("nasal", "nasal_n", "bell", "Stretched Treble", "nasal zeros against inharmonic bell partials"),
    ("vowel", "eh", "electronic", "telephone", "mid vowel to bandpass — spectral narrowing"),
    ("bell", "Berlin Freedom Bell", "vowel", "er", "low bell partials to colored vowel — mass to throat"),
    # Vowel-to-vowel — classic formant traverse
    ("vowel", "ee", "vowel", "oo", "front-to-back vowel — maximum F2 migration"),
    ("vowel", "ae", "vowel", "oo", "open-to-closed — F1 drops, F2 shifts"),
    ("vowel", "ah", "vowel", "ee", "open-back to closed-front — full vowel space diagonal"),
    ("vowel", "schwa", "vowel", "ih", "neutral to bright — subtle formant tightening"),
    # Nasal transitions
    ("vowel", "ah", "nasal", "nasal_n", "open vowel to uvular nasal — anti-formant carves the spectrum"),
    ("nasal", "nasal_m", "nasal", "nasal_n", "bilabial to uvular — anti-formant frequencies shift"),
    # Bell combinations
    ("bell", "Freiburg Hosanna", "bell", "St Mary le Tower", "two real bells — partial spacing differs"),
    ("electronic", "acid", "bell", "Stretched Treble", "sweep meets inharmonic partials"),
]


@app.post("/suggest")
def suggest_endpoint():
    """Suggest an interesting morph combination.

    Uses OpenRouter LLM if OPENROUTER_API_KEY is set, otherwise picks
    from curated cross-source pairs with acoustic rationale.
    """
    openrouter_key = os.environ.get("OPENROUTER_API_KEY")
    if openrouter_key:
        return _suggest_via_openrouter(openrouter_key)
    # Smart random from curated pairs
    sa, ka, sb, kb, rationale = random.choice(_SUGGEST_PAIRS)
    return {"source_a": sa, "key_a": ka, "source_b": sb, "key_b": kb, "rationale": rationale}


def _suggest_via_openrouter(api_key: str) -> dict:
    import httpx

    available = {
        "vowels": get_vowel_keys(),
        "nasals": get_nasal_keys(),
        "bells": get_bell_names(),
        "electronic": get_electronic_keys(),
    }
    prompt = (
        "You are a Z-plane filter body designer. Pick two sources to morph between.\n"
        f"Available: {json.dumps(available)}\n"
        "Source types: vowel, nasal, bell, electronic.\n"
        "Pick a pair that creates interesting spectral motion — "
        "formant transitions, pole-zero crossings, timbral contrast.\n"
        "Respond ONLY with JSON: "
        '{"source_a":"...","key_a":"...","source_b":"...","key_b":"...","rationale":"one sentence"}'
    )

    resp = httpx.post(
        "https://openrouter.ai/api/v1/chat/completions",
        headers={"Authorization": f"Bearer {api_key}"},
        json={
            "model": "openai/gpt-4o",
            "max_tokens": 200,
            "messages": [{"role": "user", "content": prompt}],
        },
        timeout=15.0,
    )
    if resp.status_code != 200:
        raise HTTPException(502, f"OpenRouter error: {resp.status_code}")
    text = resp.json()["choices"][0]["message"]["content"].strip()
    start = text.find("{")
    end = text.rfind("}") + 1
    if start < 0 or end <= start:
        raise HTTPException(502, f"LLM returned unparseable response: {text}")
    return json.loads(text[start:end])


@app.get("/vault")
def vault_listing():
    return list_vault_bodies()


@app.get("/vault/catalog")
def vault_catalog():
    children = sorted(VAULT_ROOT.iterdir(), key=lambda p: (not p.is_dir(), p.name.lower())) if VAULT_ROOT.exists() else []
    return {
        "root": str(VAULT_ROOT),
        "children": [_vault_node(child) for child in children if child.is_dir() or child.suffix.lower() in {".json", ".wav"}],
    }


@app.get("/vault/file")
def vault_file(relative_path: str):
    path = _vault_resolve(relative_path)
    if path.is_dir():
        return _vault_node(path)
    if path.suffix.lower() == ".json":
        with open(path) as f:
            data = json.load(f)
        return {
            "type": "json",
            "relative_path": relative_path,
            "preview": _json_preview(data),
            "data": data,
        }
    return {
        "type": "file",
        "relative_path": relative_path,
        "name": path.name,
        "size": path.stat().st_size,
    }


@app.post("/vault/push-live")
def vault_push_live(req: VaultFileRequest):
    path = _vault_resolve(req.relative_path)
    if path.suffix.lower() != ".json":
        raise HTTPException(400, "Only JSON bodies can be pushed live")
    raw_bytes = path.read_bytes()
    try:
        data = json.loads(raw_bytes.decode("utf-8"))
    except Exception as e:
        raise HTTPException(400, f"Vault JSON unreadable: {e}") from e

    # Truth-preserving path: compiled-v1 stays compiled-v1 (no StageParams round-trip).
    if isinstance(data, dict) and data.get("format") == "compiled-v1":
        Path(LIVE_PATH).write_bytes(raw_bytes)
        return {
            "status": "pushed",
            "name": str(data.get("name", "")),
            "relative_path": req.relative_path,
            "format": "compiled-v1",
        }

    body = Body.from_dict(data)
    compiled = body.to_compiled_json(provenance=f"vault-live:{req.relative_path}")
    with open(LIVE_PATH, "w", encoding="utf-8") as f:
        f.write(compiled)
    return {"status": "pushed", "name": body.name, "relative_path": req.relative_path, "format": "compiled-v1"}


@app.post("/splice")
def splice_endpoint(req: SpliceRequest):
    try:
        body_a = load_vault_body(req.body_a)
        body_b = load_vault_body(req.body_b)
    except FileNotFoundError as e:
        raise HTTPException(404, str(e)) from e

    try:
        mode = SpliceMode(req.mode)
    except ValueError as e:
        modes = [member.value for member in SpliceMode]
        raise HTTPException(400, f"Unknown mode '{req.mode}'. Use one of: {modes}") from e

    spliced = splice_corners(body_a.corners, body_b.corners, mode)
    body = Body(name=req.name, corners=spliced, boost=body_a.boost)
    return json.loads(body.to_json())


# ---------------------------------------------------------------------------
# Hostile Mode — 3-room authoring prototype
# ---------------------------------------------------------------------------

class HostileStage(BaseModel):
    zone: str | None = None
    hz: float | None = None
    role: str | None = None
    pole_hz: float | None = None
    radius: float | None = None
    gain: float | None = 0.0
    zero_hz: float | None = None
    zero_radius: float | None = None
    pole_hz_b: float | None = None
    zero_hz_b: float | None = None

class HostileCluster(BaseModel):
    name: str
    stages: list[int]
    morph_direction: str
    morph_target_zone: str
    internal_spread: str

class HostileCrossing(BaseModel):
    pair: list[int]

class HostileQLaw(BaseModel):
    compression: str
    center_hz: float
    dominance: list[dict] = Field(default_factory=list)

class HostileCompileRequest(BaseModel):
    name: str
    sentence: str
    stages: list[HostileStage]
    clusters: list[HostileCluster] = Field(default_factory=list)
    crossings: list[HostileCrossing] = Field(default_factory=list)
    q_law: HostileQLaw
    seed: int = 0

class HostileGateRequest(BaseModel):
    body: dict

class HostileSweepRequest(BaseModel):
    body: dict

class HostileTriageStartRequest(BaseModel):
    session_name: str = "triage"
    world_id: str = "bass_pressure"
    body_family: str = "speaker_knockerz"
    target_candidates: int = 20
    baseline_median_time_to_first_keep_s: float | None = None
    baseline_keeps_per_20: float | None = None

class HostileTriageDecisionRequest(BaseModel):
    session_id: str
    decision: str
    candidate_name: str = ""
    gate_passed: bool
    compile_ms: int = 0
    decision_latency_ms: int = 0
    gate_failures: list[str] = Field(default_factory=list)
    debug_open_count: int = 0
    export_path: str | None = None
    notes: str = ""


@app.get("/hostile")
def serve_hostile():
    return FileResponse(os.path.join(STATIC_DIR, "hostile.html"))


_hostile_triage_sessions: dict[str, dict] = {}


def _triage_safe_name(name: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_]+", "_", (name or "").strip())
    cleaned = re.sub(r"_+", "_", cleaned).strip("_")
    return cleaned or "triage"


def _triage_session_dir(session_id: str) -> Path:
    return (Path(TRIAGE_DIR) / session_id).resolve()


def _triage_session_path(session_id: str) -> Path:
    return _triage_session_dir(session_id) / "session.json"


def _triage_summary_path(session_id: str) -> Path:
    return _triage_session_dir(session_id) / "summary.json"


def _triage_store_session(session: dict) -> None:
    sdir = _triage_session_dir(session["id"])
    sdir.mkdir(parents=True, exist_ok=True)
    with open(_triage_session_path(session["id"]), "w", encoding="utf-8") as f:
        json.dump(session, f, indent=2)
    summary = _triage_build_summary(session)
    with open(_triage_summary_path(session["id"]), "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)


def _triage_load_session(session_id: str) -> dict | None:
    if session_id in _hostile_triage_sessions:
        return _hostile_triage_sessions[session_id]
    path = _triage_session_path(session_id)
    if not path.exists():
        return None
    with open(path, encoding="utf-8") as f:
        session = json.load(f)
    _hostile_triage_sessions[session_id] = session
    return session


def _triage_build_summary(session: dict) -> dict:
    events = session.get("events", [])
    keep_events = [e for e in events if e.get("decision") == "keep"]
    reject_events = [e for e in events if e.get("decision") == "reject"]
    again_events = [e for e in events if e.get("decision") == "again"]
    gate_pass_count = sum(1 for e in events if e.get("gate_passed"))
    latencies = [max(0, int(e.get("decision_latency_ms", 0))) for e in events]
    latencies = [x for x in latencies if x > 0]
    first_keep_s = keep_events[0]["elapsed_s"] if keep_events else None
    elapsed_s = max(0.0, time.time() - float(session.get("started_at_epoch", time.time())))
    total = len(events)
    evaluated = len(keep_events) + len(reject_events)
    keeps_per_20 = (len(keep_events) / evaluated * 20.0) if evaluated else 0.0
    gate_pass_rate = (gate_pass_count / total) if total else 0.0
    median_decision_latency_ms = statistics.median(latencies) if latencies else None
    return {
        "session_id": session["id"],
        "session_name": session.get("name", ""),
        "world_id": session.get("world_id", ""),
        "body_family": session.get("body_family", ""),
        "target_candidates": int(session.get("target_candidates", 20)),
        "started_at": session.get("started_at", ""),
        "elapsed_s": round(elapsed_s, 3),
        "total_decisions": total,
        "keep_count": len(keep_events),
        "reject_count": len(reject_events),
        "again_count": len(again_events),
        "gate_pass_count": gate_pass_count,
        "gate_pass_rate": round(gate_pass_rate, 4),
        "keeps_per_20": round(keeps_per_20, 3),
        "time_to_first_keep_s": round(first_keep_s, 3) if first_keep_s is not None else None,
        "median_decision_latency_ms": round(median_decision_latency_ms, 3) if median_decision_latency_ms is not None else None,
        "debug_open_count": int(session.get("debug_open_count", 0)),
        "baseline": {
            "median_time_to_first_keep_s": session.get("baseline_median_time_to_first_keep_s"),
            "keeps_per_20": session.get("baseline_keeps_per_20"),
        },
        "paths": {
            "session": str(_triage_session_path(session["id"])),
            "summary": str(_triage_summary_path(session["id"])),
        },
    }


@app.post("/hostile/triage/start")
def hostile_triage_start(req: HostileTriageStartRequest):
    session_name = _triage_safe_name(req.session_name)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    session_id = f"{session_name}_{stamp}"
    started_at_epoch = time.time()
    session = {
        "id": session_id,
        "name": session_name,
        "world_id": req.world_id,
        "body_family": req.body_family,
        "target_candidates": max(1, int(req.target_candidates)),
        "started_at_epoch": started_at_epoch,
        "started_at": datetime.fromtimestamp(started_at_epoch, tz=timezone.utc).isoformat(),
        "baseline_median_time_to_first_keep_s": req.baseline_median_time_to_first_keep_s,
        "baseline_keeps_per_20": req.baseline_keeps_per_20,
        "debug_open_count": 0,
        "events": [],
    }
    _hostile_triage_sessions[session_id] = session
    _triage_store_session(session)
    return _triage_build_summary(session)


@app.post("/hostile/triage/decision")
def hostile_triage_decision(req: HostileTriageDecisionRequest):
    session = _triage_load_session(req.session_id)
    if session is None:
        raise HTTPException(404, f"Unknown triage session '{req.session_id}'")

    decision = req.decision.strip().lower()
    if decision not in {"keep", "reject", "again"}:
        raise HTTPException(400, "decision must be keep, reject, or again")

    now_epoch = time.time()
    started_at_epoch = float(session.get("started_at_epoch", now_epoch))
    event = {
        "ts": datetime.fromtimestamp(now_epoch, tz=timezone.utc).isoformat(),
        "elapsed_s": max(0.0, now_epoch - started_at_epoch),
        "decision": decision,
        "candidate_name": req.candidate_name,
        "gate_passed": bool(req.gate_passed),
        "compile_ms": max(0, int(req.compile_ms)),
        "decision_latency_ms": max(0, int(req.decision_latency_ms)),
        "gate_failures": list(req.gate_failures),
        "debug_open_count": max(0, int(req.debug_open_count)),
        "export_path": req.export_path,
        "notes": req.notes,
    }
    session.setdefault("events", []).append(event)
    session["debug_open_count"] = max(int(session.get("debug_open_count", 0)), int(req.debug_open_count))
    _hostile_triage_sessions[session["id"]] = session
    _triage_store_session(session)
    return _triage_build_summary(session)


@app.get("/hostile/triage/{session_id}")
def hostile_triage_summary(session_id: str):
    session = _triage_load_session(session_id)
    if session is None:
        raise HTTPException(404, f"Unknown triage session '{session_id}'")
    return _triage_build_summary(session)


def _clip(v: float, lo: float, hi: float) -> float:
    return max(lo, min(hi, v))


def _q_pressure_params(compression: str) -> tuple[float, float]:
    key = (compression or "").strip().lower()
    if key == "extreme":
        return 0.55, 0.03
    if key == "moderate":
        return 0.75, 0.015
    return 1.0, 0.0


def _apply_q_pressure(freq_hz: float, radius: float, center_hz: float, compression: str) -> tuple[float, float]:
    squeeze, radius_boost = _q_pressure_params(compression)
    freq_q = center_hz + (freq_hz - center_hz) * squeeze
    freq_q = _clip(freq_q, 20.0, 18000.0)
    radius_q = _clip(radius + radius_boost, 0.0, 0.9995)
    return freq_q, radius_q


def _hostile_stage_params(stage: HostileStage, use_b_frame: bool = False, q_press: bool = False, q_law: HostileQLaw | None = None) -> StageParams:
    from pyruntime.stage_math import resonator, resonator_with_zero

    pole_hz_a = stage.pole_hz if stage.pole_hz is not None else (stage.hz if stage.hz is not None else 400.0)
    pole_hz_b = stage.pole_hz_b if stage.pole_hz_b is not None else pole_hz_a
    pole_hz = pole_hz_b if use_b_frame else pole_hz_a
    radius = stage.radius if stage.radius is not None else 0.95
    gain = stage.gain if stage.gain is not None else 0.0
    zero_hz_a = stage.zero_hz
    zero_hz_b = stage.zero_hz_b if stage.zero_hz_b is not None else zero_hz_a
    zero_hz = zero_hz_b if use_b_frame else zero_hz_a
    zero_radius = stage.zero_radius

    if q_press and q_law is not None:
        pole_hz, radius = _apply_q_pressure(
            freq_hz=pole_hz,
            radius=radius,
            center_hz=float(q_law.center_hz),
            compression=q_law.compression,
        )
        if zero_hz is not None and zero_radius is not None:
            zero_hz, _ = _apply_q_pressure(
                freq_hz=zero_hz,
                radius=zero_radius,
                center_hz=float(q_law.center_hz),
                compression=q_law.compression,
            )

    if zero_hz is not None and zero_radius is not None:
        return resonator_with_zero(
            freq_hz=float(_clip(pole_hz, 20.0, 18000.0)),
            radius=float(_clip(radius, 0.0, 0.9995)),
            val1=float(_clip(gain, -1.0, 1.0)),
            zero_freq_hz=float(_clip(zero_hz, 20.0, 18000.0)),
            zero_radius=float(_clip(zero_radius, 0.0, 1.2)),
        )
    return resonator(
        freq_hz=float(_clip(pole_hz, 20.0, 18000.0)),
        radius=float(_clip(radius, 0.0, 0.9995)),
        val1=float(_clip(gain, -1.0, 1.0)),
    )


def _build_hostile_direct_body(req: HostileCompileRequest) -> Body:
    stages_a = [_hostile_stage_params(s, use_b_frame=False, q_press=False, q_law=req.q_law) for s in req.stages[:6]]
    stages_c = [_hostile_stage_params(s, use_b_frame=True, q_press=False, q_law=req.q_law) for s in req.stages[:6]]
    stages_b = [_hostile_stage_params(s, use_b_frame=False, q_press=True, q_law=req.q_law) for s in req.stages[:6]]
    stages_d = [_hostile_stage_params(s, use_b_frame=True, q_press=True, q_law=req.q_law) for s in req.stages[:6]]

    while len(stages_a) < NUM_BODY_STAGES:
        stages_a.append(StageParams.passthrough())
    while len(stages_b) < NUM_BODY_STAGES:
        stages_b.append(StageParams.passthrough())
    while len(stages_c) < NUM_BODY_STAGES:
        stages_c.append(StageParams.passthrough())
    while len(stages_d) < NUM_BODY_STAGES:
        stages_d.append(StageParams.passthrough())

    corners = CornerArray(
        a=CornerState(stages=stages_a, boost=4.0),
        b=CornerState(stages=stages_b, boost=4.0),
        c=CornerState(stages=stages_c, boost=4.0),
        d=CornerState(stages=stages_d, boost=4.0),
    )
    return Body(name=req.name, corners=corners, boost=4.0)


def _build_hostile_sweep(body: Body) -> dict:
    freqs = freq_points(sr=SR)
    sweep_points = []
    for m in [0.0, 0.25, 0.5, 0.75, 1.0]:
        for q in [0.0, 0.5, 1.0]:
            enc = body.corners.interpolate(m, q)
            db = cascade_response_db(enc, freqs, SR)
            sweep_points.append({"morph": m, "q": q, "db": db.tolist()})
    return {"freqs": freqs.tolist(), "sweep": sweep_points}


def _count_peaks_in_band_local(db: list[float], freqs: list[float], lo: float, hi: float, threshold_above_mean: float = 6.0) -> int:
    xs = [(f, d) for f, d in zip(freqs, db) if lo <= f <= hi]
    if len(xs) < 3:
        return 0
    vals = [d for _, d in xs]
    mean_val = float(sum(vals) / len(vals))
    count = 0
    for i in range(1, len(vals) - 1):
        if vals[i] > vals[i - 1] and vals[i] > vals[i + 1] and vals[i] > mean_val + threshold_above_mean:
            count += 1
    return count


def _occupied_band_count_local(db: list[float], freqs: list[float]) -> int:
    bands = [(200.0, 500.0), (500.0, 1200.0), (1200.0, 2500.0), (2500.0, 5000.0)]
    if not db or not freqs:
        return 0
    peak = max(db)
    threshold = peak - 18.0
    occupied = 0
    for lo, hi in bands:
        vals = [d for f, d in zip(freqs, db) if lo <= f <= hi]
        if vals and max(vals) >= threshold:
            occupied += 1
    return occupied


def _direct_stage(s: DirectStageInput) -> StageParams:
    from pyruntime.stage_math import resonator, resonator_with_zero
    if s.zero_hz is not None and s.zero_radius is not None:
        return resonator_with_zero(s.pole_hz, s.radius, s.gain, s.zero_hz, s.zero_radius)
    return resonator(s.pole_hz, s.radius, s.gain)


def _build_direct_body(req: DirectBodyRequest) -> Body:
    from pyruntime.corner import CornerName
    fa = [_direct_stage(s) for s in req.frame_a[:6]]
    fb = [_direct_stage(s) for s in req.frame_b[:6]]
    while len(fa) < NUM_BODY_STAGES:
        fa.append(StageParams.passthrough())
    while len(fb) < NUM_BODY_STAGES:
        fb.append(StageParams.passthrough())

    # A=M0_Q0, B=M0_Q100, C=M100_Q0, D=M100_Q100
    # Heritage morph-only: A=B (Q collapsed at M0), C=D (Q collapsed at M100)
    corner_a = CornerState(stages=fa, boost=4.0)
    corner_c = CornerState(stages=fb, boost=4.0)
    corners = CornerArray(a=corner_a, b=corner_a, c=corner_c, d=corner_c)
    return Body(name=req.name, corners=corners, boost=4.0)


@app.post("/direct/response")
def direct_response(req: DirectBodyRequest):
    body = _build_direct_body(req)
    freqs = freq_points(sr=SR)
    enc = body.corners.interpolate(req.morph, req.q)
    db = cascade_response_db(enc, freqs, SR)
    return {"freqs": freqs.tolist(), "db": db.tolist()}


@app.post("/direct/body")
def direct_body(req: DirectBodyRequest):
    body = _build_direct_body(req)
    return json.loads(body.to_compiled_json(provenance="direct-author"))


@app.post("/direct/live")
def direct_live(req: DirectBodyRequest):
    body = _build_direct_body(req)
    compiled = body.to_compiled_json(provenance="direct-author")
    with open(LIVE_PATH, "w") as f:
        f.write(compiled)
    return {"status": "ok", "name": body.name}


@app.post("/direct/render")
def direct_render(req: DirectRenderRequest):
    body = _build_direct_body(req)
    return render_from_body(
        body=body, morph=req.morph, q=req.q,
        mackie_amount=req.mackie_amount, duration=req.duration,
    )


@app.post("/direct/export")
def direct_export(req: DirectBodyRequest):
    body = _build_direct_body(req)
    name = _terrain_safe_name(req.name)
    os.makedirs(VAULT_DIR, exist_ok=True)
    path = os.path.join(VAULT_DIR, f"{name}.json")
    with open(path, "w") as f:
        f.write(body.to_compiled_json(provenance="direct-author"))
    # Also write to trench_live.json for plugin hot-reload
    with open(LIVE_PATH, "w") as f:
        f.write(body.to_compiled_json(provenance="direct-author"))
    return {"path": path, "live": LIVE_PATH}


@app.get("/author")
def serve_author():
    return FileResponse(os.path.join(STATIC_DIR, "author.html"))


# ── Session endpoints ────────────────────────────────────────────────────

class CheckpointInput(BaseModel):
    name: str = "Untitled"
    morph: float = 0.0
    brightness: float = 0.0
    intensity: float = 0.0
    width: float = 0.0
    violence: float = 0.0
    anchors: list[float] = Field(default_factory=list)

class NewSessionInput(BaseModel):
    name: str
    world: str = ""
    brief: str = ""
    checkpoints: list[CheckpointInput] = Field(default_factory=list)

class UpdateCheckpointInput(BaseModel):
    session_name: str
    checkpoint_index: int
    checkpoint: CheckpointInput

class SessionMorphInput(BaseModel):
    session_name: str
    morph: float = 0.0
    q: float = 0.0


class ReplaceSessionInput(BaseModel):
    session_name: str
    world: str = ""
    brief: str = ""
    checkpoints: list[CheckpointInput] = Field(default_factory=list)


class VaultPushInput(BaseModel):
    relative_path: str

class GenerateBriefInput(BaseModel):
    brief: str
    world: str = ""

class VaultFileRequest(BaseModel):
    relative_path: str


def _vault_resolve(relative_path: str) -> Path:
    candidate = (VAULT_ROOT / relative_path).resolve()
    if not str(candidate).startswith(str(VAULT_ROOT)):
        raise HTTPException(400, "Path escapes vault")
    if not candidate.exists():
        raise HTTPException(404, f"Vault path not found: {relative_path}")
    return candidate


def _json_preview(data: dict) -> dict:
    preview: dict[str, object] = {
        "name": data.get("name", ""),
        "format": data.get("format", ""),
        "provenance": data.get("provenance", ""),
    }
    keyframes = data.get("keyframes")
    if isinstance(keyframes, list):
        preview["keyframe_count"] = len(keyframes)
        stage_counts = []
        labels = []
        for keyframe in keyframes[:4]:
            labels.append(keyframe.get("label", ""))
            stages = keyframe.get("stages", [])
            stage_counts.append(len(stages) if isinstance(stages, list) else 0)
        preview["labels"] = labels
        preview["stage_counts"] = stage_counts
    return preview


def _vault_node(path: Path) -> dict:
    rel = path.relative_to(VAULT_ROOT)
    if path.is_dir():
        summary_path = path / "_summary.json"
        summary = None
        if summary_path.exists():
            try:
                with open(summary_path) as f:
                    summary = json.load(f)
            except Exception:
                summary = {"error": "summary unreadable"}
        children = sorted(path.iterdir(), key=lambda p: (not p.is_dir(), p.name.lower()))
        return {
            "type": "dir",
            "name": path.name,
            "relative_path": str(rel).replace("\\", "/"),
            "summary": summary,
            "children": [_vault_node(child) for child in children if child.is_dir() or child.suffix.lower() in {".json", ".wav"}],
        }
    return {
        "type": "file",
        "name": path.name,
        "relative_path": str(rel).replace("\\", "/"),
        "ext": path.suffix.lower(),
        "size": path.stat().st_size,
    }


@app.get("/session")
def serve_session_page():
    return FileResponse(os.path.join(STATIC_DIR, "session.html"))


@app.get("/session/list")
def list_sessions():
    from pyruntime.session import Session
    return {"sessions": Session.list_sessions()}


@app.post("/session/new")
def new_session(req: NewSessionInput):
    from pyruntime.session import Session, Checkpoint
    cps = [Checkpoint(name=c.name, morph=c.morph, brightness=c.brightness,
                      intensity=c.intensity, width=c.width, violence=c.violence,
                      anchors=c.anchors) for c in req.checkpoints]
    s = Session(name=req.name, world=req.world, brief=req.brief,
                current_checkpoints=cps)
    s.save()
    return {"status": "created", "path": s.save()}


@app.post("/session/update")
def update_checkpoint(req: UpdateCheckpointInput):
    from pyruntime.session import Session, Checkpoint, push_live
    s = Session.load(req.session_name)
    cp = Checkpoint(name=req.checkpoint.name, morph=req.checkpoint.morph,
                    brightness=req.checkpoint.brightness,
                    intensity=req.checkpoint.intensity,
                    width=req.checkpoint.width,
                    violence=req.checkpoint.violence,
                    anchors=req.checkpoint.anchors)
    if req.checkpoint_index < len(s.current_checkpoints):
        s.current_checkpoints[req.checkpoint_index] = cp
    else:
        s.current_checkpoints.append(cp)
    s.save()
    push_live(s)
    return {"status": "updated", "pushed_live": True}


@app.post("/session/replace")
def replace_session(req: ReplaceSessionInput):
    from pyruntime.session import Session, Checkpoint, push_live
    s = Session.load(req.session_name)
    s.world = req.world
    s.brief = req.brief
    s.current_checkpoints = [
        Checkpoint(
            name=c.name,
            morph=c.morph,
            brightness=c.brightness,
            intensity=c.intensity,
            width=c.width,
            violence=c.violence,
            anchors=c.anchors,
        )
        for c in req.checkpoints
    ]
    s.save()
    push_live(s)
    return {"status": "replaced", "checkpoint_count": len(s.current_checkpoints), "pushed_live": True}


@app.post("/session/response")
def session_response(req: SessionMorphInput):
    from pyruntime.session import Session, compile_session_body
    s = Session.load(req.session_name)
    body = compile_session_body(s.current_checkpoints)
    freqs = freq_points(sr=SR)
    enc = body.corners.interpolate(req.morph, req.q)
    db = cascade_response_db(enc, freqs, SR)
    return {"freqs": freqs.tolist(), "db": db.tolist()}


@app.post("/session/sweep")
def session_sweep(req: SessionMorphInput):
    from pyruntime.session import Session, compile_session_body
    s = Session.load(req.session_name)
    body = compile_session_body(s.current_checkpoints)
    freqs = freq_points(sr=SR)
    points = []
    for m in [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]:
        enc = body.corners.interpolate(m, req.q)
        db = cascade_response_db(enc, freqs, SR)
        points.append({"morph": m, "db": db.tolist()})
    return {"freqs": freqs.tolist(), "sweep": points}


@app.post("/session/analyze")
def session_analyze(req: SessionMorphInput):
    from pyruntime.session import Session, compile_session_body, run_gates
    s = Session.load(req.session_name)
    body = compile_session_body(s.current_checkpoints)
    profile = body_profile(body)
    gates = run_gates(body, s.name)
    return {"profile": profile, "gates": gates}


@app.post("/session/save_take")
def session_save_take(req: SessionMorphInput):
    from pyruntime.session import Session, save_take
    s = Session.load(req.session_name)
    take = save_take(s, note="")
    return {"version": take.version, "gates": take.gates, "body_path": take.body_path}


@app.post("/session/push_live")
def session_push_live(req: SessionMorphInput):
    from pyruntime.session import Session, push_live
    s = Session.load(req.session_name)
    push_live(s)
    return {"status": "pushed"}


@app.get("/session/load/{name}")
def load_session(name: str):
    from pyruntime.session import Session
    from dataclasses import asdict
    s = Session.load(name)
    return asdict(s)


@app.post("/session/generate")
def session_generate(req: GenerateBriefInput):
    """Generate checkpoints from a brief using Claude API.

    Falls back to keyword matching if the API call fails.
    """
    import anthropic

    system_prompt = "You design filter morph timelines. Return ONLY valid JSON, no markdown, no explanation."
    user_prompt = f"""Generate checkpoints for this filter body.

World: {req.world}
Brief: {req.brief}

Return JSON: {{"checkpoints":[{{"name":"string","morph":0.0-1.0,"brightness":-1to1,"intensity":-1to1,"width":-1to1,"violence":-1to1,"anchors":[6 Hz values low to high]}}]}}

Rules: 3-7 checkpoints. morph=0.0 first, morph=1.0 last. 6 anchor frequencies per checkpoint.
Freq ranges: sub 30-80, chest 100-250, low-mid 250-500, mid 500-2000, presence 2000-5000, air 5000-16000.
brightness: negative=darker positive=brighter. intensity: controls resonance. width: bandwidth. violence: adds zeros/distortion."""

    try:
        import subprocess
        # prompt assembled above as full_prompt
        import shutil
        claude_path = shutil.which("claude") or "/c/Users/hooki/.local/bin/claude"
        full_prompt = system_prompt + "\n\n" + user_prompt
        result = subprocess.run(
            [claude_path, "-p", full_prompt, "--output-format", "json", "--model", "sonnet"],
            capture_output=True, text=True, timeout=90,
            cwd=os.path.dirname(__file__),
        )
        print(f"CLI rc={result.returncode} stdout_len={len(result.stdout)} stderr_len={len(result.stderr)}")
        if result.returncode == 0:
            output = result.stdout.strip()
            # Claude CLI with --output-format json wraps in {"type":"result","result":"..."}
            try:
                wrapper = json.loads(output)
                if "result" in wrapper:
                    text = wrapper["result"]
                else:
                    text = output
            except json.JSONDecodeError:
                text = output

            # Extract JSON from the text (may have markdown fences)
            if "```" in text:
                text = text.split("```json")[-1].split("```")[0] if "```json" in text else text.split("```")[1].split("```")[0]
            text = text.strip()

            data = json.loads(text)
            cps = data.get("checkpoints", data if isinstance(data, list) else [])
            if cps:
                return {"checkpoints": cps, "matched": "ai-generated"}
        print(f"CLI generation failed: rc={result.returncode} stderr={result.stderr[:200]}")
    except Exception as e:
        import traceback
        traceback.print_exc()
        print(f"AI generation failed: {e}, falling back to keyword match")
        # Fall through to keyword matching below
    """Generate checkpoints from a text brief using keyword matching.

    Maps body descriptions to checkpoint sequences using SHIPPING.md
    landmarks, Klatt formant tables, Bark scale, and instrument formants.
    """
    brief = req.brief.lower()
    world = req.world.lower()

    # ── Keyword → checkpoint library ──
    # Each entry: list of checkpoints with all fields
    LIBRARY = {
        "sub": [  # sub-dominant, bass-pressure, 808
            {"name": "Vault", "morph": 0.0, "brightness": -0.70, "intensity": 0.70, "width": 0.20, "violence": -0.20,
             "anchors": [45, 90, 150, 250, 500, 1000]},
            {"name": "Chest", "morph": 0.2, "brightness": -0.40, "intensity": 0.55, "width": 0.35, "violence": 0.00,
             "anchors": [45, 150, 220, 380, 700, 1200]},
            {"name": "Choke", "morph": 0.4, "brightness": -0.20, "intensity": 0.65, "width": -0.25, "violence": 0.25,
             "anchors": [45, 200, 320, 450, 800, 1400]},
            {"name": "Rip", "morph": 0.6, "brightness": 0.05, "intensity": 0.72, "width": -0.35, "violence": 0.60,
             "anchors": [45, 180, 400, 650, 1000, 1800]},
            {"name": "Rattle", "morph": 0.75, "brightness": 0.10, "intensity": 0.76, "width": -0.20, "violence": 0.72,
             "anchors": [45, 170, 330, 800, 1200, 2200]},
            {"name": "Cry", "morph": 0.9, "brightness": 0.35, "intensity": 0.80, "width": -0.55, "violence": 0.82,
             "anchors": [45, 160, 280, 700, 1200, 2600]},
            {"name": "Fracture", "morph": 1.0, "brightness": 0.55, "intensity": 0.85, "width": -0.60, "violence": 0.92,
             "anchors": [45, 140, 260, 520, 1200, 3200]},
        ],
        "vocal": [  # vocal formant, talking, vowel
            {"name": "Closed", "morph": 0.0, "brightness": -0.30, "intensity": 0.50, "width": 0.10, "violence": 0.00,
             "anchors": [150, 300, 870, 2240, 4000, 7000]},
            {"name": "Oo", "morph": 0.15, "brightness": -0.20, "intensity": 0.55, "width": 0.15, "violence": 0.00,
             "anchors": [150, 300, 870, 2240, 4000, 7000]},
            {"name": "Ah", "morph": 0.35, "brightness": 0.00, "intensity": 0.60, "width": 0.20, "violence": 0.10,
             "anchors": [150, 730, 1090, 2440, 4000, 7000]},
            {"name": "Ee", "morph": 0.55, "brightness": 0.20, "intensity": 0.65, "width": -0.10, "violence": 0.15,
             "anchors": [150, 270, 2290, 3010, 4800, 7000]},
            {"name": "Bite", "morph": 0.75, "brightness": 0.35, "intensity": 0.72, "width": -0.30, "violence": 0.40,
             "anchors": [150, 390, 1990, 2550, 5000, 7500]},
            {"name": "Shriek", "morph": 1.0, "brightness": 0.55, "intensity": 0.80, "width": -0.50, "violence": 0.70,
             "anchors": [150, 660, 1720, 2410, 5500, 8000]},
        ],
        "tube": [  # tube, pipe, comb, structural
            {"name": "Iron Pipe", "morph": 0.0, "brightness": -0.45, "intensity": 0.50, "width": 0.10, "violence": -0.10,
             "anchors": [90, 240, 420, 700, 1200, 2200]},
            {"name": "Bulge", "morph": 0.35, "brightness": 0.10, "intensity": 0.65, "width": 0.20, "violence": 0.30,
             "anchors": [100, 260, 500, 900, 1600, 3000]},
            {"name": "Null", "morph": 0.50, "brightness": 0.00, "intensity": 0.20, "width": -0.50, "violence": 0.80,
             "anchors": [110, 280, 520, 1000, 2200, 3600]},
            {"name": "Fracture", "morph": 0.70, "brightness": 0.30, "intensity": 0.72, "width": -0.30, "violence": 0.85,
             "anchors": [120, 320, 640, 1400, 2800, 5200]},
            {"name": "Dust", "morph": 1.0, "brightness": 0.70, "intensity": 0.40, "width": -0.60, "violence": 0.65,
             "anchors": [130, 360, 720, 1800, 4200, 8000]},
        ],
        "bright": [  # treble, air, sibilance, aluminum
            {"name": "Void", "morph": 0.0, "brightness": -0.60, "intensity": 0.30, "width": 0.30, "violence": -0.10,
             "anchors": [200, 800, 2000, 5000, 8000, 12000]},
            {"name": "Sheen", "morph": 0.25, "brightness": 0.20, "intensity": 0.50, "width": 0.20, "violence": 0.00,
             "anchors": [200, 800, 2000, 5000, 10000, 14000]},
            {"name": "Glass", "morph": 0.50, "brightness": 0.40, "intensity": 0.65, "width": -0.20, "violence": 0.35,
             "anchors": [200, 800, 2000, 5000, 7000, 12000]},
            {"name": "Tear", "morph": 0.75, "brightness": 0.60, "intensity": 0.75, "width": -0.40, "violence": 0.70,
             "anchors": [200, 800, 2000, 5000, 8000, 12000]},
            {"name": "Shatter", "morph": 1.0, "brightness": 0.80, "intensity": 0.85, "width": -0.55, "violence": 0.90,
             "anchors": [200, 800, 2000, 6000, 10000, 16000]},
        ],
    }

    # Match brief keywords to library
    matches = []
    keywords = {
        "sub": ["sub", "bass", "808", "low", "pressure", "knockerz", "speaker", "anchor", "weight"],
        "vocal": ["vocal", "vowel", "talk", "formant", "throat", "ah", "ee", "voice", "small talk"],
        "tube": ["tube", "pipe", "comb", "cul", "sac", "structural", "fracture", "iron", "null"],
        "bright": ["bright", "air", "treble", "sibilant", "aluminum", "glass", "crystal", "sheen", "high"],
    }
    for key, words in keywords.items():
        score = sum(1 for w in words if w in brief or w in world)
        if score > 0:
            matches.append((score, key))

    if matches:
        matches.sort(reverse=True)
        best = matches[0][1]
        return {"checkpoints": LIBRARY[best], "matched": best}
    else:
        # Default: generic sweep
        return {"checkpoints": LIBRARY["sub"], "matched": "sub (default)"}


@app.post("/hostile/compile")
def hostile_compile(req: HostileCompileRequest):
    """Compile hostile body and return body + multi-point response.

    Direct mode (preferred): stages include pole_hz/radius/gain(+optional zero)
    with optional B-frame pole_hz_b/zero_hz_b.
    Legacy mode (fallback): semantic zone/role stages.
    """
    direct_mode = any(s.pole_hz is not None for s in req.stages)
    if direct_mode:
        body = _build_hostile_direct_body(req)
    else:
        from pyruntime.forge_v2 import (
            SemanticBody, StageDef, ClusterDef, CrossingDef,
            QLawDef, QDominance, SemanticCompiler,
        )
        from pyruntime.forge_generator import solve_c4_surface, apply_gain_budget

        stages = [StageDef(zone=(s.zone or "boxy"), hz=(s.hz or 500.0), role=(s.role or "cluster")) for s in req.stages]
        clusters = [
            ClusterDef(
                name=c.name, stages=c.stages,
                morph_direction=c.morph_direction,
                morph_target_zone=c.morph_target_zone,
                internal_spread=c.internal_spread,
            ) for c in req.clusters
        ]
        crossings = [CrossingDef(pair=c.pair, morph_region="mid") for c in req.crossings]
        dominance = [QDominance(stage=d["stage"], role=d["role"]) for d in req.q_law.dominance]
        q_law = QLawDef(
            compression=req.q_law.compression,
            center_hz=req.q_law.center_hz,
            dominance=dominance,
        )
        sem = SemanticBody(
            name=req.name, sentence=req.sentence,
            stages=stages, clusters=clusters,
            crossings=crossings, q_law=q_law,
        )
        compiler = SemanticCompiler(sem)
        body = compiler.compile()
        body = solve_c4_surface(body, target_db=36.0)
        body = apply_gain_budget(body)

    body_dict = json.loads(body.to_compiled_json(provenance="hostile-mode"))
    sweep = _build_hostile_sweep(body)
    return {
        "body": body_dict,
        "freqs": sweep["freqs"],
        "sweep": sweep["sweep"],
    }


@app.get("/gallery")
def serve_gallery():
    return FileResponse(os.path.join(STATIC_DIR, "gallery.html"))


@app.get("/gallery/scan")
def gallery_scan(directory: str = "", limit: int = 200):
    """Scan a vault directory and return body metadata + sweep data for each."""
    import numpy as np

    if not directory:
        # List available directories
        vault = Path(VAULT_DIR).resolve()
        dirs = []
        for d in sorted(vault.iterdir()):
            if d.is_dir() and not d.name.startswith("_"):
                jsons = list(d.glob("*.json"))
                if jsons:
                    dirs.append({"name": d.name, "count": len(jsons)})
        return {"directories": dirs}

    # Resolve and validate path
    scan_dir = (Path(VAULT_DIR) / directory).resolve()
    if not str(scan_dir).startswith(str(VAULT_ROOT)):
        raise HTTPException(400, "Invalid directory")
    if not scan_dir.is_dir():
        raise HTTPException(404, f"Directory not found: {directory}")

    files = sorted(scan_dir.glob("*.json"))[:limit]
    freqs = freq_points(sr=SR)
    results = []

    for fpath in files:
        try:
            body = Body.from_json(str(fpath))
            audit = midpoint_audit(body)

            # Lightweight sweep: M0, M50, M100 at Q=0 only
            sweep = []
            for m in [0.0, 0.5, 1.0]:
                enc = body.corners.interpolate(m, 0.0)
                db = cascade_response_db(enc, freqs, SR)
                sweep.append({"morph": m, "q": 0.0, "db": db.tolist()})

            # Morph distance
            enc0 = body.corners.interpolate(0.0, 0.0)
            enc1 = body.corners.interpolate(1.0, 0.0)
            db0 = cascade_response_db(enc0, freqs, SR)
            db1 = cascade_response_db(enc1, freqs, SR)
            morph_rms = float(np.sqrt(np.mean((db0 - db1) ** 2)))

            results.append({
                "name": body.name,
                "file": fpath.name,
                "passed": audit["passed"],
                "worst_peak_db": round(audit["worst_peak_db"], 1),
                "morph_rms_db": round(morph_rms, 1),
                "sweep": sweep,
            })
        except Exception as e:
            results.append({
                "name": fpath.stem,
                "file": fpath.name,
                "passed": False,
                "error": str(e),
            })

    return {
        "directory": directory,
        "count": len(results),
        "freqs": freqs.tolist(),
        "bodies": results,
    }


@app.post("/gallery/keep")
def gallery_keep(body: dict):
    """Copy a kept body to the _triage/keeps directory."""
    name = body.get("directory", "unknown")
    filename = body.get("file", "unknown.json")
    src = (Path(VAULT_DIR) / name / filename).resolve()
    if not str(src).startswith(str(VAULT_ROOT)) or not src.exists():
        raise HTTPException(404, "Body not found")
    keeps_dir = Path(VAULT_DIR) / "_triage" / "keeps"
    keeps_dir.mkdir(parents=True, exist_ok=True)
    dst = keeps_dir / filename
    import shutil
    shutil.copy2(str(src), str(dst))
    return {"kept": str(dst)}


@app.post("/hostile/sweep")
def hostile_sweep(req: HostileSweepRequest):
    body = _load_body(req.body)
    return _build_hostile_sweep(body)


@app.post("/hostile/gate")
def hostile_gate(req: HostileGateRequest):
    """Run all validation gates on a candidate body."""
    body = _load_body(req.body)

    # 1. Midpoint audit
    audit = midpoint_audit(body)

    # 2. Structural validation
    issues = validate_corners(body.corners)
    struct_issues = [{"kind": i.kind, "severity": i.severity, "message": i.message} for i in issues]

    # 3. Morph trajectory distance
    morph_dist = morph_trajectory_distance(body)

    # 4. Q differentiation + occupied bands + sub anchor metrics
    freqs = freq_points(sr=SR)
    enc_m50_q0 = body.corners.interpolate(0.5, 0.0)
    enc_m50_q1 = body.corners.interpolate(0.5, 1.0)
    db_q0 = cascade_response_db(enc_m50_q0, freqs, SR)
    db_q1 = cascade_response_db(enc_m50_q1, freqs, SR)
    import numpy as np
    q_diff_rms = float(np.sqrt(np.mean((db_q0 - db_q1) ** 2)))
    occupied_band_count = _occupied_band_count_local(db_q1.tolist(), freqs.tolist())
    occupied_peak_count = _count_peaks_in_band_local(db_q1.tolist(), freqs.tolist(), 200.0, 5000.0, threshold_above_mean=6.0)

    test_morphs = [0.0, 0.25, 0.5, 0.75, 1.0]
    test_qs = [0.0, 0.5, 1.0]
    enc_ref = body.corners.interpolate(0.0, 0.0)
    db_ref = cascade_response_db(enc_ref, freqs, SR)
    ref_sub = max(float(v) for f, v in zip(freqs.tolist(), db_ref.tolist()) if 20.0 <= f <= 60.0)
    worst_sub_drop = 0.0
    worst_sub_point = (0.0, 0.0)
    for m in test_morphs:
        for q in test_qs:
            enc = body.corners.interpolate(m, q)
            db = cascade_response_db(enc, freqs, SR)
            sub_peak = max(float(v) for f, v in zip(freqs.tolist(), db.tolist()) if 20.0 <= f <= 60.0)
            drop = ref_sub - sub_peak
            if drop > worst_sub_drop:
                worst_sub_drop = drop
                worst_sub_point = (m, q)

    # Build gate results with plain-language verdicts
    gates = []

    # Midpoint gate
    if audit["passed"]:
        gates.append({
            "name": "Midpoint Stability",
            "passed": True,
            "detail": f"Peak {audit['worst_peak_db']:+.1f} dB — no resonance blowup in the sweep.",
        })
    else:
        gates.append({
            "name": "Midpoint Stability",
            "passed": False,
            "detail": f"Peak {audit['worst_peak_db']:+.1f} dB at morph={audit['worst_point']['morph']}, Q={audit['worst_point']['q']}. "
                       f"Poles are tearing across the unit circle mid-sweep. The interpolated coefficients produce a resonance spike "
                       f"that will distort or clip. Reduce radius on the most aggressive stage, or spread frequencies further apart.",
        })

    # Morph motion gate
    mean_morph = morph_dist["mean_rms_db"]
    if mean_morph >= 3.0:
        gates.append({
            "name": "Morph Motion",
            "passed": True,
            "detail": f"{mean_morph:.1f} dB RMS spectral change across the morph sweep. The knob does something.",
        })
    else:
        gates.append({
            "name": "Morph Motion",
            "passed": False,
            "detail": f"Only {mean_morph:.1f} dB RMS change across the morph sweep. "
                       f"The morph knob is nearly static — endpoints are too similar. "
                       f"Make Frame A and Frame B more different: change zones, swap roles, or add a cluster with real frequency movement.",
        })

    # Sub anchor gate
    if worst_sub_drop <= 6.0:
        gates.append({
            "name": "Sub Anchor",
            "passed": True,
            "detail": f"Sub held within {worst_sub_drop:.1f} dB across sweep.",
        })
    else:
        gates.append({
            "name": "Sub Anchor",
            "passed": False,
            "detail": f"Sub dropped {worst_sub_drop:.1f} dB at morph={worst_sub_point[0]}, Q={worst_sub_point[1]} (limit 6 dB).",
        })

    # Q differentiation gate
    if q_diff_rms >= 2.0:
        gates.append({
            "name": "Q Differentiation",
            "passed": True,
            "detail": f"{q_diff_rms:.1f} dB RMS difference between Q=0 and Q=1 at morph center. Pressure axis is alive.",
        })
    else:
        gates.append({
            "name": "Q Differentiation",
            "passed": False,
            "detail": f"Only {q_diff_rms:.1f} dB RMS difference between Q=0 and Q=1. "
                       f"The Q knob barely changes the sound. Add dominance entries to the Q law, "
                       f"increase compression, or assign different roles to stages so Q has something to compress.",
        })

    # Occupied bands / peak density at Q max center
    if occupied_band_count >= 2 and occupied_peak_count >= 2:
        gates.append({
            "name": "Occupied Bands",
            "passed": True,
            "detail": f"{occupied_band_count} occupied bands, {occupied_peak_count} peaks at morph=0.5, Q=1.0.",
        })
    else:
        gates.append({
            "name": "Occupied Bands",
            "passed": False,
            "detail": f"Only {occupied_band_count} occupied bands / {occupied_peak_count} peaks at morph=0.5, Q=1.0 (need >=2 each).",
        })

    # Structural issues
    errors = [i for i in struct_issues if i["severity"] == "error"]
    warnings = [i for i in struct_issues if i["severity"] == "warning"]
    if errors:
        for e in errors:
            gates.append({
                "name": f"Structure: {e['kind']}",
                "passed": False,
                "detail": e["message"],
            })
    if warnings:
        for w in warnings:
            gates.append({
                "name": f"Structure: {w['kind']}",
                "passed": True,
                "detail": f"Warning: {w['message']}",
            })

    all_passed = all(g["passed"] for g in gates)

    return {
        "passed": all_passed,
        "gates": gates,
        "audit": audit,
        "morph_distance": morph_dist,
        "q_diff_rms": q_diff_rms,
        "metrics": {
            "midpoint_peak_db": float(audit["worst_peak_db"]),
            "sub_anchor_delta_db": float(round(worst_sub_drop, 3)),
            "q_diff_rms_db": float(round(q_diff_rms, 3)),
            "occupied_band_count": int(occupied_band_count),
            "occupied_peak_count": int(occupied_peak_count),
        },
    }
