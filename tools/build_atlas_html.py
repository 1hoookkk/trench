"""Build a normalized HTML atlas for Trench bodies plus a clip-bin of frozen shapes.

Outputs (default):
  - vault/_atlas/index.html          (one page: Shapes + Bodies)
  - vault/_atlas/img/...             (normalized thumbnails)
  - vault/_clip_bin/shapes/...       (frozen static-state bodies, compiled-v1)
  - vault/_clip_bin/manifest.json    (machine-readable shape manifest)

Layer separation (do not mix or rank in one pool):
  XML       = reference only
  P2K       = complete-behavior oracle (dev/calibration)
  CUBE      = scaffold geometry
  VAULT     = benchmark taste (current product)
  GENERATED = candidates (watch drift)

Usage:
  PYTHONPATH=. python tools/build_atlas_html.py
  PYTHONPATH=. python tools/build_atlas_html.py --max-per-provenance 0   (full)
"""

from __future__ import annotations

import argparse
import dataclasses
import hashlib
import html
import json
import math
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks

ROOT = Path(__file__).resolve().parents[1]
import sys

sys.path.insert(0, str(ROOT))

from pyruntime.analysis import midpoint_audit
from pyruntime.body import Body
from pyruntime.constants import SR, TWO_PI, FREQ_MAX, FREQ_MIN
from pyruntime.corner import CornerArray, CornerName, CornerState
from pyruntime.designer_compile import compile_designer_to_body, parse_xml
from pyruntime.encode import EncodedCoeffs
from pyruntime.freq_response import cascade_response_db, freq_points
from pyruntime.stage_params import StageParams


matplotlib.rcParams.update(
    {
        "figure.facecolor": "#0a0a14",
        "axes.facecolor": "#0a0a14",
        "axes.edgecolor": "#333355",
        "axes.labelcolor": "#b5b5d6",
        "xtick.color": "#8c8cb0",
        "ytick.color": "#8c8cb0",
        "text.color": "#d7d7f5",
        "grid.color": "#232340",
        "font.family": "monospace",
        "font.size": 9,
    }
)


PROVENANCE_BADGES = {
    "XML": {"bg": "#2d6cdf", "fg": "#0a0a14"},
    "P2K": {"bg": "#df6c2d", "fg": "#0a0a14"},
    "CUBE": {"bg": "#2ddf8f", "fg": "#0a0a14"},
    "VAULT": {"bg": "#b6b6c8", "fg": "#0a0a14"},
    "GENERATED": {"bg": "#df2d6c", "fg": "#0a0a14"},
    "SHAPE": {"bg": "#8f2ddf", "fg": "#0a0a14"},
}


CORNER_LABELS: list[tuple[str, CornerName, float, float]] = [
    ("M0_Q0", CornerName.A, 0.0, 0.0),
    ("M0_Q100", CornerName.B, 0.0, 1.0),
    ("M100_Q0", CornerName.C, 1.0, 0.0),
    ("M100_Q100", CornerName.D, 1.0, 1.0),
]


def _slug(text: str) -> str:
    s = "".join(ch.lower() if ch.isalnum() else "_" for ch in text).strip("_")
    s = re.sub(r"_+", "_", s)
    return s or "untitled"


def _stable_id(*parts: str) -> str:
    h = hashlib.sha1()
    for p in parts:
        h.update(p.encode("utf-8", errors="replace"))
        h.update(b"\0")
    return h.hexdigest()[:16]


def _write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def _write_json(path: Path, obj: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2), encoding="utf-8")


def _safe_read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _hz_to_st(a_hz: float, b_hz: float) -> float:
    a = max(1.0, float(a_hz))
    b = max(1.0, float(b_hz))
    return abs(12.0 * math.log(a / b, 2.0))


def _mouth_band(centroid_hz: float) -> str:
    c = float(centroid_hz)
    if c < 400.0:
        return "chest"
    if c < 1000.0:
        return "throat"
    if c < 2500.0:
        return "mouth"
    if c < 7000.0:
        return "teeth"
    return "air"


def _stress_tag(dynamic_range_db: float, peak_db: float) -> str:
    # Corpus-derived thresholds (P2K corners) from tools/build_phoneme_atlas.py
    score = float(dynamic_range_db) + max(0.0, float(peak_db)) * 0.35
    if score < 110.0:
        return "calm"
    if score < 170.0:
        return "strained"
    return "violent"


def _pole_zero_from_kernel(enc: EncodedCoeffs) -> dict[str, float]:
    c0 = float(enc.c0)
    c1 = float(enc.c1)
    c2 = float(enc.c2)
    c3 = float(enc.c3)
    c4 = float(enc.c4)

    # Pole from denominator: 1 + c3 z^-1 + c4 z^-2
    pole_r = math.sqrt(abs(c4))
    if pole_r > 1e-6:
        cos_t = max(-1.0, min(1.0, -c3 / (2.0 * pole_r)))
        pole_hz = math.acos(cos_t) * SR / TWO_PI
    else:
        pole_hz = 0.0

    # Zero from numerator: c0 + c1 z^-1 + c2 z^-2
    if abs(c0) > 1e-9 and abs(c2) > 1e-12:
        zero_r = math.sqrt(abs(c2 / c0))
        if zero_r > 1e-6:
            zcos = max(-1.0, min(1.0, -(c1 / c0) / (2.0 * zero_r)))
            zero_hz = math.acos(zcos) * SR / TWO_PI
        else:
            zero_hz = 0.0
            zero_r = 0.0
    else:
        zero_hz = 0.0
        zero_r = 0.0

    return {
        "pole_hz": round(float(pole_hz), 1),
        "pole_r": round(float(pole_r), 6),
        "zero_hz": round(float(zero_hz), 1),
        "zero_r": round(float(zero_r), 6),
    }


def _vowel_hint_from_corner(
    enc_corner: list[EncodedCoeffs],
    vowels: dict[str, Any],
) -> dict[str, Any] | None:
    # Consider poles in [250, 3500] Hz (speech zone), pick 2 highest radius.
    stages = []
    for sidx, enc in enumerate(enc_corner[:6], start=1):
        pz = _pole_zero_from_kernel(enc)
        pz["stage"] = sidx
        stages.append(pz)

    candidates = [s for s in stages if 250.0 <= float(s.get("pole_hz", 0.0)) <= 3500.0]
    if len(candidates) < 2:
        return None

    candidates.sort(key=lambda s: (-float(s.get("pole_r", 0.0)), float(s.get("pole_hz", 0.0))))
    a = float(candidates[0]["pole_hz"])
    b = float(candidates[1]["pole_hz"])
    f1_hz = min(a, b)
    f2_hz = max(a, b)

    best_key = ""
    best_dist = float("inf")
    for key, v in vowels.items():
        d1 = _hz_to_st(f1_hz, float(v["f1"]))
        d2 = _hz_to_st(f2_hz, float(v["f2"]))
        dist = float(math.sqrt(d1 * d1 + d2 * d2))
        if dist < best_dist:
            best_dist = dist
            best_key = str(key)

    if best_key == "":
        return None

    if best_dist <= 4.0:
        conf = "high"
    elif best_dist <= 8.0:
        conf = "medium"
    elif best_dist <= 14.0:
        conf = "low"
    else:
        conf = "none"

    return {
        "key": best_key,
        "label": str(vowels.get(best_key, {}).get("label", best_key)),
        "dist_st": round(float(best_dist), 2),
        "f1_hz": round(float(f1_hz), 1),
        "f2_hz": round(float(f2_hz), 1),
        "confidence": conf,
    }


def _response_stats(body: Body, morph: float, q: float, freqs: np.ndarray) -> dict[str, Any]:
    enc = body.corners.interpolate(morph, q)
    db = cascade_response_db(enc, freqs, SR)

    peak_db = float(np.max(db)) if len(db) else -120.0
    notch_db = float(np.min(db)) if len(db) else -120.0
    dyn = float(peak_db - notch_db)

    linear = 10 ** (db / 20.0)
    total = float(np.sum(linear))
    centroid_hz = float(np.sum(freqs * linear) / total) if total > 0.0 else 0.0

    peaks, _ = find_peaks(db, prominence=3)
    notches, _ = find_peaks(-db, prominence=3)

    return {
        "peak_db": round(peak_db, 1),
        "notch_db": round(notch_db, 1),
        "dynamic_range_db": round(dyn, 1),
        "centroid_hz": round(centroid_hz, 1),
        "num_peaks": int(len(peaks)),
        "num_notches": int(len(notches)),
    }


def _stability_summary(body: Body) -> dict[str, Any]:
    # Cheap stability scan for atlas metadata (not a full trench-mcp scan).
    freqs = freq_points(n=256, sr=SR)
    grid_morphs = [0.0, 0.25, 0.5, 0.75, 1.0]
    grid_qs = [0.0, 0.5, 1.0]
    worst_peak = -999.0
    cliff = 0
    for m in grid_morphs:
        for q in grid_qs:
            enc = body.corners.interpolate(m, q)
            db = cascade_response_db(enc, freqs, SR)
            peak = float(np.max(db))
            worst_peak = max(worst_peak, peak)
            if peak > 55.0:
                cliff += 1

    audit = midpoint_audit(body, peak_limit_db=40.0)
    return {
        "worst_peak_db": round(float(worst_peak), 1),
        "cliff_count": int(cliff),
        "midpoint_audit": audit,
    }


def _apply_axes(ax: plt.Axes, *, y_min: float, y_max: float) -> None:
    ax.set_xlim(float(FREQ_MIN), float(FREQ_MAX))
    ax.set_ylim(float(y_min), float(y_max))
    ax.set_xlabel("Hz")
    ax.set_ylabel("dB")
    ax.axhline(0.0, color="#444455", linewidth=0.7, linestyle="--", alpha=0.8)
    ax.grid(True, which="both", linewidth=0.6, alpha=0.6)


def _plot_corner_card(body: Body, out_path: Path, *, y_min: float, y_max: float) -> None:
    freqs = freq_points(n=512, sr=SR)
    fig, axes = plt.subplots(2, 2, figsize=(8.6, 6.0))
    fig.suptitle(f"{body.name} / Corners", fontsize=12, fontweight="bold", color="#6fa4ff")

    for (label, _, m, q), ax in zip(CORNER_LABELS, axes.flatten()):
        enc = body.corners.interpolate(m, q)
        db = cascade_response_db(enc, freqs, SR)
        ax.semilogx(freqs, db, color="#ffffff", linewidth=1.5, alpha=0.95)
        ax.set_title(label, fontsize=10)
        _apply_axes(ax, y_min=y_min, y_max=y_max)

    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=160, bbox_inches="tight")
    plt.close(fig)


def _plot_sweep_strips(body: Body, out_path: Path, *, y_min: float, y_max: float) -> None:
    freqs = freq_points(n=512, sr=SR)
    fig, axes = plt.subplots(3, 1, figsize=(8.6, 8.8))
    fig.suptitle(f"{body.name} / Strips", fontsize=12, fontweight="bold", color="#6fa4ff")

    morph_positions = [0.0, 0.17, 0.33, 0.5, 0.67, 0.83, 1.0]
    q_positions = [0.0, 0.25, 0.5, 0.75, 1.0]
    diag_positions = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]

    # Strip 1: morph sweep at fixed Q (center)
    ax = axes[0]
    for i, m in enumerate(morph_positions):
        db = cascade_response_db(body.corners.interpolate(m, 0.5), freqs, SR)
        ax.semilogx(
            freqs,
            db,
            linewidth=1.1,
            alpha=0.9,
            color=plt.cm.cool(i / (len(morph_positions) - 1)),
        )
    ax.set_title("Morph strip @ Q=50%", fontsize=10)
    _apply_axes(ax, y_min=y_min, y_max=y_max)

    # Strip 2: Q sweep at fixed morph (center)
    ax = axes[1]
    for i, q in enumerate(q_positions):
        db = cascade_response_db(body.corners.interpolate(0.5, q), freqs, SR)
        ax.semilogx(
            freqs,
            db,
            linewidth=1.1,
            alpha=0.9,
            color=plt.cm.autumn(i / (len(q_positions) - 1)),
        )
    ax.set_title("Q strip @ M=50%", fontsize=10)
    _apply_axes(ax, y_min=y_min, y_max=y_max)

    # Strip 3: diagonal sweep (morph=Q)
    ax = axes[2]
    for i, d in enumerate(diag_positions):
        db = cascade_response_db(body.corners.interpolate(d, d), freqs, SR)
        ax.semilogx(
            freqs,
            db,
            linewidth=1.1,
            alpha=0.9,
            color=plt.cm.spring(i / (len(diag_positions) - 1)),
        )
    ax.set_title("Diagonal strip @ M=Q", fontsize=10)
    _apply_axes(ax, y_min=y_min, y_max=y_max)

    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=160, bbox_inches="tight")
    plt.close(fig)


def _plot_state_thumb(
    body: Body,
    morph: float,
    q: float,
    out_path: Path,
    *,
    title: str,
    y_min: float,
    y_max: float,
) -> None:
    freqs = freq_points(n=512, sr=SR)
    enc = body.corners.interpolate(morph, q)
    db = cascade_response_db(enc, freqs, SR)

    fig, ax = plt.subplots(1, 1, figsize=(4.8, 3.2))
    ax.semilogx(freqs, db, color="#ffffff", linewidth=1.6, alpha=0.95)
    ax.set_title(title, fontsize=10)
    _apply_axes(ax, y_min=y_min, y_max=y_max)

    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=160, bbox_inches="tight")
    plt.close(fig)


# -----------------------------------------------------------------------------
# Asset loading / atlas cards
# -----------------------------------------------------------------------------


@dataclass(frozen=True)
class AtlasCard:
    provenance: str
    title: str
    source_path: str
    card_id: str
    corner_png: str
    strips_png: str
    metadata: dict[str, Any]


def _corner_label_to_name(label: str) -> CornerName:
    mapping = {
        "M0_Q0": CornerName.A,
        "M0_Q100": CornerName.B,
        "M100_Q0": CornerName.C,
        "M100_Q100": CornerName.D,
    }
    if label not in mapping:
        raise KeyError(f"Unknown corner label: {label}")
    return mapping[label]


def _load_morpheus_cube_library(path: Path) -> list[dict[str, Any]]:
    d = _safe_read_json(path)
    cubes = d.get("cubes", [])
    if not isinstance(cubes, list):
        return []
    return cubes


def _cube_to_body(cube: dict[str, Any], *, name: str) -> Body:
    """Convert a Morpheus cube (4 corners x 6 pole-only stages) into a Body.

    Assumption: cube stage payload includes denominator (a1, a2) only; we render
    as an all-pole biquad (numerator = 1). Scaffold geometry only.
    """

    corners_raw = cube.get("corners", [])
    if not isinstance(corners_raw, list) or len(corners_raw) != 4:
        raise ValueError("Cube must have 4 corners")

    corners: list[CornerState] = []
    for corner in corners_raw:
        if not isinstance(corner, list):
            raise ValueError("Cube corner must be a list of stages")

        stages: list[StageParams] = []
        pre: list[EncodedCoeffs] = []
        for st in corner[:6]:
            a1 = float(st.get("a1", 0.0))
            a2 = float(st.get("a2", 0.0))
            r = float(st.get("r", math.sqrt(abs(a2))))

            enc = EncodedCoeffs(c0=1.0, c1=0.0, c2=0.0, c3=a1, c4=(a2 if a2 != 0.0 else (r * r)))
            pre.append(enc)

            # All-pole StageParams that re-encodes to numerator ~= [1,0,0].
            stages.append(StageParams(a1=a1, r=r, val1=0.0, val2=-a1, val3=r * r))

        while len(stages) < 12:
            stages.append(StageParams.passthrough())

        corners.append(CornerState(stages=stages, boost=4.0, _pre_encoded=list(pre)))

    ca = CornerArray(a=corners[0], b=corners[1], c=corners[2], d=corners[3])
    return Body(name=name, corners=ca, boost=4.0)


def _iter_generated_candidates() -> list[Path]:
    out: list[Path] = []
    overnight = ROOT / "vault" / "_overnight"
    if not overnight.is_dir():
        return out
    for run_dir in sorted(overnight.iterdir()):
        cand_root = run_dir / "candidates"
        if not cand_root.is_dir():
            continue
        out.extend(sorted(cand_root.rglob("*.json")))
    return sorted(out)


def _iter_vault_benchmarks() -> list[Path]:
    out: list[Path] = []
    vp = ROOT / "vault" / "_vocal_pack_2026-04-11" / "shipping"
    if vp.is_dir():
        out.extend(sorted(vp.glob("*.json")))
    return out


def _build_card(
    *,
    provenance: str,
    body: Body,
    source_path: str,
    out_img_dir: Path,
    y_min: float,
    y_max: float,
    vowels: dict[str, Any],
) -> AtlasCard:
    card_id = _stable_id(provenance, body.name, source_path)
    base = f"{_slug(provenance)}_{_slug(body.name)}_{card_id}"

    corner_png_path = out_img_dir / f"{base}__corners.png"
    strips_png_path = out_img_dir / f"{base}__strips.png"
    _plot_corner_card(body, corner_png_path, y_min=y_min, y_max=y_max)
    _plot_sweep_strips(body, strips_png_path, y_min=y_min, y_max=y_max)

    freqs = freq_points(n=512, sr=SR)
    corners_meta: dict[str, Any] = {}
    for label, cn, m, q in CORNER_LABELS:
        stats = _response_stats(body, m, q, freqs)
        enc_corner = body.corners.corner(cn).encode()
        vowel = _vowel_hint_from_corner(enc_corner, vowels)
        corners_meta[label] = {
            "response": stats,
            "mouth_band": _mouth_band(stats["centroid_hz"]),
            "stress": _stress_tag(stats["dynamic_range_db"], stats["peak_db"]),
            "closest_vowel": vowel,
        }

    stability = _stability_summary(body)
    meta = {
        "name": body.name,
        "provenance": provenance,
        "source_path": source_path,
        "corners": corners_meta,
        "stability": stability,
    }

    # HTML lives at out_img_dir.parent (out_dir); use "img/..." relative paths.
    return AtlasCard(
        provenance=provenance,
        title=body.name,
        source_path=source_path,
        card_id=card_id,
        corner_png=str(corner_png_path.relative_to(out_img_dir.parent).as_posix()),
        strips_png=str(strips_png_path.relative_to(out_img_dir.parent).as_posix()),
        metadata=meta,
    )


def _load_inventory_assets(inventory_path: Path) -> dict[str, list[dict[str, Any]]]:
    inv = _safe_read_json(inventory_path)
    assets = inv.get("assets", [])
    grouped: dict[str, list[dict[str, Any]]] = {}
    for a in assets:
        fam = str(a.get("family", ""))
        grouped.setdefault(fam, []).append(a)
    return grouped


def _select_assets(assets: list[dict[str, Any]], *, max_n: int | None, key_fn) -> list[dict[str, Any]]:
    items = sorted(list(assets), key=key_fn)
    if max_n is None or max_n <= 0:
        return items
    return items[:max_n]


# -----------------------------------------------------------------------------
# Clip bin (frozen shapes)
# -----------------------------------------------------------------------------


@dataclass(frozen=True)
class ShapeEntry:
    phoneme: str
    id: str
    corner: str
    filter_name: str
    shape_path: str
    thumb_png: str
    metadata: dict[str, Any]


def _freeze_corner_as_body(source: Body, corner_label: str, *, name: str) -> Body:
    cn = _corner_label_to_name(corner_label)
    src_corner = source.corners.corner(cn)
    pre = list(src_corner._pre_encoded) if src_corner._pre_encoded is not None else None

    # Flatten: all four corners identical.
    a = CornerState(stages=list(src_corner.stages), boost=src_corner.boost, _pre_encoded=pre)
    b = CornerState(stages=list(src_corner.stages), boost=src_corner.boost, _pre_encoded=(list(pre) if pre else None))
    c = CornerState(stages=list(src_corner.stages), boost=src_corner.boost, _pre_encoded=(list(pre) if pre else None))
    d = CornerState(stages=list(src_corner.stages), boost=src_corner.boost, _pre_encoded=(list(pre) if pre else None))
    ca = CornerArray(a=a, b=b, c=c, d=d)
    return Body(name=name, corners=ca, boost=source.boost)


def _build_clip_bin(
    *,
    phoneme_inventory_path: Path,
    out_dir: Path,
    atlas_img_dir: Path,
    y_min: float,
    y_max: float,
    vowels: dict[str, Any],
    max_phonemes: int | None,
    reps_per_phoneme: int,
) -> list[ShapeEntry]:
    inv = _safe_read_json(phoneme_inventory_path)
    phonemes = inv.get("phonemes", {})
    if not isinstance(phonemes, dict):
        return []

    out_shapes_dir = out_dir / "shapes"
    out_shapes_dir.mkdir(parents=True, exist_ok=True)

    entries: list[ShapeEntry] = []
    keys = sorted(phonemes.keys())
    if max_phonemes is not None and max_phonemes > 0:
        keys = keys[:max_phonemes]

    freqs = freq_points(n=512, sr=SR)

    for ph in keys:
        row = phonemes[ph]
        reps = row.get("representatives", [])
        if not isinstance(reps, list):
            continue

        for rep in reps[: max(1, int(reps_per_phoneme))]:
            pid = str(rep.get("id", ""))
            corner = str(rep.get("corner", ""))
            filter_name = str(rep.get("filter_name", ""))
            if pid == "" or corner == "":
                continue

            skin_path = ROOT / "datasets" / "p2k_skins" / f"{pid}.json"
            if not skin_path.is_file():
                continue

            src = Body.from_json(str(skin_path))
            shape_name = f"{ph} / {pid} / {corner}"
            shape_body = _freeze_corner_as_body(src, corner, name=shape_name)

            shape_slug = _slug(f"{ph}_{pid}_{corner}")
            shape_path = out_shapes_dir / f"{shape_slug}.json"
            _write_text(shape_path, shape_body.to_compiled_json(provenance="clip_bin"))

            thumb_png_path = atlas_img_dir / f"shape_{shape_slug}.png"
            _plot_state_thumb(
                shape_body,
                0.0,
                0.0,
                thumb_png_path,
                title=shape_name,
                y_min=y_min,
                y_max=y_max,
            )

            stats = _response_stats(shape_body, 0.0, 0.0, freqs)
            enc_corner = shape_body.corners.corner(CornerName.A).encode()
            vowel = _vowel_hint_from_corner(enc_corner, vowels)
            meta = {
                "phoneme": ph,
                "id": pid,
                "corner": corner,
                "filter_name": filter_name,
                "response": stats,
                "mouth_band": _mouth_band(stats["centroid_hz"]),
                "stress": _stress_tag(stats["dynamic_range_db"], stats["peak_db"]),
                "closest_vowel": vowel,
                "stability": _stability_summary(shape_body),
            }

            entries.append(
                ShapeEntry(
                    phoneme=ph,
                    id=pid,
                    corner=corner,
                    filter_name=filter_name,
                    shape_path=str(shape_path.relative_to(ROOT).as_posix()),
                    thumb_png=str(thumb_png_path.relative_to(atlas_img_dir.parent).as_posix()),
                    metadata=meta,
                )
            )

    _write_json(out_dir / "manifest.json", {"shapes": [dataclasses.asdict(e) for e in entries]})
    return entries


# -----------------------------------------------------------------------------
# HTML
# -----------------------------------------------------------------------------


def _badge_html(prov: str) -> str:
    info = PROVENANCE_BADGES.get(prov, {"bg": "#444444", "fg": "#ffffff"})
    bg = html.escape(info["bg"])
    fg = html.escape(info["fg"])
    return f'<span class="badge" style="background:{bg};color:{fg}">{html.escape(prov)}</span>'


def _fmt_vowel(v: dict[str, Any] | None) -> str:
    if not v:
        return ""
    key = html.escape(str(v.get("key", "")))
    conf = html.escape(str(v.get("confidence", "")))
    dist = html.escape(str(v.get("dist_st", "")))
    return f"{key} ({conf}, {dist}st)"


def _corner_meta_lines(cm: dict[str, Any]) -> str:
    r = cm.get("response", {})
    vowel = _fmt_vowel(cm.get("closest_vowel"))
    parts: list[str] = []
    parts.append(f"centroid {r.get('centroid_hz','')} Hz")
    parts.append(f"pk/not {r.get('num_peaks','')}/{r.get('num_notches','')}")
    parts.append(f"range {r.get('dynamic_range_db','')} dB")
    parts.append(str(cm.get("mouth_band", "")))
    parts.append(str(cm.get("stress", "")))
    if vowel:
        parts.append(f"vowel {vowel}")
    return " / ".join(html.escape(str(p)) for p in parts if str(p).strip() != "")


def _render_html(
    *,
    out_path: Path,
    shapes: list[ShapeEntry],
    cards: list[AtlasCard],
    max_per_provenance: int,
) -> None:
    cards_by_prov: dict[str, list[AtlasCard]] = {}
    for c in cards:
        cards_by_prov.setdefault(c.provenance, []).append(c)

    prov_order = ["XML", "P2K", "CUBE", "VAULT", "GENERATED"]
    for p in prov_order:
        cards_by_prov.setdefault(p, [])

    def _shape_row(e: ShapeEntry) -> str:
        r = e.metadata.get("response", {})
        stability = e.metadata.get("stability", {})
        audit = stability.get("midpoint_audit", {})
        audit_pass = bool(audit.get("passed", False))
        vowel = _fmt_vowel(e.metadata.get("closest_vowel"))
        return (
            "<div class='card'>"
            "<div class='card-head'>"
            f"{_badge_html('SHAPE')}"
            f"<div class='title'>{html.escape(e.phoneme)} / {html.escape(e.id)} / {html.escape(e.corner)}</div>"
            f"<div class='sub'>{html.escape(e.filter_name)}</div>"
            "</div>"
            "<div class='card-body'>"
            f"<img class='thumb' src='{html.escape(e.thumb_png)}' loading='lazy' />"
            "<div class='meta'>"
            f"<div><span class='k'>centroid</span> {html.escape(str(r.get('centroid_hz','')))} Hz</div>"
            f"<div><span class='k'>peaks/notches</span> {html.escape(str(r.get('num_peaks','')))}"
            f"/{html.escape(str(r.get('num_notches','')))}</div>"
            f"<div><span class='k'>range</span> {html.escape(str(r.get('dynamic_range_db','')))} dB</div>"
            f"<div><span class='k'>stress</span> {html.escape(str(e.metadata.get('stress','')))}</div>"
            f"<div><span class='k'>mouth</span> {html.escape(str(e.metadata.get('mouth_band','')))}</div>"
            f"<div><span class='k'>vowel</span> {vowel}</div>"
            "<div><span class='k'>midpoint</span> "
            f"<span class='{('ok' if audit_pass else 'bad')}'>"
            f"{('PASS' if audit_pass else 'FAIL')}</span>"
            f" / worst {html.escape(str(audit.get('worst_peak_db','')))} dB</div>"
            f"<div><span class='k'>cliffs</span> {html.escape(str(stability.get('cliff_count','')))}</div>"
            f"<div><span class='k'>shape</span> <span class='mono'>{html.escape(e.shape_path)}</span></div>"
            "</div>"
            "</div>"
            "</div>"
        )

    def _body_card(c: AtlasCard) -> str:
        stability = c.metadata.get("stability", {})
        audit = stability.get("midpoint_audit", {})
        audit_pass = bool(audit.get("passed", False))
        corners = c.metadata.get("corners", {})

        corner_lines = []
        for label, _, _, _ in CORNER_LABELS:
            cm = corners.get(label, {})
            corner_lines.append(
                f"<div><span class='k'>{html.escape(label)}</span> {_corner_meta_lines(cm)}</div>"
            )

        worst_point = audit.get("worst_point", {}) if isinstance(audit.get("worst_point", {}), dict) else {}
        wp_m = html.escape(str(worst_point.get("morph", "")))
        wp_q = html.escape(str(worst_point.get("q", "")))

        return (
            "<div class='card wide'>"
            "<div class='card-head'>"
            f"{_badge_html(c.provenance)}"
            f"<div class='title'>{html.escape(c.title)}</div>"
            f"<div class='sub mono'>{html.escape(c.source_path)}</div>"
            "</div>"
            "<div class='card-body wide'>"
            "<div class='imgs'>"
            f"<img class='img' src='{html.escape(c.corner_png)}' loading='lazy' />"
            f"<img class='img' src='{html.escape(c.strips_png)}' loading='lazy' />"
            "</div>"
            "<div class='meta'>"
            "<div><span class='k'>midpoint</span> "
            f"<span class='{('ok' if audit_pass else 'bad')}'>"
            f"{('PASS' if audit_pass else 'FAIL')}</span>"
            f" / worst {html.escape(str(audit.get('worst_peak_db','')))} dB @ m={wp_m} q={wp_q}</div>"
            "<div><span class='k'>cliffs</span> "
            f"{html.escape(str(stability.get('cliff_count','')))} (grid 5x3, pk>55dB)</div>"
            "<div class='spacer'></div>"
            + "".join(corner_lines)
            + "</div>"
            "</div>"
            "</div>"
        )

    max_note = (
        "<div class='note'>max-per-provenance = "
        f"{int(max_per_provenance)} (rerun with --max-per-provenance 0 for full)</div>"
        if max_per_provenance > 0
        else ""
    )

    shape_html = "\n".join(_shape_row(s) for s in shapes)

    body_sections = []
    for prov in prov_order:
        group = cards_by_prov.get(prov, [])
        body_sections.append(
            "<section>"
            f"<h2>{html.escape(prov)}</h2>"
            "<div class='grid wide-grid'>"
            + "\n".join(_body_card(c) for c in group)
            + "</div></section>"
        )

    page = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width,initial-scale=1" />
  <title>TRENCH Atlas</title>
  <style>
    :root {{
      --bg: #070711;
      --panel: #0f0f1d;
      --line: #232340;
      --text: #d7d7f5;
      --muted: #9a9ab8;
      --mono: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas,
              "Liberation Mono", "Courier New", monospace;
      --sans: ui-sans-serif, system-ui, -apple-system, Segoe UI, Roboto,
              "Helvetica Neue", Arial, "Noto Sans", "Liberation Sans", sans-serif;
    }}
    * {{ box-sizing: border-box; }}
    body {{
      margin: 0;
      background: var(--bg);
      color: var(--text);
      font: 13px/1.35 var(--sans);
    }}
    header {{
      position: sticky;
      top: 0;
      z-index: 10;
      background: rgba(7,7,17,0.92);
      backdrop-filter: blur(8px);
      border-bottom: 1px solid var(--line);
      padding: 10px 14px;
      display: flex;
      align-items: baseline;
      gap: 12px;
      flex-wrap: wrap;
    }}
    header h1 {{
      margin: 0;
      font-size: 14px;
      letter-spacing: 0.2px;
    }}
    header .sub {{
      color: var(--muted);
      font-family: var(--mono);
      font-size: 12px;
    }}
    main {{
      max-width: 1620px;
      margin: 0 auto;
      padding: 14px;
      display: grid;
      gap: 18px;
    }}
    h2 {{
      margin: 0 0 10px;
      font-size: 13px;
      text-transform: uppercase;
      letter-spacing: 0.7px;
      color: var(--muted);
    }}
    .note {{
      color: var(--muted);
      font-family: var(--mono);
      font-size: 12px;
    }}
    .grid {{
      display: grid;
      gap: 10px;
    }}
    .wide-grid {{
      grid-template-columns: 1fr;
    }}
    @media (min-width: 980px) {{
      .grid {{
        grid-template-columns: repeat(2, minmax(0, 1fr));
      }}
      .wide-grid {{
        grid-template-columns: 1fr;
      }}
    }}
    .card {{
      border: 1px solid var(--line);
      background: var(--panel);
      border-radius: 6px;
      overflow: hidden;
    }}
    .card.wide {{ grid-column: 1 / -1; }}
    .card-head {{
      padding: 10px 10px 8px;
      border-bottom: 1px solid var(--line);
      display: grid;
      grid-template-columns: auto 1fr;
      gap: 6px 10px;
      align-items: center;
    }}
    .badge {{
      display: inline-flex;
      align-items: center;
      justify-content: center;
      font-family: var(--mono);
      font-size: 11px;
      padding: 2px 8px;
      border-radius: 999px;
      height: 20px;
      letter-spacing: 0.4px;
    }}
    .title {{
      font-family: var(--mono);
      font-size: 12px;
      color: var(--text);
      white-space: nowrap;
      overflow: hidden;
      text-overflow: ellipsis;
    }}
    .sub {{
      grid-column: 2 / -1;
      color: var(--muted);
      font-family: var(--mono);
      font-size: 11px;
      white-space: nowrap;
      overflow: hidden;
      text-overflow: ellipsis;
    }}
    .card-body {{
      padding: 10px;
      display: grid;
      gap: 10px;
    }}
    .card-body.wide {{
      grid-template-columns: minmax(0, 1fr);
    }}
    @media (min-width: 980px) {{
      .card-body.wide {{
        grid-template-columns: minmax(0, 1.5fr) minmax(0, 1fr);
      }}
    }}
    .imgs {{
      display: grid;
      gap: 10px;
      grid-template-columns: minmax(0, 1fr);
    }}
    @media (min-width: 980px) {{
      .imgs {{
        grid-template-columns: repeat(2, minmax(0, 1fr));
      }}
    }}
    img {{
      display: block;
      width: 100%;
      height: auto;
      border: 1px solid var(--line);
      border-radius: 4px;
      background: #0a0a14;
    }}
    .thumb {{ max-width: 520px; }}
    .meta {{
      font-family: var(--mono);
      font-size: 12px;
      color: var(--text);
      display: grid;
      gap: 6px;
      align-content: start;
    }}
    .k {{
      color: var(--muted);
      display: inline-block;
      min-width: 94px;
    }}
    .mono {{ font-family: var(--mono); }}
    .ok {{ color: #2ddf8f; }}
    .bad {{ color: #df2d6c; }}
    .spacer {{ height: 6px; }}
  </style>
</head>
<body>
  <header>
    <h1>TRENCH Atlas</h1>
    <div class="sub">normalized plots / fixed axes / separated provenance</div>
  </header>
  <main>
    {max_note}
    <section>
      <h2>Shapes (Clip Bin)</h2>
      <div class="grid">
        {shape_html}
      </div>
    </section>
    <section>
      <h2>Bodies</h2>
      {''.join(body_sections)}
    </section>
  </main>
</body>
</html>
"""

    _write_text(out_path, page)


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------


def main() -> None:
    ap = argparse.ArgumentParser(description="Build a normalized HTML atlas plus clip bin.")
    ap.add_argument("--out-dir", type=Path, default=ROOT / "vault" / "_atlas")
    ap.add_argument("--img-dir", type=Path, default=None, help="default: <out-dir>/img")
    ap.add_argument(
        "--inventory",
        type=Path,
        default=ROOT
        / "datasets"
        / "filter_inventory"
        / "2026-03-10"
        / "filter_asset_inventory_v1.json",
    )
    ap.add_argument(
        "--phoneme-inventory",
        type=Path,
        default=ROOT / "vault" / "_phonemes" / "p2k_phoneme_inventory.json",
    )
    ap.add_argument(
        "--vowel-table",
        type=Path,
        default=ROOT / "docs" / "sonic_tables" / "tables.json",
    )
    ap.add_argument("--max-per-provenance", type=int, default=48)
    ap.add_argument("--max-phonemes", type=int, default=32)
    ap.add_argument("--reps-per-phoneme", type=int, default=2)
    ap.add_argument("--y-min", type=float, default=-80.0)
    ap.add_argument("--y-max", type=float, default=60.0)
    ap.add_argument(
        "--generated",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="include generated candidate bodies",
    )
    ap.add_argument(
        "--vault-benchmarks",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="include vault/_vocal_pack_2026-04-11/shipping benchmark bodies",
    )
    args = ap.parse_args()

    out_dir: Path = args.out_dir
    img_dir: Path = args.img_dir or (out_dir / "img")
    clip_bin_dir = ROOT / "vault" / "_clip_bin"

    tables = _safe_read_json(args.vowel_table)
    vowels = tables.get("vowels", {})
    if not isinstance(vowels, dict):
        vowels = {}

    img_dir.mkdir(parents=True, exist_ok=True)
    out_dir.mkdir(parents=True, exist_ok=True)
    clip_bin_dir.mkdir(parents=True, exist_ok=True)

    y_min = float(args.y_min)
    y_max = float(args.y_max)

    # 1) Clip bin shapes
    shapes = _build_clip_bin(
        phoneme_inventory_path=args.phoneme_inventory,
        out_dir=clip_bin_dir,
        atlas_img_dir=img_dir,
        y_min=y_min,
        y_max=y_max,
        vowels=vowels,
        max_phonemes=(None if int(args.max_phonemes) <= 0 else int(args.max_phonemes)),
        reps_per_phoneme=int(args.reps_per_phoneme),
    )

    # 2) Bodies atlas cards
    cards: list[AtlasCard] = []

    max_n = None if int(args.max_per_provenance) <= 0 else int(args.max_per_provenance)

    if args.inventory.is_file():
        grouped = _load_inventory_assets(args.inventory)

        # XML templates (reference)
        xml_assets = _select_assets(
            grouped.get("xml_template", []),
            max_n=max_n,
            key_fn=lambda a: str(a.get("name", "")),
        )
        for a in xml_assets:
            p = Path(str(a.get("path", "")))
            if not p.is_file():
                continue
            try:
                tmpl = parse_xml(str(p))
                body = compile_designer_to_body(tmpl, boost=4.0)
            except Exception:
                continue
            cards.append(
                _build_card(
                    provenance="XML",
                    body=body,
                    source_path=str(p),
                    out_img_dir=img_dir,
                    y_min=y_min,
                    y_max=y_max,
                    vowels=vowels,
                )
            )

        # P2K skins (oracle)
        p2k_assets = _select_assets(
            grouped.get("p2k_migrated_cartridge", []),
            max_n=max_n,
            key_fn=lambda a: str(a.get("name", "")),
        )
        for a in p2k_assets:
            rel = str(a.get("path", ""))
            p = (ROOT / rel) if not os.path.isabs(rel) else Path(rel)
            if not p.is_file():
                continue
            try:
                body = Body.from_json(str(p))
            except Exception:
                continue
            cards.append(
                _build_card(
                    provenance="P2K",
                    body=body,
                    source_path=str(p),
                    out_img_dir=img_dir,
                    y_min=y_min,
                    y_max=y_max,
                    vowels=vowels,
                )
            )

        # Morpheus cubes (scaffold)
        cube_assets = _select_assets(
            grouped.get("morpheus_cube", []),
            max_n=max_n,
            key_fn=lambda a: int(a.get("index", 0) or 0),
        )
        morpheus_path = ROOT / "datasets" / "morpheus_zplane_library.json"
        cubes = _load_morpheus_cube_library(morpheus_path) if morpheus_path.is_file() else []
        for a in cube_assets:
            idx: int | None
            if "index" in a:
                try:
                    idx = int(a["index"])
                except Exception:
                    idx = None
            else:
                asset_id = str(a.get("asset_id", ""))
                m = re.search(r":(\d+)$", asset_id)
                idx = int(m.group(1)) if m else None
            if idx is None:
                continue
            if idx < 0 or idx >= len(cubes):
                continue
            try:
                body = _cube_to_body(cubes[idx], name=f"Cube {idx:03d}")
            except Exception:
                continue
            cards.append(
                _build_card(
                    provenance="CUBE",
                    body=body,
                    source_path=str(morpheus_path) + f"#cube={idx}",
                    out_img_dir=img_dir,
                    y_min=y_min,
                    y_max=y_max,
                    vowels=vowels,
                )
            )

        # Promoted cartridges (benchmark)
        vault_assets = _select_assets(
            grouped.get("promoted_cartridge", []),
            max_n=max_n,
            key_fn=lambda a: str(a.get("name", "")),
        )
        for a in vault_assets:
            rel = str(a.get("path", ""))
            p = (ROOT / rel) if not os.path.isabs(rel) else Path(rel)
            if not p.is_file():
                continue
            try:
                body = Body.from_json(str(p))
            except Exception:
                continue
            cards.append(
                _build_card(
                    provenance="VAULT",
                    body=body,
                    source_path=str(p),
                    out_img_dir=img_dir,
                    y_min=y_min,
                    y_max=y_max,
                    vowels=vowels,
                )
            )

    # Extra vault benchmarks (shipping picks)
    if bool(args.vault_benchmarks):
        for p in _iter_vault_benchmarks():
            try:
                body = Body.from_json(str(p))
            except Exception:
                continue
            cards.append(
                _build_card(
                    provenance="VAULT",
                    body=body,
                    source_path=str(p),
                    out_img_dir=img_dir,
                    y_min=y_min,
                    y_max=y_max,
                    vowels=vowels,
                )
            )

    # Generated candidates (drift watch)
    if bool(args.generated):
        gen_paths = _iter_generated_candidates()
        if max_n is not None and max_n > 0:
            gen_paths = gen_paths[:max_n]
        for p in gen_paths:
            try:
                body = Body.from_json(str(p))
            except Exception:
                continue
            cards.append(
                _build_card(
                    provenance="GENERATED",
                    body=body,
                    source_path=str(p),
                    out_img_dir=img_dir,
                    y_min=y_min,
                    y_max=y_max,
                    vowels=vowels,
                )
            )

    cards.sort(key=lambda c: (c.provenance, c.title.lower(), c.source_path))

    _write_json(out_dir / "index.json", {"cards": [dataclasses.asdict(c) for c in cards]})
    _render_html(
        out_path=out_dir / "index.html",
        shapes=shapes,
        cards=cards,
        max_per_provenance=int(args.max_per_provenance),
    )

    print(f"wrote {out_dir / 'index.html'}")
    print(f"wrote {out_dir / 'index.json'}")
    print(f"wrote {clip_bin_dir / 'manifest.json'}")


if __name__ == "__main__":
    main()
