"""forge/iconic_author.py — emit RawStage cartridges for the 4 ship bodies.

Authoring touches the M0 and M100 morph corners only. The cartridge writer
auto-duplicates each morph corner into both Q slots — Q is runtime infrastructure,
not authored data, and the runtime applies its Q law on top of the morph corners.

Talking Hedz is calibration evidence for three structural laws: uniform per-cartridge
c0 attenuation (cascade gain discipline), interior-zero behavior (zero_r < pole_r,
parked away from poles), and multi-stage cancellation (zeros across stages combine).
Those laws are kept here. Heritage's specific 4-corner pattern is not a template.

Per-stage math (matches cascade.rs DF2T with positive c3, c4 in denominator):

    a1   = -2 * pole_r * cos(2π * pole_hz / sr)
    a2   = pole_r²
    b1   = -2 * zero_r * cos(2π * zero_hz / sr)
    b2   = zero_r²
    val1 = c0 - 1               (per-stage attenuation, uniform across cartridge)
    val2 = b1 - a1              (gives c1 = b1 in heritage convention)
    val3 = a2 - b2              (gives c2 = b2)
    flag = 1 active | 0 passthrough

Output:
    cartridges/factory/generated/iconic/<body_id>.json   (RawStage shipping format)
    cartridges/factory/generated/iconic/<body_id>.png    (M0 + M100 cascade plots)

Usage:
    python forge/iconic_author.py speaker_knockerz
    python forge/iconic_author.py aluminum_siding
    python forge/iconic_author.py small_talk_ah_ee
    python forge/iconic_author.py cul_de_sac
    python forge/iconic_author.py all
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from boost_normalizer import (  # noqa: E402 — after matplotlib backend set
    _cascade_peak,
    EVAL_FREQS,
    PEAK_AMPLITUDE_TARGET,
)

REPO = Path(__file__).resolve().parent.parent
OUT_DIR = REPO / "cartridges" / "factory" / "generated" / "iconic"
AUTHORING_SR = 39062.5
N_SLOTS = 12  # runtime cascade has 12 slots; unused stages are passthroughs


# ─── Stage authoring ─────────────────────────────────────────────────────────

@dataclass
class Stage:
    pole_hz: float
    pole_r: float
    zero_hz: float
    zero_r: float
    flag: int = 1


def stage_to_rawstage(s: Stage, c0: float, sr: float = AUTHORING_SR) -> dict:
    """Author stage → RawStage shipping dict {a1, r, val1, val2, val3, flag}."""
    if s.flag == 0:
        return {"a1": 0.0, "r": 0.0, "val1": 0.0, "val2": 0.0, "val3": 0.0, "flag": 0}
    pole_r = min(s.pole_r, 0.999999)
    zero_hz = min(s.zero_hz, sr * 0.49)
    a1 = -2.0 * pole_r * math.cos(2.0 * math.pi * s.pole_hz / sr)
    a2 = pole_r * pole_r
    b1 = -2.0 * s.zero_r * math.cos(2.0 * math.pi * zero_hz / sr)
    b2 = s.zero_r * s.zero_r
    return {
        "a1": a1,
        "r": pole_r,
        "val1": c0 - 1.0,
        "val2": b1 - a1,
        "val3": a2 - b2,
        "flag": 1,
    }


def rawstage_to_coefs(rs: dict) -> dict:
    """RawStage → cascade coefficients (matches runtime/trench-core/src/cascade.rs)."""
    if rs["flag"] == 0:
        return {"c0": 1.0, "c1": 0.0, "c2": 0.0, "c3": 0.0, "c4": 0.0, "_pass": True}
    a1 = rs["a1"]
    r = rs["r"]
    a2 = r * r
    return {
        "c0": 1.0 + rs["val1"],
        "c1": a1 + rs["val2"],
        "c2": a2 - rs["val3"],
        "c3": a1,
        "c4": a2,
        "_pass": False,
    }


def pad_passthrough(stages: list[dict]) -> list[dict]:
    out = list(stages)
    while len(out) < N_SLOTS:
        out.append({"a1": 0.0, "r": 0.0, "val1": 0.0, "val2": 0.0, "val3": 0.0, "flag": 0})
    return out


# ─── Cartridge assembly (Q auto-duplicated) ──────────────────────────────────

def assemble_cartridge(name: str, c0: float, m0: list[Stage], m100: list[Stage]) -> dict:
    """Emit P2k_NNN ROM-format cartridge.
    filterType 13 = MorphDesigner-class. Q corners are auto-duplicated from M
    corners (runtime infrastructure per project doctrine). boost 4.0 = heritage
    default. Top-level field order matches juce-shell/assets/sonic_tables.json.
    """
    m0_baked = pad_passthrough([stage_to_rawstage(s, c0) for s in m0])
    m100_baked = pad_passthrough([stage_to_rawstage(s, c0) for s in m100])
    n_stages = max(len(m0), len(m100))
    return {
        "name": name,
        "filterType": 13,
        "sampleRate": AUTHORING_SR,
        "stages": n_stages,
        "boost": 4.0,
        "provenance": {
            "author": "iconic_author.py",
            "c0_uniform": c0,
            "q_axis": "duplicated — runtime infrastructure, not authored data",
        },
        "keyframes": [
            {"morph": 0.0, "q": 0.0, "label": "M0_Q0",     "boost": 4.0, "stages": m0_baked},
            {"morph": 0.0, "q": 1.0, "label": "M0_Q100",   "boost": 4.0, "stages": m0_baked},
            {"morph": 1.0, "q": 0.0, "label": "M100_Q0",   "boost": 4.0, "stages": m100_baked},
            {"morph": 1.0, "q": 1.0, "label": "M100_Q100", "boost": 4.0, "stages": m100_baked},
        ],
    }


# ─── Validation cascade math + plot ──────────────────────────────────────────

def cascade_db(coefs: list[dict], freqs: np.ndarray, sr: float = AUTHORING_SR) -> np.ndarray:
    w = 2.0 * np.pi * freqs / sr
    e1 = np.exp(-1j * w)
    e2 = np.exp(-2j * w)
    mag = np.ones_like(freqs)
    for s in coefs:
        if s.get("_pass"):
            continue
        num = s["c0"] + s["c1"] * e1 + s["c2"] * e2
        den = 1.0 + s["c3"] * e1 + s["c4"] * e2
        mag *= np.abs(num / den)
    return 20.0 * np.log10(np.maximum(mag, 1e-12))


def stage_db(coef: dict, freqs: np.ndarray, sr: float = AUTHORING_SR) -> np.ndarray:
    if coef.get("_pass"):
        return np.zeros_like(freqs)
    w = 2.0 * np.pi * freqs / sr
    e1 = np.exp(-1j * w)
    e2 = np.exp(-2j * w)
    num = coef["c0"] + coef["c1"] * e1 + coef["c2"] * e2
    den = 1.0 + coef["c3"] * e1 + coef["c4"] * e2
    return 20.0 * np.log10(np.maximum(np.abs(num / den), 1e-12))


def plot_cartridge(name: str, c0: float, m0: list[Stage], m100: list[Stage],
                   cart: dict, out_path: Path, body_brief: str) -> None:
    freqs = np.geomspace(40.0, 18000.0, 4096)
    m0_baked = cart["keyframes"][0]["stages"][:len(m0)]
    m100_baked = cart["keyframes"][2]["stages"][:len(m100)]
    m0_coefs = [rawstage_to_coefs(rs) for rs in m0_baked]
    m100_coefs = [rawstage_to_coefs(rs) for rs in m100_baked]
    db_m0 = cascade_db(m0_coefs, freqs)
    db_m100 = cascade_db(m100_coefs, freqs)

    fig, axes = plt.subplots(2, 2, figsize=(16, 9), dpi=140, sharex=True)
    fig.patch.set_facecolor("#16191a")
    fig.suptitle(f"ICONIC :: {name}   (c0={c0:.3f} uniform)\n{body_brief}",
                 color="#ffba00", fontweight="bold", fontsize=12)

    stage_colors = ["#22ddff", "#88ee99", "#ffeb3b", "#ff9933", "#ff4488", "#cc66ff"]

    panels = [
        (axes[0, 0], f"M0 cascade  peak {db_m0.max():+.1f} dB", db_m0, m0, m0_coefs, "#88ee99"),
        (axes[1, 0], f"M100 cascade peak {db_m100.max():+.1f} dB", db_m100, m100, m100_coefs, "#ff5555"),
    ]
    breakdowns = [
        (axes[0, 1], "M0 per-stage", m0, m0_coefs),
        (axes[1, 1], "M100 per-stage", m100, m100_coefs),
    ]

    all_db = np.concatenate([db_m0, db_m100])
    y_lo = float(np.percentile(all_db, 0.5)) - 5.0
    y_hi = float(all_db.max()) + 5.0

    for ax, title, db, stages_meta, _coefs, color in panels:
        ax.set_facecolor("#0d0f10")
        ax.semilogx(freqs, db, color=color, lw=2.4)
        ax.set_title(title, color=color, fontsize=10, fontweight="bold")
        ax.set_ylim(y_lo, y_hi)
        ax.set_xlim(40, 18000)
        ax.grid(True, which="both", alpha=0.18, color="white")
        ax.tick_params(colors="white", labelsize=8)
        for spine in ax.spines.values():
            spine.set_color("#3a3e42")
        ax.set_ylabel("dB", color="white", fontsize=9)
        ax.axhline(0.0, color="#444", lw=0.5, alpha=0.5)
        for s in stages_meta:
            ax.axvline(s.pole_hz, color=color, lw=0.5, alpha=0.45, linestyle=":")
        # Label the most prominent feature freqs from the body brief
        for fhz, lbl in extract_landmarks(name):
            ax.axvline(fhz, color="#ffba00", lw=0.8, alpha=0.55, linestyle="--")
            ax.text(fhz, y_hi - 2.0, lbl, color="#ffba00", fontsize=7,
                    rotation=90, ha="right", va="top", alpha=0.8)

    for ax, title, stages_meta, coefs in breakdowns:
        ax.set_facecolor("#0d0f10")
        for i, (s, c) in enumerate(zip(stages_meta, coefs)):
            db = stage_db(c, freqs)
            ax.semilogx(freqs, db, color=stage_colors[i % len(stage_colors)], lw=1.4,
                        label=f"S{i}: {s.pole_hz:>5.0f} Hz r={s.pole_r:.3f}")
        ax.set_title(title, color="white", fontsize=10, fontweight="bold")
        ax.set_ylim(y_lo, y_hi)
        ax.set_xlim(40, 18000)
        ax.grid(True, which="both", alpha=0.18, color="white")
        ax.tick_params(colors="white", labelsize=8)
        for spine in ax.spines.values():
            spine.set_color("#3a3e42")
        ax.set_ylabel("dB", color="white", fontsize=9)
        ax.axhline(0.0, color="#444", lw=0.5, alpha=0.5)
        ax.legend(facecolor="#22262a", edgecolor="#444", labelcolor="white",
                  fontsize=7, loc="lower left")

    for ax in axes[-1, :]:
        ax.set_xlabel("Hz   (gold dashed = body landmark, color dotted = stage pole)",
                      color="white", fontsize=8)

    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, facecolor=fig.get_facecolor())
    plt.close(fig)


def extract_landmarks(body: str) -> list[tuple[float, str]]:
    return {
        "speaker_knockerz": [(60.0, "<60 Hz anchor"), (200.0, "Choke"),
                             (400.0, "Cardboard"), (1200.0, "Cone Cry")],
        "aluminum_siding": [(1000.0, "1k scoop"), (5000.0, "Sibilant"),
                            (8000.0, "Sibilant"), (12000.0, "Tear"),
                            (18000.0, "Dog Whistle")],
        "small_talk_ah_ee": [(300.0, "Yawn F1"), (850.0, "Yawn F2"),
                             (800.0, "Shriek F1"), (2800.0, "Shriek F2")],
        "cul_de_sac": [(250.0, "Root hum"), (800.0, "Comb"), (1900.0, "Comb"),
                       (4500.0, "Comb"), (8000.0, "Comb"), (14000.0, "Comb")],
    }.get(body, [])


# ─── Body designs ────────────────────────────────────────────────────────────
#
# Each body returns (name, c0_uniform, m0_stages, m100_stages, body_brief).
#
# Design discipline:
#   - Same number of active stages at M0 and M100 (linear morph interpolates them).
#   - Anchor/landmark stages keep their pole_hz roughly constant where the body
#     invariant requires it (sub anchor for SK, root hum for CDS).
#   - c0 is uniform across all stages of one cartridge (cascade gain discipline).
#   - Zeros sit interior (zero_r < pole_r) and away from poles; one cartridge =
#     one zero "basket" frequency that all stages share, so zeros multiply
#     constructively into a stable cascade-wide carve.
#

def body_speaker_knockerz():
    name = "speaker_knockerz"
    brief = ("Pressurized sub-harmonic resonator. Sub <60 Hz anchored. "
             "M0=Vault → M100=Fracture (sub + Choke 200 + Cone Cry 1.2k + brutal top). "
             "All radii heritage range 0.998-0.999.")
    c0 = 0.5619  # Talking_Hedz heritage uniform attenuation
    # Zero placement at 4485 Hz, zero_r=0.70 makes single-stage numerator-at-DC
    # = c0 + b1 + b2 ≈ 0 for c0=0.5619 — neutralizes DC contribution from
    # high-freq poles so only the anchor stage dominates the sub region.
    Z_HZ = 4485.0
    Z_R = 0.70
    # M0 — Vault: weighty sub anchor, mid-spread poles, gentle highs.
    m0 = [
        Stage(pole_hz=   50.0, pole_r=0.999, zero_hz=Z_HZ, zero_r=Z_R),  # Sub Anchor
        Stage(pole_hz=  180.0, pole_r=0.998, zero_hz=Z_HZ, zero_r=Z_R),  # Chest
        Stage(pole_hz=  600.0, pole_r=0.998, zero_hz=Z_HZ, zero_r=Z_R),  # Low-Mid
        Stage(pole_hz= 1800.0, pole_r=0.998, zero_hz=Z_HZ, zero_r=Z_R),  # Mid
        Stage(pole_hz= 4500.0, pole_r=0.998, zero_hz=Z_HZ, zero_r=Z_R),  # Air
        Stage(pole_hz=11000.0, pole_r=0.998, zero_hz=Z_HZ, zero_r=Z_R),  # Top
    ]
    # M100 — Total Fracture. Sub anchor identical (invariant). Sub edge tightens.
    # NO pole between 100 and 1200 → 200 Hz Choke is a natural valley.
    # Cone Cry sharp at 1200; rattle/tear/brutal top above.
    m100 = [
        Stage(pole_hz=   50.0, pole_r=0.999, zero_hz=Z_HZ, zero_r=Z_R),  # Sub Anchor (LOCKED)
        Stage(pole_hz=  100.0, pole_r=0.999, zero_hz=Z_HZ, zero_r=Z_R),  # Sub Edge
        Stage(pole_hz= 1200.0, pole_r=0.999, zero_hz=Z_HZ, zero_r=Z_R),  # CONE CRY
        Stage(pole_hz= 2800.0, pole_r=0.998, zero_hz=Z_HZ, zero_r=Z_R),  # Rattle
        Stage(pole_hz= 6000.0, pole_r=0.998, zero_hz=Z_HZ, zero_r=Z_R),  # Tear
        Stage(pole_hz=14000.0, pole_r=0.998, zero_hz=Z_HZ, zero_r=Z_R),  # Brutal top
    ]
    return name, c0, m0, m100, brief


def body_aluminum_siding():
    name = "aluminum_siding"
    brief = ("Brittle high-tension treble stressor. 1 kHz scooped permanently. "
             "M0=Dull Silver (muted phasey top) → M100=Shatter Point (Tear 12k + Dog Whistle 18k).")
    c0 = 0.45
    # All stages share zero at 1 kHz, zero_r=0.85 — cumulative numerator carves
    # a permanent ~30 dB scoop at 1 kHz that doesn't move with morph.
    # Body floor is gentle (single low pole, low r) so DC gain stays bounded.
    m0 = [
        Stage(pole_hz=  150.0, pole_r=0.85, zero_hz=1000.0, zero_r=0.85),  # gentle body floor
        Stage(pole_hz=  500.0, pole_r=0.65, zero_hz=1000.0, zero_r=0.85),  # below scoop, broad
        Stage(pole_hz= 2500.0, pole_r=0.65, zero_hz=1000.0, zero_r=0.85),  # above scoop, broad
        Stage(pole_hz= 5000.0, pole_r=0.62, zero_hz=1000.0, zero_r=0.85),  # mid HF, broad
        Stage(pole_hz= 9000.0, pole_r=0.60, zero_hz=1000.0, zero_r=0.85),  # silvery, broad
        Stage(pole_hz=14000.0, pole_r=0.58, zero_hz=1000.0, zero_r=0.85),  # ceiling
    ]
    # M100 — Shatter Point: body floor kept, mass migration up. Dog Whistle and
    # Aluminium Tear become sharp peaks; sibilance fold notch from cumulative
    # zero at 1 kHz still in place.
    m100 = [
        Stage(pole_hz=  150.0, pole_r=0.85, zero_hz=1000.0, zero_r=0.85),  # floor kept
        Stage(pole_hz=  500.0, pole_r=0.65, zero_hz=1000.0, zero_r=0.85),  # below scoop kept
        Stage(pole_hz= 5000.0, pole_r=0.985, zero_hz=1000.0, zero_r=0.85),  # Sibilant lower
        Stage(pole_hz= 8000.0, pole_r=0.99,  zero_hz=1000.0, zero_r=0.85),  # Sibilant upper
        Stage(pole_hz=12000.0, pole_r=0.992, zero_hz=1000.0, zero_r=0.85),  # ALUMINIUM TEAR
        Stage(pole_hz=18000.0, pole_r=0.992, zero_hz=1000.0, zero_r=0.85),  # DOG WHISTLE
    ]
    return name, c0, m0, m100, brief


def body_small_talk_ah_ee():
    name = "small_talk_ah_ee"
    brief = ("Vocal cavity. Two dominant formants in 200..4000 Hz at all morph. "
             "M0=Yawn (Oo, F1=300 F2=850) → M100=Shriek (Aaa, F1=800 F2=2800).")
    c0 = 0.55
    # Two formants are the body identity. F1, F2 are sharp; the rest are gentle.
    # Following Talking_Hedz: 1 anchor (low) + 2 formants + 3 air/ceiling stages.
    m0 = [
        Stage(pole_hz=  80.0, pole_r=0.85, zero_hz=5500.0, zero_r=0.70),  # gentle sub bed
        Stage(pole_hz= 300.0, pole_r=0.992, zero_hz=5500.0, zero_r=0.70),  # F1 (Oo)
        Stage(pole_hz= 850.0, pole_r=0.985, zero_hz=5500.0, zero_r=0.70),  # F2 (Oo)
        Stage(pole_hz=2500.0, pole_r=0.65, zero_hz=5500.0, zero_r=0.70),  # F3 weak
        Stage(pole_hz=4500.0, pole_r=0.60, zero_hz=5500.0, zero_r=0.70),  # soft brilliance
        Stage(pole_hz=8000.0, pole_r=0.55, zero_hz=5500.0, zero_r=0.70),  # gentle air
    ]
    # M100 — Shriek: same anchor sub bed; F1/F2 climb into Aaa positions;
    # F3 sharpens into rasp; air also tightens.
    m100 = [
        Stage(pole_hz=  80.0, pole_r=0.85, zero_hz=5500.0, zero_r=0.70),  # sub bed kept
        Stage(pole_hz= 800.0, pole_r=0.997, zero_hz=5500.0, zero_r=0.70),  # F1 throated
        Stage(pole_hz=2800.0, pole_r=0.997, zero_hz=5500.0, zero_r=0.70),  # F2 biting
        Stage(pole_hz=4500.0, pole_r=0.92, zero_hz=5500.0, zero_r=0.70),  # F3 rasp
        Stage(pole_hz=6500.0, pole_r=0.88, zero_hz=5500.0, zero_r=0.70),  # vocal fry
        Stage(pole_hz=9500.0, pole_r=0.85, zero_hz=5500.0, zero_r=0.70),  # squeal air
    ]
    return name, c0, m0, m100, brief


def body_cul_de_sac():
    name = "cul_de_sac"
    brief = ("Iron Pipe → Crystal Dust. Root hum at 250 Hz LOCKED across morph. "
             "M0 = single fat resonance, M100 = comb of 5 narrow peaks above the root.")
    c0 = 0.55  # Millennium uses 0.78 but its poles are all 3.6-14.7 kHz; ours
               # has a strong 250 Hz pole, so c0 needs to be lower to keep the
               # cumulative low-freq cascade gain bounded.
    # M0 — Iron Pipe: one dominant peak at 250 Hz; other stages broad/weak.
    m0 = [
        Stage(pole_hz=  250.0, pole_r=0.99,  zero_hz=5500.0, zero_r=0.70),  # IRON PIPE
        Stage(pole_hz=  500.0, pole_r=0.55,  zero_hz=5500.0, zero_r=0.70),  # broad
        Stage(pole_hz= 1200.0, pole_r=0.50,  zero_hz=5500.0, zero_r=0.70),  # weak
        Stage(pole_hz= 3000.0, pole_r=0.48,  zero_hz=5500.0, zero_r=0.70),  # weaker
        Stage(pole_hz= 6000.0, pole_r=0.45,  zero_hz=5500.0, zero_r=0.70),  # barely
        Stage(pole_hz=12000.0, pole_r=0.40,  zero_hz=5500.0, zero_r=0.70),  # ceiling
    ]
    # M100 — Crystal Dust: root locked + 5 sharp comb peaks log-spaced above.
    m100 = [
        Stage(pole_hz=  250.0, pole_r=0.99,  zero_hz=5500.0, zero_r=0.70),  # ROOT LOCKED
        Stage(pole_hz=  800.0, pole_r=0.985, zero_hz=5500.0, zero_r=0.70),  # comb 1
        Stage(pole_hz= 1900.0, pole_r=0.985, zero_hz=5500.0, zero_r=0.70),  # comb 2
        Stage(pole_hz= 4500.0, pole_r=0.985, zero_hz=5500.0, zero_r=0.70),  # comb 3
        Stage(pole_hz= 8000.0, pole_r=0.97,  zero_hz=5500.0, zero_r=0.70),  # comb 4
        Stage(pole_hz=14000.0, pole_r=0.96,  zero_hz=5500.0, zero_r=0.70),  # comb 5
    ]
    return name, c0, m0, m100, brief


BODIES = {
    "speaker_knockerz":  body_speaker_knockerz,
    "aluminum_siding":   body_aluminum_siding,
    "small_talk_ah_ee":  body_small_talk_ah_ee,
    "cul_de_sac":        body_cul_de_sac,
}


# ─── CLI ─────────────────────────────────────────────────────────────────────

def _compute_corner_boost(raw_stages: list[dict], target: float = PEAK_AMPLITUDE_TARGET) -> float:
    """Return boost = target / cascade_peak for a list of RawStage dicts.

    Uses emu_resonator.rs compile formula via rawstage_to_coefs, then evaluates
    |H(jw)| over the full audible band at authoring SR.  Asserts the no-fake-success
    invariant before returning.
    """
    coefs = [rawstage_to_coefs(rs) for rs in raw_stages]
    active = [c for c in coefs if not c.get("_pass")]
    peak = _cascade_peak(active, EVAL_FREQS)
    boost = target / peak if peak > 0.0 else 1.0
    check = peak * boost
    assert check <= target * 1.001, (
        f"NO-FAKE-SUCCESS: peak={peak:.6f} boost={boost:.8f} "
        f"check={check:.6f} > {target * 1.001:.6f}"
    )
    return boost


def emit(body_id: str, out_dir: Path) -> Path:
    if body_id not in BODIES:
        raise SystemExit(f"unknown body: {body_id}. choices: {list(BODIES)} or 'all'")
    name, c0, m0, m100, brief = BODIES[body_id]()
    cart = assemble_cartridge(name, c0, m0, m100)

    # ── Analytical boost normalization ────────────────────────────────────────
    # After stage construction, before JSON write: compute per-corner boost so
    # cascade_peak * boost == PEAK_AMPLITUDE_TARGET (2.5 linear / +7.96 dB).
    # This keeps agc_gain * abs_sample < 5.0 (cliff idx 4→5) even with 2× attack
    # transient overshoot headroom.  Overwrites the placeholder boost=4.0.
    for kf in cart["keyframes"]:
        boost = _compute_corner_boost(kf["stages"])
        kf["boost"] = boost
        print(f"  [boost] {name}  {kf['label']:<12}  boost={boost:.6f}")

    # Update top-level boost to the mean of the 4 corners (informational only;
    # runtime uses per-keyframe boost via interpolate_boost).
    all_boosts = [kf["boost"] for kf in cart["keyframes"]]
    cart["boost"] = sum(all_boosts) / len(all_boosts)
    # ─────────────────────────────────────────────────────────────────────────

    out_dir.mkdir(parents=True, exist_ok=True)
    json_path = out_dir / f"{body_id}.json"
    png_path = out_dir / f"{body_id}.png"
    json_path.write_text(json.dumps(cart, indent=2), encoding="utf-8")
    plot_cartridge(name, c0, m0, m100, cart, png_path, brief)
    return json_path


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("body", help="body id, or 'all'", choices=list(BODIES) + ["all"])
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR,
                    help=f"output directory (default {OUT_DIR})")
    args = ap.parse_args(argv)

    targets = list(BODIES) if args.body == "all" else [args.body]
    for body in targets:
        json_path = emit(body, args.out_dir)
        png_path = json_path.with_suffix(".png")
        print(f"  {body:20s} -> {json_path.name}   plot -> {png_path.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
