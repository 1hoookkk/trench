#!/usr/bin/env python3
"""
TRENCH Phoneme Author Studio
Standalone pole/zero cartridge authoring bench with audition playback,
loaded WAV input, virtual MIDI keyboard, and Mackie-style desk slam audition.

Run:
    python forge/phoneme_author_studio.py

Dependencies:
    Python standard library, numpy, matplotlib, tkinter.

This is not the plugin UI and does not write plugin/runtime code.
"""

from __future__ import annotations

import copy
import json
import math
import os
import sys
import tempfile
import wave
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import numpy as np

import tkinter as tk
from tkinter import ttk, filedialog, messagebox

import matplotlib
matplotlib.use("TkAgg", force=True)
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

try:
    import winsound  # Windows-only, standard library.
except Exception:  # pragma: no cover - non-Windows fallback.
    winsound = None


DEFAULT_SAMPLE_RATE = 39062.5
STAGE_COUNT = 6
FREQ_MIN = 20.0
FREQ_MAX = 18000.0
RESPONSE_POINTS = 1024
EPS = 1e-12

AUTHORING_DIR = Path("cartridges/factory/generated/authoring")
AUDITION_DIR = Path("audition")

SOURCE_SAW = "Saw sweep"
SOURCE_NOISE = "White noise"
SOURCE_IMPULSE = "Click train"
SOURCE_SUB = "808 sub"
SOURCE_MIDI = "Virtual MIDI note"
SOURCE_LOADED = "Loaded WAV input"

MORPH_SWEEP = "Sweep M0→M100"
MORPH_M0 = "Static M0"
MORPH_M50 = "Static M50"
MORPH_M100 = "Static M100"


# Frequency-band recipe knowledge base. Hover the graph; matches show as
# "recipes" in the readout, with citations to BODIES.md targets where they apply.
HOVER_RECIPES: List[Tuple[float, float, str]] = [
    (35, 80, "Sub anchor — Speaker Knockerz Vault floor / 808 fundamental"),
    (90, 170, "Chest bloom — ~150 Hz BODIES.md Chest Resonance station"),
    (180, 230, "Choke band — Speaker Knockerz invariant: ≥12 dB notch @ 200 Hz"),
    (240, 270, "Cul-de-sac root pipe — locked 250 Hz tether through full morph"),
    (270, 360, "Close-vowel F1 — EE/OO F1 (~300 Hz)"),
    (400, 580, "Mid-vowel F1 — schwa, EH; close-mid F1"),
    (600, 800, "Open-vowel F1 — AH (~700), Shriek (~800); OO F2 (~870)"),
    (820, 940, "AH F2 (~1100) approach"),
    (950, 1080, "Aluminum Siding 1 kHz scoop — invariant ≤ −12 dB vs flanking"),
    (1080, 1300, "Cone Cry — Speaker Knockerz: peak ≥10 dB above flanking @ 1.2 kHz"),
    (1300, 1700, "Nasal / lateral anti-formant region"),
    (1700, 2300, "F2 climb — Bite formant rising"),
    (2300, 2900, "EE F2 (~2400) / Shriek F2 (~2800) / Bite presence"),
    (2900, 3400, "F3 — vowel color, retroflex; F3 of EE (~3200)"),
    (3500, 5500, "Sibilant Fold lower (~5 kHz Aluminum Siding)"),
    (6000, 7500, "Glass Stress — ~7 kHz narrow piercing peak"),
    (7500, 9500, "Sibilant Fold upper (~8 kHz)"),
    (9500, 11000, "Sheen — Aluminum Siding ~10 kHz wide air boost"),
    (11000, 14000, "Aluminium Tear — ~12 kHz inharmonic ring"),
    (14000, 19000, "Dog Whistle — Aluminum Siding 18 kHz extreme push"),
]


def lookup_recipes(freq: float) -> List[str]:
    matches = []
    for lo, hi, text in HOVER_RECIPES:
        if lo <= freq <= hi:
            matches.append(text)
    return matches


@dataclass
class Stage:
    stage_index: int
    role_label: str
    pole_freq_hz: float
    pole_radius: float
    zero_freq_hz: float
    zero_radius: float
    c0: float

    @staticmethod
    def from_dict(d: dict) -> "Stage":
        return Stage(
            stage_index=int(d.get("stage_index", 0)),
            role_label=str(d.get("role_label", "Stage")),
            pole_freq_hz=float(d.get("pole_freq_hz", 1000.0)),
            pole_radius=float(d.get("pole_radius", 0.95)),
            zero_freq_hz=float(d.get("zero_freq_hz", 1000.0)),
            zero_radius=float(d.get("zero_radius", 0.6)),
            c0=float(d.get("c0", 0.7)),
        )


def clamp(v: float, lo: float, hi: float) -> float:
    return max(lo, min(hi, v))


def cartridge_slug(name: str) -> str:
    keep = []
    for ch in name.lower():
        if ch.isalnum():
            keep.append(ch)
        elif ch in (" ", "-", "_", "→", "/"):
            keep.append("_")
    slug = "".join(keep)
    while "__" in slug:
        slug = slug.replace("__", "_")
    return slug.strip("_")


def default_cartridges() -> Dict[str, dict]:
    """Formant-targeted starter cartridges.

    Six stage roles, fixed across all cartridges so the bench has anatomy:
        S0 F1 BODY        — first formant (vowel-class anchor: open vs close)
        S1 F2 TONGUE      — second formant (front vs back)
        S2 F3 TEETH       — third formant (vowel color, retroflex tilt)
        S3 LOW HOLLOW     — low-frequency anti-formant pair (cavity / nasality)
        S4 HIGH BITE      — high-frequency anti-formant pair (sibilance / edge)
        S5 AIR            — HF air / breath / stability tail

    Each cartridge carries `formant_targets` for M0 and M100 — those are the
    vertical guide lines the renderer draws, so the eye can immediately verify
    whether F1/F2/F3 cascade peaks actually land where the vowel claims.
    """

    def s(i, role, pf, pr, zf, zr, c0):
        return Stage(i, role, pf, pr, zf, zr, c0)

    carts = {
        "THROAT AH→EE": {
            "name": "THROAT AH→EE",
            "sample_rate": DEFAULT_SAMPLE_RATE,
            "formant_targets": {
                "M0":   {"F1":  750, "F2": 1150, "F3": 2550},  # AH (open back)
                "M100": {"F1":  300, "F2": 2400, "F3": 3250},  # EE (close front)
            },
            "zero_targets": {
                "M0":   {"Z1":  520, "Z2":  900, "Z3": 2100, "ZHOLLOW":  300, "ZBITE": 3600},
                "M100": {"Z1":  220, "Z2": 1850, "Z3": 2850, "ZHOLLOW":  400, "ZBITE": 4650},
            },
            "morph": {
                "M0": [
                    s(0, "F1 BODY",     750, 0.985,  520, 0.55, 0.60),
                    s(1, "F2 TONGUE",  1150, 0.975,  900, 0.55, 0.62),
                    s(2, "F3 TEETH",   2550, 0.965, 2100, 0.58, 0.66),
                    s(3, "LOW HOLLOW",  430, 0.940,  300, 0.82, 0.82),
                    s(4, "HIGH BITE",  4100, 0.945, 3600, 0.70, 0.78),
                    s(5, "AIR",        7600, 0.930, 6800, 0.65, 0.84),
                ],
                "M100": [
                    s(0, "F1 BODY",     300, 0.985,  220, 0.55, 0.58),
                    s(1, "F2 TONGUE",  2400, 0.985, 1850, 0.58, 0.60),
                    s(2, "F3 TEETH",   3250, 0.972, 2850, 0.60, 0.64),
                    s(3, "LOW HOLLOW",  520, 0.935,  400, 0.82, 0.84),
                    s(4, "HIGH BITE",  5200, 0.948, 4650, 0.76, 0.78),
                    s(5, "AIR",        8900, 0.925, 7800, 0.68, 0.86),
                ],
            },
        },
        "DIPHTHONG OO→EE": {
            # The "Oui" specialty articulation: F1 stays low and locked; F2 climbs.
            "name": "DIPHTHONG OO→EE",
            "sample_rate": DEFAULT_SAMPLE_RATE,
            "formant_targets": {
                "M0":   {"F1":  300, "F2":  870, "F3": 2300},  # OO (close back round)
                "M100": {"F1":  300, "F2": 2400, "F3": 3250},  # EE (close front)
            },
            "zero_targets": {
                "M0":   {"Z1":  220, "Z2":  680, "Z3": 1900, "ZHOLLOW":  320, "ZBITE": 3800},
                "M100": {"Z1":  220, "Z2": 1850, "Z3": 2850, "ZHOLLOW":  400, "ZBITE": 4650},
            },
            "morph": {
                "M0": [
                    s(0, "F1 BODY",     300, 0.985,  220, 0.55, 0.58),
                    s(1, "F2 TONGUE",   870, 0.978,  680, 0.55, 0.62),
                    s(2, "F3 TEETH",   2300, 0.962, 1900, 0.58, 0.66),
                    s(3, "LOW HOLLOW",  450, 0.942,  320, 0.84, 0.82),
                    s(4, "HIGH BITE",  4400, 0.940, 3800, 0.70, 0.80),
                    s(5, "AIR",        7800, 0.928, 6900, 0.65, 0.85),
                ],
                "M100": [
                    s(0, "F1 BODY",     300, 0.985,  220, 0.55, 0.58),
                    s(1, "F2 TONGUE",  2400, 0.985, 1850, 0.58, 0.60),
                    s(2, "F3 TEETH",   3250, 0.972, 2850, 0.60, 0.64),
                    s(3, "LOW HOLLOW",  520, 0.935,  400, 0.82, 0.84),
                    s(4, "HIGH BITE",  5200, 0.948, 4650, 0.76, 0.78),
                    s(5, "AIR",        8900, 0.925, 7800, 0.68, 0.86),
                ],
            },
        },
        "OPEN AH→AAA": {
            # Open vowel pushed into Shriek throat tension. F2 is the climber.
            # Aligns with the Small Talk Ah-Ee body's Open Ah → Shriek arc.
            "name": "OPEN AH→AAA",
            "sample_rate": DEFAULT_SAMPLE_RATE,
            "formant_targets": {
                "M0":   {"F1":  720, "F2": 1200, "F3": 2450},  # Open Ah
                "M100": {"F1":  800, "F2": 2800, "F3": 3200},  # Aaa (Shriek)
            },
            "zero_targets": {
                "M0":   {"Z1":  510, "Z2":  960, "Z3": 2050, "ZHOLLOW":  320, "ZBITE": 3500},
                "M100": {"Z1":  540, "Z2": 2200, "Z3": 2750, "ZHOLLOW":  420, "ZBITE": 4800},
            },
            "morph": {
                "M0": [
                    s(0, "F1 BODY",     720, 0.985,  510, 0.55, 0.60),
                    s(1, "F2 TONGUE",  1200, 0.978,  960, 0.55, 0.62),
                    s(2, "F3 TEETH",   2450, 0.962, 2050, 0.58, 0.66),
                    s(3, "LOW HOLLOW",  460, 0.938,  320, 0.80, 0.82),
                    s(4, "HIGH BITE",  4000, 0.942, 3500, 0.68, 0.78),
                    s(5, "AIR",        7400, 0.928, 6600, 0.62, 0.84),
                ],
                "M100": [
                    s(0, "F1 BODY",     800, 0.988,  540, 0.58, 0.60),
                    s(1, "F2 TONGUE",  2800, 0.988, 2200, 0.60, 0.62),
                    s(2, "F3 TEETH",   3200, 0.975, 2750, 0.62, 0.66),
                    s(3, "LOW HOLLOW",  600, 0.945,  420, 0.78, 0.78),
                    s(4, "HIGH BITE",  5400, 0.952, 4800, 0.74, 0.76),
                    s(5, "AIR",        9200, 0.930, 8100, 0.68, 0.84),
                ],
            },
        },
        "ROUND OO→AH": {
            # Back-rounded to open-back. F1 climbs; F2 climbs modestly.
            "name": "ROUND OO→AH",
            "sample_rate": DEFAULT_SAMPLE_RATE,
            "formant_targets": {
                "M0":   {"F1":  300, "F2":  870, "F3": 2300},  # OO
                "M100": {"F1":  750, "F2": 1150, "F3": 2550},  # AH
            },
            "zero_targets": {
                "M0":   {"Z1":  220, "Z2":  680, "Z3": 1900, "ZHOLLOW":  320, "ZBITE": 3800},
                "M100": {"Z1":  520, "Z2":  900, "Z3": 2100, "ZHOLLOW":  300, "ZBITE": 3600},
            },
            "morph": {
                "M0": [
                    s(0, "F1 BODY",     300, 0.985,  220, 0.55, 0.58),
                    s(1, "F2 TONGUE",   870, 0.978,  680, 0.55, 0.62),
                    s(2, "F3 TEETH",   2300, 0.962, 1900, 0.58, 0.66),
                    s(3, "LOW HOLLOW",  450, 0.942,  320, 0.84, 0.82),
                    s(4, "HIGH BITE",  4400, 0.940, 3800, 0.70, 0.80),
                    s(5, "AIR",        7800, 0.928, 6900, 0.65, 0.85),
                ],
                "M100": [
                    s(0, "F1 BODY",     750, 0.985,  520, 0.55, 0.60),
                    s(1, "F2 TONGUE",  1150, 0.975,  900, 0.55, 0.62),
                    s(2, "F3 TEETH",   2550, 0.965, 2100, 0.58, 0.66),
                    s(3, "LOW HOLLOW",  430, 0.940,  300, 0.82, 0.82),
                    s(4, "HIGH BITE",  4100, 0.945, 3600, 0.70, 0.78),
                    s(5, "AIR",        7600, 0.930, 6800, 0.65, 0.84),
                ],
            },
        },
    }

    out = {}
    for label, cart in carts.items():
        out[label] = {
            "name": cart["name"],
            "sample_rate": cart["sample_rate"],
            "formant_targets": copy.deepcopy(cart["formant_targets"]),
            "zero_targets": copy.deepcopy(cart.get("zero_targets", {})),
            "morph": {
                "M0":   [asdict(x) for x in cart["morph"]["M0"]],
                "M100": [asdict(x) for x in cart["morph"]["M100"]],
            },
        }
    return out


def compile_stage(stage: Stage, sample_rate: float) -> Tuple[float, float, float, float, float]:
    """
    Compile pole/zero authoring parameters to the runtime DF2T coefficients.

    Runtime:
        y  = c0*x + w1
        w1 = c1*x - c3*y + w2
        w2 = c2*x - c4*y

    Transfer function:
        H(z) = (c0 + c1 z^-1 + c2 z^-2) / (1 + c3 z^-1 + c4 z^-2)
    """
    theta_p = 2.0 * math.pi * stage.pole_freq_hz / sample_rate
    theta_z = 2.0 * math.pi * stage.zero_freq_hz / sample_rate

    c0 = stage.c0
    c3 = -2.0 * stage.pole_radius * math.cos(theta_p)
    c4 = stage.pole_radius ** 2

    b1 = -2.0 * stage.zero_radius * math.cos(theta_z)
    b2 = stage.zero_radius ** 2

    # c1/c2 scale with c0 so attenuation does not move zero locations.
    c1 = c0 * b1
    c2 = c0 * b2
    return c0, c1, c2, c3, c4


def stage_response(stage: Stage, freqs: np.ndarray, sample_rate: float) -> np.ndarray:
    c0, c1, c2, c3, c4 = compile_stage(stage, sample_rate)
    omega = 2.0 * np.pi * freqs / sample_rate
    z1 = np.exp(-1j * omega)
    z2 = z1 * z1
    num = c0 + c1 * z1 + c2 * z2
    den = 1.0 + c3 * z1 + c4 * z2
    return num / (den + EPS)


def cascade_response(stages: List[Stage], freqs: np.ndarray, sample_rate: float) -> Tuple[np.ndarray, np.ndarray]:
    stage_h = np.array([stage_response(stage, freqs, sample_rate) for stage in stages])
    cascade_h = np.prod(stage_h, axis=0)
    return cascade_h, stage_h


def db(x: np.ndarray) -> np.ndarray:
    return 20.0 * np.log10(np.abs(x) + EPS)


def linear_interp(a: float, b: float, t: float) -> float:
    return (1.0 - t) * a + t * b


def interpolate_stages(a: List[Stage], b: List[Stage], t: float) -> List[Stage]:
    out = []
    for sa, sb in zip(a, b):
        role = sa.role_label if t < 0.5 else sb.role_label
        out.append(
            Stage(
                stage_index=sa.stage_index,
                role_label=role,
                pole_freq_hz=linear_interp(sa.pole_freq_hz, sb.pole_freq_hz, t),
                pole_radius=linear_interp(sa.pole_radius, sb.pole_radius, t),
                zero_freq_hz=linear_interp(sa.zero_freq_hz, sb.zero_freq_hz, t),
                zero_radius=linear_interp(sa.zero_radius, sb.zero_radius, t),
                c0=linear_interp(sa.c0, sb.c0, t),
            )
        )
    return out


def dicts_to_stages(items: List[dict]) -> List[Stage]:
    return [Stage.from_dict(d) for d in items]


def stages_to_dicts(stages: List[Stage]) -> List[dict]:
    return [asdict(s) for s in stages]


def response_metrics(stages: List[Stage], sample_rate: float) -> dict:
    freqs = np.geomspace(FREQ_MIN, FREQ_MAX, RESPONSE_POINTS)
    cascade_h, stage_h = cascade_response(stages, freqs, sample_rate)
    cdb = db(cascade_h)
    sdb = db(stage_h)
    return {
        "freqs": freqs,
        "cascade_db": cdb,
        "stage_db": sdb,
        "peak_db": float(np.max(cdb)),
        "min_db": float(np.min(cdb)),
        "range_db": float(np.max(cdb) - np.min(cdb)),
        "std_db": float(np.std(cdb)),
        "dc_db": float(db(cascade_response(stages, np.array([1.0]), sample_rate)[0])[0]),
    }


def validate_endpoint(label: str, stages: List[Stage], sample_rate: float) -> List[str]:
    warnings = []
    nyquist = sample_rate * 0.5
    for st in stages:
        tag = f"{label} S{st.stage_index} {st.role_label}"
        if st.pole_radius >= 1.0:
            warnings.append(f"{tag}: unstable pole radius >= 1.0")
        if st.zero_radius > 1.0:
            warnings.append(f"{tag}: zero radius > 1.0")
        if not (0.0 < st.c0 <= 1.0):
            warnings.append(f"{tag}: c0 must satisfy 0 < c0 <= 1.0")
        if st.pole_freq_hz <= 0.0 or st.pole_freq_hz >= nyquist:
            warnings.append(f"{tag}: pole frequency outside usable range")
        if st.zero_freq_hz <= 0.0 or st.zero_freq_hz >= nyquist:
            warnings.append(f"{tag}: zero frequency outside usable range")

    metrics = response_metrics(stages, sample_rate)
    peak = metrics["peak_db"]
    rng = metrics["range_db"]
    std = metrics["std_db"]
    dc = metrics["dc_db"]

    if peak > 36.0:
        warnings.append(f"{label}: cascade peak above +36 dB ({peak:.1f} dB)")
    if dc > 9.0 and dc > peak - 3.0:
        warnings.append(f"{label}: DC shelf dominance ({dc:.1f} dB near peak)")
    if rng < 6.0:
        warnings.append(f"{label}: no meaningful peaks/notches (range {rng:.1f} dB)")
    if rng < 3.0 or std < 1.0:
        warnings.append(f"{label}: suspicious flat response")
    return warnings


def make_json_payload(cart: dict) -> dict:
    m0 = copy.deepcopy(cart["morph"]["M0"])
    m100 = copy.deepcopy(cart["morph"]["M100"])
    return {
        "name": cart["name"],
        "sample_rate": float(cart.get("sample_rate", DEFAULT_SAMPLE_RATE)),
        "stage_count": STAGE_COUNT,
        "authoring_type": "phoneme_pole_zero",
        "formant_targets": copy.deepcopy(cart.get("formant_targets", {})),
        "zero_targets": copy.deepcopy(cart.get("zero_targets", {})),
        "morph": {"M0": m0, "M100": m100},
        "corners": {
            "M0_Q0": copy.deepcopy(m0),
            "M0_Q100": copy.deepcopy(m0),
            "M100_Q0": copy.deepcopy(m100),
            "M100_Q100": copy.deepcopy(m100),
        },
    }


def read_wav_mono(path: Path, target_sample_rate: float) -> np.ndarray:
    """Read PCM WAV to mono float. Resamples with simple linear interpolation when needed."""
    with wave.open(str(path), "rb") as wf:
        channels = wf.getnchannels()
        sampwidth = wf.getsampwidth()
        framerate = wf.getframerate()
        frames = wf.getnframes()
        raw = wf.readframes(frames)

    if sampwidth == 1:
        arr = np.frombuffer(raw, dtype=np.uint8).astype(np.float64)
        arr = (arr - 128.0) / 128.0
    elif sampwidth == 2:
        arr = np.frombuffer(raw, dtype="<i2").astype(np.float64) / 32768.0
    elif sampwidth == 3:
        b = np.frombuffer(raw, dtype=np.uint8).reshape(-1, 3)
        vals = b[:, 0].astype(np.int32) | (b[:, 1].astype(np.int32) << 8) | (b[:, 2].astype(np.int32) << 16)
        vals = np.where(vals & 0x800000, vals - 0x1000000, vals)
        arr = vals.astype(np.float64) / 8388608.0
    elif sampwidth == 4:
        arr = np.frombuffer(raw, dtype="<i4").astype(np.float64) / 2147483648.0
    else:
        raise ValueError(f"Unsupported WAV sample width: {sampwidth} bytes")

    if channels > 1:
        arr = arr.reshape(-1, channels).mean(axis=1)

    if framerate != int(round(target_sample_rate)):
        old_x = np.linspace(0.0, 1.0, arr.size, endpoint=False)
        new_n = max(1, int(round(arr.size * target_sample_rate / framerate)))
        new_x = np.linspace(0.0, 1.0, new_n, endpoint=False)
        arr = np.interp(new_x, old_x, arr)

    peak = float(np.max(np.abs(arr))) if arr.size else 0.0
    if peak > 1.0:
        arr = arr / peak
    return arr.astype(np.float64)


def write_wav_mono(path: Path, audio: np.ndarray, sample_rate: float) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    audio = np.asarray(audio, dtype=np.float64)
    audio = np.nan_to_num(audio, nan=0.0, posinf=0.0, neginf=0.0)
    peak = float(np.max(np.abs(audio))) if audio.size else 0.0
    if peak > 0.98:
        audio = audio * (0.98 / peak)
    audio_i16 = np.clip(audio, -1.0, 1.0)
    audio_i16 = (audio_i16 * 32767.0).astype(np.int16)

    with wave.open(str(path), "wb") as wf:
        wf.setnchannels(1)
        wf.setsampwidth(2)
        wf.setframerate(int(round(sample_rate)))
        wf.writeframes(audio_i16.tobytes())


def note_to_freq(midi_note: int) -> float:
    return 440.0 * (2.0 ** ((midi_note - 69) / 12.0))


def source_saw_sweep(n: int, sample_rate: float) -> np.ndarray:
    progress = np.linspace(0.0, 1.0, n)
    freq = 45.0 * (9.0 ** progress)
    phase = np.cumsum(freq / sample_rate)
    return 0.22 * (2.0 * (phase % 1.0) - 1.0)


def source_noise(n: int) -> np.ndarray:
    rng = np.random.default_rng(1337)
    return 0.20 * rng.standard_normal(n)


def source_impulse_train(n: int, sample_rate: float) -> np.ndarray:
    x = np.zeros(n, dtype=np.float64)
    step = max(1, int(round(sample_rate / 9.0)))
    x[::step] = 0.85
    return x


def source_sub(n: int, sample_rate: float) -> np.ndarray:
    x = np.zeros(n, dtype=np.float64)
    hit_len = int(round(sample_rate * 0.75))
    for start in range(0, n, hit_len):
        end = min(n, start + hit_len)
        m = end - start
        tt = np.arange(m) / sample_rate
        env = np.exp(-5.0 * tt)
        freq = 72.0 - 25.0 * (1.0 - np.exp(-18.0 * tt))
        phase = 2.0 * np.pi * np.cumsum(freq / sample_rate)
        x[start:end] += 0.42 * env * np.sin(phase)
    return x


def source_midi_note(n: int, sample_rate: float, midi_note: int) -> np.ndarray:
    f = note_to_freq(midi_note)
    t = np.arange(n) / sample_rate
    env = np.ones(n)
    attack = max(1, int(0.006 * sample_rate))
    release = max(1, int(0.180 * sample_rate))
    env[:attack] = np.linspace(0.0, 1.0, attack)
    env[-release:] *= np.linspace(1.0, 0.0, release)
    phase = np.cumsum(np.full(n, f / sample_rate))
    saw = 2.0 * (phase % 1.0) - 1.0
    square = np.where((phase % 1.0) < 0.5, 1.0, -1.0)
    sub = np.sin(2.0 * np.pi * f * 0.5 * t)
    return 0.18 * env * (0.62 * saw + 0.23 * square + 0.15 * sub)


def fit_audio_length(x: np.ndarray, n: int) -> np.ndarray:
    if x.size == n:
        return x.copy()
    if x.size > n:
        return x[:n].copy()
    reps = int(math.ceil(n / max(1, x.size)))
    return np.tile(x, reps)[:n].copy()


def apply_mackie_desk_slam(audio: np.ndarray, drive_db: float, output_trim_db: float) -> np.ndarray:
    """
    Audition-only harmonic saturation.

    This is not a model of a specific circuit. It is a fast desk-slam style curve:
    gain into asymmetric tanh, mild compression, DC removal, output trim.
    """
    x = np.asarray(audio, dtype=np.float64)
    drive = 10.0 ** (drive_db / 20.0)
    trim = 10.0 ** (output_trim_db / 20.0)
    d = x * drive
    shaped = np.tanh(d + 0.065 * d * d)
    shaped = shaped / (1.0 + 0.22 * np.abs(shaped))
    shaped = shaped - float(np.mean(shaped))
    return shaped * trim


def render_df2t_morph_audio(
    source: np.ndarray,
    m0_stages: List[Stage],
    m100_stages: List[Stage],
    sample_rate: float,
    morph_mode: str,
) -> np.ndarray:
    """Offline path using the exact DF2T runtime equation."""
    n = int(source.size)
    y_out = np.zeros(n, dtype=np.float64)
    state_w1 = np.zeros(STAGE_COUNT, dtype=np.float64)
    state_w2 = np.zeros(STAGE_COUNT, dtype=np.float64)

    if morph_mode == MORPH_M0:
        fixed_t = 0.0
    elif morph_mode == MORPH_M50:
        fixed_t = 0.5
    elif morph_mode == MORPH_M100:
        fixed_t = 1.0
    else:
        fixed_t = None

    for i in range(n):
        morph_t = i / max(1, n - 1) if fixed_t is None else fixed_t
        x = float(source[i])
        stages = interpolate_stages(m0_stages, m100_stages, morph_t)
        for si, st in enumerate(stages):
            c0, c1, c2, c3, c4 = compile_stage(st, sample_rate)
            yy = c0 * x + state_w1[si]
            new_w1 = c1 * x - c3 * yy + state_w2[si]
            new_w2 = c2 * x - c4 * yy
            state_w1[si] = new_w1
            state_w2[si] = new_w2
            x = yy
        if math.isfinite(x):
            y_out[i] = x
        else:
            y_out[i] = 0.0
            state_w1[:] = 0.0
            state_w2[:] = 0.0
    return y_out


class PhonemeAuthorStudio(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("TRENCH Phoneme Author Studio")
        self.geometry("1480x900")
        self.minsize(1260, 760)

        self.bg = "#090b0d"
        self.panel = "#11161b"
        self.panel2 = "#161d23"
        self.panel3 = "#0d1115"
        self.text = "#e6edf2"
        self.subtext = "#8b96a1"
        self.line = "#26313a"
        self.accent = "#8fd3ff"
        self.m50 = "#d7ff9a"
        self.m100 = "#ffb36b"
        self.warn = "#ffc36b"
        self.danger = "#ff7474"

        self.configure(bg=self.bg)
        self._alive = True
        self.protocol("WM_DELETE_WINDOW", self._on_close)
        self.style = ttk.Style(self)
        self._configure_style()

        self.default_state = default_cartridges()
        self.cartridges = copy.deepcopy(self.default_state)
        self.current_cart_label = tk.StringVar(value=list(self.cartridges.keys())[0])
        self.last_cart_label = self.current_cart_label.get()
        self.endpoint_mode = tk.StringVar(value="Edit M0")
        self.last_endpoint_mode = self.endpoint_mode.get()
        self.sample_rate_var = tk.StringVar(value=str(DEFAULT_SAMPLE_RATE))

        self.source_var = tk.StringVar(value=SOURCE_MIDI)
        self.morph_mode_var = tk.StringVar(value=MORPH_SWEEP)
        self.duration_var = tk.StringVar(value="3.0")
        self.slam_enabled_var = tk.BooleanVar(value=False)
        self.slam_drive_var = tk.DoubleVar(value=14.0)
        self.slam_trim_var = tk.DoubleVar(value=-9.0)
        self.last_midi_note = 48  # C3
        self.loaded_wav_path: Optional[Path] = None
        self.loaded_audio: Optional[np.ndarray] = None
        self.preview_path = Path(tempfile.gettempdir()) / "trench_phoneme_author_preview.wav"

        self.formant_shift_var = tk.DoubleVar(value=0.0)
        self.resonance_var = tk.DoubleVar(value=0.0)
        self.zero_bite_var = tk.DoubleVar(value=0.0)
        self.c0_trim_var = tk.DoubleVar(value=0.0)

        self.stage_vars: List[Dict[str, tk.StringVar]] = []
        self.hover_data: Optional[dict] = None
        self._handle_index: List[dict] = []
        self._drag_state: Optional[dict] = None
        self.status_var = tk.StringVar(value="ready")

        self._build_ui()
        self._load_rows_from_state()
        self.bind("<KeyPress>", self._on_key_press)
        self.after(80, self.render_plot)

    def _configure_style(self) -> None:
        try:
            self.style.theme_use("clam")
        except tk.TclError:
            pass
        self.style.configure("TFrame", background=self.bg)
        self.style.configure("Panel.TFrame", background=self.panel)
        self.style.configure("Soft.TFrame", background=self.panel2)
        self.style.configure("Card.TFrame", background=self.panel2)
        self.style.configure("TLabel", background=self.bg, foreground=self.text, font=("Segoe UI", 10))
        self.style.configure("Muted.TLabel", background=self.bg, foreground=self.subtext, font=("Segoe UI", 9))
        self.style.configure("Panel.TLabel", background=self.panel, foreground=self.text, font=("Segoe UI", 10))
        self.style.configure("Soft.TLabel", background=self.panel2, foreground=self.text, font=("Segoe UI", 10))
        self.style.configure("MutedSoft.TLabel", background=self.panel2, foreground=self.subtext, font=("Segoe UI", 9))
        self.style.configure("Header.TLabel", background=self.bg, foreground=self.text, font=("Segoe UI", 15, "bold"))
        self.style.configure("SmallHeader.TLabel", background=self.panel, foreground=self.text, font=("Segoe UI", 10, "bold"))
        self.style.configure("TButton", background=self.panel2, foreground=self.text, padding=(8, 5), relief="flat")
        self.style.map("TButton", background=[("active", "#25313b")])
        self.style.configure("TCombobox", fieldbackground=self.panel3, background=self.panel2, foreground=self.text)
        self.style.configure("TRadiobutton", background=self.panel, foreground=self.text)
        self.style.configure("TCheckbutton", background=self.panel, foreground=self.text)
        self.style.configure("Horizontal.TScale", background=self.panel, troughcolor=self.panel3)
        self.style.configure("TEntry", fieldbackground=self.panel3, foreground=self.text, insertcolor=self.text)

    def _entry(self, parent, var, width=9):
        return tk.Entry(
            parent,
            textvariable=var,
            width=width,
            bg=self.panel3,
            fg=self.text,
            insertbackground=self.text,
            relief="flat",
            highlightthickness=1,
            highlightbackground=self.line,
            highlightcolor=self.accent,
            font=("Consolas", 9),
        )

    def _build_ui(self) -> None:
        header = ttk.Frame(self, style="TFrame")
        header.pack(fill="x", padx=14, pady=(10, 6))
        ttk.Label(header, text="TRENCH PHONEME AUTHOR", style="Header.TLabel").pack(side="left")
        ttk.Label(header, textvariable=self.status_var, style="Muted.TLabel").pack(side="right", pady=(4, 0))

        top = ttk.Frame(self, style="Panel.TFrame")
        top.pack(fill="x", padx=14, pady=(0, 8))
        self._build_top_bar(top)

        body = ttk.Frame(self, style="TFrame")
        body.pack(fill="both", expand=True, padx=14, pady=(0, 8))

        left = ttk.Frame(body, style="Panel.TFrame")
        left.pack(side="left", fill="y", padx=(0, 8))
        left.configure(width=430)
        left.pack_propagate(False)

        center = ttk.Frame(body, style="Panel.TFrame")
        center.pack(side="left", fill="both", expand=True, padx=(0, 8))

        right = ttk.Frame(body, style="Panel.TFrame")
        right.pack(side="left", fill="y")
        right.configure(width=330)
        right.pack_propagate(False)

        self._build_stage_editor(left)
        self._build_plot(center)
        self._build_audio_panel(right)

    def _build_top_bar(self, parent: ttk.Frame) -> None:
        inner = ttk.Frame(parent, style="Panel.TFrame")
        inner.pack(fill="x", padx=10, pady=10)

        ttk.Label(inner, text="cartridge", style="Panel.TLabel").grid(row=0, column=0, sticky="w")
        cart_menu = ttk.Combobox(
            inner,
            textvariable=self.current_cart_label,
            values=list(self.cartridges.keys()),
            state="readonly",
            width=25,
        )
        cart_menu.grid(row=1, column=0, sticky="ew", pady=(3, 0))
        cart_menu.bind("<<ComboboxSelected>>", self._on_cartridge_changed)

        ttk.Label(inner, text="endpoint", style="Panel.TLabel").grid(row=0, column=1, sticky="w", padx=(16, 0))
        mode_frame = ttk.Frame(inner, style="Panel.TFrame")
        mode_frame.grid(row=1, column=1, sticky="w", padx=(16, 0), pady=(3, 0))
        for i, label in enumerate(("Edit M0", "Edit M100", "View both")):
            rb = ttk.Radiobutton(mode_frame, text=label, variable=self.endpoint_mode, value=label, command=self._on_endpoint_mode_changed)
            rb.grid(row=0, column=i, sticky="w", padx=(0, 10))

        ttk.Label(inner, text="sample rate", style="Panel.TLabel").grid(row=0, column=2, sticky="w", padx=(16, 0))
        sr = self._entry(inner, self.sample_rate_var, width=12)
        sr.grid(row=1, column=2, sticky="w", padx=(16, 0), pady=(3, 0), ipady=3)

        action_frame = ttk.Frame(inner, style="Panel.TFrame")
        action_frame.grid(row=1, column=3, sticky="e", padx=(20, 0), pady=(3, 0))
        for label, cmd in (
            ("Render", self.render_plot),
            ("Calibrate", self.calibrate_c0),
            ("Save JSON", self.save_cartridge_json),
            ("Load JSON", self.load_cartridge_json),
            ("Reset", self.reset_current_cartridge),
        ):
            ttk.Button(action_frame, text=label, command=cmd).pack(side="left", padx=3)

        inner.columnconfigure(0, weight=1)
        inner.columnconfigure(3, weight=1)

    def _build_stage_editor(self, parent: ttk.Frame) -> None:
        head = ttk.Frame(parent, style="Panel.TFrame")
        head.pack(fill="x", padx=10, pady=(10, 6))
        ttk.Label(head, text="stage cards", style="SmallHeader.TLabel").pack(side="left")
        ttk.Label(head, text="numbers are editable", style="Panel.TLabel").pack(side="right")

        card_area = ttk.Frame(parent, style="Panel.TFrame")
        card_area.pack(fill="both", expand=True, padx=10, pady=(0, 8))

        for i in range(STAGE_COUNT):
            row_vars = {
                "stage_index": tk.StringVar(value=str(i)),
                "role_label": tk.StringVar(),
                "pole_freq_hz": tk.StringVar(),
                "pole_radius": tk.StringVar(),
                "zero_freq_hz": tk.StringVar(),
                "zero_radius": tk.StringVar(),
                "c0": tk.StringVar(),
            }
            self.stage_vars.append(row_vars)
            card = tk.Frame(card_area, bg=self.panel2, highlightthickness=1, highlightbackground=self.line)
            card.grid(row=i, column=0, sticky="ew", pady=4)
            card.columnconfigure(1, weight=1)
            card.columnconfigure(3, weight=1)

            title = tk.Label(card, text=f"S{i}", bg=self.panel2, fg=self.accent, font=("Segoe UI", 10, "bold"))
            title.grid(row=0, column=0, sticky="w", padx=8, pady=(7, 2))
            role = self._entry(card, row_vars["role_label"], width=16)
            role.grid(row=0, column=1, columnspan=3, sticky="ew", padx=(4, 8), pady=(7, 2), ipady=2)

            labels = [("pole", "pole_freq_hz"), ("r", "pole_radius"), ("zero", "zero_freq_hz"), ("zr", "zero_radius"), ("c0", "c0")]
            for idx, (lab, key) in enumerate(labels):
                r = 1 + idx // 2
                c = (idx % 2) * 2
                tk.Label(card, text=lab, bg=self.panel2, fg=self.subtext, font=("Segoe UI", 8)).grid(row=r, column=c, sticky="w", padx=(8, 2), pady=2)
                ent = self._entry(card, row_vars[key], width=9)
                ent.grid(row=r, column=c + 1, sticky="ew", padx=(2, 8), pady=2, ipady=2)

        macro = ttk.Frame(parent, style="Panel.TFrame")
        macro.pack(fill="x", padx=10, pady=(0, 10))
        ttk.Label(macro, text="quick sculpt selected endpoint", style="SmallHeader.TLabel").pack(anchor="w", pady=(0, 4))
        self._macro_slider(macro, "formant shift", self.formant_shift_var, -12, 12, " semitones")
        self._macro_slider(macro, "resonance", self.resonance_var, -50, 50, "")
        self._macro_slider(macro, "zero bite", self.zero_bite_var, -50, 50, "")
        self._macro_slider(macro, "c0 trim", self.c0_trim_var, -18, 6, " dB")
        btns = ttk.Frame(macro, style="Panel.TFrame")
        btns.pack(fill="x", pady=(6, 0))
        ttk.Button(btns, text="Apply sculpt", command=self.apply_sculpt).pack(side="left", expand=True, fill="x", padx=(0, 3))
        ttk.Button(btns, text="Copy M0→M100", command=lambda: self.copy_endpoint("M0", "M100")).pack(side="left", expand=True, fill="x", padx=3)
        ttk.Button(btns, text="Swap", command=self.swap_endpoints).pack(side="left", expand=True, fill="x", padx=(3, 0))

    def _macro_slider(self, parent, label, var, lo, hi, suffix):
        row = ttk.Frame(parent, style="Panel.TFrame")
        row.pack(fill="x", pady=1)
        txt = ttk.Label(row, text=label, style="Panel.TLabel", width=13)
        txt.pack(side="left")
        scale = ttk.Scale(row, from_=lo, to=hi, variable=var, orient="horizontal")
        scale.pack(side="left", fill="x", expand=True, padx=5)
        val_label = ttk.Label(row, text="", style="Panel.TLabel", width=9)
        val_label.pack(side="right")

        def update(*_):
            val_label.configure(text=f"{var.get():+.1f}{suffix}")
        var.trace_add("write", update)
        update()

    def _build_plot(self, parent: ttk.Frame) -> None:
        self.fig = Figure(figsize=(7.5, 5.3), dpi=100, facecolor=self.panel)
        self.ax = self.fig.add_subplot(111)
        self.fig.subplots_adjust(left=0.085, right=0.985, top=0.94, bottom=0.12)
        self.canvas = FigureCanvasTkAgg(self.fig, master=parent)
        widget = self.canvas.get_tk_widget()
        widget.configure(bg=self.panel, highlightthickness=0)
        widget.pack(fill="both", expand=True, padx=8, pady=(8, 0))
        self._hover_cid = self.canvas.mpl_connect("motion_notify_event", self._on_plot_hover)
        self._press_cid = self.canvas.mpl_connect("button_press_event", self._on_plot_press)
        self._release_cid = self.canvas.mpl_connect("button_release_event", self._on_plot_release)

        bottom = ttk.Frame(parent, style="Panel.TFrame")
        bottom.pack(fill="x", padx=8, pady=8)
        self.info_text = tk.Text(
            bottom,
            height=9,
            wrap="word",
            bg=self.panel3,
            fg=self.text,
            insertbackground=self.text,
            relief="flat",
            highlightthickness=1,
            highlightbackground=self.line,
            font=("Consolas", 9),
        )
        self.info_text.pack(fill="x")
        self.info_text.configure(state="disabled")

    def _build_audio_panel(self, parent: ttk.Frame) -> None:
        wrap = ttk.Frame(parent, style="Panel.TFrame")
        wrap.pack(fill="both", expand=True, padx=10, pady=10)
        ttk.Label(wrap, text="audition", style="SmallHeader.TLabel").pack(anchor="w")

        self._combo_row(wrap, "source", self.source_var, [SOURCE_MIDI, SOURCE_SAW, SOURCE_NOISE, SOURCE_IMPULSE, SOURCE_SUB, SOURCE_LOADED])
        self._combo_row(wrap, "morph", self.morph_mode_var, [MORPH_SWEEP, MORPH_M0, MORPH_M50, MORPH_M100])
        self._entry_row(wrap, "duration", self.duration_var, "seconds")

        load_row = ttk.Frame(wrap, style="Panel.TFrame")
        load_row.pack(fill="x", pady=(6, 2))
        ttk.Button(load_row, text="Load input WAV", command=self.load_input_wav).pack(side="left", fill="x", expand=True)
        self.loaded_label = ttk.Label(load_row, text="none", style="Panel.TLabel")
        self.loaded_label.pack(side="right", padx=(8, 0))

        slam_box = tk.Frame(wrap, bg=self.panel, highlightthickness=1, highlightbackground=self.line)
        slam_box.pack(fill="x", pady=(10, 8))
        ttk.Checkbutton(slam_box, text="Mackie desk slam", variable=self.slam_enabled_var, command=self._update_slam_status).pack(anchor="w", padx=8, pady=(7, 2))
        self._slam_slider(slam_box, "drive", self.slam_drive_var, 0, 36, " dB")
        self._slam_slider(slam_box, "trim", self.slam_trim_var, -24, 0, " dB")
        ttk.Label(slam_box, text="audition-only saturation before the DF2T cascade", style="Panel.TLabel").pack(anchor="w", padx=8, pady=(0, 7))

        play_row = ttk.Frame(wrap, style="Panel.TFrame")
        play_row.pack(fill="x", pady=4)
        ttk.Button(play_row, text="Play preview", command=self.play_preview).pack(side="left", expand=True, fill="x", padx=(0, 3))
        ttk.Button(play_row, text="Stop", command=self.stop_audio).pack(side="left", expand=True, fill="x", padx=(3, 0))

        export_row = ttk.Frame(wrap, style="Panel.TFrame")
        export_row.pack(fill="x", pady=4)
        ttk.Button(export_row, text="Export current WAV", command=self.export_current_wav).pack(side="left", expand=True, fill="x", padx=(0, 3))
        ttk.Button(export_row, text="Export all", command=self.export_wav_audition).pack(side="left", expand=True, fill="x", padx=(3, 0))

        ttk.Label(wrap, text="virtual MIDI keyboard", style="SmallHeader.TLabel").pack(anchor="w", pady=(14, 4))
        ttk.Label(wrap, text="click keys or use computer keys: z/s/x/d/c/v/g/b/h/n/j/m", style="Panel.TLabel").pack(anchor="w", pady=(0, 5))
        self._build_midi_keyboard(wrap)

    def _combo_row(self, parent, label, var, values):
        row = ttk.Frame(parent, style="Panel.TFrame")
        row.pack(fill="x", pady=(7, 0))
        ttk.Label(row, text=label, style="Panel.TLabel", width=9).pack(side="left")
        cb = ttk.Combobox(row, textvariable=var, values=values, state="readonly")
        cb.pack(side="left", fill="x", expand=True)

    def _entry_row(self, parent, label, var, suffix):
        row = ttk.Frame(parent, style="Panel.TFrame")
        row.pack(fill="x", pady=(7, 0))
        ttk.Label(row, text=label, style="Panel.TLabel", width=9).pack(side="left")
        ent = self._entry(row, var, width=8)
        ent.pack(side="left", fill="x", expand=True)
        ttk.Label(row, text=suffix, style="Panel.TLabel").pack(side="right", padx=(6, 0))

    def _slam_slider(self, parent, label, var, lo, hi, suffix):
        row = ttk.Frame(parent, style="Panel.TFrame")
        row.pack(fill="x", padx=8, pady=2)
        ttk.Label(row, text=label, style="Panel.TLabel", width=6).pack(side="left")
        scale = ttk.Scale(row, from_=lo, to=hi, variable=var, orient="horizontal")
        scale.pack(side="left", fill="x", expand=True, padx=5)
        val = ttk.Label(row, text="", style="Panel.TLabel", width=8)
        val.pack(side="right")
        def update(*_):
            val.configure(text=f"{var.get():.1f}{suffix}")
        var.trace_add("write", update)
        update()

    def _build_midi_keyboard(self, parent):
        keys = [
            (48, "C2", False), (49, "C#", True), (50, "D", False), (51, "D#", True),
            (52, "E", False), (53, "F", False), (54, "F#", True), (55, "G", False),
            (56, "G#", True), (57, "A", False), (58, "A#", True), (59, "B", False),
            (60, "C3", False), (61, "C#", True), (62, "D", False), (63, "D#", True),
            (64, "E", False), (65, "F", False), (66, "F#", True), (67, "G", False),
            (68, "G#", True), (69, "A", False), (70, "A#", True), (71, "B", False),
            (72, "C4", False),
        ]
        grid = tk.Frame(parent, bg=self.panel)
        grid.pack(fill="x")
        for i, (note, label, is_black) in enumerate(keys):
            bg = "#050607" if is_black else "#e8edf0"
            fg = "#e8edf0" if is_black else "#101418"
            btn = tk.Button(
                grid,
                text=label,
                bg=bg,
                fg=fg,
                activebackground=self.accent,
                activeforeground="#050607",
                relief="flat",
                width=3,
                height=3 if not is_black else 2,
                command=lambda n=note: self.play_midi_note(n),
            )
            btn.grid(row=0, column=i, padx=1, sticky="nsew")
            grid.columnconfigure(i, weight=1)

    def _safe_canvas_draw(self) -> None:
        if not getattr(self, "_alive", False):
            return
        try:
            self.canvas.draw()
            self.canvas.flush_events()
        except tk.TclError:
            return

    def _on_close(self) -> None:
        self._alive = False
        self.stop_audio()
        try:
            if hasattr(self, "canvas") and hasattr(self, "_hover_cid"):
                self.canvas.mpl_disconnect(self._hover_cid)
            if hasattr(self, "canvas"):
                self.canvas.get_tk_widget().destroy()
        except Exception:
            pass
        self.quit()
        self.destroy()

    def _parse_sample_rate(self) -> float:
        sr = float(self.sample_rate_var.get())
        if sr <= 1000:
            raise ValueError("Sample rate must be above 1000 Hz.")
        return sr

    def _parse_duration(self) -> float:
        dur = float(self.duration_var.get())
        return clamp(dur, 0.15, 12.0)

    def _rows_to_stages(self) -> List[Stage]:
        stages = []
        errors = []
        for i, row in enumerate(self.stage_vars):
            try:
                stages.append(Stage(
                    stage_index=i,
                    role_label=row["role_label"].get().strip() or f"Stage {i}",
                    pole_freq_hz=float(row["pole_freq_hz"].get()),
                    pole_radius=float(row["pole_radius"].get()),
                    zero_freq_hz=float(row["zero_freq_hz"].get()),
                    zero_radius=float(row["zero_radius"].get()),
                    c0=float(row["c0"].get()),
                ))
            except ValueError:
                errors.append(f"S{i}: could not parse numeric values")
        if errors:
            raise ValueError("\n".join(errors))
        return stages

    def _stages_to_rows(self, stages: List[Stage]) -> None:
        for row, st in zip(self.stage_vars, stages):
            row["stage_index"].set(str(st.stage_index))
            row["role_label"].set(st.role_label)
            row["pole_freq_hz"].set(f"{st.pole_freq_hz:.3f}".rstrip("0").rstrip("."))
            row["pole_radius"].set(f"{st.pole_radius:.5f}".rstrip("0").rstrip("."))
            row["zero_freq_hz"].set(f"{st.zero_freq_hz:.3f}".rstrip("0").rstrip("."))
            row["zero_radius"].set(f"{st.zero_radius:.5f}".rstrip("0").rstrip("."))
            row["c0"].set(f"{st.c0:.5f}".rstrip("0").rstrip("."))

    def _current_cart(self, label: Optional[str] = None) -> dict:
        return self.cartridges[label or self.current_cart_label.get()]

    def _save_rows_to_state(self, cart_label: Optional[str] = None, mode: Optional[str] = None) -> None:
        mode = mode or self.endpoint_mode.get()
        if mode not in ("Edit M0", "Edit M100"):
            return
        endpoint = "M0" if mode == "Edit M0" else "M100"
        cart = self._current_cart(cart_label)
        cart["sample_rate"] = self._parse_sample_rate()
        cart["morph"][endpoint] = stages_to_dicts(self._rows_to_stages())

    def _load_rows_from_state(self) -> None:
        cart = self._current_cart()
        self.sample_rate_var.set(str(cart.get("sample_rate", DEFAULT_SAMPLE_RATE)))
        endpoint = "M100" if self.endpoint_mode.get() == "Edit M100" else "M0"
        self._stages_to_rows(dicts_to_stages(cart["morph"][endpoint]))
        self._set_row_editing(self.endpoint_mode.get() != "View both")

    def _set_row_editing(self, enabled: bool) -> None:
        def visit(w):
            for child in w.winfo_children():
                if isinstance(child, tk.Entry):
                    child.configure(state="normal" if enabled else "disabled")
                visit(child)
        visit(self)
        # Keep always-edit fields active.
        for var in (self.sample_rate_var, self.duration_var):
            pass

    def _on_cartridge_changed(self, _event=None) -> None:
        old_label = self.last_cart_label
        new_label = self.current_cart_label.get()
        try:
            if old_label in self.cartridges and self.endpoint_mode.get() in ("Edit M0", "Edit M100"):
                self._save_rows_to_state(old_label)
        except Exception as exc:
            self.current_cart_label.set(old_label)
            messagebox.showerror("Cartridge switch blocked", str(exc))
            return
        self.last_cart_label = new_label
        self._load_rows_from_state()
        self.render_plot()

    def _on_endpoint_mode_changed(self) -> None:
        old_mode = self.last_endpoint_mode
        new_mode = self.endpoint_mode.get()
        try:
            if old_mode in ("Edit M0", "Edit M100"):
                self._save_rows_to_state(mode=old_mode)
            self.last_endpoint_mode = new_mode
            if new_mode in ("Edit M0", "Edit M100"):
                self._load_rows_from_state()
            else:
                self._set_row_editing(False)
            self.render_plot()
        except Exception as exc:
            self.endpoint_mode.set(old_mode)
            messagebox.showerror("Endpoint change failed", str(exc))

    def _set_info_text(self, text: str) -> None:
        self.info_text.configure(state="normal")
        self.info_text.delete("1.0", "end")
        self.info_text.insert("1.0", text)
        self.info_text.configure(state="disabled")

    def _collect_validation_text(self, sample_rate: float, m0: List[Stage], m100: List[Stage]) -> str:
        warnings = []
        warnings.extend(validate_endpoint("M0", m0, sample_rate))
        warnings.extend(validate_endpoint("M100", m100, sample_rate))
        m50 = interpolate_stages(m0, m100, 0.5)
        warnings.extend(validate_endpoint("M50", m50, sample_rate))
        m0m = response_metrics(m0, sample_rate)
        m100m = response_metrics(m100, sample_rate)
        m50m = response_metrics(m50, sample_rate)
        lines = [
            f"SR {sample_rate:.1f} Hz | {self.source_var.get()} | {self.morph_mode_var.get()}",
            f"M0   peak {m0m['peak_db']:+.1f} dB | range {m0m['range_db']:.1f} dB",
            f"M50  peak {m50m['peak_db']:+.1f} dB | range {m50m['range_db']:.1f} dB",
            f"M100 peak {m100m['peak_db']:+.1f} dB | range {m100m['range_db']:.1f} dB",
            "",
        ]
        if warnings:
            lines.append("WARNINGS")
            lines.extend([f"- {w}" for w in warnings])
        else:
            lines.append("safe: no validation warnings")
        return "\n".join(lines)

    def render_plot(self) -> None:
        try:
            # Skip syncing rows-to-state during a drag — the cart is the
            # source of truth and the rows haven't caught up yet.
            if self._drag_state is None and self.endpoint_mode.get() in ("Edit M0", "Edit M100"):
                self._save_rows_to_state()
            sample_rate = self._parse_sample_rate()
            cart = self._current_cart()
            m0 = dicts_to_stages(cart["morph"]["M0"])
            m100 = dicts_to_stages(cart["morph"]["M100"])
        except Exception as exc:
            messagebox.showerror("Render failed", str(exc))
            return

        freqs = np.geomspace(FREQ_MIN, FREQ_MAX, RESPONSE_POINTS)
        self.ax.clear()
        self.ax.set_facecolor(self.panel3)
        self.ax.set_xscale("log")
        self.ax.set_xlim(FREQ_MIN, FREQ_MAX)
        self.ax.set_ylim(-72, 42)
        self.ax.set_xlabel("frequency Hz", color=self.subtext)
        self.ax.set_ylabel("dB", color=self.subtext)
        self.ax.tick_params(axis="both", colors=self.subtext, labelsize=9)
        for spine in self.ax.spines.values():
            spine.set_color(self.line)
        self.ax.grid(True, which="major", color="#27313a", alpha=0.75, linewidth=0.8)
        self.ax.grid(True, which="minor", color="#1d252c", alpha=0.55, linewidth=0.5)
        self.ax.axhline(0.0, color="#4b5964", linewidth=0.8)
        self.ax.set_title(cart["name"], color=self.text, fontsize=12, fontweight="bold")

        mode = self.endpoint_mode.get()
        if mode == "Edit M0":
            line_sets = [("M0", m0, self.accent, "-", 2.6)]
        elif mode == "Edit M100":
            line_sets = [("M100", m100, self.m100, "-", 2.6)]
        else:
            line_sets = [
                ("M0", m0, self.accent, "-", 2.3),
                ("M50", interpolate_stages(m0, m100, 0.5), self.m50, "-.", 2.0),
                ("M100", m100, self.m100, "--", 2.3),
            ]

        # Formant target guide lines — vertical markers at the vowel's
        # claimed F1/F2/F3 so the eye can verify cascade peaks land there.
        formant_targets = cart.get("formant_targets") or {}
        target_color = {"M0": self.accent, "M100": self.m100}
        if mode == "Edit M0":
            visible_targets = {"M0": formant_targets.get("M0", {})}
        elif mode == "Edit M100":
            visible_targets = {"M100": formant_targets.get("M100", {})}
        else:
            visible_targets = formant_targets
        for endpoint, formants in visible_targets.items():
            color = target_color.get(endpoint, self.subtext)
            for fname, freq in (formants or {}).items():
                try:
                    f = float(freq)
                except (TypeError, ValueError):
                    continue
                if f <= 0:
                    continue
                self.ax.axvline(f, color=color, linewidth=1.0, alpha=0.45, linestyle=":")
                y_anchor = 38.0 if endpoint == "M0" else 30.0
                self.ax.text(
                    f, y_anchor, f"{endpoint} {fname}",
                    color=color, fontsize=8, alpha=0.85,
                    rotation=90, ha="right", va="top",
                )

        # Zero target guide lines — dashed dimmed cousins of formant guides.
        # These mark intended anti-formant / notch placements per endpoint.
        zero_targets = cart.get("zero_targets") or {}
        zero_color = {"M0": "#5a8fc4", "M100": "#c4885a"}
        if mode == "Edit M0":
            visible_zeros = {"M0": zero_targets.get("M0", {})}
        elif mode == "Edit M100":
            visible_zeros = {"M100": zero_targets.get("M100", {})}
        else:
            visible_zeros = zero_targets
        for endpoint, zeros in visible_zeros.items():
            color = zero_color.get(endpoint, self.subtext)
            for zname, freq in (zeros or {}).items():
                try:
                    f = float(freq)
                except (TypeError, ValueError):
                    continue
                if f <= 0:
                    continue
                self.ax.axvline(f, color=color, linewidth=0.8, alpha=0.30, linestyle="--")
                y_anchor = -52.0 if endpoint == "M0" else -60.0
                self.ax.text(
                    f, y_anchor, f"{endpoint} {zname}",
                    color=color, fontsize=7, alpha=0.7,
                    rotation=90, ha="right", va="bottom",
                )

        hover_choice = None
        for label, stages, color, linestyle, width in line_sets:
            cascade_h, stage_h = cascade_response(stages, freqs, sample_rate)
            cdb = np.clip(db(cascade_h), -120, 120)
            sdb = np.clip(db(stage_h), -120, 120)
            for si in range(stage_h.shape[0]):
                self.ax.plot(freqs, sdb[si], color=color, linestyle=":", linewidth=0.75, alpha=0.23)
            self.ax.plot(freqs, cdb, color=color, linestyle=linestyle, linewidth=width, label=label)
            if label == "M50" or hover_choice is None:
                hover_choice = {"label": label, "freqs": freqs, "cascade_db": cdb, "stage_db": sdb, "stages": stages}
        # Draggable stage handles — Pole = X (top), Zero = O (bottom).
        # Only shown for the active edit endpoint so drag target is unambiguous.
        # In "View both" mode no handles render and dragging is disabled.
        self._handle_index = []
        if mode == "Edit M0":
            handle_sets = [("M0", m0, self.accent)]
        elif mode == "Edit M100":
            handle_sets = [("M100", m100, self.m100)]
        else:
            handle_sets = []
        for endpoint, stages, color in handle_sets:
            for st in stages:
                self.ax.scatter(
                    [st.pole_freq_hz], [6.0], marker="X",
                    s=120, c=color, edgecolors="white",
                    linewidth=1.2, alpha=0.95, zorder=5,
                )
                self.ax.scatter(
                    [st.zero_freq_hz], [-26.0], marker="o",
                    s=80, c="none", edgecolors=color,
                    linewidth=1.4, alpha=0.95, zorder=5,
                )
                self.ax.annotate(
                    f"S{st.stage_index}", (st.pole_freq_hz, 6.0),
                    textcoords="offset points", xytext=(7, 5),
                    color=color, fontsize=8, alpha=0.9,
                )
                self._handle_index.append({"stage_index": st.stage_index, "kind": "pole",
                                           "endpoint": endpoint, "freq": st.pole_freq_hz})
                self._handle_index.append({"stage_index": st.stage_index, "kind": "zero",
                                           "endpoint": endpoint, "freq": st.zero_freq_hz})

        self.ax.legend(loc="upper right", frameon=False, labelcolor=self.text)
        self.hover_data = hover_choice
        self._safe_canvas_draw()
        self._set_info_text(self._collect_validation_text(sample_rate, m0, m100))
        self.status_var.set("rendered")

    def _on_plot_hover(self, event) -> None:
        if event.inaxes != self.ax or event.xdata is None:
            return
        # Drag mode: update the active stage's pole or zero frequency.
        if self._drag_state is not None:
            self._apply_drag(event.xdata)
            return
        if self.hover_data is None:
            return
        freqs = self.hover_data["freqs"]
        idx = int(np.argmin(np.abs(np.log(freqs) - math.log(max(FREQ_MIN, event.xdata)))))
        freq = freqs[idx]
        cdb = self.hover_data["cascade_db"][idx]
        sdb = self.hover_data["stage_db"][:, idx]
        stages = self.hover_data["stages"]
        lines = [
            f"hover {self.hover_data['label']} | {freq:8.1f} Hz | cascade {cdb:+7.2f} dB",
            "",
        ]
        for st, val in zip(stages, sdb):
            lines.append(f"S{st.stage_index} {st.role_label[:12]:12s} {val:+8.2f} dB")
        recipes = lookup_recipes(freq)
        if recipes:
            lines.append("")
            lines.append("recipes:")
            for r in recipes:
                lines.append(f"  - {r}")
        nearest = self._nearest_handle(event.xdata)
        if nearest is not None and nearest["distance"] < 0.06:
            lines.append("")
            lines.append(f"drag handle: S{nearest['stage_index']} "
                         f"{nearest['kind']} ({nearest['endpoint']}) "
                         f"@ {nearest['freq']:.1f} Hz  -- click+drag to move")
        self._set_info_text("\n".join(lines))

    def _nearest_handle(self, x_data: float) -> Optional[dict]:
        """Find the handle whose log-frequency is closest to x_data."""
        if not self._handle_index or x_data is None or x_data <= 0:
            return None
        target = math.log10(max(1.0, x_data))
        best = None
        best_dist = math.inf
        for h in self._handle_index:
            d = abs(math.log10(max(1.0, h["freq"])) - target)
            if d < best_dist:
                best_dist = d
                best = dict(h, distance=d)
        return best

    def _on_plot_press(self, event) -> None:
        if event.inaxes != self.ax or event.button != 1 or event.xdata is None:
            return
        nearest = self._nearest_handle(event.xdata)
        if nearest is None or nearest["distance"] >= 0.06:
            return
        # Capture which stage handle is being dragged.
        self._drag_state = {
            "stage_index": nearest["stage_index"],
            "kind": nearest["kind"],
            "endpoint": nearest["endpoint"],
        }
        self.status_var.set(
            f"dragging S{nearest['stage_index']} {nearest['kind']} ({nearest['endpoint']})"
        )

    def _on_plot_release(self, event) -> None:
        if self._drag_state is None:
            return
        if event.xdata is not None:
            self._apply_drag(event.xdata)
        self._drag_state = None
        # Sync row entries to the freshly mutated cart.
        self._load_rows_from_state()
        self.status_var.set("drag committed")

    def _apply_drag(self, x_data: float) -> None:
        if self._drag_state is None or x_data is None:
            return
        try:
            sample_rate = self._parse_sample_rate()
        except Exception:
            return
        nyq = sample_rate * 0.5
        new_freq = clamp(float(x_data), FREQ_MIN, nyq - 100.0)
        sidx = self._drag_state["stage_index"]
        kind = self._drag_state["kind"]
        endpoint = self._drag_state["endpoint"]
        cart = self._current_cart()
        stages_data = cart["morph"][endpoint]
        if not (0 <= sidx < len(stages_data)):
            return
        key = "pole_freq_hz" if kind == "pole" else "zero_freq_hz"
        stages_data[sidx][key] = new_freq
        self.render_plot()

    def _get_morph_stages(self) -> Tuple[float, List[Stage], List[Stage]]:
        if self.endpoint_mode.get() in ("Edit M0", "Edit M100"):
            self._save_rows_to_state()
        sample_rate = self._parse_sample_rate()
        cart = self._current_cart()
        return sample_rate, dicts_to_stages(cart["morph"]["M0"]), dicts_to_stages(cart["morph"]["M100"])

    def _make_source(self, source_name: str, sample_rate: float, duration: float, midi_note: Optional[int] = None) -> np.ndarray:
        n = max(1, int(round(sample_rate * duration)))
        if source_name == SOURCE_SAW:
            return source_saw_sweep(n, sample_rate)
        if source_name == SOURCE_NOISE:
            return source_noise(n)
        if source_name == SOURCE_IMPULSE:
            return source_impulse_train(n, sample_rate)
        if source_name == SOURCE_SUB:
            return source_sub(n, sample_rate)
        if source_name == SOURCE_LOADED:
            if self.loaded_audio is None:
                raise ValueError("No WAV input loaded.")
            return fit_audio_length(self.loaded_audio, n)
        return source_midi_note(n, sample_rate, midi_note if midi_note is not None else self.last_midi_note)

    def _render_current_audio(self, source_name: Optional[str] = None, midi_note: Optional[int] = None) -> np.ndarray:
        source_name = source_name or self.source_var.get()
        sample_rate, m0, m100 = self._get_morph_stages()
        duration = self._parse_duration()
        x = self._make_source(source_name, sample_rate, duration, midi_note)
        if self.slam_enabled_var.get():
            x = apply_mackie_desk_slam(x, float(self.slam_drive_var.get()), float(self.slam_trim_var.get()))
        y = render_df2t_morph_audio(x, m0, m100, sample_rate, self.morph_mode_var.get())
        return y

    def play_preview(self) -> None:
        try:
            sample_rate = self._parse_sample_rate()
            y = self._render_current_audio()
            write_wav_mono(self.preview_path, y, sample_rate)
            self._play_wav_file(self.preview_path)
            self.status_var.set(f"playing {self.preview_path}")
        except Exception as exc:
            messagebox.showerror("Play failed", str(exc))

    def play_midi_note(self, midi_note: int) -> None:
        self.last_midi_note = midi_note
        self.source_var.set(SOURCE_MIDI)
        try:
            sample_rate = self._parse_sample_rate()
            old_duration = self.duration_var.get()
            if float(old_duration) > 2.5:
                self.duration_var.set("1.8")
            y = self._render_current_audio(SOURCE_MIDI, midi_note=midi_note)
            self.duration_var.set(old_duration)
            write_wav_mono(self.preview_path, y, sample_rate)
            self._play_wav_file(self.preview_path)
            self.status_var.set(f"MIDI {midi_note} / {note_to_freq(midi_note):.1f} Hz")
        except Exception as exc:
            messagebox.showerror("MIDI note failed", str(exc))

    def _on_key_press(self, event):
        keymap = {
            "z": 48, "s": 49, "x": 50, "d": 51, "c": 52, "v": 53,
            "g": 54, "b": 55, "h": 56, "n": 57, "j": 58, "m": 59,
            ",": 60, "l": 61, ".": 62, ";": 63, "/": 64,
        }
        ch = event.char.lower() if event.char else ""
        if ch in keymap:
            self.play_midi_note(keymap[ch])

    def _play_wav_file(self, path: Path) -> None:
        if winsound is not None:
            winsound.PlaySound(str(path), winsound.SND_FILENAME | winsound.SND_ASYNC)
        else:
            messagebox.showinfo("WAV rendered", f"Playback is automatic on Windows. WAV rendered here:\n{path}")

    def stop_audio(self) -> None:
        if winsound is not None:
            try:
                winsound.PlaySound(None, winsound.SND_PURGE)
            except Exception:
                pass
        self.status_var.set("stopped")

    def load_input_wav(self) -> None:
        try:
            sample_rate = self._parse_sample_rate()
            path = filedialog.askopenfilename(title="Load input WAV", filetypes=[("WAV files", "*.wav"), ("All files", "*.*")])
            if not path:
                return
            self.loaded_wav_path = Path(path)
            self.loaded_audio = read_wav_mono(self.loaded_wav_path, sample_rate)
            self.source_var.set(SOURCE_LOADED)
            self.loaded_label.configure(text=self.loaded_wav_path.name[:20])
            self.status_var.set(f"loaded {self.loaded_wav_path.name}")
        except Exception as exc:
            messagebox.showerror("Load WAV failed", str(exc))

    def export_current_wav(self) -> None:
        try:
            sample_rate = self._parse_sample_rate()
            y = self._render_current_audio()
            cart = self._current_cart()
            source_slug = cartridge_slug(self.source_var.get())
            morph_slug = cartridge_slug(self.morph_mode_var.get())
            slam = "_slam" if self.slam_enabled_var.get() else ""
            AUDITION_DIR.mkdir(parents=True, exist_ok=True)
            path = AUDITION_DIR / f"{cartridge_slug(cart['name'])}_{source_slug}_{morph_slug}{slam}.wav"
            write_wav_mono(path, y, sample_rate)
            self.status_var.set(f"saved {path}")
            messagebox.showinfo("Exported", f"Saved WAV:\n{path}")
        except Exception as exc:
            messagebox.showerror("Export failed", str(exc))

    def export_wav_audition(self) -> None:
        try:
            sample_rate = self._parse_sample_rate()
            cart = self._current_cart()
            sources = [SOURCE_SAW, SOURCE_NOISE, SOURCE_IMPULSE, SOURCE_SUB, SOURCE_MIDI]
            if self.loaded_audio is not None:
                sources.append(SOURCE_LOADED)
            AUDITION_DIR.mkdir(parents=True, exist_ok=True)
            paths = []
            old_source = self.source_var.get()
            for source_name in sources:
                self.source_var.set(source_name)
                y = self._render_current_audio(source_name)
                source_slug = cartridge_slug(source_name)
                morph_slug = cartridge_slug(self.morph_mode_var.get())
                slam = "_slam" if self.slam_enabled_var.get() else ""
                path = AUDITION_DIR / f"{cartridge_slug(cart['name'])}_{source_slug}_{morph_slug}{slam}.wav"
                write_wav_mono(path, y, sample_rate)
                paths.append(str(path))
            self.source_var.set(old_source)
            self.status_var.set("exported audition set")
            messagebox.showinfo("Audition exported", "Saved WAVs:\n" + "\n".join(paths))
        except Exception as exc:
            messagebox.showerror("WAV export failed", str(exc))

    def _update_slam_status(self):
        state = "on" if self.slam_enabled_var.get() else "off"
        self.status_var.set(f"Mackie desk slam {state}")

    def apply_sculpt(self) -> None:
        try:
            mode = self.endpoint_mode.get()
            if mode not in ("Edit M0", "Edit M100"):
                messagebox.showwarning("Sculpt", "Select Edit M0 or Edit M100 before applying sculpt.")
                return
            stages = self._rows_to_stages()
            sr = self._parse_sample_rate()
            nyq = sr * 0.5
            freq_mul = 2.0 ** (self.formant_shift_var.get() / 12.0)
            res_delta = self.resonance_var.get() * 0.00035
            zero_delta = self.zero_bite_var.get() * 0.0014
            c0_mul = 10.0 ** (self.c0_trim_var.get() / 20.0)
            for st in stages:
                st.pole_freq_hz = clamp(st.pole_freq_hz * freq_mul, 20.0, nyq - 100.0)
                st.zero_freq_hz = clamp(st.zero_freq_hz * freq_mul, 20.0, nyq - 100.0)
                st.pole_radius = clamp(st.pole_radius + res_delta, 0.10, 0.999)
                st.zero_radius = clamp(st.zero_radius + zero_delta, 0.0, 1.0)
                st.c0 = clamp(st.c0 * c0_mul, 0.000001, 1.0)
            self._stages_to_rows(stages)
            self._save_rows_to_state()
            self.render_plot()
            self.status_var.set("sculpt applied")
        except Exception as exc:
            messagebox.showerror("Sculpt failed", str(exc))

    def copy_endpoint(self, src: str, dst: str) -> None:
        try:
            if self.endpoint_mode.get() in ("Edit M0", "Edit M100"):
                self._save_rows_to_state()
            cart = self._current_cart()
            cart["morph"][dst] = copy.deepcopy(cart["morph"][src])
            self.endpoint_mode.set(f"Edit {dst}")
            self.last_endpoint_mode = f"Edit {dst}"
            self._load_rows_from_state()
            self.render_plot()
            self.status_var.set(f"copied {src} to {dst}")
        except Exception as exc:
            messagebox.showerror("Copy failed", str(exc))

    def swap_endpoints(self) -> None:
        try:
            if self.endpoint_mode.get() in ("Edit M0", "Edit M100"):
                self._save_rows_to_state()
            cart = self._current_cart()
            cart["morph"]["M0"], cart["morph"]["M100"] = cart["morph"]["M100"], cart["morph"]["M0"]
            self._load_rows_from_state()
            self.render_plot()
            self.status_var.set("endpoints swapped")
        except Exception as exc:
            messagebox.showerror("Swap failed", str(exc))

    def calibrate_c0(self) -> None:
        try:
            if self.endpoint_mode.get() in ("Edit M0", "Edit M100"):
                self._save_rows_to_state()
            sample_rate = self._parse_sample_rate()
            cart = self._current_cart()
            m0 = dicts_to_stages(cart["morph"]["M0"])
            m100 = dicts_to_stages(cart["morph"]["M100"])
            m50 = interpolate_stages(m0, m100, 0.5)
            peaks = [response_metrics(x, sample_rate)["peak_db"] for x in (m0, m100, m50)]
            worst_peak = max(peaks)
            target_high = 28.0
            if worst_peak <= target_high:
                messagebox.showinfo("Calibrate c0", f"No c0 reduction needed. Worst peak is {worst_peak:+.1f} dB.")
                return
            scale = 10.0 ** ((target_high - worst_peak) / (20.0 * STAGE_COUNT))
            scale = min(1.0, max(0.001, scale))
            for endpoint in ("M0", "M100"):
                stages = dicts_to_stages(cart["morph"][endpoint])
                for st in stages:
                    st.c0 = max(1e-6, min(1.0, st.c0 * scale))
                cart["morph"][endpoint] = stages_to_dicts(stages)
            self._load_rows_from_state()
            self.render_plot()
            self.status_var.set(f"c0 scaled {scale:.4f}")
        except Exception as exc:
            messagebox.showerror("Calibrate failed", str(exc))

    def save_cartridge_json(self) -> None:
        try:
            if self.endpoint_mode.get() in ("Edit M0", "Edit M100"):
                self._save_rows_to_state()
            cart = self._current_cart()
            payload = make_json_payload(cart)
            AUTHORING_DIR.mkdir(parents=True, exist_ok=True)
            default_name = cartridge_slug(cart["name"]) + ".json"
            path = filedialog.asksaveasfilename(
                title="Save cartridge JSON",
                initialdir=str(AUTHORING_DIR),
                initialfile=default_name,
                defaultextension=".json",
                filetypes=[("JSON files", "*.json"), ("All files", "*.*")],
            )
            if not path:
                return
            with open(path, "w", encoding="utf-8") as f:
                json.dump(payload, f, indent=2)
            self.status_var.set(f"saved {path}")
        except Exception as exc:
            messagebox.showerror("Save failed", str(exc))

    def load_cartridge_json(self) -> None:
        try:
            path = filedialog.askopenfilename(
                title="Load cartridge JSON",
                initialdir=str(AUTHORING_DIR),
                filetypes=[("JSON files", "*.json"), ("All files", "*.*")],
            )
            if not path:
                return
            with open(path, "r", encoding="utf-8") as f:
                payload = json.load(f)
            if "morph" not in payload or "M0" not in payload["morph"] or "M100" not in payload["morph"]:
                raise ValueError("JSON must contain morph.M0 and morph.M100 stage arrays.")
            name = payload.get("name", Path(path).stem)
            if name in self.cartridges:
                name = name + " (loaded)"
            self.cartridges[name] = {
                "name": name,
                "sample_rate": float(payload.get("sample_rate", DEFAULT_SAMPLE_RATE)),
                "formant_targets": copy.deepcopy(payload.get("formant_targets", {})),
                "zero_targets": copy.deepcopy(payload.get("zero_targets", {})),
                "morph": {
                    "M0": stages_to_dicts(dicts_to_stages(payload["morph"]["M0"])),
                    "M100": stages_to_dicts(dicts_to_stages(payload["morph"]["M100"])),
                },
            }
            self.current_cart_label.set(name)
            self.last_cart_label = name
            self._refresh_cartridge_menu_values()
            self.endpoint_mode.set("Edit M0")
            self.last_endpoint_mode = "Edit M0"
            self._load_rows_from_state()
            self.render_plot()
            self.status_var.set(f"loaded {name}")
        except Exception as exc:
            messagebox.showerror("Load failed", str(exc))

    def _refresh_cartridge_menu_values(self) -> None:
        def visit(w):
            for child in w.winfo_children():
                if isinstance(child, ttk.Combobox) and str(child.cget("textvariable")) == str(self.current_cart_label):
                    child.configure(values=list(self.cartridges.keys()))
                visit(child)
        visit(self)

    def reset_current_cartridge(self) -> None:
        label = self.current_cart_label.get()
        if label not in self.default_state:
            messagebox.showwarning("Reset", "Loaded custom cartridges do not have a built-in default to reset to.")
            return
        self.cartridges[label] = copy.deepcopy(self.default_state[label])
        self.endpoint_mode.set("Edit M0")
        self.last_endpoint_mode = "Edit M0"
        self._load_rows_from_state()
        self.render_plot()
        self.status_var.set("reset to default")


def main() -> None:
    app = PhonemeAuthorStudio()
    app.mainloop()


if __name__ == "__main__":
    main()
