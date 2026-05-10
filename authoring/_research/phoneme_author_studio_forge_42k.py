#!/usr/bin/env python3
"""
TRENCH Phoneme Author Studio
Standalone pole/zero cartridge authoring bench.

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
import wave
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import numpy as np

import tkinter as tk
from tkinter import ttk, filedialog, messagebox

# Force TkAgg early and avoid asynchronous idle draws on Windows.
# This tool does not need animated blitting; it redraws only after edits.
import matplotlib
matplotlib.use("TkAgg", force=True)

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure


DEFAULT_SAMPLE_RATE = 39062.5
STAGE_COUNT = 6
FREQ_MIN = 20.0
FREQ_MAX = 18000.0
RESPONSE_POINTS = 1024
EPS = 1e-12

AUTHORING_DIR = Path("cartridges/factory/generated/authoring")
AUDITION_DIR = Path("audition")


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
    """Starter cartridges. These are deliberately sane seeds, not doctrine."""

    def s(i, role, pf, pr, zf, zr, c0):
        return Stage(i, role, pf, pr, zf, zr, c0)

    carts = {
        "THROAT AH→EE": {
            "name": "THROAT AH→EE",
            "sample_rate": DEFAULT_SAMPLE_RATE,
            "morph": {
                "M0": [
                    s(0, "Root AH", 720, 0.982, 520, 0.64, 0.62),
                    s(1, "Chest", 1120, 0.975, 930, 0.60, 0.64),
                    s(2, "Tongue", 2450, 0.968, 1850, 0.58, 0.66),
                    s(3, "Nasal", 3350, 0.960, 3020, 0.54, 0.69),
                    s(4, "Teeth", 5200, 0.955, 4550, 0.50, 0.72),
                    s(5, "Air", 8250, 0.945, 7600, 0.48, 0.74),
                ],
                "M100": [
                    s(0, "Root EE", 310, 0.986, 260, 0.58, 0.58),
                    s(1, "Pinch", 2260, 0.990, 1780, 0.62, 0.58),
                    s(2, "Needle", 2960, 0.985, 2520, 0.64, 0.60),
                    s(3, "Glass", 3880, 0.974, 3380, 0.60, 0.63),
                    s(4, "Fric", 6100, 0.962, 5580, 0.56, 0.68),
                    s(5, "Breath", 9800, 0.948, 8900, 0.52, 0.72),
                ],
            },
        },
        "HOLLOW OH→OO": {
            "name": "HOLLOW OH→OO",
            "sample_rate": DEFAULT_SAMPLE_RATE,
            "morph": {
                "M0": [
                    s(0, "Cavity", 450, 0.956, 300, 0.78, 0.82),
                    s(1, "Bell", 760, 0.950, 620, 0.76, 0.82),
                    s(2, "Tube", 1180, 0.942, 980, 0.72, 0.84),
                    s(3, "Void", 2150, 0.936, 1850, 0.70, 0.86),
                    s(4, "Box", 3650, 0.928, 3300, 0.68, 0.88),
                    s(5, "Dust", 7200, 0.920, 6550, 0.66, 0.90),
                ],
                "M100": [
                    s(0, "Mouth", 290, 0.968, 240, 0.84, 0.76),
                    s(1, "Tunnel", 610, 0.960, 480, 0.82, 0.78),
                    s(2, "Round", 930, 0.952, 770, 0.80, 0.80),
                    s(3, "Deep", 1550, 0.944, 1280, 0.78, 0.83),
                    s(4, "Shell", 2750, 0.936, 2360, 0.76, 0.86),
                    s(5, "Cap", 6100, 0.926, 5450, 0.72, 0.89),
                ],
            },
        },
        "GRINDER UH→AE": {
            "name": "GRINDER UH→AE",
            "sample_rate": DEFAULT_SAMPLE_RATE,
            "morph": {
                "M0": [
                    s(0, "Sub Jaw", 145, 0.990, 95, 0.55, 0.15),
                    s(1, "Pressure", 310, 0.985, 220, 0.58, 0.16),
                    s(2, "Grind", 720, 0.972, 610, 0.62, 0.18),
                    s(3, "Knuckle", 1280, 0.965, 1080, 0.60, 0.20),
                    s(4, "Scrape", 3220, 0.952, 2900, 0.66, 0.23),
                    s(5, "Ash", 7900, 0.940, 7000, 0.70, 0.26),
                ],
                "M100": [
                    s(0, "Bite Low", 240, 0.982, 180, 0.60, 0.36),
                    s(1, "Throat", 690, 0.978, 520, 0.64, 0.38),
                    s(2, "Edge", 1720, 0.970, 1420, 0.68, 0.40),
                    s(3, "Rip", 2850, 0.965, 2440, 0.72, 0.42),
                    s(4, "Steel", 5150, 0.955, 4620, 0.76, 0.45),
                    s(5, "Shard", 10800, 0.944, 9700, 0.78, 0.48),
                ],
            },
        },
        "RAZOR EE→EH": {
            "name": "RAZOR EE→EH",
            "sample_rate": DEFAULT_SAMPLE_RATE,
            "morph": {
                "M0": [
                    s(0, "Needle", 360, 0.984, 300, 0.92, 0.76),
                    s(1, "Slot", 2100, 0.990, 1700, 0.94, 0.72),
                    s(2, "Cut", 2950, 0.986, 2580, 0.95, 0.72),
                    s(3, "Blade", 4200, 0.976, 3800, 0.94, 0.76),
                    s(4, "Wire", 6750, 0.962, 6100, 0.92, 0.82),
                    s(5, "Dust", 11800, 0.945, 10600, 0.90, 0.88),
                ],
                "M100": [
                    s(0, "Jaw", 540, 0.980, 430, 0.90, 0.78),
                    s(1, "Flat", 1780, 0.986, 1460, 0.93, 0.74),
                    s(2, "Hook", 2480, 0.982, 2120, 0.94, 0.74),
                    s(3, "Saw", 3580, 0.970, 3180, 0.94, 0.78),
                    s(4, "Fang", 5600, 0.958, 5050, 0.92, 0.83),
                    s(5, "Bright", 9900, 0.942, 9000, 0.90, 0.88),
                ],
            },
        },
    }

    # Convert Stage objects to plain dictionaries for JSON-safe app state.
    out = {}
    for label, cart in carts.items():
        out[label] = {
            "name": cart["name"],
            "sample_rate": cart["sample_rate"],
            "morph": {
                "M0": [asdict(x) for x in cart["morph"]["M0"]],
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

    # Critical: c1/c2 scale with c0 so attenuation does not move zero locations.
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
        "morph": {
            "M0": m0,
            "M100": m100,
        },
        "corners": {
            "M0_Q0": copy.deepcopy(m0),
            "M0_Q100": copy.deepcopy(m0),
            "M100_Q0": copy.deepcopy(m100),
            "M100_Q100": copy.deepcopy(m100),
        },
    }


def write_wav_mono(path: Path, audio: np.ndarray, sample_rate: float) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    audio = np.asarray(audio, dtype=np.float64)
    audio = np.nan_to_num(audio, nan=0.0, posinf=0.0, neginf=0.0)
    peak = float(np.max(np.abs(audio))) if audio.size else 0.0
    if peak > 0.98:
        audio = audio * (0.98 / peak)
    audio_i16 = np.clip(audio, -1.0, 1.0)
    audio_i16 = (audio_i16 * 32767.0).astype(np.int16)

    # RIFF/WAV stores sample rate as an integer. DSP still uses DEFAULT_SAMPLE_RATE.
    wav_rate = int(round(sample_rate))
    with wave.open(str(path), "wb") as wf:
        wf.setnchannels(1)
        wf.setsampwidth(2)
        wf.setframerate(wav_rate)
        wf.writeframes(audio_i16.tobytes())


def source_saw_sweep(n: int, sample_rate: float) -> np.ndarray:
    t = np.arange(n) / sample_rate
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
    starts = range(0, n, hit_len)
    for start in starts:
        end = min(n, start + hit_len)
        m = end - start
        if m <= 0:
            continue
        tt = np.arange(m) / sample_rate
        env = np.exp(-5.0 * tt)
        freq = 72.0 - 25.0 * (1.0 - np.exp(-18.0 * tt))
        phase = 2.0 * np.pi * np.cumsum(freq / sample_rate)
        x[start:end] += 0.42 * env * np.sin(phase)
    return x


def render_df2t_morph_audio(
    source: np.ndarray,
    m0_stages: List[Stage],
    m100_stages: List[Stage],
    sample_rate: float,
) -> np.ndarray:
    """Offline audition path using the exact DF2T runtime equation."""
    n = int(source.size)
    y_out = np.zeros(n, dtype=np.float64)
    state_w1 = np.zeros(STAGE_COUNT, dtype=np.float64)
    state_w2 = np.zeros(STAGE_COUNT, dtype=np.float64)

    # Coefficient interpolation is done in pole/zero authoring space, then compiled.
    # This is intentionally simple for the authoring bench.
    for i in range(n):
        morph_t = i / max(1, n - 1)
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
            # Keep rendering instead of creating a corrupt WAV.
            y_out[i] = 0.0
            state_w1[:] = 0.0
            state_w2[:] = 0.0
    return y_out


class PhonemeAuthorStudio(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("TRENCH Phoneme Author Studio")
        self.geometry("1320x850")
        self.minsize(1180, 720)

        self.bg = "#101316"
        self.panel = "#151a1f"
        self.panel2 = "#1b2229"
        self.text = "#d6dde3"
        self.subtext = "#8f9aa6"
        self.line = "#2b343d"
        self.accent = "#9fd0ff"
        self.warn = "#ffc36b"

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

        self.stage_vars: List[Dict[str, tk.StringVar]] = []
        self.hover_data: Optional[dict] = None

        self._build_ui()
        self._load_rows_from_state()
        # Defer the first render until Tk has fully created the window image.
        self.after(50, self.render_plot)

    def _configure_style(self) -> None:
        try:
            self.style.theme_use("clam")
        except tk.TclError:
            pass
        self.style.configure("TFrame", background=self.bg)
        self.style.configure("Panel.TFrame", background=self.panel)
        self.style.configure("TLabel", background=self.bg, foreground=self.text, font=("TkDefaultFont", 10))
        self.style.configure("Muted.TLabel", background=self.bg, foreground=self.subtext, font=("TkDefaultFont", 9))
        self.style.configure("Panel.TLabel", background=self.panel, foreground=self.text)
        self.style.configure("Header.TLabel", background=self.bg, foreground=self.text, font=("TkDefaultFont", 16, "bold"))
        self.style.configure("SmallHeader.TLabel", background=self.panel, foreground=self.text, font=("TkDefaultFont", 10, "bold"))
        self.style.configure("TButton", background=self.panel2, foreground=self.text, padding=(9, 6), relief="flat")
        self.style.map("TButton", background=[("active", "#26313a")])
        self.style.configure("TCombobox", fieldbackground=self.panel2, background=self.panel2, foreground=self.text)
        self.style.configure("TRadiobutton", background=self.panel, foreground=self.text)
        self.style.configure("TEntry", fieldbackground="#0d1013", foreground=self.text, insertcolor=self.text)

    def _build_ui(self) -> None:
        header = ttk.Frame(self, style="TFrame")
        header.pack(fill="x", padx=14, pady=(12, 8))
        ttk.Label(header, text="TRENCH / PHONEME AUTHOR STUDIO", style="Header.TLabel").pack(side="left")
        ttk.Label(
            header,
            text="standalone pole-zero bench · DF2T cascade · no plugin writes",
            style="Muted.TLabel",
        ).pack(side="left", padx=(14, 0), pady=(4, 0))

        root = ttk.Frame(self, style="TFrame")
        root.pack(fill="both", expand=True, padx=14, pady=(0, 14))

        left = ttk.Frame(root, style="Panel.TFrame")
        left.pack(side="left", fill="y", padx=(0, 10))
        left.configure(width=500)
        left.pack_propagate(False)

        right = ttk.Frame(root, style="Panel.TFrame")
        right.pack(side="left", fill="both", expand=True)

        self._build_controls(left)
        self._build_plot(right)

    def _safe_canvas_draw(self) -> None:
        """Draw the Matplotlib canvas without scheduling Tk idle blits."""
        if not getattr(self, "_alive", False):
            return
        try:
            self.canvas.draw()
            self.canvas.flush_events()
        except tk.TclError:
            # Stale Tk image/canvas during shutdown or rapid redraw.
            return

    def _on_close(self) -> None:
        self._alive = False
        try:
            if hasattr(self, "canvas") and hasattr(self, "_hover_cid"):
                self.canvas.mpl_disconnect(self._hover_cid)
            if hasattr(self, "canvas"):
                self.canvas.get_tk_widget().destroy()
        except Exception:
            pass
        self.quit()
        self.destroy()

    def _build_controls(self, parent: ttk.Frame) -> None:
        top = ttk.Frame(parent, style="Panel.TFrame")
        top.pack(fill="x", padx=12, pady=12)

        ttk.Label(top, text="Cartridge", style="Panel.TLabel").grid(row=0, column=0, sticky="w")
        cart_menu = ttk.Combobox(
            top,
            textvariable=self.current_cart_label,
            values=list(self.cartridges.keys()),
            state="readonly",
            width=28,
        )
        cart_menu.grid(row=1, column=0, sticky="ew", pady=(4, 10))
        cart_menu.bind("<<ComboboxSelected>>", self._on_cartridge_changed)

        ttk.Label(top, text="Sample rate", style="Panel.TLabel").grid(row=0, column=1, sticky="w", padx=(10, 0))
        sr = ttk.Entry(top, textvariable=self.sample_rate_var, width=12)
        sr.grid(row=1, column=1, sticky="ew", padx=(10, 0), pady=(4, 10))

        mode_frame = ttk.Frame(top, style="Panel.TFrame")
        mode_frame.grid(row=2, column=0, columnspan=2, sticky="ew")
        for i, label in enumerate(("Edit M0", "Edit M100", "View both")):
            rb = ttk.Radiobutton(
                mode_frame,
                text=label,
                variable=self.endpoint_mode,
                value=label,
                command=self._on_endpoint_mode_changed,
            )
            rb.grid(row=0, column=i, sticky="w", padx=(0, 12))
        top.columnconfigure(0, weight=1)

        table_outer = ttk.Frame(parent, style="Panel.TFrame")
        table_outer.pack(fill="x", padx=12, pady=(4, 8))
        ttk.Label(table_outer, text="Six stage endpoint editor", style="SmallHeader.TLabel").pack(anchor="w", pady=(0, 8))

        table = tk.Frame(table_outer, bg=self.panel)
        table.pack(fill="x")
        headings = ["S", "role", "pole Hz", "pole r", "zero Hz", "zero r", "c0"]
        widths = [3, 12, 9, 7, 9, 7, 6]
        for col, (heading, width) in enumerate(zip(headings, widths)):
            lbl = tk.Label(
                table,
                text=heading,
                bg=self.panel,
                fg=self.subtext,
                font=("TkDefaultFont", 8, "bold"),
                anchor="w",
            )
            lbl.grid(row=0, column=col, sticky="ew", padx=2, pady=(0, 4))
            table.columnconfigure(col, weight=1 if col == 1 else 0)

        for row in range(STAGE_COUNT):
            row_vars = {
                "stage_index": tk.StringVar(value=str(row)),
                "role_label": tk.StringVar(),
                "pole_freq_hz": tk.StringVar(),
                "pole_radius": tk.StringVar(),
                "zero_freq_hz": tk.StringVar(),
                "zero_radius": tk.StringVar(),
                "c0": tk.StringVar(),
            }
            self.stage_vars.append(row_vars)
            keys = ["stage_index", "role_label", "pole_freq_hz", "pole_radius", "zero_freq_hz", "zero_radius", "c0"]
            for col, key in enumerate(keys):
                ent = tk.Entry(
                    table,
                    textvariable=row_vars[key],
                    width=widths[col],
                    bg="#0d1013",
                    fg=self.text,
                    insertbackground=self.text,
                    relief="flat",
                    highlightthickness=1,
                    highlightbackground=self.line,
                    highlightcolor=self.accent,
                    font=("TkDefaultFont", 9),
                )
                if key == "stage_index":
                    ent.configure(state="readonly", readonlybackground="#11161b", fg=self.subtext)
                ent.grid(row=row + 1, column=col, sticky="ew", padx=2, pady=2, ipady=3)

        btns = ttk.Frame(parent, style="Panel.TFrame")
        btns.pack(fill="x", padx=12, pady=8)
        button_specs = [
            ("Render plot", self.render_plot),
            ("Calibrate c0", self.calibrate_c0),
            ("Save cartridge JSON", self.save_cartridge_json),
            ("Load cartridge JSON", self.load_cartridge_json),
            ("Export WAV audition", self.export_wav_audition),
            ("Reset to default cartridge", self.reset_current_cartridge),
        ]
        for i, (label, cmd) in enumerate(button_specs):
            b = ttk.Button(btns, text=label, command=cmd)
            b.grid(row=i // 2, column=i % 2, sticky="ew", padx=3, pady=3)
        btns.columnconfigure(0, weight=1)
        btns.columnconfigure(1, weight=1)

        readout_frame = ttk.Frame(parent, style="Panel.TFrame")
        readout_frame.pack(fill="both", expand=True, padx=12, pady=(8, 12))
        ttk.Label(readout_frame, text="Hover readout / warnings", style="SmallHeader.TLabel").pack(anchor="w", pady=(0, 6))
        self.info_text = tk.Text(
            readout_frame,
            height=16,
            wrap="word",
            bg="#0d1013",
            fg=self.text,
            insertbackground=self.text,
            relief="flat",
            highlightthickness=1,
            highlightbackground=self.line,
            font=("Menlo", 10) if os.name != "nt" else ("Consolas", 10),
        )
        self.info_text.pack(fill="both", expand=True)
        self.info_text.configure(state="disabled")

    def _build_plot(self, parent: ttk.Frame) -> None:
        self.fig = Figure(figsize=(7.7, 5.3), dpi=100, facecolor=self.panel)
        self.ax = self.fig.add_subplot(111)
        self.fig.subplots_adjust(left=0.08, right=0.98, top=0.94, bottom=0.12)
        self.canvas = FigureCanvasTkAgg(self.fig, master=parent)
        widget = self.canvas.get_tk_widget()
        widget.configure(bg=self.panel, highlightthickness=0)
        widget.pack(fill="both", expand=True, padx=8, pady=(8, 0))

        footer = ttk.Frame(parent, style="Panel.TFrame")
        footer.pack(fill="x", padx=8, pady=(0, 8))
        ttk.Label(footer, text="hover graph for stage readout · render after edits", style="Panel.TLabel").pack(side="left")

        self._hover_cid = self.canvas.mpl_connect("motion_notify_event", self._on_plot_hover)

    def _parse_sample_rate(self) -> float:
        try:
            sr = float(self.sample_rate_var.get())
            if sr <= 1000:
                raise ValueError
            return sr
        except ValueError:
            raise ValueError("Sample rate must be a number above 1000 Hz.")

    def _rows_to_stages(self) -> List[Stage]:
        stages = []
        errors = []
        for i, row in enumerate(self.stage_vars):
            try:
                st = Stage(
                    stage_index=i,
                    role_label=row["role_label"].get().strip() or f"Stage {i}",
                    pole_freq_hz=float(row["pole_freq_hz"].get()),
                    pole_radius=float(row["pole_radius"].get()),
                    zero_freq_hz=float(row["zero_freq_hz"].get()),
                    zero_radius=float(row["zero_radius"].get()),
                    c0=float(row["c0"].get()),
                )
                stages.append(st)
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
        stages = self._rows_to_stages()
        cart = self._current_cart(cart_label)
        cart["sample_rate"] = self._parse_sample_rate()
        cart["morph"][endpoint] = stages_to_dicts(stages)

    def _load_rows_from_state(self) -> None:
        cart = self._current_cart()
        self.sample_rate_var.set(str(cart.get("sample_rate", DEFAULT_SAMPLE_RATE)))
        mode = self.endpoint_mode.get()
        endpoint = "M0" if mode != "Edit M100" else "M100"
        stages = dicts_to_stages(cart["morph"][endpoint])
        self._stages_to_rows(stages)
        self._set_row_editing(mode != "View both")

    def _set_row_editing(self, enabled: bool) -> None:
        # Preserve stage index as readonly. Other entries become disabled in view mode.
        table_widgets = []
        for child in self.winfo_children():
            table_widgets.extend(child.winfo_children())
        for row in self.stage_vars:
            # Entry state is handled by searching all widgets with matching textvariables.
            pass
        # Simpler: inspect descendants for Entry widgets. Keep index entries readonly.
        def visit(w):
            for child in w.winfo_children():
                if isinstance(child, tk.Entry):
                    tv = str(child.cget("textvariable"))
                    is_index = any(str(row["stage_index"]) == tv for row in self.stage_vars)
                    if is_index:
                        child.configure(state="readonly")
                    else:
                        child.configure(state="normal" if enabled else "disabled")
                visit(child)
        visit(self)

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
            f"Sample rate: {sample_rate:.3f} Hz",
            "",
            "Peaks:",
            f"  M0   peak {m0m['peak_db']:+.1f} dB | range {m0m['range_db']:.1f} dB",
            f"  M50  peak {m50m['peak_db']:+.1f} dB | range {m50m['range_db']:.1f} dB",
            f"  M100 peak {m100m['peak_db']:+.1f} dB | range {m100m['range_db']:.1f} dB",
            "",
            "Warnings:",
        ]
        if warnings:
            lines.extend([f"  - {w}" for w in warnings])
        else:
            lines.append("  none")
        lines.append("")
        lines.append("Hover the plot for per-stage contribution.")
        return "\n".join(lines)

    def render_plot(self) -> None:
        try:
            if self.endpoint_mode.get() in ("Edit M0", "Edit M100"):
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
        self.ax.set_facecolor("#0d1013")
        self.ax.set_xscale("log")
        self.ax.set_xlim(FREQ_MIN, FREQ_MAX)
        self.ax.set_ylim(-72, 42)
        self.ax.set_xlabel("frequency Hz", color=self.subtext)
        self.ax.set_ylabel("dB", color=self.subtext)
        self.ax.tick_params(axis="both", colors=self.subtext, labelsize=9)
        for spine in self.ax.spines.values():
            spine.set_color(self.line)
        self.ax.grid(True, which="major", color="#27313a", alpha=0.7, linewidth=0.8)
        self.ax.grid(True, which="minor", color="#1d252c", alpha=0.5, linewidth=0.5)
        self.ax.axhline(0.0, color="#46515c", linewidth=0.8)

        mode = self.endpoint_mode.get()
        title = cart["name"]
        self.ax.set_title(title, color=self.text, fontsize=12, fontweight="bold")

        line_sets = []
        if mode == "Edit M0":
            line_sets.append(("M0", m0, "#9fd0ff", "-", 2.5))
        elif mode == "Edit M100":
            line_sets.append(("M100", m100, "#ffb86b", "-", 2.5))
        else:
            m50 = interpolate_stages(m0, m100, 0.5)
            line_sets.extend([
                ("M0", m0, "#9fd0ff", "-", 2.3),
                ("M50", m50, "#d8ff9f", "-.", 2.0),
                ("M100", m100, "#ffb86b", "--", 2.3),
            ])

        hover_choice = None
        for label, stages, color, linestyle, width in line_sets:
            cascade_h, stage_h = cascade_response(stages, freqs, sample_rate)
            cdb = np.clip(db(cascade_h), -120, 120)
            sdb = np.clip(db(stage_h), -120, 120)
            for si in range(stage_h.shape[0]):
                self.ax.plot(freqs, sdb[si], color=color, linestyle=":", linewidth=0.7, alpha=0.22)
            self.ax.plot(freqs, cdb, color=color, linestyle=linestyle, linewidth=width, label=label)
            if label == "M50" or hover_choice is None:
                hover_choice = {
                    "label": label,
                    "freqs": freqs,
                    "cascade_db": cdb,
                    "stage_db": sdb,
                    "stages": stages,
                }
        self.ax.legend(loc="upper right", frameon=False, labelcolor=self.text)
        self.hover_data = hover_choice
        self._safe_canvas_draw()
        self._set_info_text(self._collect_validation_text(sample_rate, m0, m100))

    def _on_plot_hover(self, event) -> None:
        if event.inaxes != self.ax or self.hover_data is None or event.xdata is None:
            return
        freqs = self.hover_data["freqs"]
        idx = int(np.argmin(np.abs(np.log(freqs) - math.log(max(FREQ_MIN, event.xdata)))))
        freq = freqs[idx]
        cdb = self.hover_data["cascade_db"][idx]
        sdb = self.hover_data["stage_db"][:, idx]
        stages = self.hover_data["stages"]
        label = self.hover_data["label"]

        lines = [
            f"Hover readout: {label}",
            f"freq:    {freq:8.1f} Hz",
            f"cascade: {cdb:+8.2f} dB",
            "",
            "stage contributions:",
        ]
        for st, val in zip(stages, sdb):
            lines.append(f"  S{st.stage_index} {st.role_label[:10]:10s} {val:+8.2f} dB")
        lines.append("")
        lines.append("Render plot to refresh validation warnings.")
        self._set_info_text("\n".join(lines))

    def calibrate_c0(self) -> None:
        try:
            if self.endpoint_mode.get() in ("Edit M0", "Edit M100"):
                self._save_rows_to_state()
            sample_rate = self._parse_sample_rate()
            cart = self._current_cart()
            m0 = dicts_to_stages(cart["morph"]["M0"])
            m100 = dicts_to_stages(cart["morph"]["M100"])
            m50 = interpolate_stages(m0, m100, 0.5)
            peaks = [
                response_metrics(m0, sample_rate)["peak_db"],
                response_metrics(m100, sample_rate)["peak_db"],
                response_metrics(m50, sample_rate)["peak_db"],
            ]
            worst_peak = max(peaks)
            target_high = 28.0
            if worst_peak <= target_high:
                messagebox.showinfo(
                    "Calibrate c0",
                    f"No c0 reduction needed. Worst peak is {worst_peak:+.1f} dB.",
                )
                return
            # Scaling all six stage c0 values by s changes cascade level by 20*6*log10(s).
            scale = 10.0 ** ((target_high - worst_peak) / (20.0 * STAGE_COUNT))
            scale = min(1.0, max(0.001, scale))
            for endpoint in ("M0", "M100"):
                stages = dicts_to_stages(cart["morph"][endpoint])
                for st in stages:
                    st.c0 = max(1e-6, min(1.0, st.c0 * scale))
                cart["morph"][endpoint] = stages_to_dicts(stages)
            self._load_rows_from_state()
            self.render_plot()
            messagebox.showinfo(
                "Calibrate c0",
                f"Scaled all c0 values by {scale:.4f}. Worst peak was {worst_peak:+.1f} dB; target high is +28 dB.",
            )
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
            messagebox.showinfo("Saved", f"Saved cartridge JSON:\n{path}")
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
            messagebox.showinfo("Loaded", f"Loaded cartridge:\n{name}")
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

    def export_wav_audition(self) -> None:
        try:
            if self.endpoint_mode.get() in ("Edit M0", "Edit M100"):
                self._save_rows_to_state()
            sample_rate = self._parse_sample_rate()
            cart = self._current_cart()
            m0 = dicts_to_stages(cart["morph"]["M0"])
            m100 = dicts_to_stages(cart["morph"]["M100"])
            n = int(round(sample_rate * 4.0))
            slug = cartridge_slug(cart["name"])
            sources = {
                "saw": source_saw_sweep(n, sample_rate),
                "noise": source_noise(n),
                "impulse": source_impulse_train(n, sample_rate),
                "sub": source_sub(n, sample_rate),
            }
            AUDITION_DIR.mkdir(parents=True, exist_ok=True)
            paths = []
            for name, x in sources.items():
                y = render_df2t_morph_audio(x, m0, m100, sample_rate)
                path = AUDITION_DIR / f"{slug}_{name}.wav"
                write_wav_mono(path, y, sample_rate)
                paths.append(str(path))
            messagebox.showinfo("Audition exported", "Saved WAVs:\n" + "\n".join(paths))
        except Exception as exc:
            messagebox.showerror("WAV export failed", str(exc))


def main() -> None:
    app = PhonemeAuthorStudio()
    app.mainloop()


if __name__ == "__main__":
    main()
