"""forge/sift.py — TRENCH cartridge sculptor.

Live Z-plane authoring with direct manipulation. Replaces the previous static
sift bench; merges the studio's drag handles with sift's real-time audio path.

  - Drag pole (X) and zero (O) handles on the plot:
      horizontal motion = frequency
      vertical motion   = radius (sensitivity tuned to heritage range 0.95..0.999)
  - Pinknoise loops through the cascade at the current MORPH position
  - One uniform c0 spinbox (cascade gain discipline — per-cartridge attenuation)
  - KEEP / REJECT / AGAIN

Markers are color-coded: blue = M0 endpoint, orange = M100 endpoint.
Both endpoints' markers are visible at once so you can see the full morph map.

Usage:
    python forge/sift.py <cartridge.json>
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import threading
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional

import numpy as np

import tkinter as tk

import matplotlib
matplotlib.use("TkAgg", force=True)
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

try:
    import sounddevice as sd
except ImportError:
    sd = None
try:
    from scipy.signal import sosfilt
except ImportError:
    sosfilt = None


REPO = Path(__file__).resolve().parent.parent
KEPT_DIR = REPO / "cartridges" / "factory" / "kept"
# Plugin hot-reload slot — TrenchAuthoringSlot watches this path and reloads
# the live VST3 in your DAW within ~500 ms of every save.
PLUGIN_SLOT = Path.home() / "Documents" / "TRENCH" / "authoring_slot.json"
AUTHORING_SR = 39062.5
PLAYBACK_SR = 44100
BLOCK = 512
N_STAGES = 6

# Drag tuning
RADIUS_SENS = 0.0025      # radius units per dB of vertical drag
POLE_Y = 36.0             # marker visual y for poles (above response)
ZERO_Y = -32.0            # marker visual y for zeros (below response)
PICK_THRESH = 0.18        # nearest-handle log-distance threshold for pick

# Heritage-bounded ranges
POLE_R_MIN, POLE_R_MAX = 0.5, 0.99999
ZERO_R_MIN, ZERO_R_MAX = 0.05, 0.99999
FREQ_MIN_AUTH = 20.0
FREQ_MAX_AUTH = AUTHORING_SR * 0.49


@dataclass
class Stage:
    pole_freq_hz: float
    pole_r: float
    zero_freq_hz: float
    zero_r: float
    flag: int = 1


# ─── Cartridge format conversion ─────────────────────────────────────────────

def rawstage_to_authoring(rs: dict, sr: float = AUTHORING_SR) -> Stage:
    """RawStage {a1, r, val1, val2, val3, flag} → authoring Stage (heritage form)."""
    if rs.get("flag", 1) == 0:
        return Stage(0.0, 0.0, 0.0, 0.0, flag=0)
    r = float(rs["r"])
    a1 = float(rs["a1"])
    if r > 0:
        cos_t = max(-1.0, min(1.0, -a1 / (2.0 * r)))
        pole_freq = sr / (2.0 * math.pi) * math.acos(cos_t)
    else:
        pole_freq = 0.0
    b1 = a1 + float(rs["val2"])
    b2 = r * r - float(rs["val3"])
    zero_r = math.sqrt(max(0.0, b2))
    if zero_r > 0:
        cos_tz = max(-1.0, min(1.0, -b1 / (2.0 * zero_r)))
        zero_freq = sr / (2.0 * math.pi) * math.acos(cos_tz)
    else:
        zero_freq = 0.0
    return Stage(pole_freq, r, zero_freq, zero_r, flag=1)


def authoring_to_rawstage(st: Stage, c0: float, sr: float = AUTHORING_SR) -> dict:
    if st.flag == 0:
        return {"a1": 0.0, "r": 0.0, "val1": 0.0, "val2": 0.0, "val3": 0.0, "flag": 0}
    pole_r = max(0.0, min(0.999999, st.pole_r))
    zero_freq = min(st.zero_freq_hz, sr * 0.49)
    a1 = -2.0 * pole_r * math.cos(2.0 * math.pi * st.pole_freq_hz / sr)
    a2 = pole_r * pole_r
    b1 = -2.0 * st.zero_r * math.cos(2.0 * math.pi * zero_freq / sr)
    b2 = st.zero_r * st.zero_r
    return {
        "a1": a1, "r": pole_r,
        "val1": c0 - 1.0,
        "val2": b1 - a1,
        "val3": a2 - b2,
        "flag": 1,
    }


def stage_to_coefs(st: Stage, c0: float, sr: float) -> dict:
    """Authoring Stage → DF2T coefs (heritage convention: c1=b1, c2=b2)."""
    if st.flag == 0:
        return {"c0": 1.0, "c1": 0.0, "c2": 0.0, "c3": 0.0, "c4": 0.0}
    a1 = -2.0 * st.pole_r * math.cos(2.0 * math.pi * st.pole_freq_hz / sr)
    b1 = -2.0 * st.zero_r * math.cos(2.0 * math.pi * st.zero_freq_hz / sr)
    return {
        "c0": c0,
        "c1": b1,
        "c2": st.zero_r * st.zero_r,
        "c3": a1,
        "c4": st.pole_r * st.pole_r,
    }


def lerp_stage(a: Stage, b: Stage, t: float) -> Stage:
    return Stage(
        pole_freq_hz=(1 - t) * a.pole_freq_hz + t * b.pole_freq_hz,
        pole_r=(1 - t) * a.pole_r + t * b.pole_r,
        zero_freq_hz=(1 - t) * a.zero_freq_hz + t * b.zero_freq_hz,
        zero_r=(1 - t) * a.zero_r + t * b.zero_r,
        flag=1 if (a.flag or b.flag) else 0,
    )


def cascade_db(stages: List[Stage], c0: float, sr: float, freqs: np.ndarray) -> np.ndarray:
    w = 2.0 * np.pi * freqs / sr
    e1 = np.exp(-1j * w)
    e2 = np.exp(-2j * w)
    mag = np.ones_like(freqs)
    for st in stages:
        c = stage_to_coefs(st, c0, sr)
        if c["c3"] == 0.0 and c["c4"] == 0.0:
            continue
        mag *= np.abs((c["c0"] + c["c1"] * e1 + c["c2"] * e2) /
                      (1.0 + c["c3"] * e1 + c["c4"] * e2))
    return 20.0 * np.log10(np.maximum(mag, 1e-12))


def synth_pinknoise(n: int, seed: int = 42) -> np.ndarray:
    rng = np.random.default_rng(seed)
    rows = 16
    cols = np.zeros((rows, n), dtype=np.float64)
    cols[0] = rng.standard_normal(n)
    for k in range(1, rows):
        period = 1 << k
        idx = np.arange(n) // period
        unique = rng.standard_normal(idx.max() + 1)
        cols[k] = unique[idx]
    pink = cols.sum(axis=0)
    rms = float(np.sqrt(np.mean(pink ** 2)))
    if rms > 0:
        pink *= 0.18 / rms
    return pink


# ─── Sculptor UI ─────────────────────────────────────────────────────────────

class Sculptor(tk.Tk):
    def __init__(self, cart_path: Path):
        super().__init__()
        self.title(f"TRENCH SCULPT :: {cart_path.stem}")
        self.geometry("1280x800")
        self.configure(bg="#0a0c0e")

        self.cart_path = cart_path
        cart = json.loads(cart_path.read_text(encoding="utf-8"))
        self.cart_name = cart.get("name", cart_path.stem)
        m0, m100, c0 = self._load_cart(cart)
        while len(m0) < N_STAGES: m0.append(Stage(0, 0, 0, 0, flag=0))
        while len(m100) < N_STAGES: m100.append(Stage(0, 0, 0, 0, flag=0))
        self.m0 = m0[:N_STAGES]
        self.m100 = m100[:N_STAGES]
        self.c0 = c0

        self.morph = 0.0
        self.lock = threading.Lock()
        self.playing = False
        self.paused = False
        self.stream = None
        self.state = None
        self.source = synth_pinknoise(int(PLAYBACK_SR * 6.0))
        self.source_pos = 0
        self._drag: Optional[dict] = None
        self._handle_index: List[dict] = []

        self._build_ui()
        self._render_plot()
        self.protocol("WM_DELETE_WINDOW", self._on_close)

    def _load_cart(self, cart: dict):
        if "keyframes" in cart:
            kf = {k.get("label"): k for k in cart["keyframes"]}
            need = ["M0_Q0", "M100_Q0"]
            for n in need:
                if n not in kf:
                    raise ValueError(f"cartridge missing keyframe: {n}")
            m0 = [rawstage_to_authoring(rs) for rs in kf["M0_Q0"]["stages"][:N_STAGES]]
            m100 = [rawstage_to_authoring(rs) for rs in kf["M100_Q0"]["stages"][:N_STAGES]]
            c0_uniform = float(cart.get("provenance", {}).get("c0_uniform", 0.5619))
            return m0, m100, c0_uniform
        if "morph" in cart:
            def conv(d):
                return Stage(
                    pole_freq_hz=float(d["pole_freq_hz"]),
                    pole_r=float(d["pole_radius"]),
                    zero_freq_hz=float(d["zero_freq_hz"]),
                    zero_r=float(d["zero_radius"]),
                )
            m0 = [conv(s) for s in cart["morph"]["M0"][:N_STAGES]]
            m100 = [conv(s) for s in cart["morph"]["M100"][:N_STAGES]]
            c0_first = float(cart["morph"]["M0"][0].get("c0", 0.5619))
            return m0, m100, c0_first
        raise ValueError("unknown cartridge format (no 'keyframes' or 'morph' key)")

    # ─── UI ──────────────────────────────────────────────────────────────

    def _build_ui(self):
        top = tk.Frame(self, bg="#0a0c0e")
        top.pack(fill="x", padx=14, pady=(12, 4))
        tk.Label(top, text=self.cart_name, fg="#ffba00", bg="#0a0c0e",
                 font=("Segoe UI", 14, "bold")).pack(side="left")
        tk.Label(top,
                 text="  drag X (poles, top) / O (zeros, bottom)  -  horizontal=freq, vertical=radius",
                 fg="#8b96a1", bg="#0a0c0e", font=("Segoe UI", 9)).pack(side="left", padx=10)

        self.fig = Figure(figsize=(10, 4.6), dpi=100, facecolor="#11161b")
        self.ax = self.fig.add_subplot(111)
        self.fig.subplots_adjust(left=0.06, right=0.985, top=0.92, bottom=0.10)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self)
        self.canvas.get_tk_widget().pack(fill="both", expand=True, padx=14, pady=4)
        self.canvas.mpl_connect("button_press_event", self._on_press)
        self.canvas.mpl_connect("motion_notify_event", self._on_motion)
        self.canvas.mpl_connect("button_release_event", self._on_release)

        morph_row = tk.Frame(self, bg="#0a0c0e")
        morph_row.pack(fill="x", padx=14, pady=(8, 4))
        tk.Label(morph_row, text="MORPH", fg="#e6edf2", bg="#0a0c0e",
                 font=("Segoe UI", 11, "bold"), width=8, anchor="w").pack(side="left")
        self.morph_var = tk.DoubleVar(value=0.0)
        self.morph_label = tk.Label(morph_row, text="0.000", fg="#d8ff9a",
                                    bg="#0a0c0e", font=("Consolas", 11), width=8)
        self.morph_label.pack(side="right")
        tk.Scale(morph_row, from_=0.0, to=1.0, resolution=0.001,
                 orient="horizontal", command=self._on_morph,
                 variable=self.morph_var,
                 bg="#0a0c0e", fg="#e6edf2",
                 troughcolor="#1b2229", highlightthickness=0,
                 showvalue=0, length=600, sliderrelief="flat").pack(
            side="left", fill="x", expand=True, padx=10)

        ctrl_row = tk.Frame(self, bg="#0a0c0e")
        ctrl_row.pack(fill="x", padx=14, pady=(4, 4))
        tk.Label(ctrl_row, text="c0", fg="#e6edf2", bg="#0a0c0e",
                 font=("Segoe UI", 10, "bold"), width=4).pack(side="left")
        self.c0_var = tk.StringVar(value=f"{self.c0:.4f}")
        c0_entry = tk.Entry(ctrl_row, textvariable=self.c0_var, width=8,
                            bg="#1b2229", fg="#e6edf2", insertbackground="#e6edf2",
                            relief="flat", font=("Consolas", 10),
                            highlightthickness=1, highlightbackground="#26313a")
        c0_entry.pack(side="left", padx=(0, 12), ipady=3)
        c0_entry.bind("<Return>", self._on_c0_change)
        c0_entry.bind("<FocusOut>", self._on_c0_change)

        self.play_btn = tk.Button(ctrl_row, text="PLAY", fg="white", bg="#22aa44",
                                  font=("Segoe UI", 10, "bold"), width=8,
                                  command=self._toggle_play, relief="flat")
        self.play_btn.pack(side="left", padx=2)

        self.status = tk.Label(ctrl_row, text="silent — press PLAY",
                               fg="#8b96a1", bg="#0a0c0e", font=("Segoe UI", 10))
        self.status.pack(side="left", padx=12)

        sift_row = tk.Frame(self, bg="#0a0c0e")
        sift_row.pack(fill="x", padx=14, pady=(8, 12))
        for label, color, cmd in (
            ("KEEP",   "#22aa44", self._on_keep),
            ("REJECT", "#cc3333", self._on_reject),
            ("AGAIN",  "#3366cc", self._on_again),
        ):
            tk.Button(sift_row, text=label, fg="white", bg=color,
                      font=("Segoe UI", 12, "bold"), height=2,
                      command=cmd, relief="flat").pack(
                side="left", padx=4, expand=True, fill="x")

    # ─── Plot ────────────────────────────────────────────────────────────

    def _render_plot(self):
        freqs = np.geomspace(20.0, 18000.0, 1024)
        with self.lock:
            t = self.morph
            c0 = self.c0
        m0_db = cascade_db(self.m0, c0, AUTHORING_SR, freqs)
        m100_db = cascade_db(self.m100, c0, AUTHORING_SR, freqs)
        cur_stages = [lerp_stage(a, b, t) for a, b in zip(self.m0, self.m100)]
        cur_db = cascade_db(cur_stages, c0, AUTHORING_SR, freqs)

        self.ax.clear()
        self.ax.set_facecolor("#0d1013")
        self.ax.set_xscale("log")
        self.ax.set_xlim(20, 18000)
        self.ax.set_ylim(-60, 45)
        self.ax.tick_params(colors="#8b96a1", labelsize=8)
        for spine in self.ax.spines.values():
            spine.set_color("#26313a")
        self.ax.grid(True, which="major", alpha=0.18, color="#27313a")
        self.ax.axhline(0.0, color="#46515c", lw=0.5)
        self.ax.axhline(POLE_Y, color="#26313a", lw=0.4, ls="--", alpha=0.5)
        self.ax.axhline(ZERO_Y, color="#26313a", lw=0.4, ls="--", alpha=0.5)
        self.ax.text(22, POLE_Y + 1, "POLES", color="#46515c", fontsize=7)
        self.ax.text(22, ZERO_Y - 4, "ZEROS", color="#46515c", fontsize=7)

        self.ax.plot(freqs, m0_db, color="#5f95c4", lw=1.0, alpha=0.5, label="M0")
        self.ax.plot(freqs, m100_db, color="#c4885a", lw=1.0, alpha=0.5, label="M100")
        self.ax.plot(freqs, cur_db, color="#d8ff9a", lw=2.4, label=f"M={t:.2f}", zorder=4)

        self._handle_index = []
        for endpoint, stages, p_color, z_color in (
            ("m0",   self.m0,   "#8fd3ff", "#5f95c4"),
            ("m100", self.m100, "#ffb36b", "#c4885a"),
        ):
            for i, st in enumerate(stages):
                if st.flag == 0:
                    continue
                self.ax.scatter([st.pole_freq_hz], [POLE_Y], marker="X",
                                s=130, c=p_color, edgecolors="white",
                                linewidth=1.2, zorder=6)
                self.ax.text(st.pole_freq_hz, POLE_Y + 1.5, f"S{i}",
                             color=p_color, fontsize=7, ha="center")
                self._handle_index.append({"endpoint": endpoint, "kind": "pole",
                                           "stage_idx": i,
                                           "freq": st.pole_freq_hz,
                                           "radius": st.pole_r})
                self.ax.scatter([st.zero_freq_hz], [ZERO_Y], marker="o",
                                s=85, c="none", edgecolors=z_color,
                                linewidth=1.4, zorder=6)
                self._handle_index.append({"endpoint": endpoint, "kind": "zero",
                                           "stage_idx": i,
                                           "freq": st.zero_freq_hz,
                                           "radius": st.zero_r})

        self.ax.set_title(
            f"{self.cart_name}   c0={c0:.4f}   "
            f"M0_peak={m0_db.max():+.1f}   M100_peak={m100_db.max():+.1f}   "
            f"now={cur_db.max():+.1f} dB",
            color="#e6edf2", fontsize=9)
        self.ax.legend(loc="lower left", facecolor="#11161b",
                       labelcolor="#e6edf2", fontsize=8, frameon=False)
        try:
            self.canvas.draw_idle()
        except Exception:
            pass

    # ─── Drag ────────────────────────────────────────────────────────────

    def _nearest_handle(self, x_data, y_data):
        if not self._handle_index or x_data is None or y_data is None:
            return None
        best = None
        best_d = math.inf
        log_x = math.log10(max(20.0, x_data))
        for h in self._handle_index:
            target_y = POLE_Y if h["kind"] == "pole" else ZERO_Y
            log_f = math.log10(max(20.0, h["freq"]))
            dx = abs(log_x - log_f)
            dy = abs(y_data - target_y) / 80.0
            d = math.hypot(dx, dy)
            if d < best_d:
                best_d = d
                best = dict(h, distance=d)
        return best

    def _on_press(self, event):
        if event.inaxes != self.ax or event.button != 1:
            return
        nearest = self._nearest_handle(event.xdata, event.ydata)
        if nearest is None or nearest["distance"] > PICK_THRESH:
            return
        self._drag = {
            "endpoint": nearest["endpoint"],
            "kind": nearest["kind"],
            "stage_idx": nearest["stage_idx"],
            "anchor_y": event.ydata,
            "anchor_radius": nearest["radius"],
        }
        self.status.config(text=f"dragging {nearest['endpoint']} S{nearest['stage_idx']} "
                                f"{nearest['kind']}  freq={nearest['freq']:.0f}  "
                                f"r={nearest['radius']:.4f}")

    def _on_motion(self, event):
        if self._drag is None or event.inaxes != self.ax:
            return
        if event.xdata is None or event.ydata is None:
            return
        new_freq = max(FREQ_MIN_AUTH, min(FREQ_MAX_AUTH, float(event.xdata)))
        y_delta = event.ydata - self._drag["anchor_y"]
        new_radius = self._drag["anchor_radius"] + y_delta * RADIUS_SENS
        if self._drag["kind"] == "pole":
            new_radius = max(POLE_R_MIN, min(POLE_R_MAX, new_radius))
        else:
            new_radius = max(ZERO_R_MIN, min(ZERO_R_MAX, new_radius))
        ep_list = self.m0 if self._drag["endpoint"] == "m0" else self.m100
        with self.lock:
            st = ep_list[self._drag["stage_idx"]]
            if self._drag["kind"] == "pole":
                st.pole_freq_hz = new_freq
                st.pole_r = new_radius
            else:
                st.zero_freq_hz = new_freq
                st.zero_r = new_radius
        self.status.config(text=f"{self._drag['endpoint']} S{self._drag['stage_idx']} "
                                f"{self._drag['kind']}  freq={new_freq:.0f}  "
                                f"r={new_radius:.4f}")
        self._render_plot()

    def _on_release(self, event):
        if self._drag is not None:
            self._drag = None
            self.status.config(text="committed")

    # ─── c0 ──────────────────────────────────────────────────────────────

    def _on_c0_change(self, _event=None):
        try:
            v = float(self.c0_var.get())
        except ValueError:
            self.c0_var.set(f"{self.c0:.4f}")
            return
        v = max(0.001, min(1.0, v))
        with self.lock:
            self.c0 = v
        self.c0_var.set(f"{v:.4f}")
        self._render_plot()

    # ─── Morph ───────────────────────────────────────────────────────────

    def _on_morph(self, val):
        with self.lock:
            self.morph = float(val)
        self.morph_label.config(text=f"{self.morph:.3f}")
        self._render_plot()

    # ─── Audio ───────────────────────────────────────────────────────────

    def _audio_callback(self, outdata, frames, time_info, status):
        if status:
            print(status, file=sys.stderr)
        with self.lock:
            if self.paused:
                outdata[:] = 0.0
                return
            t = self.morph
            c0 = self.c0
            m0 = list(self.m0)
            m100 = list(self.m100)
        sos_rows = []
        for sa, sb in zip(m0, m100):
            st = lerp_stage(sa, sb, t)
            c = stage_to_coefs(st, c0, PLAYBACK_SR)
            if c["c3"] == 0.0 and c["c4"] == 0.0:
                continue
            sos_rows.append([c["c0"], c["c1"], c["c2"], 1.0, c["c3"], c["c4"]])
        end = self.source_pos + frames
        if end <= len(self.source):
            block = self.source[self.source_pos:end]
            self.source_pos = end
        else:
            n_remain = len(self.source) - self.source_pos
            block = np.concatenate([
                self.source[self.source_pos:],
                self.source[: frames - n_remain],
            ])
            self.source_pos = frames - n_remain
        if not sos_rows:
            outdata[:, 0] = block.astype(np.float32) * 0.5
            return
        sos = np.asarray(sos_rows, dtype=np.float64)
        if self.state is None or self.state.shape[0] != len(sos_rows):
            self.state = np.zeros((len(sos_rows), 2), dtype=np.float64)
        out, self.state = sosfilt(sos, block, zi=self.state)
        out = np.tanh(out * 0.45) * 0.6
        outdata[:, 0] = out.astype(np.float32)

    def _toggle_play(self):
        if sd is None:
            self.status.config(text="ERROR: pip install sounddevice"); return
        if sosfilt is None:
            self.status.config(text="ERROR: pip install scipy"); return
        if not self.playing:
            self._start(); return
        with self.lock:
            self.paused = not self.paused
        if self.paused:
            self.play_btn.config(text="PLAY", bg="#22aa44")
            self.status.config(text="paused")
        else:
            self.play_btn.config(text="PAUSE", bg="#cc8833")
            self.status.config(text="streaming — drag handles live")

    def _start(self):
        try:
            self.stream = sd.OutputStream(
                samplerate=PLAYBACK_SR, channels=1, blocksize=BLOCK,
                callback=self._audio_callback, dtype="float32",
            )
            self.stream.start()
            self.playing = True
            self.paused = False
            self.play_btn.config(text="PAUSE", bg="#cc8833")
            self.status.config(text="streaming — drag handles live")
        except Exception as exc:
            self.status.config(text=f"audio error: {exc}")

    def _stop(self):
        if self.stream is not None:
            try:
                self.stream.stop()
                self.stream.close()
            except Exception:
                pass
            self.stream = None
        self.state = None
        self.playing = False
        self.paused = False

    # ─── KEEP / REJECT / AGAIN ───────────────────────────────────────────

    def _on_keep(self):
        KEPT_DIR.mkdir(parents=True, exist_ok=True)
        dest = KEPT_DIR / self.cart_path.name
        m0_baked = [authoring_to_rawstage(s, self.c0) for s in self.m0]
        m100_baked = [authoring_to_rawstage(s, self.c0) for s in self.m100]
        passthrough = {"a1": 0.0, "r": 0.0, "val1": 0.0, "val2": 0.0, "val3": 0.0, "flag": 0}
        while len(m0_baked) < 12:
            m0_baked.append(dict(passthrough))
        while len(m100_baked) < 12:
            m100_baked.append(dict(passthrough))
        payload = {
            "name": self.cart_name,
            "filterType": 13,
            "sampleRate": AUTHORING_SR,
            "stages": N_STAGES,
            "boost": 4.0,
            "provenance": {"author": "sift-sculptor", "c0_uniform": self.c0},
            "keyframes": [
                {"morph": 0.0, "q": 0.0, "label": "M0_Q0",     "boost": 4.0, "stages": m0_baked},
                {"morph": 0.0, "q": 1.0, "label": "M0_Q100",   "boost": 4.0, "stages": m0_baked},
                {"morph": 1.0, "q": 0.0, "label": "M100_Q0",   "boost": 4.0, "stages": m100_baked},
                {"morph": 1.0, "q": 1.0, "label": "M100_Q100", "boost": 4.0, "stages": m100_baked},
            ],
        }
        cart_text = json.dumps(payload, indent=2)
        dest.write_text(cart_text, encoding="utf-8")
        # Also publish to the plugin's hot-reload slot for instant DAW audition.
        try:
            PLUGIN_SLOT.parent.mkdir(parents=True, exist_ok=True)
            PLUGIN_SLOT.write_text(cart_text, encoding="utf-8")
            print(f"KEEP -> {dest}\n     -> {PLUGIN_SLOT} (plugin slot)")
            self.status.config(text=f"KEPT -> {dest.name}  +  plugin slot")
        except Exception as exc:
            print(f"KEEP -> {dest}  (plugin slot write failed: {exc})")
            self.status.config(text=f"KEPT -> {dest.name}  (plugin slot fail)")

    def _on_reject(self):
        print(f"REJECT {self.cart_path.name}")
        self.status.config(text="REJECTED")

    def _on_again(self):
        rng = np.random.default_rng()
        with self.lock:
            for stages in (self.m0, self.m100):
                for st in stages:
                    if st.flag == 0:
                        continue
                    st.pole_freq_hz *= 1.0 + rng.uniform(-0.05, 0.05)
                    st.pole_freq_hz = max(FREQ_MIN_AUTH, min(FREQ_MAX_AUTH, st.pole_freq_hz))
                    st.zero_freq_hz *= 1.0 + rng.uniform(-0.05, 0.05)
                    st.zero_freq_hz = max(FREQ_MIN_AUTH, min(FREQ_MAX_AUTH, st.zero_freq_hz))
                    st.pole_r = max(POLE_R_MIN, min(POLE_R_MAX,
                                    st.pole_r * (1.0 + rng.uniform(-0.01, 0.01))))
                    st.zero_r = max(ZERO_R_MIN, min(ZERO_R_MAX,
                                    st.zero_r * (1.0 + rng.uniform(-0.05, 0.05))))
        self._render_plot()
        self.status.config(text="AGAIN — perturbed")

    def _on_close(self):
        self._stop()
        self.quit()
        self.destroy()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("cartridge", type=Path)
    args = ap.parse_args()
    if not args.cartridge.exists():
        print(f"cartridge not found: {args.cartridge}", file=sys.stderr); return 1
    Sculptor(args.cartridge).mainloop()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
