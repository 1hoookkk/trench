"""forge/authoring_studio.py — minimal hand-authoring GUI.

Scope:
  - corner tabs (M0 / M100)
  - stage list with spinboxes (pole_hz, pole_r, zero_hz, zero_r, c0)
  - z-plane editor for selected stage (drag X = pole, drag O = zero)
  - cascade plot (M0 + M100)
  - + Add stage / ✕ Delete
  - Audition (renders WAVs through forge/audition.py)
  - Save / Load cartridge JSON

State auto-persists in cartridges/factory/generated/qlaw/_studio_state.json.

Run:  python forge/authoring_studio.py
"""
from __future__ import annotations

import json
import math
import os
import subprocess
import sys
import threading
import tkinter as tk
from pathlib import Path
from tkinter import ttk, messagebox, filedialog

import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
import numpy as np


REPO = Path(__file__).resolve().parent.parent
OUT_DIR = REPO / "cartridges" / "factory" / "generated" / "qlaw"
STATE_PATH = OUT_DIR / "_studio_state.json"
DEFAULT_CART_PATH = OUT_DIR / "_studio_cartridge.json"
AUDITION_SCRIPT = REPO / "forge" / "audition.py"

AUTHORING_SR = 39062.5
NYQUIST = AUTHORING_SR * 0.5
STATE_SCHEMA = "studio.v3"

PARAMS = [
    # name,      lo,        hi,        step,     fmt
    ("pole_hz", 20.0,     19500.0,    1.0,    "%.0f"),
    ("pole_r",  0.50,     0.9999,     0.001,  "%.4f"),
    ("zero_hz", 1.0,      19500.0,    1.0,    "%.0f"),
    ("zero_r",  0.0,      1.0,        0.01,   "%.3f"),
    ("c0",      0.001,    1.0,        0.01,   "%.3f"),
]
STAGE_COLORS = ["#22ddff", "#88ee99", "#ffeb3b", "#ff9933", "#ff4488",
                "#cc66ff", "#ff8855", "#66ddaa", "#ddaaff", "#aabb22",
                "#ff77bb", "#44ccdd"]


# ─── DSP ──────────────────────────────────────────────────────────────────────

def make_stage(pole_hz: float = 1000.0) -> dict:
    pole_hz = float(max(20.0, min(NYQUIST - 100, pole_hz)))
    zero_hz = max(1.0, min(NYQUIST - 100, pole_hz * (2.0 ** (-3.0/12.0))))
    return {"pole_hz": pole_hz, "pole_r": 0.99, "zero_hz": zero_hz,
            "zero_r": 0.92, "c0": 0.85}


def stage_to_compiled(s: dict, sr: float = AUTHORING_SR) -> dict:
    """Studio convention: c0 multiplies the unit-gain numerator polynomial.
    c1 = c0·b1, c2 = c0·b2 — keeps zero locations stable under c0 attenuation."""
    pole_hz = float(s["pole_hz"]); pole_r = float(s["pole_r"])
    zero_hz = float(s["zero_hz"]); zero_r = float(s["zero_r"])
    c0 = float(s["c0"])
    if zero_hz >= sr * 0.5: zero_hz = sr * 0.49
    theta_p = 2 * math.pi * pole_hz / sr
    theta_z = 2 * math.pi * zero_hz / sr
    b1 = -2.0 * zero_r * math.cos(theta_z)
    b2 = zero_r * zero_r
    return {"c0": c0, "c1": c0 * b1, "c2": c0 * b2,
            "c3": -2.0 * pole_r * math.cos(theta_p), "c4": pole_r * pole_r}


def cascade_db(compiled: list[dict], freqs: np.ndarray, sr: float = AUTHORING_SR) -> np.ndarray:
    w = 2 * np.pi * freqs / sr
    e1 = np.exp(-1j * w); e2 = np.exp(-2j * w)
    mag = np.ones_like(freqs)
    for s in compiled:
        if s["c3"] == 0.0 and s["c4"] == 0.0: continue
        num = s["c0"] + s["c1"]*e1 + s["c2"]*e2
        den = 1.0 + s["c3"]*e1 + s["c4"]*e2
        mag *= np.abs(num/den)
    return 20*np.log10(np.maximum(mag, 1e-12))


def build_cartridge_json(state: dict, name: str) -> dict:
    def compile_corner(stages):
        out = [stage_to_compiled(s) for s in stages]
        while len(out) < 12:
            out.append({"c0":1.,"c1":0.,"c2":0.,"c3":0.,"c4":0.})
        return out
    return {
        "format": "compiled-v1",
        "name": name,
        "sampleRate": AUTHORING_SR,
        "recipe": "studio_hand_authored",
        "provenance": {"mode": "studio", "schema": STATE_SCHEMA,
                       "stages_m0": state["M0"], "stages_m100": state["M100"]},
        "keyframes": [
            {"label": "M0_Q0",    "morph": 0.0, "q": 0.0, "boost": 1.0, "stages": compile_corner(state["M0"])},
            {"label": "M0_Q100",  "morph": 0.0, "q": 1.0, "boost": 1.0, "stages": compile_corner(state["M0"])},
            {"label": "M100_Q0",  "morph": 1.0, "q": 0.0, "boost": 1.0, "stages": compile_corner(state["M100"])},
            {"label": "M100_Q100","morph": 1.0, "q": 1.0, "boost": 1.0, "stages": compile_corner(state["M100"])},
        ],
    }


# ─── App ──────────────────────────────────────────────────────────────────────

def empty_state() -> dict:
    return {"schema": STATE_SCHEMA, "M0": [], "M100": []}


class StudioApp:
    def __init__(self, root: tk.Tk):
        self.root = root
        self.root.title("TRENCH Authoring Studio")
        self.root.geometry("1500x900")
        self.root.configure(bg="#16191a")

        self.state = self._load_state()
        self.active_corner = "M0"
        self.selected_idx: int | None = None
        self.stage_vars: dict[tuple[str,int,str], tk.DoubleVar] = {}
        self._drag_state: dict | None = None
        self._zp_drag_state: dict | None = None
        self._save_after_id = None

        self._build_layout()
        self._refresh_all()

    def _load_state(self) -> dict:
        if not STATE_PATH.exists(): return empty_state()
        try:
            d = json.loads(STATE_PATH.read_text())
            if d.get("schema") != STATE_SCHEMA: return empty_state()
            return {"schema": STATE_SCHEMA, "M0": list(d.get("M0", [])), "M100": list(d.get("M100", []))}
        except Exception:
            return empty_state()

    def _save_state(self):
        if self._save_after_id is not None:
            try: self.root.after_cancel(self._save_after_id)
            except Exception: pass
        self._save_after_id = self.root.after(250, self._save_state_now)

    def _save_state_now(self):
        self._save_after_id = None
        try:
            OUT_DIR.mkdir(parents=True, exist_ok=True)
            STATE_PATH.write_text(json.dumps(self.state, indent=2))
        except Exception as e:
            self.status_var.set(f"save failed: {e}")

    def _build_layout(self):
        main = tk.Frame(self.root, bg="#16191a")
        main.pack(fill=tk.BOTH, expand=True)

        # Left: corner selector + stage list + buttons
        ctl = tk.Frame(main, bg="#16191a", width=620)
        ctl.pack(side=tk.LEFT, fill=tk.Y, padx=10, pady=10)
        ctl.pack_propagate(False)

        ts = tk.Frame(ctl, bg="#16191a"); ts.pack(fill=tk.X, pady=4)
        tk.Label(ts, text="Corner:", bg="#16191a", fg="white",
                 font=("Arial",10,"bold")).pack(side=tk.LEFT, padx=4)
        self.corner_var = tk.StringVar(value="M0")
        for c in ("M0","M100"):
            tk.Radiobutton(ts, text=c, variable=self.corner_var, value=c,
                           command=self._on_corner_change,
                           bg="#16191a", fg="white", selectcolor="#2a2d2e",
                           activebackground="#16191a", activeforeground="white",
                           font=("Arial",10,"bold")).pack(side=tk.LEFT, padx=8)

        ab = tk.Frame(ctl, bg="#16191a"); ab.pack(fill=tk.X, pady=4)
        tk.Button(ab, text="+ Add stage", command=self._add_stage,
                  bg="#22ddff", fg="black", font=("Arial",10,"bold"),
                  padx=12, pady=4).pack(side=tk.LEFT, padx=4)
        tk.Button(ab, text="Clear corner", command=self._clear_corner,
                  bg="#444", fg="white", padx=8).pack(side=tk.LEFT, padx=4)

        sf = tk.Frame(ctl, bg="#16191a"); sf.pack(fill=tk.BOTH, expand=True, pady=4)
        self.stage_canvas = tk.Canvas(sf, bg="#16191a", highlightthickness=0)
        sb = ttk.Scrollbar(sf, orient=tk.VERTICAL, command=self.stage_canvas.yview)
        self.stage_canvas.configure(yscrollcommand=sb.set)
        self.stage_canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        sb.pack(side=tk.RIGHT, fill=tk.Y)
        self.stage_inner = tk.Frame(self.stage_canvas, bg="#16191a")
        self.stage_canvas.create_window((0,0), window=self.stage_inner, anchor="nw")
        self.stage_inner.bind("<Configure>",
            lambda e: self.stage_canvas.configure(scrollregion=self.stage_canvas.bbox("all")))

        bf = tk.Frame(ctl, bg="#16191a", pady=8); bf.pack(fill=tk.X, side=tk.BOTTOM)
        tk.Button(bf, text="Audition", command=self._audition,
                  bg="#ff5555", fg="white", font=("Arial",10,"bold"),
                  padx=12, pady=4).pack(side=tk.LEFT, padx=4)
        tk.Button(bf, text="Save…", command=self._save_cartridge,
                  bg="#ffba00", fg="black", font=("Arial",10,"bold"),
                  padx=10).pack(side=tk.LEFT, padx=4)
        tk.Button(bf, text="Load…", command=self._load_cartridge,
                  bg="#88ee99", fg="black", padx=10).pack(side=tk.LEFT, padx=4)
        self.status_var = tk.StringVar(value="ready. + Add stage to begin.")
        tk.Label(bf, textvariable=self.status_var, bg="#16191a", fg="#ffba00",
                 font=("Consolas",9), anchor="w").pack(side=tk.LEFT, padx=14, fill=tk.X, expand=True)

        # Right: cascade + z-plane plots
        right = tk.Frame(main, bg="#16191a")
        right.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=10, pady=10)

        self.fig = plt.figure(figsize=(10, 9), dpi=100)
        self.fig.patch.set_facecolor("#16191a")
        gs = self.fig.add_gridspec(3, 2, height_ratios=[3, 3, 4],
                                   hspace=0.35, wspace=0.18,
                                   left=0.07, right=0.97, top=0.96, bottom=0.06)
        self.ax_cas_m0   = self.fig.add_subplot(gs[0, :])
        self.ax_cas_m100 = self.fig.add_subplot(gs[1, :], sharex=self.ax_cas_m0)
        self.ax_zp_m0    = self.fig.add_subplot(gs[2, 0])
        self.ax_zp_m100  = self.fig.add_subplot(gs[2, 1])

        self.canvas = FigureCanvasTkAgg(self.fig, master=right)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        self.freqs = np.geomspace(20.0, 18000.0, 4096)
        self.freqs_fast = np.geomspace(20.0, 18000.0, 768)

        tk.Label(right,
            text="Cascade: click empty=add, drag handle=freq+c0   |   "
                 "Z-plane: click stage to select, drag X (pole) or O (zero)",
            bg="#16191a", fg="#888", font=("Consolas",9)).pack(side=tk.BOTTOM, pady=4)

        self._handle_radius_pts = 14
        self._zp_radius_pts = 12
        self.canvas.mpl_connect("button_press_event",   self._on_press)
        self.canvas.mpl_connect("motion_notify_event",  self._on_motion)
        self.canvas.mpl_connect("button_release_event", self._on_release)
        self.canvas.mpl_connect("scroll_event",         self._on_scroll)

    # ─── Stage rows ────────────────────────────────────────────────────────

    def _render_stages(self):
        for w in self.stage_inner.winfo_children(): w.destroy()
        self.stage_vars.clear()

        corner = self.active_corner
        stages = self.state[corner]
        if not stages:
            tk.Label(self.stage_inner,
                     text=f"\n{corner} empty.\nClick + Add stage.\n",
                     bg="#16191a", fg="#888", font=("Arial",11)).pack(pady=40)
            return

        # Header
        hdr = tk.Frame(self.stage_inner, bg="#16191a"); hdr.pack(fill=tk.X, padx=4, pady=2)
        tk.Label(hdr, text="#", bg="#16191a", fg="white", width=3,
                 font=("Arial",8,"bold")).grid(row=0, column=0)
        for col, (pname, *_) in enumerate(PARAMS):
            tk.Label(hdr, text=pname, bg="#16191a", fg="white", width=9,
                     font=("Arial",8,"bold")).grid(row=0, column=col+1, padx=2)
        tk.Label(hdr, text="", width=3).grid(row=0, column=len(PARAMS)+1)

        for i, s in enumerate(stages):
            color = STAGE_COLORS[i % len(STAGE_COLORS)]
            is_sel = (self.selected_idx == i)
            border = "#ffffff" if is_sel else color
            thickness = 3 if is_sel else 1
            row = tk.Frame(self.stage_inner, bg="#0d0f10",
                           highlightbackground=border, highlightthickness=thickness)
            row.pack(fill=tk.X, padx=4, pady=2)

            sbtn = tk.Label(row, text=f"S{i}", bg="#0d0f10", fg=color, width=3,
                            font=("Arial",10,"bold"), cursor="hand2")
            sbtn.grid(row=0, column=0, padx=2, pady=4)
            sbtn.bind("<Button-1>", lambda e, idx=i: self._select_stage(idx))

            for col, (pname, lo, hi, step, fmt) in enumerate(PARAMS):
                var = tk.DoubleVar(value=float(s[pname]))
                self.stage_vars[(corner, i, pname)] = var
                spin = tk.Spinbox(row, from_=lo, to=hi, increment=step, format=fmt,
                                  textvariable=var, width=8,
                                  bg="#1a1d1e", fg="white", insertbackground="white",
                                  buttonbackground="#333", font=("Consolas",9))
                spin.grid(row=0, column=col+1, padx=2, pady=4)
                var.trace_add("write",
                              lambda *_, c=corner, idx=i, p=pname: self._on_var(c, idx, p))

            tk.Button(row, text="✕", command=lambda idx=i: self._remove_stage(idx),
                      bg="#660000", fg="white", width=3,
                      font=("Arial",9,"bold")).grid(row=0, column=len(PARAMS)+1, padx=2)

    def _on_var(self, corner, idx, pname):
        try:
            val = float(self.stage_vars[(corner, idx, pname)].get())
        except (tk.TclError, ValueError):
            return
        if idx >= len(self.state[corner]): return
        self.state[corner][idx][pname] = val
        self._save_state()
        self._refresh_plot()

    def _on_corner_change(self):
        self.active_corner = self.corner_var.get()
        self.selected_idx = None
        self._refresh_all()

    def _select_stage(self, idx):
        if idx is not None and idx >= len(self.state[self.active_corner]): idx = None
        self.selected_idx = idx
        self._refresh_all()

    def _add_stage(self):
        self.state[self.active_corner].append(make_stage(1000.0))
        self.selected_idx = len(self.state[self.active_corner]) - 1
        self._save_state(); self._refresh_all()
        self.status_var.set(f"+ S{self.selected_idx} @ 1000 Hz")

    def _remove_stage(self, idx):
        if idx >= len(self.state[self.active_corner]): return
        self.state[self.active_corner].pop(idx)
        if self.selected_idx == idx: self.selected_idx = None
        elif self.selected_idx is not None and self.selected_idx > idx:
            self.selected_idx -= 1
        self._save_state(); self._refresh_all()
        self.status_var.set(f"removed S{idx}")

    def _clear_corner(self):
        if not messagebox.askyesno("Clear", f"Remove all stages from {self.active_corner}?"): return
        self.state[self.active_corner].clear()
        self.selected_idx = None
        self._save_state(); self._refresh_all()

    # ─── Plots ─────────────────────────────────────────────────────────────

    def _refresh_all(self):
        self._render_stages(); self._refresh_plot()

    def _draw_zplane(self, ax, corner: str):
        ax.clear()
        ax.set_facecolor("#0d0f10"); ax.set_aspect("equal")
        ax.set_xlim(-1.15, 1.15); ax.set_ylim(-1.15, 1.15)
        theta = np.linspace(0, 2*np.pi, 200)
        ax.plot(np.cos(theta), np.sin(theta), color="#666", lw=1.0)
        ax.axhline(0, color="#444", lw=0.5); ax.axvline(0, color="#444", lw=0.5)
        ax.tick_params(colors="white", labelsize=7)
        for spine in ax.spines.values(): spine.set_color("#3a3e42")

        idx = self.selected_idx
        title_color = "#ff5555" if corner == "M0" else "#ff7777"
        if idx is None or idx >= len(self.state[corner]):
            ax.set_title(f"{corner} z-plane — click stage to select",
                         color="#888", fontsize=10)
            return
        s = self.state[corner][idx]
        stage_color = STAGE_COLORS[idx % len(STAGE_COLORS)]

        pole_r = float(s["pole_r"])
        theta_p = 2*math.pi * float(s["pole_hz"]) / AUTHORING_SR
        zero_r = float(s["zero_r"])
        theta_z = 2*math.pi * float(s["zero_hz"]) / AUTHORING_SR

        ax.plot([pole_r*math.cos(theta_p)], [pole_r*math.sin(theta_p)], "x",
                color=stage_color, markersize=14, markeredgewidth=3, zorder=5)
        ax.plot([pole_r*math.cos(theta_p)], [-pole_r*math.sin(theta_p)], "x",
                color=stage_color, markersize=10, markeredgewidth=2, alpha=0.55, zorder=4)
        ax.plot([zero_r*math.cos(theta_z)], [zero_r*math.sin(theta_z)], "o",
                color="none", markeredgecolor=stage_color, markersize=14,
                markeredgewidth=3, zorder=5)
        ax.plot([zero_r*math.cos(theta_z)], [-zero_r*math.sin(theta_z)], "o",
                color="none", markeredgecolor=stage_color, markersize=10,
                markeredgewidth=2, alpha=0.55, zorder=4)

        ax.set_title(f"{corner} z-plane — S{idx}",
                     color=title_color, fontsize=10, fontweight="bold")
        ax.text(0.02, 0.98,
                f"pole {s['pole_hz']:.0f}Hz r={pole_r:.3f}\n"
                f"zero {s['zero_hz']:.0f}Hz r={zero_r:.3f}",
                transform=ax.transAxes, color=stage_color, fontsize=8,
                family="monospace", verticalalignment="top",
                bbox=dict(facecolor="#0d0f10", edgecolor="#3a3e42", alpha=0.8))

    def _refresh_plot(self):
        dragging = (self._drag_state is not None) or (self._zp_drag_state is not None)
        freqs = self.freqs_fast if dragging else self.freqs

        for ax in (self.ax_cas_m0, self.ax_cas_m100):
            ax.clear(); ax.set_facecolor("#0d0f10")

        for corner, ax in (("M0", self.ax_cas_m0), ("M100", self.ax_cas_m100)):
            stages = self.state[corner]
            if not stages:
                ax.text(0.5, 0.5, f"{corner} empty — click to add a stage",
                        transform=ax.transAxes, ha="center", va="center",
                        color="#666", fontsize=14, fontweight="bold")
                ax.set_xscale("log"); ax.set_xlim(20, 18000); ax.set_ylim(-40, 30)
                continue
            compiled = [stage_to_compiled(s) for s in stages]
            cas = cascade_db(compiled, freqs)
            peak = float(np.max(cas))
            color = "#ff5555" if corner == "M0" else "#ff7777"
            ax.semilogx(freqs, cas, color=color, lw=2.4)
            ax.set_title(f"{corner} cascade ({peak:+.1f} dB,  {len(stages)} stages)",
                         color=color, fontsize=10, fontweight="bold")

            for i, s in enumerate(stages):
                hx = s["pole_hz"]
                hy = 20.0 * math.log10(max(float(s["c0"]), 1e-6))
                stage_color = STAGE_COLORS[i % len(STAGE_COLORS)]
                is_sel = (corner == self.active_corner) and (self.selected_idx == i)
                ec = "#ffffff" if is_sel else "#000000"
                lw = 3.0 if is_sel else 1.0
                size = 200 if is_sel else 140
                ax.scatter([hx], [hy], s=size, c=[stage_color], edgecolors=ec,
                           linewidths=lw, zorder=5)
                ax.text(hx, hy, str(i), ha="center", va="center",
                        color="black", fontsize=8, fontweight="bold", zorder=6)

        for ax in (self.ax_cas_m0, self.ax_cas_m100):
            ax.set_xlim(20, 18000); ax.set_ylim(-50, 40)
            ax.grid(True, which="both", alpha=0.18, color="white")
            ax.tick_params(colors="white", labelsize=8)
            for spine in ax.spines.values(): spine.set_color("#3a3e42")
            ax.set_ylabel("dB", color="white", fontsize=9)
            ax.axhline(0.0, color="#444", lw=0.5, alpha=0.5)
        self.ax_cas_m100.set_xlabel("Hz", color="white", fontsize=9)

        self._draw_zplane(self.ax_zp_m0,   "M0")
        self._draw_zplane(self.ax_zp_m100, "M100")
        self.canvas.draw_idle()

    # ─── Mouse interaction ─────────────────────────────────────────────────

    def _corner_for_axis(self, ax):
        if ax is self.ax_cas_m0 or ax is self.ax_zp_m0: return "M0"
        if ax is self.ax_cas_m100 or ax is self.ax_zp_m100: return "M100"
        return None

    def _is_zp(self, ax):
        return ax is self.ax_zp_m0 or ax is self.ax_zp_m100

    def _hit_cas(self, event):
        if event.inaxes is None or event.x is None: return None
        if self._is_zp(event.inaxes): return None
        corner = self._corner_for_axis(event.inaxes)
        if corner is None: return None
        for i, s in enumerate(self.state[corner]):
            hx = s["pole_hz"]
            hy = 20.0 * math.log10(max(float(s["c0"]), 1e-6))
            px, py = event.inaxes.transData.transform((hx, hy))
            if math.hypot(px - event.x, py - event.y) <= self._handle_radius_pts:
                return (corner, i)
        return None

    def _hit_zp(self, event):
        if event.inaxes is None or event.x is None: return None
        if not self._is_zp(event.inaxes): return None
        corner = self._corner_for_axis(event.inaxes)
        if corner is None or self.selected_idx is None: return None
        idx = self.selected_idx
        if idx >= len(self.state[corner]): return None
        s = self.state[corner][idx]
        ax = event.inaxes
        pr = float(s["pole_r"]); thp = 2*math.pi * float(s["pole_hz"]) / AUTHORING_SR
        for sign in (+1, -1):
            px, py = ax.transData.transform((pr*math.cos(thp), sign*pr*math.sin(thp)))
            if math.hypot(px - event.x, py - event.y) <= self._zp_radius_pts:
                return (corner, idx, "pole")
        zr = float(s["zero_r"]); thz = 2*math.pi * float(s["zero_hz"]) / AUTHORING_SR
        for sign in (+1, -1):
            zx, zy = ax.transData.transform((zr*math.cos(thz), sign*zr*math.sin(thz)))
            if math.hypot(zx - event.x, zy - event.y) <= self._zp_radius_pts:
                return (corner, idx, "zero")
        return None

    def _on_press(self, event):
        if event.inaxes is None: return
        if self._is_zp(event.inaxes):
            hit = self._hit_zp(event)
            if hit is None: return
            corner, idx, kind = hit
            if corner != self.active_corner:
                self.active_corner = corner; self.corner_var.set(corner)
            self.selected_idx = idx
            self._zp_drag_state = {"corner": corner, "idx": idx, "kind": kind}
            self.status_var.set(f"dragging {kind} of S{idx}")
            return

        # Cascade plot
        corner = self._corner_for_axis(event.inaxes)
        hit = self._hit_cas(event)
        if event.button == 3 and hit is not None:
            self._remove_stage(hit[1]); return
        if event.button == 1:
            if hit is not None:
                hcorner, idx = hit
                if hcorner != self.active_corner:
                    self.active_corner = hcorner; self.corner_var.set(hcorner)
                self.selected_idx = idx
                self._render_stages()
                self._drag_state = {"corner": hcorner, "idx": idx}
                self.status_var.set(f"dragging S{idx}")
            else:
                if event.xdata is None: return
                hz = float(event.xdata)
                if hz < 20 or hz > NYQUIST: return
                if corner != self.active_corner:
                    self.active_corner = corner; self.corner_var.set(corner)
                self.state[corner].append(make_stage(hz))
                self.selected_idx = len(self.state[corner]) - 1
                self._save_state(); self._refresh_all()
                self.status_var.set(f"+ S{self.selected_idx} @ {hz:.0f} Hz")

    def _on_motion(self, event):
        if self._zp_drag_state is not None:
            if event.inaxes is None or event.xdata is None: return
            if not self._is_zp(event.inaxes): return
            ds = self._zp_drag_state
            if ds["corner"] != self._corner_for_axis(event.inaxes): return
            x, y = float(event.xdata), float(event.ydata)
            r = max(0.001, min(0.9999, math.hypot(x, y)))
            theta = max(0.0, min(math.pi - 1e-4, math.atan2(abs(y), x)))
            hz = max(20.0, min(NYQUIST - 100, theta * AUTHORING_SR / (2*math.pi)))
            s = self.state[ds["corner"]][ds["idx"]]
            if ds["kind"] == "pole":
                s["pole_r"] = r; s["pole_hz"] = hz
            else:
                s["zero_r"] = r; s["zero_hz"] = hz
            self._refresh_plot()
            self.status_var.set(f"S{ds['idx']} {ds['kind']}: r={r:.4f} hz={hz:.0f}")
            return

        if self._drag_state is not None:
            if event.inaxes is None or event.xdata is None: return
            ds = self._drag_state
            if ds["idx"] >= len(self.state[ds["corner"]]): return
            s = self.state[ds["corner"]][ds["idx"]]
            hz = float(max(20.0, min(NYQUIST - 100, event.xdata)))
            c0 = max(0.001, min(1.0, 10.0 ** (float(event.ydata) / 20.0)))
            s["pole_hz"] = hz; s["c0"] = c0
            self._refresh_plot()
            self.status_var.set(f"S{ds['idx']}: {hz:.0f} Hz  c0={c0:.3f}")

    def _on_release(self, event):
        if self._zp_drag_state is not None:
            self._save_state(); self._render_stages(); self._zp_drag_state = None
            return
        if self._drag_state is not None:
            self._save_state(); self._render_stages(); self._drag_state = None

    def _on_scroll(self, event):
        hit = self._hit_cas(event)
        if hit is None: return
        corner, idx = hit
        s = self.state[corner][idx]
        delta = 0.002 * (1 if event.button == "up" else -1)
        s["pole_r"] = max(0.50, min(0.9999, float(s["pole_r"]) + delta))
        self._save_state(); self._refresh_plot()
        if corner == self.active_corner:
            try: self.stage_vars[(corner, idx, "pole_r")].set(s["pole_r"])
            except KeyError: pass
        self.status_var.set(f"S{idx} pole_r={s['pole_r']:.4f}")

    # ─── Audition / save / load ────────────────────────────────────────────

    def _audition(self):
        if not self.state["M0"] and not self.state["M100"]:
            self.status_var.set("nothing to audition — add stages first"); return
        cart = build_cartridge_json(self.state, "_studio_cartridge")
        DEFAULT_CART_PATH.write_text(json.dumps(cart, indent=2))
        self.status_var.set("rendering…"); self.root.update_idletasks()
        def run():
            try:
                env = {**os.environ, "PYTHONIOENCODING":"utf-8", "PYTHONUTF8":"1"}
                p = subprocess.run(
                    [sys.executable, str(AUDITION_SCRIPT), str(DEFAULT_CART_PATH)],
                    capture_output=True, text=True, encoding="utf-8", env=env,
                )
                if p.returncode == 0:
                    self.root.after(0, lambda: self.status_var.set("WAVs → audition/_studio_cartridge_*.wav"))
                else:
                    err = (p.stderr or p.stdout).strip().split("\n")[-1][:80]
                    self.root.after(0, lambda: self.status_var.set(f"audition failed: {err}"))
            except Exception as e:
                self.root.after(0, lambda: self.status_var.set(f"error: {e}"))
        threading.Thread(target=run, daemon=True).start()

    def _save_cartridge(self):
        path = filedialog.asksaveasfilename(
            initialdir=str(OUT_DIR), initialfile="my_cartridge.json",
            defaultextension=".json", filetypes=[("Cartridge JSON","*.json")])
        if not path: return
        cart = build_cartridge_json(self.state, name=Path(path).stem)
        Path(path).write_text(json.dumps(cart, indent=2))
        self.status_var.set(f"saved → {Path(path).name}")

    def _load_cartridge(self):
        path = filedialog.askopenfilename(
            initialdir=str(OUT_DIR),
            filetypes=[("Cartridge JSON","*.json")])
        if not path: return
        try:
            cart = json.loads(Path(path).read_text())
            prov = cart.get("provenance", {})
            if "stages_m0" in prov and "stages_m100" in prov:
                self.state["M0"]   = list(prov["stages_m0"])
                self.state["M100"] = list(prov["stages_m100"])
                self.selected_idx = None
                self._save_state(); self._refresh_all()
                self.status_var.set(f"loaded ← {Path(path).name}")
            else:
                messagebox.showwarning("Load",
                    "No studio provenance — only studio-saved cartridges editable here.")
        except Exception as e:
            messagebox.showerror("Load failed", str(e))


def main():
    root = tk.Tk()
    StudioApp(root)
    root.mainloop()


if __name__ == "__main__":
    main()
