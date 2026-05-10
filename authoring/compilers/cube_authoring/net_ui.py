"""Flat pill-net / cube-net UI — interactive spectrum tiles with gate feedback.

Pill mode shows a 2×2 grid (morph × Q). Cube mode shows a 2×4 grid (XY
plane at Z=0 left, Z=1 right). Gate-pass edges are green; failures are red.
Click any tile to audition that corner/keyframe through the selected probe.
Press 'r' to reroll all tiles.

Usage:
  python tools/cube_authoring/net_ui.py --pill --scope hf_resonance
  python tools/cube_authoring/net_ui.py --cube --scope hf_resonance
  python tools/cube_authoring/net_ui.py --input cartridges/pills/x.pill.json
"""
from __future__ import annotations

import argparse
import json
import os
import struct
import subprocess
import sys
import tempfile
import wave
from pathlib import Path

import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.widgets import Button, RadioButtons
import numpy as np

TOOLS_DIR = Path(__file__).parent
REPO_DIR   = TOOLS_DIR.parent.parent

SR          = 44100
IMPULSE_LEN = 512

PROBE_IDS = ["sine_sweep", "reese", "drum_loop", "noise"]

# ── Pill layout ────────────────────────────────────────────────────────────────
# Row 0 = high-Q, Row 1 = low-Q. Col 0 = lo-morph, Col 1 = hi-morph.
PILL_LAYOUT: dict[str, tuple[int, int]] = {
    "M0_Q100":    (0, 0),
    "M100_Q100":  (0, 1),
    "M0_Q0":      (1, 0),
    "M100_Q0":    (1, 1),
}
PILL_EDGES = [
    ("M0_Q0",    "M100_Q0",    "morph"),
    ("M0_Q100",  "M100_Q100",  "morph"),
    ("M0_Q0",    "M0_Q100",    "q"),
    ("M100_Q0",  "M100_Q100",  "q"),
]

# ── Cube layout ────────────────────────────────────────────────────────────────
# Cols 0-1 = Z=0 face; cols 2-3 = Z=1 face. Row 0 = Y=1, Row 1 = Y=0.
CUBE_LAYOUT: dict[str, tuple[int, int]] = {
    "c010": (0, 0), "c110": (0, 1), "c011": (0, 2), "c111": (0, 3),
    "c000": (1, 0), "c100": (1, 1), "c001": (1, 2), "c101": (1, 3),
}
CUBE_EDGES = [
    ("c000", "c100", "x"), ("c010", "c110", "x"),
    ("c001", "c101", "x"), ("c011", "c111", "x"),
    ("c000", "c010", "y"), ("c100", "c110", "y"),
    ("c001", "c011", "y"), ("c101", "c111", "y"),
    ("c000", "c001", "z"), ("c100", "c101", "z"),
    ("c010", "c011", "z"), ("c110", "c111", "z"),
]


# ── DSP helpers ───────────────────────────────────────────────────────────────

def _f32(x: float) -> float:
    return struct.unpack("<f", struct.pack("<f", x))[0]


def _stage_tuples(stages: list[dict]) -> list[tuple]:
    return [(s["c0"], s["c1"], s["c2"], s["c3"], s["c4"]) for s in stages]


def _run_filter(signal: np.ndarray, stages: list[tuple], boost: float) -> np.ndarray:
    state = [(0.0, 0.0)] * len(stages)
    out = np.empty(len(signal), dtype=np.float32)
    for n, x in enumerate(signal):
        sig = float(x)
        new_state = []
        for c0, c1, c2, c3, c4 in stages:
            w1, w2 = state[len(new_state)]
            y = c0 * sig + w1
            new_state.append((c1 * sig - c3 * y + w2, c2 * sig - c4 * y))
            sig = y
        state = new_state
        out[n] = _f32(sig * boost)
    return out


def _spectrum_db(stages: list[dict], boost: float = 1.0) -> tuple[np.ndarray, np.ndarray]:
    impulse = np.zeros(IMPULSE_LEN, dtype=np.float32)
    impulse[0] = 1.0
    ir = _run_filter(impulse, _stage_tuples(stages), boost)
    spec = np.fft.rfft(ir, n=4096)
    db   = 20.0 * np.log10(np.maximum(np.abs(spec), 1e-12))
    freqs = np.fft.rfftfreq(4096, d=1.0 / SR)
    return freqs, db


# ── Subprocess helpers ────────────────────────────────────────────────────────

def _run(cmd: list[str]) -> str:
    r = subprocess.run([sys.executable] + cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"subprocess failed:\n{r.stderr}")
    return r.stdout


def _reroll_pill(scope: str, seed: int | None) -> dict:
    cmd = [str(TOOLS_DIR / "reroll.py"), "--scope", scope,
           "--pill", "--pill-id", f"pill.{scope}.ui"]
    if seed is not None:
        cmd += ["--seed", str(seed)]
    return json.loads(_run(cmd))


def _reroll_cube(scope: str, seed: int | None) -> dict:
    cmd = [str(TOOLS_DIR / "reroll.py"), "--scope", scope,
           "--cube", "--cube-id", f"cube.{scope}.ui"]
    if seed is not None:
        cmd += ["--seed", str(seed)]
    return json.loads(_run(cmd))


def _compile_pill(auth_doc: dict) -> dict:
    with tempfile.NamedTemporaryFile(suffix=".json", mode="w", delete=False) as f:
        json.dump(auth_doc, f)
        fname = f.name
    out = tempfile.mktemp(suffix=".json")
    subprocess.run(
        [sys.executable, str(REPO_DIR / "tools" / "compile_pill.py"),
         fname, "--out", out],
        check=True, capture_output=True,
    )
    result = json.loads(Path(out).read_text())
    os.unlink(fname)
    os.unlink(out)
    return result


def _compile_cube(auth_doc: dict) -> dict:
    with tempfile.NamedTemporaryFile(suffix=".json", mode="w", delete=False) as f:
        json.dump(auth_doc, f)
        fname = f.name
    out = tempfile.mktemp(suffix=".json")
    subprocess.run(
        [sys.executable, str(REPO_DIR / "tools" / "compile_cube.py"),
         fname, "--out", out],
        check=True, capture_output=True,
    )
    result = json.loads(Path(out).read_text())
    os.unlink(fname)
    os.unlink(out)
    return result


def _gate_pill(compiled: dict) -> dict | None:
    with tempfile.NamedTemporaryFile(suffix=".json", mode="w", delete=False) as f:
        json.dump(compiled, f)
        fname = f.name
    r = subprocess.run(
        [sys.executable, str(TOOLS_DIR / "pill_gate.py"), fname],
        capture_output=True, text=True,
    )
    os.unlink(fname)
    if r.returncode not in (0, 1):
        return None
    return json.loads(r.stdout)


def _gate_cube(compiled: dict) -> dict | None:
    with tempfile.NamedTemporaryFile(suffix=".json", mode="w", delete=False) as f:
        json.dump(compiled, f)
        fname = f.name
    r = subprocess.run(
        [sys.executable, str(TOOLS_DIR / "edge_gate.py"), fname],
        capture_output=True, text=True,
    )
    os.unlink(fname)
    if r.returncode not in (0, 1):
        return None
    return json.loads(r.stdout)


# ── Edge verdict helpers ──────────────────────────────────────────────────────

def _pill_edge_passed(verdict: dict | None, la: str, lb: str) -> bool | None:
    """True=pass, False=fail, None=unknown."""
    if verdict is None:
        return None
    keys = {f"{la}-{lb}", f"{lb}-{la}"}
    result = True
    for probe in verdict.get("per_probe", []):
        for edge in probe.get("edges", []):
            if edge["edge"] in keys:
                if not edge["passed"]:
                    return False
    return result


def _cube_edge_passed(verdict: dict | None, la: str, lb: str) -> bool | None:
    if verdict is None:
        return None
    keys = {f"{la}-{lb}", f"{lb}-{la}"}
    for entry in verdict.get("edges", []):
        if entry.get("edge") in keys:
            return entry.get("passed", None)
    return None


# ── Net UI ────────────────────────────────────────────────────────────────────

class NetUI:
    def __init__(self, mode: str, scope: str, seed: int | None,
                 input_file: str | None) -> None:
        self.mode     = mode
        self.scope    = scope
        self.seed     = seed
        self.compiled: dict | None = None
        self.verdict:  dict | None = None
        self.probe     = "sine_sweep"
        self._tmp_wav: str | None = None
        self._setup_figure()
        self._load(input_file)
        self._draw()
        plt.show()

    # ── Figure setup ──────────────────────────────────────────────────────────

    def _setup_figure(self) -> None:
        layout = PILL_LAYOUT if self.mode == "pill" else CUBE_LAYOUT
        nrows  = 2
        ncols  = 2 if self.mode == "pill" else 4
        fw     = 10 if self.mode == "pill" else 16
        self.layout = layout
        self.edges  = PILL_EDGES if self.mode == "pill" else CUBE_EDGES

        self.fig = plt.figure(figsize=(fw, 7), facecolor="#1a1a1a")
        gs = self.fig.add_gridspec(
            nrows, ncols, top=0.90, bottom=0.22, hspace=0.04, wspace=0.04,
        )
        self.axes: dict[str, plt.Axes] = {}
        for label, (r, c) in layout.items():
            ax = self.fig.add_subplot(gs[r, c])
            ax.set_facecolor("#111111")
            for sp in ax.spines.values():
                sp.set_color("#444444")
                sp.set_linewidth(1)
            ax.tick_params(left=False, bottom=False,
                           labelleft=False, labelbottom=False)
            self.axes[label] = ax

        # Probe selector
        ax_radio = self.fig.add_axes([0.02, 0.02, 0.18, 0.16],
                                      facecolor="#1a1a1a")
        self.radio = RadioButtons(ax_radio, PROBE_IDS, active=0,
                                   activecolor="#44ff88")
        for lbl in self.radio.labels:
            lbl.set_color("white")
            lbl.set_fontsize(8)
        self.radio.on_clicked(self._on_probe)

        # Reroll all button
        ax_btn = self.fig.add_axes([0.82, 0.02, 0.15, 0.08])
        self.btn_reroll = Button(ax_btn, "↺  Reroll All",
                                  color="#333333", hovercolor="#555555")
        self.btn_reroll.label.set_color("white")
        self.btn_reroll.on_clicked(self._on_reroll_all)

        # Status
        self.status_txt = self.fig.text(
            0.28, 0.08, "", color="white", fontsize=11,
            fontweight="bold", va="center", ha="left",
        )

        self.fig.suptitle(
            f"TRENCH {self.mode.upper()} NET — {self.scope}",
            color="white", fontsize=13,
        )
        self.fig.canvas.mpl_connect("button_press_event", self._on_click)
        self.fig.canvas.mpl_connect("key_press_event", self._on_key)

    # ── Load / gate cycle ────────────────────────────────────────────────────

    def _load(self, input_file: str | None = None) -> None:
        self.status_txt.set_text("loading…")
        self.fig.canvas.draw_idle()
        try:
            if input_file:
                auth_doc = json.loads(Path(input_file).read_text())
                if self.mode == "pill":
                    self.compiled = _compile_pill(auth_doc)
                    self.verdict  = _gate_pill(self.compiled)
                else:
                    self.compiled = _compile_cube(auth_doc)
                    self.verdict  = _gate_cube(self.compiled)
            elif self.mode == "pill":
                auth_doc      = _reroll_pill(self.scope, self.seed)
                self.compiled = _compile_pill(auth_doc)
                self.verdict  = _gate_pill(self.compiled)
            else:
                auth_doc      = _reroll_cube(self.scope, self.seed)
                self.compiled = _compile_cube(auth_doc)
                self.verdict  = _gate_cube(self.compiled)
        except Exception as exc:
            self.status_txt.set_text(f"error: {exc}")
            return
        self._update_status()

    def _update_status(self) -> None:
        if self.verdict is None:
            self.status_txt.set_text("verdict unavailable")
            self.status_txt.set_color("white")
            return
        admitted = self.verdict.get("admitted", False)
        self.status_txt.set_text("ADMITTED" if admitted else "REJECTED")
        self.status_txt.set_color("#44ff88" if admitted else "#ff4444")

    # ── Drawing ───────────────────────────────────────────────────────────────

    def _draw(self) -> None:
        if self.compiled is None:
            return
        tiles = self._tile_map()
        for label, ax in self.axes.items():
            tile = tiles.get(label)
            ax.cla()
            ax.set_facecolor("#111111")
            if tile:
                freqs, db = _spectrum_db(tile["stages"], tile.get("boost", 1.0))
                mask = (freqs >= 40) & (freqs <= 22000)
                ax.semilogx(freqs[mask], db[mask], color="#44aaff", lw=1.0)
            ax.set_xlim(40, 22000)
            ax.set_ylim(-60, 40)
            ax.set_title(label, color="#aaaaaa", fontsize=8, pad=2)
            for sp in ax.spines.values():
                sp.set_color("#444444")
                sp.set_linewidth(1)
            ax.tick_params(left=False, bottom=False,
                           labelleft=False, labelbottom=False)

        self._draw_edges()
        self.fig.canvas.draw_idle()

    def _tile_map(self) -> dict[str, dict]:
        if self.mode == "pill":
            return {kf["label"]: kf for kf in self.compiled["keyframes"]}
        else:
            return {k: v for k, v in self.compiled["corners"].items()}

    def _draw_edges(self) -> None:
        get_passed = (_pill_edge_passed if self.mode == "pill"
                      else _cube_edge_passed)
        for la, lb, _ in self.edges:
            passed = get_passed(self.verdict, la, lb)
            color  = ("#44ff88" if passed is True
                      else "#ff4444" if passed is False
                      else "#666666")
            width  = 3 if passed is not None else 1
            ra, ca = self.layout[la]
            rb, cb = self.layout[lb]
            ax_a = self.axes[la]
            ax_b = self.axes[lb]
            # Vertical neighbors (same col, adjacent rows)
            if ca == cb and abs(ra - rb) == 1:
                top_ax  = ax_a if ra < rb else ax_b
                bot_ax  = ax_b if ra < rb else ax_a
                top_ax.spines["bottom"].set_color(color)
                top_ax.spines["bottom"].set_linewidth(width)
                bot_ax.spines["top"].set_color(color)
                bot_ax.spines["top"].set_linewidth(width)
            # Horizontal neighbors (same row, adjacent cols)
            elif ra == rb and abs(ca - cb) == 1:
                left_ax  = ax_a if ca < cb else ax_b
                right_ax = ax_b if ca < cb else ax_a
                left_ax.spines["right"].set_color(color)
                left_ax.spines["right"].set_linewidth(width)
                right_ax.spines["left"].set_color(color)
                right_ax.spines["left"].set_linewidth(width)

    # ── Audition ──────────────────────────────────────────────────────────────

    def _audition(self, label: str) -> None:
        tiles = self._tile_map()
        tile  = tiles.get(label)
        if tile is None:
            return
        sys.path.insert(0, str(TOOLS_DIR))
        from probes import get_probe  # type: ignore
        probe_audio = get_probe(self.probe, sr=SR)
        filtered = _run_filter(probe_audio, _stage_tuples(tile["stages"]),
                               tile.get("boost", 1.0))
        rms = np.sqrt(np.mean(filtered ** 2))
        if rms > 1e-9:
            filtered = filtered * (0.2 / rms)
        filtered = np.clip(filtered, -1.0, 1.0)
        pcm = (filtered * 32767).astype(np.int16)

        if self._tmp_wav and os.path.exists(self._tmp_wav):
            try:
                os.unlink(self._tmp_wav)
            except OSError:
                pass
        tmp = tempfile.mktemp(suffix=".wav")
        with wave.open(tmp, "wb") as w:
            w.setnchannels(1)
            w.setsampwidth(2)
            w.setframerate(SR)
            w.writeframes(pcm.tobytes())
        self._tmp_wav = tmp
        import winsound
        winsound.PlaySound(tmp, winsound.SND_FILENAME | winsound.SND_ASYNC)

    # ── Event callbacks ───────────────────────────────────────────────────────

    def _on_probe(self, label: str) -> None:
        self.probe = label

    def _on_reroll_all(self, _event) -> None:
        self._load()
        self._draw()

    def _on_click(self, event) -> None:
        for label, ax in self.axes.items():
            if ax == event.inaxes:
                self._audition(label)
                break

    def _on_key(self, event) -> None:
        if event.key == "r":
            self._on_reroll_all(event)


# ── CLI ───────────────────────────────────────────────────────────────────────

def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(
        description="Interactive pill/cube net UI with gate feedback.",
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--pill", action="store_true",
                        help="Show 2×2 pill net (morph × Q).")
    group.add_argument("--cube", action="store_true",
                        help="Show 2×4 cube net (XY plane × Z).")
    parser.add_argument("--scope",  default="hf_resonance")
    parser.add_argument("--seed",   type=int, default=None)
    parser.add_argument("--input",  default=None,
                         help="Load existing pill/cube JSON (skip reroll).")
    args = parser.parse_args(argv)
    mode = "pill" if args.pill else "cube"
    NetUI(mode, args.scope, seed=args.seed, input_file=args.input)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
