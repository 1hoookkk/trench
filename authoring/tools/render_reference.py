"""forge/render_reference.py — heritage calibration vs forge output, same renderer.

Walks reference/calibration/*.json (the silicon-extracted heritage cartridges:
Talking_Hedz, Ear_Bender, Razor_Blades, BassBox_303, Freak_Shifta, Fuzzi_Face,
Meaty_Gizmo, Millennium, Radio_Craze, Early_Rizer, Ooh_to_Eee_(approx)),
converts each from the heritage (a1, r, val1, val2, val3) format into our
compiled-v1 (c0..c4) form via the Rosetta Stone, then plots all of them in a
single composite image alongside our most recent forge output.

Output: cartridges/factory/generated/qlaw/_heritage_vs_forge.png

Use this to eyeball where the synthesized forge sits relative to the
silicon-extracted ground truth. If our shape doesn't match what's in
reference/calibration/, our doctrine is off.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


REPO = Path(__file__).resolve().parent.parent
CALIB_DIR = REPO / "reference" / "calibration"
FORGE_OUT_DIR = REPO / "cartridges" / "factory" / "generated" / "qlaw"
GENERATED = FORGE_OUT_DIR / "ah_to_ee_pathA.json"
OUT_PATH = FORGE_OUT_DIR / "_heritage_vs_forge.png"

AUTHORING_SR_DEFAULT = 39062.5


def _heritage_stage_to_compiled(s: dict, sr: float = AUTHORING_SR_DEFAULT) -> dict[str, float]:
    """Convert one heritage stage to compiled-v1 c0..c4.

    Two heritage formats exist in the calibration tree:
      A) reference/calibration/*.json — verbose RE format with `radius`,
         `pole_freq_hz`, `val1/2/3`, and a `zeros: {b1, b2, ...}` sub-dict.
      B) vault/*.json — compact format with `a1, r, val1, val2, val3`.
    """
    if "zeros" in s and "radius" in s:
        # Format A — verbose calibration RE
        radius = float(s.get("radius", 0.0))
        pole_hz = float(s.get("pole_freq_hz", 0.0))
        if radius == 0.0 or pole_hz == 0.0:
            return {"c0": 1.0, "c1": 0.0, "c2": 0.0, "c3": 0.0, "c4": 0.0}
        theta = 2.0 * math.pi * pole_hz / sr
        a1 = -2.0 * radius * math.cos(theta)
        z = s["zeros"]
        b1 = float(z.get("b1", 0.0))
        b2 = float(z.get("b2", 0.0))
        val1 = float(s.get("val1", 0.0))
        return {
            "c0": 1.0 + val1,
            "c1": b1,
            "c2": b2,
            "c3": a1,
            "c4": radius * radius,
        }

    # Format B — compact vault-style
    a1 = float(s.get("a1", 0.0))
    r = float(s.get("r", 0.0))
    val1 = float(s.get("val1", 0.0))
    val2 = float(s.get("val2", 0.0))
    val3 = float(s.get("val3", 0.0))
    if r == 0.0:
        return {"c0": 1.0, "c1": 0.0, "c2": 0.0, "c3": 0.0, "c4": 0.0}
    r2 = r * r
    return {
        "c0": 1.0 + val1,
        "c1": a1 + val2,
        "c2": r2 - val3,
        "c3": a1,
        "c4": r2,
    }


def heritage_to_compiled(ref: dict) -> dict | None:
    if "corners" not in ref:
        return None
    sample_rate = float(ref.get("sample_rate", AUTHORING_SR_DEFAULT))
    keyframes = []
    for label, corner in ref["corners"].items():
        stages_in = corner.get("stages", [])
        stages_out = [_heritage_stage_to_compiled(s, sample_rate) for s in stages_in]
        while len(stages_out) < 12:
            stages_out.append({"c0": 1.0, "c1": 0.0, "c2": 0.0, "c3": 0.0, "c4": 0.0})
        morph = 0.0 if "M0" in label else 1.0
        q = 0.0 if "Q0" in label else 1.0
        keyframes.append({
            "label": label,
            "morph": corner.get("morph", morph),
            "q": corner.get("q", q),
            "boost": float(corner.get("boost", ref.get("boost", 1.0))),
            "stages": stages_out,
        })
    return {
        "format": "compiled-v1",
        "name": ref.get("name", "unknown"),
        "sampleRate": sample_rate,
        "keyframes": keyframes,
    }


def composite_db(cart: dict, label: str, freqs: np.ndarray) -> np.ndarray | None:
    kf = next((k for k in cart["keyframes"] if k["label"] == label), None)
    if kf is None:
        return None
    sr = float(cart.get("sampleRate", AUTHORING_SR_DEFAULT))
    boost = float(kf.get("boost", 1.0))
    w = 2.0 * np.pi * freqs / sr
    e1 = np.exp(-1j * w)
    e2 = np.exp(-2j * w)
    mag = np.ones_like(freqs)
    for s in kf["stages"]:
        c0 = float(s["c0"]); c1 = float(s["c1"]); c2 = float(s["c2"])
        c3 = float(s["c3"]); c4 = float(s["c4"])
        if c3 == 0.0 and c4 == 0.0:
            continue
        num = c0 + c1 * e1 + c2 * e2
        den = 1.0 + c3 * e1 + c4 * e2
        mag *= np.abs(num / den)
    return 20.0 * np.log10(np.maximum(mag * boost, 1e-9))


def main() -> int:
    panels: list[tuple[str, dict, str]] = []  # (display_name, cartridge, accent)

    # 1) the user's most recent forge output, if it exists
    if GENERATED.exists():
        gen = json.loads(GENERATED.read_text(encoding="utf-8"))
        panels.append((f"FORGE → {gen.get('name','?')}", gen, "ours"))

    # 2) every heritage reference
    for jf in sorted(CALIB_DIR.glob("*.json")):
        try:
            data = json.loads(jf.read_text(encoding="utf-8"))
        except Exception:
            continue
        compiled = heritage_to_compiled(data)
        if compiled is None:
            continue
        panels.append((f"HERITAGE :: {jf.stem}", compiled, "ref"))

    if not panels:
        print("nothing to render")
        return 1

    n = len(panels)
    cols = 3
    rows = (n + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(16.0, 3.2 * rows), dpi=140, sharex=True)
    fig.patch.set_facecolor("#16191a")
    fig.suptitle(
        "TRENCH FORGE   ::   heritage calibration  vs.  forge output\n"
        "M0_Q0 (bold amber)   M100_Q0 (cyan)",
        color="#ffba00", fontweight="bold", fontsize=13,
    )

    freqs = np.geomspace(40.0, 18000.0, 1024)
    flat = axes.flatten() if hasattr(axes, "flatten") else [axes]

    for ax, (name, cart, kind) in zip(flat, panels):
        ax.set_facecolor("#0d0f10")
        m0 = composite_db(cart, "M0_Q0", freqs)
        m100 = composite_db(cart, "M100_Q0", freqs)
        if m0 is not None:
            ax.semilogx(freqs, m0, color="#ffba00", lw=2.0, label="M0")
        if m100 is not None:
            ax.semilogx(freqs, m100, color="#22ddff", lw=1.0, alpha=0.85, label="M100")

        title_color = "#ff5555" if kind == "ours" else "#cccccc"
        ax.set_title(name, color=title_color, fontsize=10, fontweight="bold")
        ax.grid(True, which="both", alpha=0.18, color="white")
        ax.tick_params(colors="white", labelsize=8)
        for spine in ax.spines.values():
            spine.set_color("#3a3e42")
        ax.set_ylabel("dB", color="white", fontsize=8)
        ax.legend(fontsize=7, loc="lower left", framealpha=0.4,
                  facecolor="#22262a", edgecolor="#444444", labelcolor="white")

    for j in range(len(panels), len(flat)):
        flat[j].axis("off")

    for ax in flat[-cols:]:
        ax.set_xlabel("Hz", color="white", fontsize=8)

    fig.tight_layout()
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PATH, facecolor=fig.get_facecolor())
    plt.close(fig)
    print(OUT_PATH)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
