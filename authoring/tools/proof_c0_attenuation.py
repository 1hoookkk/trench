"""forge/proof_c0_attenuation.py — visual proof.

Claim under test: in a 6-stage cascade with heritage-shaped poles and zeros,
c0 attenuation is the load-bearing variable. Strict Rule 5 (c0 = 1.0) destroys
formant character regardless of pole/zero placement.

Method: take two heritage M0_Q0 frames verbatim (Talking_Hedz, 6 active;
Ooh_to_Eee_(approx), 3 active). Render each twice through the SAME composite_db
renderer the rest of the forge uses:

    LEFT column   — heritage c0 (silicon truth)
    RIGHT column  — same b1/b2/c3/c4, but c0 forced to 1.0 (Rule 5 strict)

If LEFT shows multi-formant peaks and RIGHT shows shelves, the proof holds.
Same poles. Same zeros. Only c0 differs.

Output: cartridges/factory/generated/qlaw/_proof_c0_attenuation.png
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
OUT_PATH = REPO / "cartridges" / "factory" / "generated" / "qlaw" / "_proof_c0_attenuation.png"
AUTHORING_SR = 39062.5


def heritage_stage_to_coefs(s: dict, sr: float = AUTHORING_SR) -> dict[str, float]:
    """Match render_reference.py exactly."""
    radius = float(s.get("radius", 0.0))
    pole_hz = float(s.get("pole_freq_hz", 0.0))
    if radius == 0.0 or pole_hz == 0.0:
        return {"c0": 1.0, "c1": 0.0, "c2": 0.0, "c3": 0.0, "c4": 0.0}
    theta = 2.0 * math.pi * pole_hz / sr
    a1 = -2.0 * radius * math.cos(theta)
    z = s["zeros"]
    return {
        "c0": 1.0 + float(s.get("val1", 0.0)),
        "c1": float(z.get("b1", 0.0)),
        "c2": float(z.get("b2", 0.0)),
        "c3": a1,
        "c4": radius * radius,
    }


def composite_db(stages: list[dict], freqs: np.ndarray, sr: float = AUTHORING_SR,
                 boost: float = 1.0, c0_override: float | None = None) -> np.ndarray:
    """Identical math to render_reference.composite_db.

    If c0_override is set, replaces c0 on every active stage. Used to isolate
    Rule 5 enforcement from heritage attenuation.
    """
    w = 2.0 * np.pi * freqs / sr
    e1 = np.exp(-1j * w)
    e2 = np.exp(-2j * w)
    mag = np.ones_like(freqs)
    for s in stages:
        c0 = float(s["c0"]) if c0_override is None else float(c0_override)
        c1 = float(s["c1"]); c2 = float(s["c2"])
        c3 = float(s["c3"]); c4 = float(s["c4"])
        if c3 == 0.0 and c4 == 0.0:
            continue
        num = c0 + c1 * e1 + c2 * e2
        den = 1.0 + c3 * e1 + c4 * e2
        mag *= np.abs(num / den)
    return 20.0 * np.log10(np.maximum(mag * boost, 1e-9))


def load_heritage_m0(name: str) -> tuple[list[dict], float, float]:
    fp = CALIB_DIR / f"{name}.json"
    d = json.loads(fp.read_text(encoding="utf-8"))
    sr = float(d.get("sample_rate", AUTHORING_SR))
    boost = float(d.get("boost", 1.0))
    stages_in = d["corners"]["M0_Q0"]["stages"]
    stages_out = [heritage_stage_to_coefs(s, sr) for s in stages_in]
    while len(stages_out) < 12:
        stages_out.append({"c0": 1.0, "c1": 0.0, "c2": 0.0, "c3": 0.0, "c4": 0.0})
    return stages_out, sr, boost


def main() -> int:
    cases = [
        ("Talking_Hedz",        "c0 = 0.5619 (heritage)",  "c0 = 1.0 (Rule 5 strict)"),
        ("Ooh_to_Eee_(approx)", "c0 = 0.0146 (heritage)",  "c0 = 1.0 (Rule 5 strict)"),
    ]

    freqs = np.geomspace(40.0, 18000.0, 2048)

    fig, axes = plt.subplots(2, 2, figsize=(14.0, 8.0), dpi=140, sharex=True)
    fig.patch.set_facecolor("#16191a")
    fig.suptitle(
        "PATH OF TRUTH  ::  c0 attenuation is the load-bearing variable\n"
        "same poles, same zeros — only c0 differs",
        color="#ffba00", fontweight="bold", fontsize=14,
    )

    for row, (name, left_label, right_label) in enumerate(cases):
        stages, sr, boost = load_heritage_m0(name)

        # heritage c0 (left)
        db_heritage = composite_db(stages, freqs, sr=sr, boost=boost)
        # rule-5 strict (right)
        db_rule5 = composite_db(stages, freqs, sr=sr, boost=boost, c0_override=1.0)

        # heritage c0 actually used (for label)
        active_c0 = [s["c0"] for s in stages if not (s["c3"] == 0.0 and s["c4"] == 0.0)]
        n_active = len(active_c0)

        ax_l = axes[row, 0]
        ax_l.set_facecolor("#0d0f10")
        ax_l.semilogx(freqs, db_heritage, color="#ffba00", lw=2.2, label="heritage")
        ax_l.set_title(f"{name}  —  {left_label}  —  {n_active} active stages",
                       color="#ffba00", fontsize=10, fontweight="bold")

        ax_r = axes[row, 1]
        ax_r.set_facecolor("#0d0f10")
        ax_r.semilogx(freqs, db_rule5, color="#ff5555", lw=2.2, label="rule 5 strict")
        ax_r.set_title(f"{name}  —  {right_label}",
                       color="#ff5555", fontsize=10, fontweight="bold")

        # share y-axis range across the row so the visual comparison is honest
        all_vals = np.concatenate([db_heritage, db_rule5])
        y_lo, y_hi = float(np.min(all_vals)) - 5.0, float(np.max(all_vals)) + 5.0
        for ax in (ax_l, ax_r):
            ax.set_ylim(y_lo, y_hi)
            ax.grid(True, which="both", alpha=0.18, color="white")
            ax.tick_params(colors="white", labelsize=8)
            for spine in ax.spines.values():
                spine.set_color("#3a3e42")
            ax.set_ylabel("dB", color="white", fontsize=8)

    for ax in axes[-1, :]:
        ax.set_xlabel("Hz", color="white", fontsize=8)

    fig.tight_layout()
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PATH, facecolor=fig.get_facecolor())
    plt.close(fig)
    print(OUT_PATH)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
