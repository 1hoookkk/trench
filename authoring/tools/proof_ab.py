"""forge/proof_ab.py — direct A/B overlay of extreme_carve vs Talking_Hedz."""
from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


REPO = Path(__file__).resolve().parent.parent
EC_PATH = REPO / "cartridges" / "factory" / "generated" / "qlaw" / "_extreme_carve_proof.json"
TH_PATH = REPO / "reference" / "calibration" / "Talking_Hedz.json"
OUT_PATH = REPO / "cartridges" / "factory" / "generated" / "qlaw" / "_proof_ab_overlay.png"
SR = 39062.5


def composite_db(stages, freqs, sr=SR, boost=1.0):
    w = 2 * np.pi * freqs / sr
    e1 = np.exp(-1j * w); e2 = np.exp(-2j * w)
    mag = np.ones_like(freqs)
    for s in stages:
        c0, c1, c2, c3, c4 = s["c0"], s["c1"], s["c2"], s["c3"], s["c4"]
        if c3 == 0.0 and c4 == 0.0: continue
        num = c0 + c1*e1 + c2*e2
        den = 1.0 + c3*e1 + c4*e2
        mag *= np.abs(num/den)
    return 20*np.log10(np.maximum(mag * boost, 1e-12))


def load_compiled(path):
    return json.loads(path.read_text(encoding="utf-8"))


def heritage_to_compiled(path):
    d = json.loads(path.read_text(encoding="utf-8"))
    sr = float(d.get("sample_rate", SR))
    boost = float(d.get("boost", 1.0))
    keyframes = []
    for label, corner in d["corners"].items():
        stages_out = []
        for s in corner["stages"]:
            r = float(s.get("radius", 0))
            f = float(s.get("pole_freq_hz", 0))
            if r == 0 or f == 0:
                stages_out.append({"c0":1.,"c1":0.,"c2":0.,"c3":0.,"c4":0.}); continue
            a1 = -2*r*math.cos(2*math.pi*f/sr); a2 = r*r
            v1 = float(s.get("val1",0)); v2 = float(s.get("val2",0)); v3 = float(s.get("val3",0))
            stages_out.append({"c0":1.+v1,"c1":a1+v2,"c2":a2-v3,"c3":a1,"c4":a2})
        while len(stages_out) < 12:
            stages_out.append({"c0":1.,"c1":0.,"c2":0.,"c3":0.,"c4":0.})
        keyframes.append({"label":label,"boost":boost,"stages":stages_out})
    return {"keyframes":keyframes, "boost":boost, "sampleRate":sr}


def get_corner(cart, label):
    for kf in cart["keyframes"]:
        if kf["label"] == label: return kf
    return None


def main():
    ec = load_compiled(EC_PATH)
    th = heritage_to_compiled(TH_PATH)
    freqs = np.geomspace(40.0, 18000.0, 4096)

    fig, axes = plt.subplots(2, 1, figsize=(13.0, 8.0), dpi=140, sharex=True)
    fig.patch.set_facecolor("#16191a")
    fig.suptitle("EXTREME_CARVE  vs  TALKING_HEDZ   ::  same axes, same renderer",
                 color="#ffba00", fontweight="bold", fontsize=14)

    pairs = [("M0_Q0", axes[0], "M0"), ("M100_Q0", axes[1], "M100")]
    for label, ax, title in pairs:
        ec_kf = get_corner(ec, label)
        th_kf = get_corner(th, label)
        ec_db = composite_db(ec_kf["stages"], freqs, boost=ec_kf.get("boost",1.0))
        th_db = composite_db(th_kf["stages"], freqs, boost=th_kf.get("boost",1.0))

        ax.set_facecolor("#0d0f10")
        ax.semilogx(freqs, th_db, color="#ffba00", lw=2.0, label="Talking_Hedz (silicon)")
        ax.semilogx(freqs, ec_db, color="#ff5555", lw=2.2, label="extreme_carve (forge)")
        ax.set_title(title, color="white", fontsize=12, fontweight="bold")
        ax.grid(True, which="both", alpha=0.18, color="white")
        ax.tick_params(colors="white", labelsize=9)
        for spine in ax.spines.values(): spine.set_color("#3a3e42")
        ax.set_ylabel("dB", color="white")
        ax.legend(facecolor="#22262a", edgecolor="#444", labelcolor="white", fontsize=10)
        ax.axhline(0, color="#444", lw=0.5)

    axes[1].set_xlabel("Hz", color="white")
    fig.tight_layout()
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PATH, facecolor=fig.get_facecolor())
    plt.close(fig)
    print(OUT_PATH)


if __name__ == "__main__":
    main()
