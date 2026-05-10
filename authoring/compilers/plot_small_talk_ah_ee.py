"""Plot the 4-corner magnitude response of Small Talk Ah-Ee.

Evaluates at sr=44100 (per memory rule: compiled-v1 sampleRate field is
stale authoring metadata; magnitude math uses the actual execution rate
of 44100). Marks the adult-male formant targets the body was authored
for so the bake's structure is verifiable against intent.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parents[1]
CART = REPO / "cartridges" / "engine" / "character" / "small_talk_ah_ee.json"
OUT = REPO / "parity_plots" / "small_talk_ah_ee_response.png"

NFFT = 8192
SR_EVAL = 44100.0
DB_FLOOR = 1e-12
PASSTHROUGH = (1.0, 0.0, 0.0, 0.0, 0.0)

# Adult-male formant targets the body was authored for.
AH_FORMANTS = (730.0, 1090.0, 2440.0)   # Ah (/ɑ/) → expected at M0
EE_FORMANTS = (270.0, 2290.0, 3010.0)   # Ee (/i/) → expected at M100


def cascade_response(stages: list[dict], sr: float) -> tuple[np.ndarray, np.ndarray, int]:
    w = np.linspace(0.0, np.pi, NFFT)
    z_inv = np.exp(-1j * w)
    h = np.ones_like(w, dtype=complex)
    active = 0
    for s in stages:
        c0, c1, c2, c3, c4 = s["c0"], s["c1"], s["c2"], s["c3"], s["c4"]
        if (c0, c1, c2, c3, c4) == PASSTHROUGH:
            continue
        active += 1
        num = c0 + c1 * z_inv + c2 * (z_inv * z_inv)
        den = 1.0 + c3 * z_inv + c4 * (z_inv * z_inv)
        h *= num / den
    freqs = w * sr / (2 * np.pi)
    mag_db = 20.0 * np.log10(np.maximum(np.abs(h), DB_FLOOR))
    return freqs, mag_db, active


def main() -> int:
    cart = json.loads(CART.read_text())
    boost = cart["keyframes"][0].get("boost", 1.0)
    boost_db = 20.0 * np.log10(boost)

    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True, sharey=True)
    fig.suptitle(
        f'Small Talk Ah-Ee — magnitude response @ sr={SR_EVAL:.0f} Hz '
        f'(boost {boost:.1f}× = {boost_db:+.1f} dB applied)',
        fontsize=12,
    )

    layout = {
        "M0_Q0":     (0, 0, AH_FORMANTS, "Ah-like (M=0)"),
        "M0_Q100":   (0, 1, AH_FORMANTS, "Ah-like, Q=1"),
        "M100_Q0":   (1, 0, EE_FORMANTS, "Ee-like (M=1)"),
        "M100_Q100": (1, 1, EE_FORMANTS, "Ee-like, Q=1"),
    }

    for kf in cart["keyframes"]:
        label = kf["label"]
        if label not in layout:
            continue
        row, col, targets, title = layout[label]
        ax = axes[row][col]
        freqs, mag_db, active = cascade_response(kf["stages"], SR_EVAL)
        mag_db = mag_db + boost_db
        ax.semilogx(freqs[1:], mag_db[1:], color="#222222", linewidth=1.4)
        for f_hz in targets:
            ax.axvline(f_hz, color="#cc3333", linestyle="--", linewidth=0.8, alpha=0.6)
            ax.text(f_hz, ax.get_ylim()[1] - 3 if ax.get_ylim()[1] > -10 else -10,
                    f' {f_hz:.0f}', color="#cc3333", fontsize=8, va="top")
        ax.set_title(f'{label} — {title} ({active} active)', fontsize=10)
        ax.set_xlim(20, SR_EVAL / 2)
        ax.grid(True, which="both", alpha=0.3)
        ax.set_xlabel("Hz")
        ax.set_ylabel("dB")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(OUT, dpi=120)
    print(f"wrote {OUT}")

    # Summary: scan each corner for the strongest peaks.
    print("\nDominant peaks (top 4 per corner, prominence > 6 dB):")
    for kf in cart["keyframes"]:
        freqs, mag_db, _ = cascade_response(kf["stages"], SR_EVAL)
        peaks = []
        for i in range(2, len(mag_db) - 2):
            if mag_db[i] > mag_db[i - 1] and mag_db[i] > mag_db[i + 1]:
                local = mag_db[max(0, i - 50): i + 50]
                prom = mag_db[i] - np.min(local)
                if prom > 6 and freqs[i] > 50:
                    peaks.append((freqs[i], mag_db[i] + boost_db, prom))
        peaks.sort(key=lambda t: -t[2])
        head = peaks[:4]
        head_str = ", ".join(f"{f:.0f}Hz @ {m:+.1f}dB" for f, m, _ in head) or "(none)"
        print(f"  {kf['label']:10s} {head_str}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
