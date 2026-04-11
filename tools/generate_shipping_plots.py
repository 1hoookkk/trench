"""Generate stage and frequency plots for current shipping winners.

Defaults:
- Winners source: vault/_scorecards/shipping_release_gate.json
- Candidate source: vault/_shipping_finalists/manifest.json
- Output dir: vault/_plots/shipping
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
import numpy as np
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
import sys
sys.path.insert(0, str(ROOT))

from pyruntime.body import Body
from pyruntime.constants import SR
from pyruntime.corner import CornerName
from pyruntime.freq_response import cascade_response_db, freq_points, stage_response

matplotlib.rcParams.update({
    "figure.facecolor": "#0a0a14",
    "axes.facecolor": "#0a0a14",
    "axes.edgecolor": "#333355",
    "axes.labelcolor": "#b5b5d6",
    "xtick.color": "#8c8cb0",
    "ytick.color": "#8c8cb0",
    "text.color": "#d7d7f5",
    "grid.color": "#232340",
    "font.family": "monospace",
    "font.size": 9,
})

BODY_KEY_BY_PUBLIC = {
    "Speaker Knockerz": "speaker_knockerz",
    "Aluminum Siding": "aluminum_siding",
    "Small Talk Ah-Ee": "small_talk",
    "Cul-De-Sac": "cul_de_sac",
}

CORNER_LABELS = [
    (CornerName.A, "M0_Q0"),
    (CornerName.B, "M0_Q100"),
    (CornerName.C, "M100_Q0"),
    (CornerName.D, "M100_Q100"),
]


def _slug(text: str) -> str:
    return "".join(ch.lower() if ch.isalnum() else "_" for ch in text).strip("_")


def _is_passthrough(enc) -> bool:
    return (
        abs(enc.c0 - 1.0) < 1e-6 and
        abs(enc.c1) < 1e-6 and
        abs(enc.c2) < 1e-6 and
        abs(enc.c3) < 1e-6 and
        abs(enc.c4) < 1e-6
    )


def _resolve_winners(gate_path: Path, manifest_path: Path) -> list[tuple[str, str]]:
    gate = json.loads(gate_path.read_text(encoding="utf-8"))
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))

    selected = gate.get("selected", {})
    out: list[tuple[str, str]] = []
    for public_name, body_key in BODY_KEY_BY_PUBLIC.items():
        row = selected.get(public_name)
        if not row:
            raise KeyError(f"Missing selected body row for {public_name} in {gate_path}")
        winner_name = row["design_candidate"]
        candidates = manifest.get(body_key, [])
        match = next((c for c in candidates if c.get("name") == winner_name), None)
        if match is None:
            raise KeyError(f"Could not resolve winner '{winner_name}' in manifest for {public_name}")
        out.append((public_name, match["path"]))
    return out


def _plot_frequency(body: Body, out_path: Path) -> None:
    freqs = freq_points(n=1024, sr=SR)

    fig, axes = plt.subplots(2, 2, figsize=(14, 9))
    fig.suptitle(f"{body.name} · Frequency Sweeps", fontsize=14, fontweight="bold", color="#6fa4ff")

    morph_positions = [0.0, 0.17, 0.33, 0.5, 0.67, 0.83, 1.0]
    q_positions = [0.0, 0.25, 0.5, 0.75, 1.0]

    ax = axes[0, 0]
    for i, m in enumerate(morph_positions):
        db = cascade_response_db(body.corners.interpolate(m, 0.0), freqs, SR)
        ax.semilogx(freqs, db, linewidth=1.2, alpha=0.9, label=f"M={m:.0%}",
                    color=plt.cm.cool(i / (len(morph_positions) - 1)))
    ax.set_title("Morph sweep @ Q=0%")

    ax = axes[0, 1]
    for i, m in enumerate(morph_positions):
        db = cascade_response_db(body.corners.interpolate(m, 1.0), freqs, SR)
        ax.semilogx(freqs, db, linewidth=1.2, alpha=0.9, label=f"M={m:.0%}",
                    color=plt.cm.winter(i / (len(morph_positions) - 1)))
    ax.set_title("Morph sweep @ Q=100%")

    ax = axes[1, 0]
    for i, q in enumerate(q_positions):
        db = cascade_response_db(body.corners.interpolate(0.0, q), freqs, SR)
        ax.semilogx(freqs, db, linewidth=1.2, alpha=0.9, label=f"Q={q:.0%}",
                    color=plt.cm.autumn(i / (len(q_positions) - 1)))
    ax.set_title("Q sweep @ M=0%")

    ax = axes[1, 1]
    diag = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
    for i, d in enumerate(diag):
        db = cascade_response_db(body.corners.interpolate(d, d), freqs, SR)
        ax.semilogx(freqs, db, linewidth=1.2, alpha=0.9, label=f"M=Q={d:.0%}",
                    color=plt.cm.spring(i / (len(diag) - 1)))
    ax.set_title("Diagonal sweep @ M=Q")

    for ax in axes.flatten():
        ax.set_xlim(20, (SR / 2.0) - 1.0)
        ax.set_ylim(-70, 50)
        ax.set_xlabel("Frequency (Hz)")
        ax.set_ylabel("Magnitude (dB)")
        ax.grid(True, which="both", linewidth=0.6, alpha=0.6)
        ax.legend(fontsize=7, framealpha=0.25, loc="upper right")

    fig.tight_layout()
    fig.savefig(out_path, dpi=170, bbox_inches="tight")
    plt.close(fig)


def _plot_stage(body: Body, out_path: Path) -> None:
    freqs = freq_points(n=512, sr=SR)
    # Stage curves at center point.
    center_enc = body.corners.interpolate(0.5, 0.5)

    active_indices = [i for i, s in enumerate(center_enc) if not _is_passthrough(s)]
    if not active_indices:
        active_indices = list(range(len(center_enc)))

    fig, axes = plt.subplots(1, 2, figsize=(15, 6))
    fig.suptitle(f"{body.name} · Stage Analysis", fontsize=14, fontweight="bold", color="#6fa4ff")

    ax = axes[0]
    for i, idx in enumerate(active_indices):
        h = stage_response(center_enc[idx], freqs, SR)
        db = 20.0 * np.log10(np.maximum(np.abs(h), 1e-20))
        ax.semilogx(freqs, db, linewidth=1.0, alpha=0.85,
                    color=plt.cm.plasma(i / max(1, len(active_indices) - 1)),
                    label=f"S{idx + 1}")
    cascade_db = cascade_response_db(center_enc, freqs, SR)
    ax.semilogx(freqs, cascade_db, color="#ffffff", linewidth=2.1, alpha=0.95, label="Cascade")
    ax.set_title("Per-stage responses @ M=50%, Q=50%")
    ax.set_xlim(20, (SR / 2.0) - 1.0)
    ax.set_ylim(-45, 45)
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Magnitude (dB)")
    ax.grid(True, which="both", linewidth=0.6, alpha=0.6)
    ax.legend(fontsize=7, ncol=2, framealpha=0.25, loc="upper right")

    # Corner stage peak maps.
    peak_freq = np.full((len(center_enc), len(CORNER_LABELS)), np.nan, dtype=float)
    peak_db = np.full((len(center_enc), len(CORNER_LABELS)), np.nan, dtype=float)
    for cidx, (cn, _) in enumerate(CORNER_LABELS):
        enc = body.corners.corner(cn).encode()
        for sidx, stage in enumerate(enc):
            if _is_passthrough(stage):
                continue
            h = stage_response(stage, freqs, SR)
            db = 20.0 * np.log10(np.maximum(np.abs(h), 1e-20))
            p = int(np.argmax(db))
            peak_freq[sidx, cidx] = freqs[p]
            peak_db[sidx, cidx] = db[p]

    ax = axes[1]
    display = np.log10(np.maximum(peak_freq, 20.0))
    im = ax.imshow(display, aspect="auto", cmap="viridis", interpolation="nearest")
    ax.set_title("Stage peak frequency map (corners)")
    ax.set_xlabel("Corner")
    ax.set_ylabel("Stage")
    ax.set_xticks(range(len(CORNER_LABELS)))
    ax.set_xticklabels([label for _, label in CORNER_LABELS], rotation=20)
    ax.set_yticks(range(len(center_enc)))
    ax.set_yticklabels([f"S{i + 1}" for i in range(len(center_enc))])
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("log10(Hz)")

    for sidx in range(len(center_enc)):
        for cidx in range(len(CORNER_LABELS)):
            if np.isnan(peak_freq[sidx, cidx]):
                txt = "PT"
            else:
                txt = f"{int(round(peak_freq[sidx, cidx]))}\n{peak_db[sidx, cidx]:.1f}dB"
            ax.text(cidx, sidx, txt, ha="center", va="center", fontsize=6, color="#e7e7ff")

    fig.tight_layout()
    fig.savefig(out_path, dpi=170, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate shipping stage + frequency plots.")
    parser.add_argument(
        "--gate",
        type=Path,
        default=ROOT / "vault" / "_scorecards" / "shipping_release_gate.json",
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=ROOT / "vault" / "_shipping_finalists" / "manifest.json",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=ROOT / "vault" / "_plots" / "shipping",
    )
    args = parser.parse_args()

    if not args.gate.exists():
        raise FileNotFoundError(f"Gate file not found: {args.gate}")
    if not args.manifest.exists():
        raise FileNotFoundError(f"Manifest file not found: {args.manifest}")

    args.out_dir.mkdir(parents=True, exist_ok=True)

    winners = _resolve_winners(args.gate, args.manifest)
    index: list[dict] = []
    for public_name, body_path in winners:
        body = Body.from_json(body_path)
        slug = _slug(public_name)
        freq_out = args.out_dir / f"{slug}_frequency.png"
        stage_out = args.out_dir / f"{slug}_stage.png"
        _plot_frequency(body, freq_out)
        _plot_stage(body, stage_out)
        index.append({
            "body": public_name,
            "body_path": body_path,
            "frequency_plot": str(freq_out),
            "stage_plot": str(stage_out),
        })
        print(f"{public_name}\n  frequency={freq_out}\n  stage={stage_out}")

    (args.out_dir / "index.json").write_text(json.dumps(index, indent=2), encoding="utf-8")
    print(args.out_dir / "index.json")


if __name__ == "__main__":
    main()
