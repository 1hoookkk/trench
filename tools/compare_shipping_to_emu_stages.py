"""Compare shipping winners to selected E-MU P2k bodies stage-by-stage.

Outputs:
- One stage-overlay PNG per (shipping winner, E-MU body) pair
- JSON summary with per-stage mean absolute dB deltas
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
import sys
sys.path.insert(0, str(ROOT))

from pyruntime.body import Body
from pyruntime.constants import SR
from pyruntime.freq_response import freq_points, stage_response

matplotlib.rcParams.update({
    "figure.facecolor": "#0a0a14",
    "axes.facecolor": "#0a0a14",
    "axes.edgecolor": "#2f2f55",
    "axes.labelcolor": "#c0c0e0",
    "xtick.color": "#8f8fb3",
    "ytick.color": "#8f8fb3",
    "text.color": "#dfdffd",
    "grid.color": "#212140",
    "font.family": "monospace",
    "font.size": 9,
})

BODY_KEY_BY_PUBLIC = {
    "Speaker Knockerz": "speaker_knockerz",
    "Aluminum Siding": "aluminum_siding",
    "Small Talk Ah-Ee": "small_talk",
    "Cul-De-Sac": "cul_de_sac",
}

DEFAULT_EMU_INDICES = [0, 1, 2, 8, 9, 10, 13, 18, 25, 26]
DEFAULT_MAPPED_BY_BODY = {
    "Speaker Knockerz": 6,
    "Aluminum Siding": 26,
    "Small Talk Ah-Ee": 13,
    "Cul-De-Sac": 10,
}


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
            raise KeyError(f"Missing selected row for {public_name}")
        winner_name = row["design_candidate"]
        candidates = manifest.get(body_key, [])
        match = next((c for c in candidates if c.get("name") == winner_name), None)
        if match is None:
            raise KeyError(f"Winner '{winner_name}' not found in manifest for {public_name}")
        out.append((public_name, match["path"]))
    return out


def _load_emu_names(names_path: Path) -> dict[int, str]:
    obj = json.loads(names_path.read_text(encoding="utf-8"))
    names = obj.get("names", {})
    out = {}
    for k, v in names.items():
        try:
            out[int(k)] = str(v)
        except ValueError:
            continue
    return out


def _plot_stage_pair(
    shipping_name: str,
    shipping: Body,
    emu_idx: int,
    emu_name: str,
    emu: Body,
    *,
    morph: float,
    q: float,
    out_path: Path,
) -> dict:
    freqs = freq_points(n=512, sr=SR)
    ship_enc = shipping.corners.interpolate(morph, q)
    emu_enc = emu.corners.interpolate(morph, q)
    n = min(len(ship_enc), len(emu_enc))

    active_indices = [i for i in range(n) if (not _is_passthrough(ship_enc[i]) or not _is_passthrough(emu_enc[i]))]
    if not active_indices:
        active_indices = list(range(n))
    # Keep the first 6 active stages for readability.
    active_indices = active_indices[:6]

    fig, axes = plt.subplots(2, 3, figsize=(16, 9))
    fig.suptitle(
        f"{shipping_name} vs P2k_{emu_idx:03d} ({emu_name}) @ M={morph:.2f}, Q={q:.2f}",
        fontsize=13,
        fontweight="bold",
        color="#74a9ff",
    )
    stage_metrics: list[dict] = []

    for k, ax in enumerate(axes.flatten()):
        if k >= len(active_indices):
            ax.set_visible(False)
            continue
        sidx = active_indices[k]
        h_ship = stage_response(ship_enc[sidx], freqs, SR)
        h_emu = stage_response(emu_enc[sidx], freqs, SR)
        db_ship = 20.0 * np.log10(np.maximum(np.abs(h_ship), 1e-20))
        db_emu = 20.0 * np.log10(np.maximum(np.abs(h_emu), 1e-20))
        mad = float(np.mean(np.abs(db_ship - db_emu)))

        stage_metrics.append({
            "stage": sidx + 1,
            "mean_abs_db_delta": mad,
            "shipping_peak_db": float(np.max(db_ship)),
            "emu_peak_db": float(np.max(db_emu)),
        })

        ax.semilogx(freqs, db_ship, color="#e0e0ff", linewidth=1.8, label=shipping_name)
        ax.semilogx(freqs, db_emu, color="#ff5b5b", linewidth=1.2, alpha=0.95, label=f"P2k_{emu_idx:03d}")
        ax.set_title(f"Stage {sidx + 1} · mean|Δ|={mad:.2f} dB")
        ax.set_xlim(20, (SR / 2.0) - 1.0)
        ax.set_ylim(-35, 35)
        ax.set_xlabel("Hz")
        ax.set_ylabel("dB")
        ax.grid(True, which="both", linewidth=0.6, alpha=0.6)
        ax.legend(fontsize=7, framealpha=0.25, loc="upper right")

    fig.tight_layout()
    fig.savefig(out_path, dpi=170, bbox_inches="tight")
    plt.close(fig)

    mean_mad = float(np.mean([m["mean_abs_db_delta"] for m in stage_metrics])) if stage_metrics else 0.0
    return {
        "shipping_body": shipping_name,
        "emu_index": emu_idx,
        "emu_name": emu_name,
        "morph": morph,
        "q": q,
        "mean_abs_db_delta": mean_mad,
        "stage_metrics": stage_metrics,
        "plot": str(out_path),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Compare shipping winners to E-MU bodies stage-by-stage.")
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
        "--emu-dir",
        type=Path,
        default=ROOT / "cartridges" / "p2k",
    )
    parser.add_argument(
        "--names",
        type=Path,
        default=ROOT / "datasets" / "p2k_filter_names.json",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=ROOT / "vault" / "_plots" / "shipping_vs_emu_stages",
    )
    parser.add_argument(
        "--indices",
        type=str,
        default="",
        help="Comma-separated P2k indices. Empty = default subset.",
    )
    parser.add_argument("--morph", type=float, default=0.5)
    parser.add_argument("--q", type=float, default=0.5)
    args = parser.parse_args()

    if not args.gate.exists():
        raise FileNotFoundError(f"Missing gate file: {args.gate}")
    if not args.manifest.exists():
        raise FileNotFoundError(f"Missing manifest file: {args.manifest}")
    if not args.emu_dir.exists():
        raise FileNotFoundError(f"Missing E-MU directory: {args.emu_dir}")
    if not args.names.exists():
        raise FileNotFoundError(f"Missing names file: {args.names}")
    if not (0.0 <= args.morph <= 1.0 and 0.0 <= args.q <= 1.0):
        raise ValueError("--morph and --q must be in [0,1].")

    names = _load_emu_names(args.names)
    winners = _resolve_winners(args.gate, args.manifest)

    if args.indices.strip():
        indices = sorted({int(x.strip()) for x in args.indices.split(",") if x.strip()})
    else:
        indices = sorted(set(DEFAULT_EMU_INDICES) | set(DEFAULT_MAPPED_BY_BODY.values()))

    args.out_dir.mkdir(parents=True, exist_ok=True)
    results: list[dict] = []

    for public_name, body_path in winners:
        shipping = Body.from_json(body_path)
        ship_slug = _slug(public_name)
        body_dir = args.out_dir / ship_slug
        body_dir.mkdir(parents=True, exist_ok=True)

        body_indices = sorted(set(indices + [DEFAULT_MAPPED_BY_BODY[public_name]]))
        for idx in body_indices:
            emu_path = args.emu_dir / f"P2k_{idx:03d}.json"
            if not emu_path.exists():
                continue
            emu = Body.from_json(str(emu_path))
            emu_name = names.get(idx, emu.name)
            out_plot = body_dir / f"{ship_slug}__vs_p2k_{idx:03d}_stages.png"
            rec = _plot_stage_pair(
                public_name,
                shipping,
                idx,
                emu_name,
                emu,
                morph=args.morph,
                q=args.q,
                out_path=out_plot,
            )
            rec["shipping_path"] = body_path
            rec["emu_path"] = str(emu_path)
            results.append(rec)
            print(f"{public_name} vs P2k_{idx:03d} -> {out_plot}")

    results.sort(key=lambda r: (r["shipping_body"], r["mean_abs_db_delta"]))
    out_json = args.out_dir / "stage_comparison_summary.json"
    out_json.write_text(
        json.dumps({
            "morph": args.morph,
            "q": args.q,
            "indices": indices,
            "results": results,
        }, indent=2),
        encoding="utf-8",
    )
    print(out_json)


if __name__ == "__main__":
    main()
