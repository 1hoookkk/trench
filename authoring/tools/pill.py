"""forge/pill.py — Identity Scoring Oracle.

Bake N candidates per Throw. Score each by Character (peak-to-notch ratio +
centroid monotonicity + slope steepness, 0-100). Survival of the elite —
keep only the top K. Label each survivor with its closest heritage relative
("92% Razor_Blades match"). Output a contact sheet of new Launch Bodies.

Not 16 ways to pass a gate. 16 candidates with identity.

Usage:
    python forge/pill.py ah_to_ee --candidates 1000 --survivors 16 --seed 42
"""
from __future__ import annotations

import argparse
import json
import random
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import forge_cli  # noqa: E402
from mine_variants import (  # noqa: E402
    VariantParams, random_params, apply_params, restore_params, trajectory_db,
)
from render_reference import heritage_to_compiled  # noqa: E402


REPO = Path(__file__).resolve().parent.parent
CALIB_DIR = REPO / "reference" / "calibration"
PILL_DIR = REPO / "cartridges" / "factory" / "generated" / "qlaw" / "pills"


# ─── score components (each 0-100) ────────────────────────────────────────

def _spectral_centroid(db: np.ndarray, freqs: np.ndarray) -> float:
    """Spectral centroid in Hz, weighted by linear magnitude."""
    lin = np.power(10.0, db / 20.0)
    total = lin.sum()
    if total < 1e-12:
        return 0.0
    return float(np.sum(freqs * lin) / total)


def _bell(value: float, target: float, sigma: float) -> float:
    """Bell-curve scoring: 100 at target, falls off both sides."""
    return float(100.0 * np.exp(-((value - target) ** 2) / (2.0 * sigma * sigma)))


def peak_to_notch_score(traces: list[np.ndarray]) -> float:
    """Reward MUSICAL dynamic range — bell curve centered at 50 dB.
    Heritage cartridges (Talking_Hedz, Razor_Blades) sit 40-60 dB. Above
    100 dB is volume-bomb territory; below 25 dB is washed out."""
    if not traces:
        return 0.0
    ranges = [float(t.max() - t.min()) for t in traces]
    avg = float(np.mean(ranges))
    return _bell(avg, target=50.0, sigma=22.0)


def centroid_monotonicity_score(traces: list[np.ndarray], freqs: np.ndarray) -> float:
    """How monotonically does the spectral centroid move M0→M100?
    Linear in |corr| × 100. Static (no motion) → 0."""
    if len(traces) < 3:
        return 0.0
    centroids = np.array([_spectral_centroid(t, freqs) for t in traces])
    if centroids.std() < 1.0:
        return 0.0
    steps = np.arange(len(traces), dtype=float)
    corr = np.corrcoef(steps, centroids)[0, 1]
    return float(abs(corr) * 100.0) if np.isfinite(corr) else 0.0


def slope_steepness_score(traces: list[np.ndarray], freqs: np.ndarray) -> float:
    """Reward MODERATE rolloff slope — bell curve centered at 50 dB/oct.
    Heritage cartridges roll off 30-70 dB/oct. Above 150 = brittle digital
    spike; below 15 = mush. Bell at 50, σ=30."""
    if not traces:
        return 0.0
    log_freqs = np.log2(freqs)
    max_slope = 0.0
    for t in traces:
        slope = np.abs(np.gradient(t, log_freqs))
        m = float(slope.max())
        if m > max_slope:
            max_slope = m
    return _bell(max_slope, target=50.0, sigma=30.0)


def composure_score(traces: list[np.ndarray]) -> float:
    """Penalize cartridges that just scream loud everywhere (mean dB >> 0)
    or fully under-driven (mean dB << -30). Bell at 0 dB mean, σ=30."""
    if not traces:
        return 0.0
    mean_db = float(np.mean([t.mean() for t in traces]))
    return _bell(mean_db, target=0.0, sigma=30.0)


def trajectory_smoothness_score(traces: list[np.ndarray]) -> float:
    """Reward smooth step-to-step transitions. Penalize jumpy morph paths."""
    if len(traces) < 2:
        return 0.0
    arr = np.stack(traces)  # (n_steps, n_freqs)
    jumps = np.linalg.norm(arr[1:] - arr[:-1], axis=1) / arr.shape[1]
    if jumps.std() < 1e-6:
        return 100.0
    cv = jumps.std() / (jumps.mean() + 1e-6)  # coefficient of variation
    # Low CV = smooth. CV=0 → 100, CV=1 → ~0.
    return float(100.0 * np.exp(-cv * 2.0))


def character_score(traces: list[np.ndarray], freqs: np.ndarray) -> dict[str, float]:
    """0-100 Character Score, weighted:
        peak-to-notch       30%   (bell @ 50 dB)
        centroid monotonicity 25% (linear)
        slope steepness     20%   (bell @ 50 dB/oct)
        composure           15%   (bell @ 0 dB mean)
        smoothness          10%   (low step-to-step variance)
    """
    p2n = peak_to_notch_score(traces)
    mono = centroid_monotonicity_score(traces, freqs)
    slope = slope_steepness_score(traces, freqs)
    comp = composure_score(traces)
    smooth = trajectory_smoothness_score(traces)
    total = (0.30 * p2n + 0.25 * mono + 0.20 * slope
             + 0.15 * comp + 0.10 * smooth)
    return {
        "total": float(total),
        "p2n": p2n, "mono": mono, "slope": slope,
        "comp": comp, "smooth": smooth,
    }


# ─── heritage matching ────────────────────────────────────────────────────

def load_heritage_traces(freqs: np.ndarray, n_steps: int) -> dict[str, list[np.ndarray]]:
    out: dict[str, list[np.ndarray]] = {}
    for jf in sorted(CALIB_DIR.glob("*.json")):
        try:
            data = json.loads(jf.read_text(encoding="utf-8"))
        except Exception:
            continue
        cart = heritage_to_compiled(data)
        if cart is None:
            continue
        try:
            out[jf.stem] = trajectory_db(cart, freqs, n_steps=n_steps)
        except Exception:
            continue
    return out


def _dc_normalize(arr: np.ndarray) -> np.ndarray:
    """Subtract mean → only the SHAPE is compared, not gross loudness."""
    return arr - arr.mean()


def heritage_match(candidate_traces: list[np.ndarray],
                   heritage: dict[str, list[np.ndarray]]) -> tuple[str, float]:
    """DC-normalized cosine similarity between flattened trajectories.
    Subtracting the mean removes the gross-loudness bias so we compare the
    SHAPE of the spectrum, not just the level."""
    cand = _dc_normalize(np.concatenate(candidate_traces))
    cand_norm = np.linalg.norm(cand)
    if cand_norm < 1e-9 or not heritage:
        return ("none", 0.0)
    best_name = "none"
    best_sim = -1.0
    for name, ref_traces in heritage.items():
        ref = _dc_normalize(np.concatenate(ref_traces))
        if ref.shape != cand.shape:
            continue
        ref_norm = np.linalg.norm(ref)
        if ref_norm < 1e-9:
            continue
        sim = float(np.dot(cand, ref) / (cand_norm * ref_norm))
        if sim > best_sim:
            best_sim = sim
            best_name = name
    return (best_name, max(0.0, best_sim) * 100.0)


# ─── plotting ─────────────────────────────────────────────────────────────

@dataclass
class Survivor:
    idx: int
    params: VariantParams
    cart: dict
    traces: list[np.ndarray]
    score: dict[str, float]
    heritage: tuple[str, float]


def render_pill_sheet(arc: forge_cli.Arc, survivors: list[Survivor],
                      n_total_baked: int, out_path: Path) -> None:
    n = len(survivors)
    cols = 4
    rows = (n + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(16.0, 2.9 * rows), dpi=140, sharex=True)
    fig.patch.set_facecolor("#16191a")
    fig.suptitle(
        f"FORGE :: PILL — {arc.name}  ({arc.recipe})\n"
        f"top {n} of {n_total_baked} candidates  ·  Character Score + Heritage Match",
        color="#ffba00", fontweight="bold", fontsize=12,
    )
    freqs = np.geomspace(40.0, 18000.0, 256)
    flat = axes.flatten() if hasattr(axes, "flatten") else [axes]

    cmap = plt.get_cmap("plasma")
    n_steps = len(survivors[0].traces) if survivors else 9
    line_colors = [cmap(i / (n_steps - 1)) for i in range(n_steps)]

    for ax, s in zip(flat, survivors):
        ax.set_facecolor("#0d0f10")
        for col, tr in zip(line_colors, s.traces):
            ax.semilogx(freqs, tr, color=col, lw=0.9, alpha=0.85)
        ax.semilogx(freqs, s.traces[0],  color="#22ddff", lw=1.6, alpha=0.95)
        ax.semilogx(freqs, s.traces[-1], color="#ffba00", lw=1.6, alpha=0.95)

        title = (f"#{s.idx:02d}  CHAR {s.score['total']:.0f}/100\n"
                 f"{s.heritage[1]:.0f}% {s.heritage[0]}")
        ax.set_title(title, color="white", fontsize=9, fontweight="bold")
        sub = (f"p2n {s.score['p2n']:.0f}  mono {s.score['mono']:.0f}  "
               f"slope {s.score['slope']:.0f}  ·  r={s.params.radius:.3f}")
        ax.text(0.02, 0.04, sub, transform=ax.transAxes,
                color="#9aa0a6", fontsize=7, ha="left", va="bottom",
                bbox=dict(facecolor="#0d0f10", edgecolor="none", pad=2, alpha=0.6))
        ax.grid(True, which="both", alpha=0.18, color="white")
        ax.tick_params(colors="white", labelsize=7)
        for spine in ax.spines.values():
            spine.set_color("#3a3e42")

    for j in range(n, len(flat)):
        flat[j].axis("off")

    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, facecolor=fig.get_facecolor())
    plt.close(fig)


# ─── main ────────────────────────────────────────────────────────────────

def main(argv: list[str]) -> int:
    p = argparse.ArgumentParser(
        description="Identity Scoring Oracle — bake N candidates, survive top K, label heritage."
    )
    p.add_argument("arc", help="arc name (run forge_cli list to see them)")
    p.add_argument("--candidates", type=int, default=1000,
                   help="number of candidates to bake (default 1000)")
    p.add_argument("--survivors", type=int, default=16,
                   help="top-K to keep (default 16)")
    p.add_argument("--seed", type=int, default=None)
    p.add_argument("--policy", default="A", choices=["A", "B"])
    args = p.parse_args(argv)

    rng = random.Random(args.seed)
    arc = forge_cli.find_arc(args.arc)

    score_freqs = np.geomspace(40.0, 18000.0, 256)
    n_steps = 9

    print(f"loading heritage references…", flush=True)
    heritage = load_heritage_traces(score_freqs, n_steps=n_steps)
    print(f"  → {len(heritage)} heritage cartridges loaded")

    print(f"baking {args.candidates} candidates for {arc.name} ({arc.recipe})…", flush=True)
    pool: list[tuple[VariantParams, dict, list[np.ndarray], dict[str, float]]] = []
    rejects = 0
    progress_every = max(1, args.candidates // 20)

    for i in range(args.candidates):
        params = random_params(rng)
        snap = apply_params(params)
        try:
            frame_map = forge_cli.arc_to_frame_map(arc)
            stages, meta = forge_cli.plan(frame_map, arc.recipe)
            solved = forge_cli.solve(stages, args.policy)
            cart = forge_cli.build_cartridge(meta, solved, args.policy)
            gate_result = forge_cli.gate(cart)
        finally:
            restore_params(snap)

        if not gate_result.passed:
            rejects += 1
            continue

        traces = trajectory_db(cart, score_freqs, n_steps=n_steps)
        score = character_score(traces, score_freqs)
        pool.append((params, cart, traces, score))

        if (i + 1) % progress_every == 0:
            print(f"  {i + 1}/{args.candidates} baked ({rejects} gate-rejected, "
                  f"{len(pool)} pool)", flush=True)

    print(f"\n{len(pool)} candidates passed gate ({rejects} rejected)")
    if not pool:
        print("no survivors — gate rejected everything", file=sys.stderr)
        return 1

    pool.sort(key=lambda t: t[3]["total"], reverse=True)
    elites = pool[:args.survivors]

    survivors: list[Survivor] = []
    print(f"\nELITE SURVIVORS — {arc.name}")
    print(f"{'#':>3}  {'CHAR':>5}  {'p2n':>4} {'mono':>4} {'slope':>5}  "
          f"{'rad':>5}  heritage")
    for new_idx, (params, cart, traces, score) in enumerate(elites):
        match = heritage_match(traces, heritage)
        survivors.append(Survivor(
            idx=new_idx, params=params, cart=cart,
            traces=traces, score=score, heritage=match,
        ))
        print(f"  {new_idx:02d}  {score['total']:5.1f}  "
              f"{score['p2n']:4.0f} {score['mono']:4.0f} {score['slope']:5.0f}  "
              f"{params.radius:5.3f}  {match[1]:5.1f}% {match[0]}")

    arc_dir = PILL_DIR / arc.name
    arc_dir.mkdir(parents=True, exist_ok=True)
    for s in survivors:
        char_int = int(round(s.score["total"]))
        match_int = int(round(s.heritage[1]))
        path = arc_dir / (
            f"{arc.name}_pill_{s.idx:02d}_"
            f"char{char_int:03d}_"
            f"{match_int:02d}pct_{s.heritage[0]}.json"
        )
        path.write_text(json.dumps(s.cart, indent=2) + "\n", encoding="utf-8")

    sheet = arc_dir / "_pill_sheet.png"
    render_pill_sheet(arc, survivors, len(pool), sheet)
    print(f"\npill sheet → {sheet}")

    if sys.platform == "win32":
        try:
            import os
            os.startfile(str(sheet))  # type: ignore[attr-defined]
        except Exception:
            pass
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
