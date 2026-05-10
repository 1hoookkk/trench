"""forge/morph_grid_scanner.py — read-only morph-grid diagnostic.

Reads the 4 bundled body cartridges from juce-shell/assets/cartridges/,
evaluates |H(jw)| over a 21x5 (morph x q) grid at each (morph, q) point,
and emits PNG heatmaps + a stdout summary.

Does NOT write to any cartridge JSON. Does NOT invoke cmake, cargo, or any
authoring script. Read-only except for PNG output to forge/diagnostics/morph_grid/.

Corner order convention (matching cartridge.rs):
    corners[0] = M0_Q0
    corners[1] = M100_Q0
    corners[2] = M0_Q100
    corners[3] = M100_Q100

Bilinear interpolation (matching cartridge.rs::interpolate):
    q_m0   = corners[0] + (corners[2] - corners[0]) * q
    q_m1   = corners[1] + (corners[3] - corners[1]) * q
    result = q_m0 + (q_m1 - q_m0) * morph
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np

# ─── Paths ────────────────────────────────────────────────────────────────────
REPO = Path(__file__).resolve().parent.parent
CARTRIDGE_DIR = REPO / "juce-shell" / "assets" / "cartridges"
OUT_DIR = Path(__file__).resolve().parent / "diagnostics" / "morph_grid"

AUTHORING_SR: float = 39062.5
N_FREQ = 1024
FREQS = np.geomspace(20.0, 20000.0, N_FREQ)
PEAK_AMPLITUDE_TARGET: float = 2.5  # linear — from boost_normalizer.py
NUM_ACTIVE_STAGES = 6

BODY_FILES = [
    "speaker_knockerz.json",
    "aluminum_siding.json",
    "small_talk_ah_ee.json",
    "cul_de_sac.json",
]

CORNER_LABELS = ["M0_Q0", "M100_Q0", "M0_Q100", "M100_Q100"]
# cartridge.rs corner index order: [M0_Q0, M100_Q0, M0_Q100, M100_Q100]
# corners[0]=M0_Q0, corners[1]=M100_Q0, corners[2]=M0_Q100, corners[3]=M100_Q100


# ─── Step 1: RawStage → DF2T compile (mirrors emu_resonator.rs) ───────────────
def compile_rawstage(rs: dict[str, Any]) -> list[float] | None:
    """Compile one RawStage dict to [c0, c1, c2, c3, c4].

    Formula (emu_resonator.rs lines 26-40):
        a1 = -2 * r * cos(2π * freq / sr)  [stored directly as rs["a1"]]
        a2 = r * r                           [stored directly as rs["r"]^2]
        c0 = 1.0 + val1
        c1 = a1  + val2
        c2 = a2  - val3
        c3 = a1
        c4 = a2

    Returns None for passthrough (r==0 and a1==0).
    """
    a1 = float(rs["a1"])
    r = float(rs["r"])
    if r == 0.0 and a1 == 0.0:
        return None  # passthrough — cartridge.rs uses identity, not this compile path
    a2 = r * r
    return [
        1.0 + float(rs["val1"]),  # c0
        a1 + float(rs["val2"]),   # c1
        a2 - float(rs["val3"]),   # c2
        a1,                        # c3
        a2,                        # c4
    ]


def passthrough_coeffs() -> list[float]:
    """Identity DF2T stage: [c0=1, c1=0, c2=0, c3=0, c4=0]."""
    return [1.0, 0.0, 0.0, 0.0, 0.0]


def parse_corner(keyframe: dict[str, Any]) -> tuple[list[list[float]], float]:
    """Parse one keyframe → (coeffs_6x5, boost).

    Takes first 6 stages (NUM_ACTIVE_STAGES). Passthrough (r=0,a1=0) → identity.
    Returns list of 6 coefficient lists [c0,c1,c2,c3,c4].
    """
    stages_raw = keyframe["stages"][:NUM_ACTIVE_STAGES]
    coeffs = []
    for rs in stages_raw:
        c = compile_rawstage(rs)
        if c is None:
            c = passthrough_coeffs()
        coeffs.append(c)
    # Pad with passthrough if fewer than NUM_ACTIVE_STAGES stages present
    while len(coeffs) < NUM_ACTIVE_STAGES:
        coeffs.append(passthrough_coeffs())
    boost = float(keyframe.get("boost", 1.0))
    return coeffs, boost


def load_cartridge(path: Path) -> tuple[list[list[list[float]]], list[float], str]:
    """Load cartridge JSON → (corners_4x6x5, boosts_4, name).

    Corner index: 0=M0_Q0, 1=M100_Q0, 2=M0_Q100, 3=M100_Q100
    (matches cartridge.rs corner array ordering).
    """
    data = json.loads(path.read_text(encoding="utf-8"))
    name = data.get("name", path.stem)
    keyframes = {kf["label"]: kf for kf in data["keyframes"]}
    corners = []
    boosts = []
    for lbl in CORNER_LABELS:
        kf = keyframes[lbl]
        c, b = parse_corner(kf)
        corners.append(c)
        boosts.append(b)
    return corners, boosts, name


# ─── Step 2: H(jw) evaluation ─────────────────────────────────────────────────
def eval_cascade_mag(coeffs_6x5: list[list[float]], freqs: np.ndarray) -> np.ndarray:
    """Evaluate |H(jw)| cascade product over freqs.

    Each stage: H_k(z) = (c0 + c1*z^-1 + c2*z^-2) / (1 + c3*z^-1 + c4*z^-2)
    Cascade = product of all stage H_k(z).

    This is equivalent to the runtime DF2T serial chain for stable LTI stages:
    output of stage k feeds input of stage k+1, which in the frequency domain
    means the total transfer function is the product of individual H_k(e^jw).
    """
    w = 2.0 * np.pi * freqs / AUTHORING_SR
    e1 = np.exp(-1j * w)
    e2 = np.exp(-2j * w)
    mag = np.ones(len(freqs), dtype=np.float64)
    for c in coeffs_6x5:
        c0, c1, c2, c3, c4 = c
        num = c0 + c1 * e1 + c2 * e2
        den = 1.0 + c3 * e1 + c4 * e2
        mag *= np.abs(num / den)
    return mag


# ─── Step 3: Bilinear interpolation (matches cartridge.rs::interpolate) ────────
def bilinear_interp_coeffs(
    corners: list[list[list[float]]],  # [4][6][5]
    morph: float,
    q: float,
) -> list[list[float]]:
    """Bilinearly interpolate 4 corner coefficient arrays.

    Exactly mirrors cartridge.rs::interpolate (lines 167-177):
        q_m0   = corners[0] + (corners[2] - corners[0]) * q
        q_m1   = corners[1] + (corners[3] - corners[1]) * q
        result = q_m0 + (q_m1 - q_m0) * morph

    corners[0]=M0_Q0, corners[1]=M100_Q0, corners[2]=M0_Q100, corners[3]=M100_Q100
    """
    result = []
    for stage in range(NUM_ACTIVE_STAGES):
        interp_stage = []
        for c_idx in range(5):
            v00 = corners[0][stage][c_idx]  # M0_Q0
            v10 = corners[1][stage][c_idx]  # M100_Q0
            v01 = corners[2][stage][c_idx]  # M0_Q100
            v11 = corners[3][stage][c_idx]  # M100_Q100
            # Q interpolation at each morph extreme
            q_m0 = v00 + (v01 - v00) * q
            q_m1 = v10 + (v11 - v10) * q
            # Morph interpolation
            val = q_m0 + (q_m1 - q_m0) * morph
            interp_stage.append(val)
        result.append(interp_stage)
    return result


def bilinear_interp_boost(boosts: list[float], morph: float, q: float) -> float:
    """Bilinearly interpolate 4 corner boosts.

    Matches cartridge.rs::interpolate_boost (lines 179-183).
    boosts[0]=M0_Q0, [1]=M100_Q0, [2]=M0_Q100, [3]=M100_Q100
    """
    q_m0 = boosts[0] + (boosts[2] - boosts[0]) * q
    q_m1 = boosts[1] + (boosts[3] - boosts[1]) * q
    return q_m0 + (q_m1 - q_m0) * morph


# ─── boost_normalizer-equivalent: compute target boost per corner ──────────────
def compute_corner_boost(coeffs_6x5: list[list[float]], target: float = PEAK_AMPLITUDE_TARGET) -> float:
    """Return target / cascade_peak, matching boost_normalizer.py logic."""
    mag = eval_cascade_mag(coeffs_6x5, FREQS)
    peak = float(np.max(mag))
    return target / peak if peak > 0.0 else 1.0


# ─── Step 3 (shape preservation): verify boost is a scalar vertical shift ──────
def shape_preservation_check(
    coeffs_6x5: list[list[float]],
    boost: float,
) -> float:
    """Compute r(w) = 20*log10(|H_boosted(w)|) - 20*log10(|H_raw(w)|).

    This should equal 20*log10(boost) at every frequency within float epsilon.
    Returns max |deviation| in dB across all freq points.
    """
    mag = eval_cascade_mag(coeffs_6x5, FREQS)
    mag = np.maximum(mag, 1e-300)
    db_raw = 20.0 * np.log10(mag)
    db_boosted = 20.0 * np.log10(mag * boost)
    expected_shift = 20.0 * math.log10(max(boost, 1e-300))
    residual = np.abs((db_boosted - db_raw) - expected_shift)
    return float(np.max(residual))


# ─── Main grid scan ────────────────────────────────────────────────────────────
def scan_body(json_path: Path) -> dict:
    """Full morph-grid scan for one body.

    Returns dict with all measurements needed for report and plots.
    """
    corners, boosts_stored, name = load_cartridge(json_path)

    morph_vals = np.linspace(0.0, 1.0, 21)
    q_vals = np.linspace(0.0, 1.0, 5)

    # Grid: peak[i_morph, i_q] and peak_times_boost[i_morph, i_q]
    peak_grid = np.zeros((21, 5))
    boosted_grid_stored = np.zeros((21, 5))   # peak * interpolated stored boost
    boosted_grid_computed = np.zeros((21, 5)) # peak * interpolated computed boost

    # Compute per-corner boosts from boost_normalizer logic
    boosts_computed = [compute_corner_boost(corners[k]) for k in range(4)]

    for i_m, morph in enumerate(morph_vals):
        for i_q, q in enumerate(q_vals):
            interp_coeffs = bilinear_interp_coeffs(corners, morph, q)
            mag = eval_cascade_mag(interp_coeffs, FREQS)
            peak = float(np.max(mag))
            peak_grid[i_m, i_q] = peak

            boost_s = bilinear_interp_boost(boosts_stored, morph, q)
            boost_c = bilinear_interp_boost(boosts_computed, morph, q)
            boosted_grid_stored[i_m, i_q] = peak * boost_s
            boosted_grid_computed[i_m, i_q] = peak * boost_c

    # Find global max and argmax
    peak_max = float(np.max(peak_grid))
    flat_idx = int(np.argmax(peak_grid))
    i_m_max, i_q_max = np.unravel_index(flat_idx, peak_grid.shape)
    morph_star = float(morph_vals[i_m_max])
    q_star = float(q_vals[i_q_max])

    required_body_boost = PEAK_AMPLITUDE_TARGET / peak_max if peak_max > 0 else 1.0

    # Shape preservation: per corner, max deviation in dB
    shape_devs: dict[str, float] = {}
    for k, lbl in enumerate(CORNER_LABELS):
        b = boosts_computed[k]
        dev = shape_preservation_check(corners[k], b)
        shape_devs[lbl] = dev

    return {
        "name": name,
        "body_id": json_path.stem,
        "corners": corners,
        "boosts_stored": boosts_stored,
        "boosts_computed": boosts_computed,
        "morph_vals": morph_vals,
        "q_vals": q_vals,
        "peak_grid": peak_grid,
        "boosted_grid_stored": boosted_grid_stored,
        "boosted_grid_computed": boosted_grid_computed,
        "peak_max": peak_max,
        "morph_star": morph_star,
        "q_star": q_star,
        "required_body_boost": required_body_boost,
        "shape_devs": shape_devs,
        "max_stored_boosted": float(np.max(boosted_grid_stored)),
        "max_computed_boosted": float(np.max(boosted_grid_computed)),
    }


# ─── Plotting ─────────────────────────────────────────────────────────────────
def plot_body(result: dict, out_dir: Path) -> Path:
    """Two-panel heatmap: left=peak_grid, right=peak*computed_boost.

    Both on log color scale. Worst point marked with X.
    Target line at 2.5 marked with dashed contour on right panel.
    """
    body_id = result["body_id"]
    name = result["name"]
    morph_vals = result["morph_vals"]
    q_vals = result["q_vals"]
    peak_grid = result["peak_grid"]
    boosted_grid = result["boosted_grid_computed"]
    peak_max = result["peak_max"]
    morph_star = result["morph_star"]
    q_star = result["q_star"]

    # Indices for worst points
    def worst_idx(grid):
        flat = int(np.argmax(grid))
        return np.unravel_index(flat, grid.shape)

    i_m_raw, i_q_raw = worst_idx(peak_grid)
    i_m_bst, i_q_bst = worst_idx(boosted_grid)

    fig, axes = plt.subplots(1, 2, figsize=(16, 6), dpi=120)
    fig.patch.set_facecolor("#16191a")
    fig.suptitle(
        f"MORPH-GRID DIAGNOSTIC :: {name}  ({body_id})\n"
        f"21 morph × 5 q points  |  SR = 39062.5 Hz  |  6 active stages",
        color="#ffba00", fontweight="bold", fontsize=11,
    )

    # ── Left: raw peak ──
    ax0 = axes[0]
    ax0.set_facecolor("#0d0f10")
    norm_raw = mcolors.LogNorm(
        vmin=max(peak_grid.min(), 1e-6),
        vmax=peak_grid.max(),
    )
    im0 = ax0.pcolormesh(
        q_vals, morph_vals, peak_grid,
        cmap="inferno", norm=norm_raw, shading="auto",
    )
    # Mark worst point
    ax0.scatter(
        q_vals[i_q_raw], morph_vals[i_m_raw],
        marker="x", color="#22ddff", s=160, linewidths=2.5,
        label=f"worst  peak={peak_grid[i_m_raw, i_q_raw]:.3f}\n"
              f"morph={morph_vals[i_m_raw]:.2f} q={q_vals[i_q_raw]:.2f}",
        zorder=10,
    )
    cb0 = fig.colorbar(im0, ax=ax0, pad=0.02)
    cb0.set_label("|H| linear", color="white", fontsize=9)
    cb0.ax.yaxis.set_tick_params(color="white", labelcolor="white")
    ax0.set_xlabel("q", color="white", fontsize=9)
    ax0.set_ylabel("morph", color="white", fontsize=9)
    ax0.set_title(
        f"peak(morph, q)   max={peak_max:.3f}  ({20*math.log10(max(peak_max,1e-12)):+.1f} dB)",
        color="#88ee99", fontsize=10, fontweight="bold",
    )
    ax0.tick_params(colors="white", labelsize=8)
    for spine in ax0.spines.values():
        spine.set_color("#3a3e42")
    ax0.legend(facecolor="#22262a", edgecolor="#444", labelcolor="white",
               fontsize=7, loc="upper left")

    # ── Right: peak * computed boost ──
    ax1 = axes[1]
    ax1.set_facecolor("#0d0f10")
    bst_max = float(np.max(boosted_grid))
    norm_bst = mcolors.LogNorm(
        vmin=max(boosted_grid.min(), 1e-6),
        vmax=max(bst_max, PEAK_AMPLITUDE_TARGET * 1.01),
    )
    im1 = ax1.pcolormesh(
        q_vals, morph_vals, boosted_grid,
        cmap="plasma", norm=norm_bst, shading="auto",
    )
    # Target contour at 2.5
    try:
        cs = ax1.contour(
            q_vals, morph_vals, boosted_grid,
            levels=[PEAK_AMPLITUDE_TARGET],
            colors=["#44ee44"], linewidths=[1.8],
        )
        ax1.clabel(cs, fmt=f"target={PEAK_AMPLITUDE_TARGET}", fontsize=7, colors="#44ee44")
    except Exception:
        pass
    # Cliff contour at 5.0
    try:
        cs2 = ax1.contour(
            q_vals, morph_vals, boosted_grid,
            levels=[5.0],
            colors=["#ff4444"], linewidths=[1.4], linestyles=["--"],
        )
        ax1.clabel(cs2, fmt="cliff=5.0", fontsize=7, colors="#ff4444")
    except Exception:
        pass
    # Mark worst boosted point
    ax1.scatter(
        q_vals[i_q_bst], morph_vals[i_m_bst],
        marker="x", color="#22ddff", s=160, linewidths=2.5,
        label=f"worst  peak×boost={boosted_grid[i_m_bst, i_q_bst]:.3f}\n"
              f"morph={morph_vals[i_m_bst]:.2f} q={q_vals[i_q_bst]:.2f}",
        zorder=10,
    )
    cb1 = fig.colorbar(im1, ax=ax1, pad=0.02)
    cb1.set_label("|H| × boost  linear", color="white", fontsize=9)
    cb1.ax.yaxis.set_tick_params(color="white", labelcolor="white")
    ax1.set_xlabel("q", color="white", fontsize=9)
    ax1.set_ylabel("morph", color="white", fontsize=9)
    ax1.set_title(
        f"peak × interp_boost   max={bst_max:.3f}  target={PEAK_AMPLITUDE_TARGET}",
        color="#ffba00", fontsize=10, fontweight="bold",
    )
    ax1.tick_params(colors="white", labelsize=8)
    for spine in ax1.spines.values():
        spine.set_color("#3a3e42")
    ax1.legend(facecolor="#22262a", edgecolor="#444", labelcolor="white",
               fontsize=7, loc="upper left")

    fig.tight_layout()
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{body_id}.png"
    fig.savefig(out_path, facecolor=fig.get_facecolor())
    plt.close(fig)
    return out_path


# ─── Entry point ──────────────────────────────────────────────────────────────
def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 72)
    print("MORPH-GRID SCANNER  |  forge/morph_grid_scanner.py")
    print("=" * 72)

    all_results = []
    for fname in BODY_FILES:
        json_path = CARTRIDGE_DIR / fname
        if not json_path.exists():
            print(f"  ERROR: not found: {json_path}")
            sys.exit(1)
        print(f"\nScanning {fname} ...")
        result = scan_body(json_path)
        all_results.append(result)
        png_path = plot_body(result, OUT_DIR)
        print(f"  -> plot written: {png_path}")

    # ── Summary table ──────────────────────────────────────────────────────────
    print()
    print("=" * 72)
    print("PER-BODY MORPH-GRID SUMMARY")
    print("=" * 72)
    hdr = f"{'body':<22}  {'peak_max':>10}  {'peak_max_dB':>11}  {'morph*':>7}  {'q*':>5}  {'req_body_boost':>15}"
    print(hdr)
    print("-" * len(hdr))
    for r in all_results:
        peak_db = 20.0 * math.log10(max(r["peak_max"], 1e-12))
        print(
            f"{r['body_id']:<22}  {r['peak_max']:>10.4f}  {peak_db:>+11.2f} dB"
            f"  {r['morph_star']:>7.2f}  {r['q_star']:>5.2f}  {r['required_body_boost']:>15.6f}"
        )

    # ── Per-corner boost adequacy ──────────────────────────────────────────────
    print()
    print("=" * 72)
    print("PER-CORNER BOOST ADEQUACY  (stored boosts from JSON vs computed target boosts)")
    print("Threshold: max(peak × interp_boost) > 5.0 → UNSAFE")
    print("=" * 72)
    for r in all_results:
        body = r["body_id"]
        max_s = r["max_stored_boosted"]
        max_c = r["max_computed_boosted"]
        verdict_s = "UNSAFE" if max_s > 5.0 else "SAFE"
        verdict_c = "UNSAFE" if max_c > 5.0 else "SAFE"
        overshoot_s = 20.0 * math.log10(max_s / 5.0) if max_s > 5.0 else 0.0
        overshoot_c = 20.0 * math.log10(max_c / 5.0) if max_c > 5.0 else 0.0
        print(f"\n  {body}")
        print(f"    stored  boosts {r['boosts_stored']}:")
        print(f"      max(peak × interp_stored_boost)   = {max_s:.4f}  [{verdict_s}]"
              + (f"  overshoot {overshoot_s:+.2f} dB above cliff" if verdict_s == "UNSAFE" else ""))
        print(f"    computed boosts {[f'{b:.4f}' for b in r['boosts_computed']]}:")
        print(f"      max(peak × interp_computed_boost) = {max_c:.4f}  [{verdict_c}]"
              + (f"  overshoot {overshoot_c:+.2f} dB above cliff" if verdict_c == "UNSAFE" else ""))

    # ── Shape preservation ─────────────────────────────────────────────────────
    print()
    print("=" * 72)
    print("SHAPE PRESERVATION  (boost = vertical shift; max |r(w)| per (body, corner))")
    print("r(w) = 20*log10(|H_boosted(w)|) - 20*log10(|H_raw(w)|)")
    print("Expected = 20*log10(boost) = constant.  Anything > 1e-6 dB is a bug.")
    print("=" * 72)
    print(f"{'body':<22}  {'M0_Q0':>12}  {'M100_Q0':>12}  {'M0_Q100':>12}  {'M100_Q100':>12}")
    print("-" * 76)
    for r in all_results:
        sd = r["shape_devs"]
        vals = [sd.get(lbl, float("nan")) for lbl in CORNER_LABELS]
        row = f"{r['body_id']:<22}"
        for v in vals:
            row += f"  {v:>12.6e}"
        print(row)

    # ── Files written ──────────────────────────────────────────────────────────
    print()
    print("=" * 72)
    print("FILES WRITTEN")
    print("=" * 72)
    for fname in BODY_FILES:
        body_id = Path(fname).stem
        print(f"  {OUT_DIR / (body_id + '.png')}")

    print()
    print("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
