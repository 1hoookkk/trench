"""forge/iconic_author_dry_run.py — read-only morph-grid diagnostic for the
patched iconic_author.py.  Does NOT invoke emit(), main(), or any function
that writes cartridge JSON or PNG to the factory output directory.

Obtains in-memory RawStage lists by calling the four body design functions
(body_speaker_knockerz, body_aluminum_siding, body_small_talk_ah_ee,
body_cul_de_sac) directly from iconic_author.py.  These functions have zero
write side-effects — they return plain Python data structures.  Only emit()
(iconic_author.py lines 428-455) performs JSON/PNG writes; we never call it.

Outputs (written only to forge/diagnostics/iconic_dry_run/):
    speaker_knockerz.png
    aluminum_siding.png
    small_talk_ah_ee.png
    cul_de_sac.png
    summary.txt
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np

# ─── Paths ────────────────────────────────────────────────────────────────────
FORGE = Path(__file__).resolve().parent
REPO  = FORGE.parent
sys.path.insert(0, str(FORGE))

OUT_DIR = FORGE / "diagnostics" / "iconic_dry_run"
SHIPPED_DIR = REPO / "juce-shell" / "assets" / "cartridges"
ONDISK_DIR  = REPO / "cartridges" / "factory" / "generated" / "iconic"

# ─── Import design functions WITHOUT triggering writes ────────────────────────
# Only body_* functions are imported.  emit(), main(), and all globals that
# reference OUT_DIR are never called.  Import is safe because iconic_author.py
# has no module-level side-effects (all writes are inside emit()).
from iconic_author import (  # noqa: E402
    body_speaker_knockerz,
    body_aluminum_siding,
    body_small_talk_ah_ee,
    body_cul_de_sac,
    stage_to_rawstage,
    rawstage_to_coefs,
    pad_passthrough,
    AUTHORING_SR,
    N_SLOTS,
)
from boost_normalizer import _cascade_peak, EVAL_FREQS, PEAK_AMPLITUDE_TARGET

# ─── Grid parameters (matching morph_grid_scanner.py) ────────────────────────
N_MORPH       = 21
N_Q           = 5
MORPH_VALS    = np.linspace(0.0, 1.0, N_MORPH)
Q_VALS        = np.linspace(0.0, 1.0, N_Q)
N_FREQ        = 1024
FREQS         = np.geomspace(20.0, 20000.0, N_FREQ)
CLIFF         = 5.0   # absolute linear amplitude — AGC cliff threshold


# ─── Step 1: compile RawStage list to coef dicts ─────────────────────────────

def compile_stages(raw_list: list[dict]) -> list[dict]:
    """Convert list of RawStage dicts to DF2T coef dicts (skipping passthroughs)."""
    return [c for rs in raw_list
            if (c := rawstage_to_coefs(rs)) and not c.get("_pass")]


# ─── Step 2: build 4 corners for one body ────────────────────────────────────
# iconic_author.py auto-duplicates Q corners from M corners (Q is runtime
# infrastructure, not authored data).  M0_Q0 == M0_Q100 and M100_Q0 == M100_Q100
# in the patched author.  We replicate that here exactly as assemble_cartridge()
# does (lines 123-124), then compile each corner to [c0..c4] float lists for
# bilinear interpolation.

def build_corners(body_fn) -> tuple[list[list[list[float]]], list[float], str, float]:
    """Call body design function, bake stages, compile to coef arrays.

    Returns:
        corners  : list[4][n_active][5]  — [M0_Q0, M100_Q0, M0_Q100, M100_Q100]
        boosts_computed : list[4] — PEAK_AMPLITUDE_TARGET / corner_peak
        name     : str
        c0_uniform : float
    """
    name, c0, m0_stages, m100_stages, _brief = body_fn()

    # Bake stages exactly as assemble_cartridge() does
    m0_baked   = pad_passthrough([stage_to_rawstage(s, c0) for s in m0_stages])
    m100_baked = pad_passthrough([stage_to_rawstage(s, c0) for s in m100_stages])

    # 4 corners: M0_Q0, M100_Q0, M0_Q100, M100_Q100
    # Q corners auto-duplicated from M corners
    raw_corners = [m0_baked, m100_baked, m0_baked, m100_baked]

    # Convert to [c0,c1,c2,c3,c4] float arrays (6 active stages; passthrough = identity)
    def to_coef_arrays(raw_list: list[dict]) -> list[list[float]]:
        result = []
        for rs in raw_list[:6]:
            if rs["r"] == 0.0 and rs["a1"] == 0.0:
                result.append([1.0, 0.0, 0.0, 0.0, 0.0])  # identity passthrough
            else:
                c = rawstage_to_coefs(rs)
                result.append([c["c0"], c["c1"], c["c2"], c["c3"], c["c4"]])
        while len(result) < 6:
            result.append([1.0, 0.0, 0.0, 0.0, 0.0])
        return result

    corners = [to_coef_arrays(raw_corners[k]) for k in range(4)]

    # Per-corner peaks and boost
    def corner_peak(corner_coefs: list[list[float]]) -> float:
        mag = eval_cascade_mag(corner_coefs, FREQS)
        return float(np.max(mag))

    corner_peaks = [corner_peak(corners[k]) for k in range(4)]
    boosts_computed = [
        PEAK_AMPLITUDE_TARGET / p if p > 0.0 else 1.0
        for p in corner_peaks
    ]

    return corners, boosts_computed, corner_peaks, name, c0


# ─── Step 3: transfer function evaluation ────────────────────────────────────

def eval_cascade_mag(coef_arrays: list[list[float]], freqs: np.ndarray) -> np.ndarray:
    """Evaluate |H(jw)| over freqs for one set of compiled biquad stages."""
    w  = 2.0 * np.pi * freqs / AUTHORING_SR
    e1 = np.exp(-1j * w)
    e2 = np.exp(-2j * w)
    mag = np.ones(len(freqs), dtype=np.float64)
    for c in coef_arrays:
        c0, c1, c2, c3, c4 = c
        num = c0 + c1 * e1 + c2 * e2
        den = 1.0 + c3 * e1 + c4 * e2
        mag *= np.abs(num / den)
    return mag


# ─── Step 4: bilinear interpolation (matches cartridge.rs::interpolate) ──────

def bilinear_interp(corners: list[list[list[float]]], morph: float, q: float) -> list[list[float]]:
    """Bilinearly interpolate 4 corner coefficient arrays.
    corners[0]=M0_Q0, [1]=M100_Q0, [2]=M0_Q100, [3]=M100_Q100
    """
    result = []
    for s in range(6):
        row = []
        for ci in range(5):
            v00 = corners[0][s][ci]
            v10 = corners[1][s][ci]
            v01 = corners[2][s][ci]
            v11 = corners[3][s][ci]
            q_m0 = v00 + (v01 - v00) * q
            q_m1 = v10 + (v11 - v10) * q
            row.append(q_m0 + (q_m1 - q_m0) * morph)
        result.append(row)
    return result


def bilinear_interp_boost(boosts: list[float], morph: float, q: float) -> float:
    """Bilinearly interpolate 4 corner boosts.
    boosts[0]=M0_Q0, [1]=M100_Q0, [2]=M0_Q100, [3]=M100_Q100
    """
    q_m0 = boosts[0] + (boosts[2] - boosts[0]) * q
    q_m1 = boosts[1] + (boosts[3] - boosts[1]) * q
    return q_m0 + (q_m1 - q_m0) * morph


# ─── Step 5: morph-grid scan ─────────────────────────────────────────────────

def scan_body(body_fn) -> dict:
    """Full morph-grid scan for one body design function."""
    corners, boosts_computed, corner_peaks, name, c0 = build_corners(body_fn)

    peak_grid   = np.zeros((N_MORPH, N_Q), dtype=np.float64)
    boosted_grid = np.zeros((N_MORPH, N_Q), dtype=np.float64)

    for i_m, morph in enumerate(MORPH_VALS):
        for i_q, q in enumerate(Q_VALS):
            interp_coefs = bilinear_interp(corners, morph, q)
            mag = eval_cascade_mag(interp_coefs, FREQS)
            peak = float(np.max(mag))
            boost = bilinear_interp_boost(boosts_computed, morph, q)
            peak_grid[i_m, i_q]    = peak
            boosted_grid[i_m, i_q] = peak * boost

    peak_max    = float(np.max(peak_grid))
    flat_idx    = int(np.argmax(peak_grid))
    i_m_max, i_q_max = np.unravel_index(flat_idx, peak_grid.shape)
    morph_star  = float(MORPH_VALS[i_m_max])
    q_star      = float(Q_VALS[i_q_max])
    required_body_boost = PEAK_AMPLITUDE_TARGET / peak_max if peak_max > 0.0 else 1.0

    max_boosted = float(np.max(boosted_grid))
    verdict     = "SAFE" if max_boosted <= CLIFF else "UNSAFE"
    db_over     = 20.0 * math.log10(max_boosted / CLIFF) if max_boosted > CLIFF else 0.0

    return {
        "name":               name,
        "c0_uniform":         c0,
        "corners":            corners,
        "boosts_computed":    boosts_computed,
        "corner_peaks":       corner_peaks,
        "peak_grid":          peak_grid,
        "boosted_grid":       boosted_grid,
        "peak_max":           peak_max,
        "morph_star":         morph_star,
        "q_star":             q_star,
        "required_body_boost": required_body_boost,
        "max_boosted":        max_boosted,
        "verdict":            verdict,
        "db_over":            db_over,
    }


# ─── Step 6: per-corner adequacy ─────────────────────────────────────────────

def per_corner_adequacy(result: dict) -> str:
    """Return formatted per-corner adequacy line for summary."""
    boosts  = result["boosts_computed"]
    c_peaks = result["corner_peaks"]
    labels  = ["M0_Q0", "M100_Q0", "M0_Q100", "M100_Q100"]
    boost_str = "  ".join(f"{lbl}={b:.6f}" for lbl, b in zip(labels, boosts))
    peak_str  = "  ".join(f"{lbl}={p:.4f}" for lbl, p in zip(labels, c_peaks))
    return (
        f"  boosts: {boost_str}\n"
        f"  corner_peaks: {peak_str}\n"
        f"  max(peak*boost) over grid = {result['max_boosted']:.5f}  "
        f"verdict={result['verdict']}"
        + (f"  dB_over_cliff={result['db_over']:+.2f} dB" if result['verdict'] == "UNSAFE" else "")
    )


# ─── Step 7: load shipped / on-disk JSON and scan ────────────────────────────

def load_and_scan_json(json_path: Path) -> dict | None:
    """Run identical morph-grid scan on an existing cartridge JSON file."""
    import json as _json
    if not json_path.exists():
        return None
    data = _json.loads(json_path.read_text(encoding="utf-8"))
    name = data.get("name", json_path.stem)

    CORNER_LABELS = ["M0_Q0", "M100_Q0", "M0_Q100", "M100_Q100"]

    def parse_corner_json(kf: dict) -> list[list[float]]:
        result = []
        for rs in kf["stages"][:6]:
            a1 = float(rs["a1"]); r = float(rs["r"])
            if r == 0.0 and a1 == 0.0:
                result.append([1.0, 0.0, 0.0, 0.0, 0.0])
            else:
                a2 = r * r
                result.append([
                    1.0 + float(rs["val1"]),
                    a1 + float(rs["val2"]),
                    a2 - float(rs["val3"]),
                    a1, a2,
                ])
        while len(result) < 6:
            result.append([1.0, 0.0, 0.0, 0.0, 0.0])
        return result

    kf_map = {kf["label"]: kf for kf in data["keyframes"]}
    corners = [parse_corner_json(kf_map[lbl]) for lbl in CORNER_LABELS]

    def corner_peak(coefs: list[list[float]]) -> float:
        return float(np.max(eval_cascade_mag(coefs, FREQS)))

    corner_peaks = [corner_peak(corners[k]) for k in range(4)]
    boosts_computed = [PEAK_AMPLITUDE_TARGET / p if p > 0.0 else 1.0 for p in corner_peaks]

    peak_grid   = np.zeros((N_MORPH, N_Q), dtype=np.float64)
    boosted_grid = np.zeros((N_MORPH, N_Q), dtype=np.float64)
    for i_m, morph in enumerate(MORPH_VALS):
        for i_q, q in enumerate(Q_VALS):
            ic = bilinear_interp(corners, morph, q)
            mag = eval_cascade_mag(ic, FREQS)
            peak = float(np.max(mag))
            boost = bilinear_interp_boost(boosts_computed, morph, q)
            peak_grid[i_m, i_q]    = peak
            boosted_grid[i_m, i_q] = peak * boost

    peak_max = float(np.max(peak_grid))
    return {"name": name, "peak_max": peak_max}


# ─── Step 8: heatmap plot ─────────────────────────────────────────────────────

def plot_body(result: dict, out_path: Path) -> None:
    name    = result["name"]
    pg      = result["peak_grid"]
    bg      = result["boosted_grid"]
    peak_max = result["peak_max"]
    morph_star = result["morph_star"]
    q_star  = result["q_star"]

    def worst_idx(grid):
        flat = int(np.argmax(grid))
        return np.unravel_index(flat, grid.shape)

    i_m_r, i_q_r = worst_idx(pg)
    i_m_b, i_q_b = worst_idx(bg)

    fig, axes = plt.subplots(1, 2, figsize=(16, 6), dpi=120)
    fig.patch.set_facecolor("#16191a")
    fig.suptitle(
        f"ICONIC DRY-RUN MORPH-GRID :: {name}   c0={result['c0_uniform']:.4f}\n"
        f"21 morph x 5 q  |  SR={AUTHORING_SR} Hz  |  6 stages",
        color="#ffba00", fontweight="bold", fontsize=11,
    )

    # Left: raw peak
    ax0 = axes[0]
    ax0.set_facecolor("#0d0f10")
    norm0 = mcolors.LogNorm(vmin=max(pg.min(), 1e-6), vmax=pg.max())
    im0 = ax0.pcolormesh(Q_VALS, MORPH_VALS, pg, cmap="inferno", norm=norm0, shading="auto")
    ax0.scatter(Q_VALS[i_q_r], MORPH_VALS[i_m_r], marker="x", color="#22ddff", s=160, lw=2.5,
                label=f"worst  peak={pg[i_m_r,i_q_r]:.3f}\nm={MORPH_VALS[i_m_r]:.2f} q={Q_VALS[i_q_r]:.2f}",
                zorder=10)
    cb0 = fig.colorbar(im0, ax=ax0, pad=0.02)
    cb0.set_label("|H| linear", color="white", fontsize=9)
    cb0.ax.yaxis.set_tick_params(color="white", labelcolor="white")
    ax0.set_xlabel("q", color="white"); ax0.set_ylabel("morph", color="white")
    ax0.set_title(
        f"peak(morph,q)   max={peak_max:.4f}  ({20*math.log10(max(peak_max,1e-12)):+.2f} dB)",
        color="#88ee99", fontsize=10, fontweight="bold",
    )
    ax0.tick_params(colors="white", labelsize=8)
    for sp in ax0.spines.values(): sp.set_color("#3a3e42")
    ax0.legend(facecolor="#22262a", edgecolor="#444", labelcolor="white", fontsize=7)

    # Right: peak * computed boost
    ax1 = axes[1]
    ax1.set_facecolor("#0d0f10")
    bst_max = float(np.max(bg))
    norm1 = mcolors.LogNorm(vmin=max(bg.min(), 1e-6), vmax=max(bst_max, PEAK_AMPLITUDE_TARGET * 1.01))
    im1 = ax1.pcolormesh(Q_VALS, MORPH_VALS, bg, cmap="plasma", norm=norm1, shading="auto")
    # Target contour at 2.5
    try:
        cs = ax1.contour(Q_VALS, MORPH_VALS, bg, levels=[PEAK_AMPLITUDE_TARGET],
                         colors=["#44ee44"], linewidths=[1.8])
        ax1.clabel(cs, fmt=f"target={PEAK_AMPLITUDE_TARGET}", fontsize=7, colors="#44ee44")
    except Exception:
        pass
    # Cliff contour at 5.0
    try:
        cs2 = ax1.contour(Q_VALS, MORPH_VALS, bg, levels=[CLIFF],
                          colors=["#ff4444"], linewidths=[1.4], linestyles=["--"])
        ax1.clabel(cs2, fmt="cliff=5.0", fontsize=7, colors="#ff4444")
    except Exception:
        pass
    ax1.scatter(Q_VALS[i_q_b], MORPH_VALS[i_m_b], marker="x", color="#22ddff", s=160, lw=2.5,
                label=f"worst  peak*boost={bg[i_m_b,i_q_b]:.3f}\nm={MORPH_VALS[i_m_b]:.2f} q={Q_VALS[i_q_b]:.2f}",
                zorder=10)
    # Mark target=2.5 reference line on colorbar
    cb1 = fig.colorbar(im1, ax=ax1, pad=0.02)
    cb1.set_label("|H| * boost  linear", color="white", fontsize=9)
    cb1.ax.yaxis.set_tick_params(color="white", labelcolor="white")
    ax1.set_xlabel("q", color="white"); ax1.set_ylabel("morph", color="white")
    ax1.set_title(
        f"peak*interp_boost   max={bst_max:.4f}  target={PEAK_AMPLITUDE_TARGET}  "
        + (f"[UNSAFE +{result['db_over']:.2f}dB]" if result['verdict']=="UNSAFE" else "[SAFE]"),
        color="#ff4444" if result['verdict']=="UNSAFE" else "#ffba00",
        fontsize=10, fontweight="bold",
    )
    ax1.tick_params(colors="white", labelsize=8)
    for sp in ax1.spines.values(): sp.set_color("#3a3e42")
    ax1.legend(facecolor="#22262a", edgecolor="#444", labelcolor="white", fontsize=7)

    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, facecolor=fig.get_facecolor())
    plt.close(fig)


# ─── Step 9: entry point ─────────────────────────────────────────────────────

BODY_FNS = [
    ("speaker_knockerz",  body_speaker_knockerz),
    ("aluminum_siding",   body_aluminum_siding),
    ("small_talk_ah_ee",  body_small_talk_ah_ee),
    ("cul_de_sac",        body_cul_de_sac),
]

CLIFF_DB = 20.0 * math.log10(CLIFF)  # +14.0 dB


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    lines: list[str] = []

    header = (
        "=" * 80 + "\n"
        "ICONIC DRY-RUN DIAGNOSTIC  |  forge/iconic_author_dry_run.py\n"
        "Stage source: iconic_author.body_* functions (NO write side-effects)\n"
        "Boost source: PEAK_AMPLITUDE_TARGET / corner_peak (re-derived analytically)\n"
        f"PEAK_AMPLITUDE_TARGET = {PEAK_AMPLITUDE_TARGET}  CLIFF = {CLIFF}  SR = {AUTHORING_SR}\n"
        "=" * 80
    )
    print(header)
    lines.append(header)

    results: list[dict] = []
    for body_id, body_fn in BODY_FNS:
        print(f"\nScanning {body_id} ...")
        result = scan_body(body_fn)
        result["body_id"] = body_id
        results.append(result)
        png_path = OUT_DIR / f"{body_id}.png"
        plot_body(result, png_path)
        print(f"  -> plot: {png_path}")

    # ── B: Per-body morph-grid table ─────────────────────────────────────────
    section = "\n" + "=" * 80 + "\nB. PER-BODY MORPH-GRID SUMMARY (21x5 grid)\n" + "=" * 80
    print(section); lines.append(section)

    hdr = f"{'body':<22}  {'peak_max lin':>14}  {'peak_max dB':>12}  {'morph*':>7}  {'q*':>5}  {'req_body_boost':>15}"
    print(hdr); lines.append(hdr)
    div = "-" * len(hdr)
    print(div); lines.append(div)
    for r in results:
        pk_db = 20.0 * math.log10(max(r["peak_max"], 1e-12))
        row = (
            f"{r['body_id']:<22}  {r['peak_max']:>14.6f}  {pk_db:>+12.4f} dB"
            f"  {r['morph_star']:>7.2f}  {r['q_star']:>5.2f}  {r['required_body_boost']:>15.8f}"
        )
        print(row); lines.append(row)

    # ── C: Per-corner adequacy ────────────────────────────────────────────────
    section = "\n" + "=" * 80 + "\nC. PER-CORNER ADEQUACY  (computed boost via PEAK_AMPLITUDE_TARGET/corner_peak)\n" + "=" * 80
    print(section); lines.append(section)

    corner_labels = ["M0_Q0", "M100_Q0", "M0_Q100", "M100_Q100"]
    for r in results:
        bid   = r["body_id"]
        boosts = r["boosts_computed"]
        cpeaks = r["corner_peaks"]
        max_b  = r["max_boosted"]
        verd   = r["verdict"]
        db_ov  = r["db_over"]

        row = f"\n  {bid}"
        boost_parts = "  ".join(f"{lbl}={b:.6f}" for lbl, b in zip(corner_labels, boosts))
        peak_parts  = "  ".join(f"{lbl}={p:.5f}" for lbl, p in zip(corner_labels, cpeaks))
        row += f"\n    corner peaks:  {peak_parts}"
        row += f"\n    corner boosts: {boost_parts}"
        row += f"\n    max(peak*interp_boost) over grid = {max_b:.6f}"
        row += f"  [{verd}]"
        if verd == "UNSAFE":
            row += f"  dB_over_cliff = {db_ov:+.4f} dB"
        print(row); lines.append(row)

    # ── D: Comparison to shipped assets ──────────────────────────────────────
    section = "\n" + "=" * 80 + "\nD. COMPARISON TO SHIPPED ASSETS  (juce-shell/assets/cartridges/)\n" + "=" * 80
    print(section); lines.append(section)

    for r in results:
        bid    = r["body_id"]
        dry_pk = r["peak_max"]
        shipped = load_and_scan_json(SHIPPED_DIR / f"{bid}.json")
        if shipped is None:
            row = f"  {bid:<22}  DRY={dry_pk:.6f}  SHIPPED=NOT FOUND"
        else:
            sh_pk = shipped["peak_max"]
            ratio  = dry_pk / sh_pk if sh_pk > 0 else float("inf")
            ratio_db = 20.0 * math.log10(max(ratio, 1e-12))
            comparison = ("patched-iconic much tamer" if ratio < 0.1
                          else "patched-iconic tamer" if ratio < 0.9
                          else "comparable" if ratio < 1.1
                          else "patched-iconic higher")
            row = (
                f"  {bid:<22}  DRY={dry_pk:.6f} ({20*math.log10(max(dry_pk,1e-12)):+.2f}dB)  "
                f"SHIPPED={sh_pk:.6f} ({20*math.log10(max(sh_pk,1e-12)):+.2f}dB)  "
                f"ratio={ratio:.4f} ({ratio_db:+.2f}dB)  -> {comparison}"
            )
        print(row); lines.append(row)

    # ── E: Comparison to on-disk generated/iconic ────────────────────────────
    section = "\n" + "=" * 80 + "\nE. COMPARISON TO ON-DISK generated/iconic/\n" + "=" * 80
    print(section); lines.append(section)

    for r in results:
        bid    = r["body_id"]
        dry_pk = r["peak_max"]
        ondisk = load_and_scan_json(ONDISK_DIR / f"{bid}.json")
        if ondisk is None:
            row = f"  {bid:<22}  DRY={dry_pk:.6f}  ON-DISK=NOT FOUND"
        else:
            od_pk = ondisk["peak_max"]
            ratio  = dry_pk / od_pk if od_pk > 0 else float("inf")
            ratio_db = 20.0 * math.log10(max(ratio, 1e-12))
            if abs(ratio - 1.0) < 0.001:
                verdict = "IDENTICAL (patched source == on-disk)"
            elif ratio < 1.0:
                verdict = f"patched tamer by {-ratio_db:.2f} dB"
            else:
                verdict = f"patched hotter by {ratio_db:.2f} dB"
            row = (
                f"  {bid:<22}  DRY={dry_pk:.6f} ({20*math.log10(max(dry_pk,1e-12)):+.2f}dB)  "
                f"ON-DISK={od_pk:.6f} ({20*math.log10(max(od_pk,1e-12)):+.2f}dB)  -> {verdict}"
            )
        print(row); lines.append(row)

    # ── F: Files written ──────────────────────────────────────────────────────
    section = "\n" + "=" * 80 + "\nF. FILES WRITTEN\n" + "=" * 80
    print(section); lines.append(section)
    for body_id, _ in BODY_FNS:
        p = OUT_DIR / f"{body_id}.png"
        print(f"  {p}"); lines.append(f"  {p}")
    summary_path = OUT_DIR / "summary.txt"
    print(f"  {summary_path}"); lines.append(f"  {summary_path}")

    # Write summary.txt
    summary_path.write_text("\n".join(lines), encoding="utf-8")
    print(f"\nDone. Summary saved -> {summary_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
