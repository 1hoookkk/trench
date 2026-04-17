"""Stdlib-only SVG plotter for compiled-v1 cartridge cascades.

Generates frequency response plots (log Hz × dB) comparing filters side
by side. No matplotlib, no numpy — just math + string templating.

Usage:
    python plot_filters.py
    # Writes plots/p2k_013_talking_hedz.svg
    #       plots/speaker_knockerz_morph_path.svg
    #       plots/comparison_p2k_vs_sk.svg
"""
from __future__ import annotations

import json
import math
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
PLOTS_DIR = REPO / "plots"
# P2k_013 (Talking Hedz) lives in git history — this script extracts it
# via `git show` at run time rather than keeping a recovered copy in-tree,
# because the file was deleted in commit f146fb3 ("drop pre-pill cartridge
# archives") and we don't want to un-delete it just to plot it.
P2K_SKIN_NAME = "P2k_013"
P2K_GIT_REF = "f146fb3^"
P2K_GIT_PATH = "cartridges/p2k/P2k_013.json"
SK_DIR = REPO / "cartridges" / "engine" / "_source" / "shapes" / "speaker_knockerz"

SR = 39062.5

# Frequency grid: log-spaced from 20 Hz to Nyquist. 2048 points gives
# ~0.14% frequency resolution — dense enough to catch narrow formant
# peaks (typical bw ~1-3% of center frequency) without aliasing the
# curves. A sparser grid was originally used and hid multi-formant
# structure in the P2k_013 plot — see this commit's companion fix.
F_MIN = 20.0
F_MAX = 18000.0
N_FREQS = 2048

# dB display range
DB_MIN = -40.0
DB_MAX = 40.0

# SVG canvas
CANVAS_W = 960
CANVAS_H = 540
MARGIN_L = 70
MARGIN_R = 200  # room for legend
MARGIN_T = 60
MARGIN_B = 60
PLOT_W = CANVAS_W - MARGIN_L - MARGIN_R
PLOT_H = CANVAS_H - MARGIN_T - MARGIN_B


# ---------------------------------------------------------------------------
# DSP
# ---------------------------------------------------------------------------

def biquad_mag(c: list, omega: float) -> float:
    cos1 = math.cos(omega)
    sin1 = math.sin(omega)
    cos2 = math.cos(2.0 * omega)
    sin2 = math.sin(2.0 * omega)
    num_re = c[0] + c[1] * cos1 + c[2] * cos2
    num_im = -c[1] * sin1 - c[2] * sin2
    den_re = 1.0 + c[3] * cos1 + c[4] * cos2
    den_im = -c[3] * sin1 - c[4] * sin2
    num_mag = math.sqrt(num_re * num_re + num_im * num_im)
    den_mag = math.sqrt(den_re * den_re + den_im * den_im)
    if den_mag < 1e-18:
        return 1e6
    return num_mag / den_mag


def cascade_response(stages: list, boost: float, freqs: list) -> list:
    out = []
    for f in freqs:
        omega = 2.0 * math.pi * f / SR
        m = boost
        for s in stages:
            m *= biquad_mag(s, omega)
        out.append(m)
    return out


def db(x: float) -> float:
    if x <= 1e-10:
        return -200.0
    return 20.0 * math.log10(x)


def stages_from_keyframe(kf: dict) -> list:
    return [[s["c0"], s["c1"], s["c2"], s["c3"], s["c4"]] for s in kf["stages"]]


def log_freqs(n: int) -> list:
    lo = math.log10(F_MIN)
    hi = math.log10(F_MAX)
    return [10 ** (lo + (hi - lo) * i / (n - 1)) for i in range(n)]


# ---------------------------------------------------------------------------
# SVG geometry helpers
# ---------------------------------------------------------------------------

def x_for(f: float) -> float:
    lo = math.log10(F_MIN)
    hi = math.log10(F_MAX)
    return MARGIN_L + (math.log10(f) - lo) / (hi - lo) * PLOT_W


def y_for(dbv: float) -> float:
    clamped = max(DB_MIN, min(DB_MAX, dbv))
    return MARGIN_T + (DB_MAX - clamped) / (DB_MAX - DB_MIN) * PLOT_H


# ---------------------------------------------------------------------------
# SVG elements
# ---------------------------------------------------------------------------

def svg_header(title: str, subtitle: str = "") -> str:
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {CANVAS_W} {CANVAS_H}" '
        f'font-family="monospace" font-size="11">',
        f'<rect width="{CANVAS_W}" height="{CANVAS_H}" fill="#1a1a20"/>',
        f'<text x="{CANVAS_W/2}" y="28" fill="#eeeef0" text-anchor="middle" '
        f'font-size="16" font-weight="bold">{title}</text>',
    ]
    if subtitle:
        parts.append(
            f'<text x="{CANVAS_W/2}" y="45" fill="#999" text-anchor="middle" '
            f'font-size="11">{subtitle}</text>'
        )
    return "\n".join(parts)


def svg_grid() -> str:
    parts = []
    # Plot area background
    parts.append(
        f'<rect x="{MARGIN_L}" y="{MARGIN_T}" width="{PLOT_W}" height="{PLOT_H}" '
        f'fill="#0e0e14" stroke="#3a3a46" stroke-width="1"/>'
    )

    # Vertical grid lines at standard frequencies
    for f in [20, 50, 100, 200, 500, 1000, 2000, 5000, 10000]:
        x = x_for(f)
        parts.append(
            f'<line x1="{x:.1f}" y1="{MARGIN_T}" x2="{x:.1f}" '
            f'y2="{MARGIN_T + PLOT_H}" stroke="#2a2a35" stroke-width="0.5"/>'
        )
        label = f"{f}" if f < 1000 else f"{f // 1000}k"
        parts.append(
            f'<text x="{x:.1f}" y="{MARGIN_T + PLOT_H + 16}" fill="#888" '
            f'text-anchor="middle">{label}</text>'
        )

    # Horizontal grid lines every 10 dB
    for dbv in range(int(DB_MIN), int(DB_MAX) + 1, 10):
        y = y_for(dbv)
        stroke = "#555" if dbv == 0 else "#2a2a35"
        width = 1.0 if dbv == 0 else 0.5
        parts.append(
            f'<line x1="{MARGIN_L}" y1="{y:.1f}" x2="{MARGIN_L + PLOT_W}" '
            f'y2="{y:.1f}" stroke="{stroke}" stroke-width="{width}"/>'
        )
        parts.append(
            f'<text x="{MARGIN_L - 8}" y="{y + 4:.1f}" fill="#888" '
            f'text-anchor="end">{dbv:+d}</text>'
        )

    # Axis labels
    parts.append(
        f'<text x="{MARGIN_L + PLOT_W / 2}" y="{CANVAS_H - 12}" fill="#aaa" '
        f'text-anchor="middle">Frequency (Hz, log)</text>'
    )
    parts.append(
        f'<text x="20" y="{MARGIN_T + PLOT_H / 2}" fill="#aaa" '
        f'text-anchor="middle" transform="rotate(-90 20 {MARGIN_T + PLOT_H / 2})">'
        f'Magnitude (dB)</text>'
    )
    return "\n".join(parts)


def svg_curve(freqs: list, mags: list, color: str, label: str = "") -> str:
    pts = []
    for f, m in zip(freqs, mags):
        x = x_for(f)
        y = y_for(db(m))
        pts.append(f"{x:.1f},{y:.1f}")
    return (
        f'<polyline points="{" ".join(pts)}" fill="none" stroke="{color}" '
        f'stroke-width="1.6" opacity="0.9"><title>{label}</title></polyline>'
    )


def svg_legend(entries: list[tuple[str, str, str]]) -> str:
    """entries: list of (label, color, subtitle)"""
    x = MARGIN_L + PLOT_W + 14
    y = MARGIN_T + 10
    parts = [
        f'<rect x="{x - 6}" y="{y - 14}" width="186" height="{14 + len(entries) * 22}" '
        f'fill="#121218" stroke="#3a3a46" stroke-width="0.5"/>',
    ]
    for i, (label, color, sub) in enumerate(entries):
        ly = y + i * 22
        parts.append(
            f'<line x1="{x}" y1="{ly}" x2="{x + 22}" y2="{ly}" '
            f'stroke="{color}" stroke-width="3"/>'
        )
        parts.append(
            f'<text x="{x + 28}" y="{ly + 4}" fill="#e0e0e8" '
            f'font-weight="bold">{label}</text>'
        )
        if sub:
            parts.append(
                f'<text x="{x + 28}" y="{ly + 16}" fill="#888" font-size="9">{sub}</text>'
            )
    return "\n".join(parts)


def svg_footer() -> str:
    return "</svg>"


def find_local_maxima(freqs: list, mags: list, min_prominence_db: float = 3.0) -> list:
    """Return (hz, db) tuples for every local maximum in the response
    that rises at least `min_prominence_db` above its surrounding
    valleys. Filters out trivial ripples so we only mark real formants."""
    n = len(mags)
    raw = []
    for i in range(1, n - 1):
        if mags[i] > mags[i - 1] and mags[i] > mags[i + 1]:
            raw.append(i)
    if not raw:
        return []

    # Compute per-peak prominence as (peak_db - max(adjacent valley floor))
    mags_db = [db(m) for m in mags]
    peaks = []
    for idx in raw:
        # Walk left until we find a higher peak or the boundary
        left_floor = mags_db[idx]
        j = idx - 1
        while j > 0 and mags_db[j] < mags_db[idx]:
            if mags_db[j] < left_floor:
                left_floor = mags_db[j]
            j -= 1
        # Walk right similarly
        right_floor = mags_db[idx]
        j = idx + 1
        while j < n - 1 and mags_db[j] < mags_db[idx]:
            if mags_db[j] < right_floor:
                right_floor = mags_db[j]
            j += 1
        prominence = mags_db[idx] - max(left_floor, right_floor)
        if prominence >= min_prominence_db:
            peaks.append((freqs[idx], mags_db[idx]))
    return peaks


def svg_peak_markers(peaks: list, color: str) -> str:
    """Draw small circles + hz labels at each peak."""
    parts = []
    for hz, dbv in peaks:
        x = x_for(hz)
        y = y_for(dbv)
        parts.append(
            f'<circle cx="{x:.1f}" cy="{y:.1f}" r="3" fill="{color}" '
            f'stroke="#ffffff" stroke-width="0.8"/>'
        )
        parts.append(
            f'<text x="{x:.1f}" y="{y - 8:.1f}" fill="#eeeeee" font-size="8" '
            f'text-anchor="middle">{hz:.0f}</text>'
        )
    return "\n".join(parts)


# ---------------------------------------------------------------------------
# Plot 1: P2k_013 Talking Hedz, all 4 corners
# ---------------------------------------------------------------------------

P2K_COLORS = {
    "M0_Q0":    "#4aa3ff",  # blue
    "M100_Q0":  "#ff6a4a",  # red-orange
    "M0_Q100":  "#6aff9a",  # green
    "M100_Q100":"#ffb84a",  # amber
}

P2K_DESCRIPTIONS = {
    "M0_Q0":    "morph=0, Q=0",
    "M100_Q0":  "morph=1, Q=0",
    "M0_Q100":  "morph=0, Q=1",
    "M100_Q100":"morph=1, Q=1",
}


def load_p2k_from_history() -> dict:
    import subprocess
    result = subprocess.run(
        ["git", "-C", str(REPO), "show", f"{P2K_GIT_REF}:{P2K_GIT_PATH}"],
        capture_output=True, text=True, check=True,
    )
    return json.loads(result.stdout)


def plot_p2k() -> None:
    data = load_p2k_from_history()
    freqs = log_freqs(N_FREQS)

    parts = [
        svg_header(
            "P2K #013 — Talking Hedz (real E-mu vocal morph filter)",
            "4-corner morph surface, shape-bank authored, "
            f"sampleRate={data.get('sampleRate', SR)} Hz",
        ),
        svg_grid(),
    ]

    legend_entries = []
    for kf in data["keyframes"]:
        label = kf["label"]
        stages = stages_from_keyframe(kf)
        boost = kf.get("boost", 1.0)
        mags = cascade_response(stages, boost, freqs)
        color = P2K_COLORS.get(label, "#ffffff")
        parts.append(svg_curve(freqs, mags, color, label))

        # Mark every significant local maximum so the multi-formant
        # structure is visible (fixes the earlier single-peak bug).
        peaks = find_local_maxima(freqs, mags, min_prominence_db=3.0)
        parts.append(svg_peak_markers(peaks, color))

        # Legend subtitle: count of peaks + top peak
        if peaks:
            top_hz, top_db = max(peaks, key=lambda p: p[1])
            legend_entries.append(
                (label, color, f"{len(peaks)} peaks, top {top_hz:.0f} Hz ({top_db:+.0f} dB)")
            )
        else:
            legend_entries.append((label, color, "no peaks (flat/rolloff)"))

    parts.append(svg_legend(legend_entries))
    parts.append(svg_footer())

    out_path = PLOTS_DIR / "p2k_013_talking_hedz.svg"
    out_path.write_text("\n".join(parts), encoding="utf-8")
    print(f"wrote {out_path.relative_to(REPO)}")


# ---------------------------------------------------------------------------
# Plot 1b: P2K per-corner subplot grid, matching the E-mu UI's single-curve
# view so direct visual comparison against the plugin screenshots is easy.
# ---------------------------------------------------------------------------

def plot_p2k_per_corner() -> None:
    data = load_p2k_from_history()
    freqs = log_freqs(N_FREQS)

    # 2×2 subplot grid: M0_Q0 | M100_Q0 / M0_Q100 | M100_Q100
    # Each subplot gets its own mini axes. Canvas grows to fit.
    sub_w = 460
    sub_h = 240
    margin_l_sub = 55
    margin_r_sub = 18
    margin_t_sub = 30
    margin_b_sub = 40
    plot_w_sub = sub_w - margin_l_sub - margin_r_sub
    plot_h_sub = sub_h - margin_t_sub - margin_b_sub
    canvas_w = 2 * sub_w + 40
    canvas_h = 2 * sub_h + 90

    def x_for_sub(f, origin_x):
        lo = math.log10(F_MIN)
        hi = math.log10(F_MAX)
        return origin_x + margin_l_sub + (math.log10(f) - lo) / (hi - lo) * plot_w_sub

    def y_for_sub(dbv, origin_y):
        clamped = max(DB_MIN, min(DB_MAX, dbv))
        return origin_y + margin_t_sub + (DB_MAX - clamped) / (DB_MAX - DB_MIN) * plot_h_sub

    def subplot_grid(origin_x, origin_y, title):
        parts = []
        parts.append(
            f'<rect x="{origin_x + margin_l_sub}" y="{origin_y + margin_t_sub}" '
            f'width="{plot_w_sub}" height="{plot_h_sub}" fill="#0e0e14" '
            f'stroke="#3a3a46" stroke-width="1"/>'
        )
        for f in [50, 100, 200, 500, 1000, 2000, 5000, 10000]:
            x = x_for_sub(f, origin_x)
            parts.append(
                f'<line x1="{x:.1f}" y1="{origin_y + margin_t_sub}" '
                f'x2="{x:.1f}" y2="{origin_y + margin_t_sub + plot_h_sub}" '
                f'stroke="#2a2a35" stroke-width="0.4"/>'
            )
            label_f = f"{f}" if f < 1000 else f"{f // 1000}k"
            parts.append(
                f'<text x="{x:.1f}" y="{origin_y + margin_t_sub + plot_h_sub + 13}" '
                f'fill="#888" font-size="9" text-anchor="middle">{label_f}</text>'
            )
        for dbv in range(int(DB_MIN), int(DB_MAX) + 1, 20):
            y = y_for_sub(dbv, origin_y)
            stroke = "#555" if dbv == 0 else "#2a2a35"
            parts.append(
                f'<line x1="{origin_x + margin_l_sub}" y1="{y:.1f}" '
                f'x2="{origin_x + margin_l_sub + plot_w_sub}" y2="{y:.1f}" '
                f'stroke="{stroke}" stroke-width="0.4"/>'
            )
            parts.append(
                f'<text x="{origin_x + margin_l_sub - 6}" y="{y + 3:.1f}" '
                f'fill="#888" font-size="9" text-anchor="end">{dbv:+d}</text>'
            )
        parts.append(
            f'<text x="{origin_x + sub_w/2}" y="{origin_y + margin_t_sub - 10}" '
            f'fill="#e0e0e8" font-size="13" font-weight="bold" text-anchor="middle">{title}</text>'
        )
        return "\n".join(parts)

    def subplot_curve(origin_x, origin_y, freqs, mags, color):
        pts = []
        for f, m in zip(freqs, mags):
            x = x_for_sub(f, origin_x)
            y = y_for_sub(db(m), origin_y)
            pts.append(f"{x:.1f},{y:.1f}")
        return (
            f'<polyline points="{" ".join(pts)}" fill="none" stroke="{color}" '
            f'stroke-width="1.6"/>'
        )

    def subplot_peaks(origin_x, origin_y, peaks, color):
        parts = []
        for hz, dbv in peaks:
            x = x_for_sub(hz, origin_x)
            y = y_for_sub(dbv, origin_y)
            parts.append(
                f'<circle cx="{x:.1f}" cy="{y:.1f}" r="2.5" fill="{color}" '
                f'stroke="#ffffff" stroke-width="0.6"/>'
            )
            parts.append(
                f'<text x="{x:.1f}" y="{y - 6:.1f}" fill="#eeeeee" font-size="8" '
                f'text-anchor="middle">{hz:.0f}</text>'
            )
        return "\n".join(parts)

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {canvas_w} {canvas_h}" '
        f'font-family="monospace" font-size="11">',
        f'<rect width="{canvas_w}" height="{canvas_h}" fill="#1a1a20"/>',
        f'<text x="{canvas_w/2}" y="28" fill="#eeeef0" text-anchor="middle" '
        f'font-size="16" font-weight="bold">P2K #013 — Talking Hedz, all 4 corners as subplots</text>',
        f'<text x="{canvas_w/2}" y="46" fill="#999" text-anchor="middle" font-size="11">'
        f'matches the E-mu UI one-curve-per-setting view; peaks marked at ≥3 dB prominence</text>',
    ]

    # 4 corners in a 2×2 grid — matching the UI screenshots' arrangement
    corner_layout = [
        ("M0_Q0",    20,           70),
        ("M100_Q0",  sub_w + 40,   70),
        ("M0_Q100",  20,           sub_h + 90),
        ("M100_Q100",sub_w + 40,   sub_h + 90),
    ]

    kf_by_label = {kf["label"]: kf for kf in data["keyframes"]}

    for label, ox, oy in corner_layout:
        kf = kf_by_label[label]
        stages = stages_from_keyframe(kf)
        boost = kf.get("boost", 1.0)
        mags = cascade_response(stages, boost, freqs)
        peaks = find_local_maxima(freqs, mags, min_prominence_db=3.0)
        color = P2K_COLORS[label]
        title = f"{label}   ({len(peaks)} formants)"
        parts.append(subplot_grid(ox, oy, title))
        parts.append(subplot_curve(ox, oy, freqs, mags, color))
        parts.append(subplot_peaks(ox, oy, peaks, color))

    parts.append("</svg>")

    out_path = PLOTS_DIR / "p2k_013_per_corner.svg"
    out_path.write_text("\n".join(parts), encoding="utf-8")
    print(f"wrote {out_path.relative_to(REPO)}")

    # Also print the formant inventory to stdout for the report
    print("\nP2K #013 formant inventory (shape-bank data):")
    for label, _, _ in corner_layout:
        kf = kf_by_label[label]
        stages = stages_from_keyframe(kf)
        boost = kf.get("boost", 1.0)
        mags = cascade_response(stages, boost, freqs)
        peaks = find_local_maxima(freqs, mags, min_prominence_db=3.0)
        if peaks:
            peak_strs = [f"{hz:.0f}Hz({dbv:+.0f}dB)" for hz, dbv in peaks]
            print(f"  {label:<11} ({len(peaks)} peaks): {', '.join(peak_strs)}")
        else:
            print(f"  {label:<11} (no peaks — flat/rolloff)")


# ---------------------------------------------------------------------------
# Plot 2: Speaker Knockerz — 7 pill morph-path stations overlaid
# ---------------------------------------------------------------------------

SK_STATIONS = [
    ("sk_vault",    "VLT — The Vault",        "clean sub, clamped highs",   "#3a78ff"),
    ("sk_chest",    "CST — Chest Resonance",  "warm 150 Hz bloom",          "#d68a3a"),
    ("sk_choke",    "CHK — The Choke",        "−28 dB notch at 200 Hz",     "#c8324a"),
    ("sk_rip",      "RIP — Cardboard Rip",    "400 Hz resonance flare",     "#ff8a3a"),
    ("sk_rattle",   "RAT — The Rattle",       "800 Hz comb interference",   "#d2c83a"),
    ("sk_cry",      "CRY — Cone Cry",         "1.2 kHz screaming peak",     "#3ad8d8"),
    ("sk_fracture", "FRC — Total Fracture",   "3.2 kHz top, mids killed",   "#b84ad8"),
]


def plot_speaker_knockerz() -> None:
    freqs = log_freqs(N_FREQS)

    parts = [
        svg_header(
            "Speaker Knockerz — 7 morph-path stations (hand-authored via v2 pill template)",
            "Each pill is a static spectral state. Sub at 50 Hz normalized to unity across all 7 (invariant satisfied).",
        ),
        svg_grid(),
    ]

    # 0 dB reference line gets a special annotation
    legend_entries = []

    for key, label, desc, color in SK_STATIONS:
        pill_path = SK_DIR / f"{key}.json"
        data = json.loads(pill_path.read_text())
        kf = data["keyframes"][0]  # morph-invariant: M0_Q0 same as others
        stages = stages_from_keyframe(kf)
        boost = kf.get("boost", 1.0)
        mags = cascade_response(stages, boost, freqs)
        parts.append(svg_curve(freqs, mags, color, label))
        legend_entries.append((label, color, desc))

    parts.append(svg_legend(legend_entries))

    # Annotate the unity sub line at 50 Hz
    x_50 = x_for(50)
    y_0 = y_for(0)
    parts.append(
        f'<circle cx="{x_50:.1f}" cy="{y_0:.1f}" r="4" fill="none" '
        f'stroke="#ffffff" stroke-width="1.5" stroke-dasharray="2,2"/>'
    )
    parts.append(
        f'<text x="{x_50 + 8:.1f}" y="{y_0 - 8:.1f}" fill="#ffffff" '
        f'font-size="10">sub anchor: 50 Hz @ 0 dB (invariant)</text>'
    )

    parts.append(svg_footer())

    out_path = PLOTS_DIR / "speaker_knockerz_morph_path.svg"
    out_path.write_text("\n".join(parts), encoding="utf-8")
    print(f"wrote {out_path.relative_to(REPO)}")


# ---------------------------------------------------------------------------
# Plot 3: side-by-side comparison
# ---------------------------------------------------------------------------

def plot_comparison() -> None:
    freqs = log_freqs(N_FREQS)

    parts = [
        svg_header(
            "Real E-mu vs hand-authored — frequency response comparison",
            "P2K #013 Talking Hedz (4 corners, dashed) vs Speaker Knockerz 7 pills (solid, hand-authored)",
        ),
        svg_grid(),
    ]

    # P2K corners as dashed grey family
    data = load_p2k_from_history()
    for i, kf in enumerate(data["keyframes"]):
        label = kf["label"]
        stages = stages_from_keyframe(kf)
        boost = kf.get("boost", 1.0)
        mags = cascade_response(stages, boost, freqs)
        pts = []
        for f, m in zip(freqs, mags):
            pts.append(f"{x_for(f):.1f},{y_for(db(m)):.1f}")
        parts.append(
            f'<polyline points="{" ".join(pts)}" fill="none" '
            f'stroke="#6a6a74" stroke-width="1.0" opacity="0.55" '
            f'stroke-dasharray="4,3"><title>P2K {label}</title></polyline>'
        )

    # Speaker Knockerz pills as solid colored family
    legend_entries = [("P2K #013 corners", "#6a6a74", "real E-mu, dashed")]
    for key, label, desc, color in SK_STATIONS:
        pill_path = SK_DIR / f"{key}.json"
        pdata = json.loads(pill_path.read_text())
        kf = pdata["keyframes"][0]
        stages = stages_from_keyframe(kf)
        boost = kf.get("boost", 1.0)
        mags = cascade_response(stages, boost, freqs)
        parts.append(svg_curve(freqs, mags, color, label))
        legend_entries.append((label, color, desc))

    parts.append(svg_legend(legend_entries))
    parts.append(svg_footer())

    out_path = PLOTS_DIR / "comparison_p2k_vs_sk.svg"
    out_path.write_text("\n".join(parts), encoding="utf-8")
    print(f"wrote {out_path.relative_to(REPO)}")


def main() -> int:
    PLOTS_DIR.mkdir(exist_ok=True)
    plot_p2k()
    plot_p2k_per_corner()
    plot_speaker_knockerz()
    plot_comparison()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
