"""forge/stage_rca.py — per-stage root-cause analysis for cascade peak disasters.

Read-only diagnostic.  Outputs to forge/diagnostics/stage_rca/.
No cartridge JSON is written; no runtime/plugin/AGC code is touched.

Targets (morph=0.00, q=1.00):
  speaker_knockerz.json  — peak 6.24e10
  cul_de_sac.json        — peak 7.53e10

Compile path (mirrors emu_resonator.rs exactly):
    c0 = 1.0 + val1
    c1 = a1  + val2
    c2 = r^2 - val3
    c3 = a1
    c4 = r^2

Interpolation (mirrors cartridge.rs::interpolate, corner order below):
    corners[0] = M0_Q0,   corners[1] = M100_Q0
    corners[2] = M0_Q100, corners[3] = M100_Q100

    q_m0   = corners[0][s][c] + (corners[2][s][c] - corners[0][s][c]) * q
    q_m1   = corners[1][s][c] + (corners[3][s][c] - corners[1][s][c]) * q
    result = q_m0 + (q_m1 - q_m0) * morph
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np

# ─── Paths ────────────────────────────────────────────────────────────────────
REPO        = Path(__file__).resolve().parent.parent
CART_DIR    = REPO / "juce-shell" / "assets" / "cartridges"
OUT_DIR     = Path(__file__).resolve().parent / "diagnostics" / "stage_rca"

AUTHORING_SR: float = 39062.5
N_FREQ = 4096
FREQS = np.geomspace(20.0, 20_000.0, N_FREQ)

CORNER_ORDER = ["M0_Q0", "M100_Q0", "M0_Q100", "M100_Q100"]
NUM_ACTIVE = 6

# ─── Compile ──────────────────────────────────────────────────────────────────

def compile_rawstage(rs: dict[str, Any]) -> list[float] | None:
    """RawStage -> [c0, c1, c2, c3, c4].  Returns None for passthrough (r=0, a1=0)."""
    a1  = float(rs["a1"])
    r   = float(rs["r"])
    if r == 0.0 and a1 == 0.0:
        return None
    a2 = r * r
    return [
        1.0 + float(rs["val1"]),   # c0
        a1  + float(rs["val2"]),   # c1
        a2  - float(rs["val3"]),   # c2
        a1,                         # c3
        a2,                         # c4
    ]


def load_corners(json_path: Path) -> dict[str, list[list[float] | None]]:
    """Load all 4 corners.  Returns dict[label] -> list of 6 compiled coeff vectors
    (None = passthrough — will be treated as identity in cascade eval)."""
    data = json.loads(json_path.read_text(encoding="utf-8"))
    out: dict[str, list[list[float] | None]] = {}
    for kf in data["keyframes"]:
        label = kf["label"]
        if label not in CORNER_ORDER:
            continue
        stages = kf["stages"][:NUM_ACTIVE]
        out[label] = [compile_rawstage(s) for s in stages]
    return out


# ─── Bilinear interpolation ───────────────────────────────────────────────────

def bilinear_interpolate(
    corners: dict[str, list[list[float] | None]],
    morph: float, q: float,
) -> list[list[float]]:
    """Interpolate all 6 active stages at (morph, q).
    Passthrough (None) treated as [1,0,0,0,0] (identity).
    Returns list of 6 [c0..c4] arrays."""
    def get(label: str, stage: int, coef: int) -> float:
        v = corners[label][stage]
        if v is None:
            return 1.0 if coef == 0 else 0.0   # identity biquad
        return v[coef]

    result = []
    for s in range(NUM_ACTIVE):
        row = []
        for c in range(5):
            q_m0 = get("M0_Q0",    s, c) + (get("M0_Q100",   s, c) - get("M0_Q0",    s, c)) * q
            q_m1 = get("M100_Q0",  s, c) + (get("M100_Q100", s, c) - get("M100_Q0",  s, c)) * q
            row.append(q_m0 + (q_m1 - q_m0) * morph)
        result.append(row)
    return result


# ─── Transfer functions ───────────────────────────────────────────────────────

def stage_response(coeffs: list[float], freqs: np.ndarray) -> np.ndarray:
    """Evaluate |H_k(jw)| over freqs for one biquad [c0,c1,c2,c3,c4]."""
    w  = 2.0 * np.pi * freqs / AUTHORING_SR
    e1 = np.exp(-1j * w)
    e2 = np.exp(-2j * w)
    num = coeffs[0] + coeffs[1]*e1 + coeffs[2]*e2
    den = 1.0           + coeffs[3]*e1 + coeffs[4]*e2
    return np.abs(num / den)


def cascade_response(stages: list[list[float]], freqs: np.ndarray) -> np.ndarray:
    """Product |H(jw)| = prod_k |H_k(jw)|."""
    mag = np.ones(len(freqs), dtype=np.float64)
    for c in stages:
        mag *= stage_response(c, freqs)
    return mag


# ─── Pole / zero analysis ─────────────────────────────────────────────────────

def pole_roots(c3: float, c4: float) -> list[complex]:
    """Roots of 1 + c3*z^{-1} + c4*z^{-2} = 0.
    Equivalently roots of z^2 + c3*z + c4 = 0."""
    disc = c3*c3 - 4.0*c4
    if disc >= 0.0:
        # Real poles
        sq = math.sqrt(disc)
        return [complex((-c3 + sq)/2.0), complex((-c3 - sq)/2.0)]
    else:
        sq = math.sqrt(-disc)
        re = -c3/2.0
        return [complex(re, sq/2.0), complex(re, -sq/2.0)]


def zero_roots(c0: float, c1: float, c2: float) -> list[complex]:
    """Roots of c0 + c1*z^{-1} + c2*z^{-2} = 0.
    Equivalently roots of c2*z^2 + c1*z + c0 = 0 (if c2 != 0)
    else linear."""
    if abs(c2) < 1e-15:
        if abs(c1) < 1e-15:
            return []
        return [complex(-c0/c1)]
    disc = c1*c1 - 4.0*c2*c0
    if disc >= 0.0:
        sq = math.sqrt(disc)
        return [complex((-c1 + sq)/(2.0*c2)), complex((-c1 - sq)/(2.0*c2))]
    else:
        sq = math.sqrt(-disc)
        re = -c1/(2.0*c2)
        return [complex(re, sq/(2.0*c2)), complex(re, -sq/(2.0*c2))]


def pole_freq_hz(p: complex) -> float:
    """Angle of pole -> Hz at authoring SR."""
    angle = abs(math.atan2(p.imag, p.real))
    return angle * AUTHORING_SR / (2.0 * math.pi)


def nearest_zero_distance(poles: list[complex], zeros: list[complex]) -> float:
    """Minimum z-plane distance from any pole to any zero."""
    if not poles or not zeros:
        return float("inf")
    dists = [abs(p - z) for p in poles for z in zeros]
    return min(dists)


def is_real_pole(c3: float, c4: float) -> bool:
    return c3*c3 >= 4.0*c4


# ─── Leave-one-out ────────────────────────────────────────────────────────────

def leave_one_out_peaks(stages: list[list[float]], freqs: np.ndarray) -> list[float]:
    """For each k in 0..5, compute cascade peak with stage k removed."""
    peaks = []
    for k in range(NUM_ACTIVE):
        reduced = [s for i, s in enumerate(stages) if i != k]
        mag = cascade_response(reduced, freqs)
        peaks.append(float(np.max(mag)))
    return peaks


# ─── Per-corner |p| table for interpolation amplification test ───────────────

def corner_pole_magnitudes(
    corners: dict[str, list[list[float] | None]],
    stage_idx: int
) -> dict[str, float]:
    """Return |p| (max of |pole| over 2 roots) for stage_idx at each corner."""
    out = {}
    for label in CORNER_ORDER:
        v = corners[label][stage_idx]
        if v is None:
            out[label] = 0.0
            continue
        c3, c4 = v[3], v[4]
        ps = pole_roots(c3, c4)
        out[label] = max(abs(p) for p in ps)
    return out


# ─── Main analysis for one body ───────────────────────────────────────────────

def analyze_body(json_path: Path, morph_target: float, q_target: float) -> str:
    """Run full RCA for one cartridge at (morph_target, q_target).
    Returns a formatted report string."""
    body_name = json.loads(json_path.read_text(encoding="utf-8")).get("name", json_path.stem)
    corners   = load_corners(json_path)
    interp    = bilinear_interpolate(corners, morph_target, q_target)

    freqs     = FREQS
    full_mag  = cascade_response(interp, freqs)
    full_peak = float(np.max(full_mag))
    full_peak_freq = float(freqs[np.argmax(full_mag)])

    # Per-stage analysis
    stage_peaks   = []
    stage_peak_hz = []
    for s_idx, c in enumerate(interp):
        mag = stage_response(c, freqs)
        pk  = float(np.max(mag))
        pkf = float(freqs[np.argmax(mag)])
        stage_peaks.append(pk)
        stage_peak_hz.append(pkf)

    # Leave-one-out
    loo_peaks = leave_one_out_peaks(interp, freqs)

    # Per-stage pole/zero diagnostics
    rows = []
    for s_idx, c in enumerate(interp):
        c0, c1, c2, c3, c4 = c
        ps  = pole_roots(c3, c4)
        zs  = zero_roots(c0, c1, c2)
        p_mag      = max(abs(p) for p in ps) if ps else 0.0
        p_freq     = pole_freq_hz(ps[0]) if ps else 0.0
        z_mag      = max(abs(z) for z in zs) if zs else 0.0
        z_freq     = pole_freq_hz(zs[0]) if zs else 0.0
        nzd        = nearest_zero_distance(ps, zs)
        real_poles = is_real_pole(c3, c4)

        # Per-corner |p| for interpolation amplification test
        corner_pm  = corner_pole_magnitudes(corners, s_idx)
        max_corner_p = max(corner_pm.values())
        interp_amplifies = p_mag > max_corner_p + 1e-9

        rows.append({
            "stage":    s_idx,
            "c0": c0, "c1": c1, "c2": c2, "c3": c3, "c4": c4,
            "p_mag":    p_mag,
            "p_freq":   p_freq,
            "real_poles": real_poles,
            "z_mag":    z_mag,
            "z_freq":   z_freq,
            "nzd":      nzd,
            "peak_H":   stage_peaks[s_idx],
            "peak_H_hz": stage_peak_hz[s_idx],
            "loo_peak": loo_peaks[s_idx],
            "corner_p_max": max_corner_p,
            "interp_amplifies": interp_amplifies,
            "corner_pm": corner_pm,
        })

    # Sort by peak |H| descending for table
    rows_sorted = sorted(rows, key=lambda r: r["peak_H"], reverse=True)

    # ── Q-axis radius runaway check ──────────────────────────────────────────
    # Compare |p| at (morph=0, q=0) vs (morph=0, q=1)
    interp_q0 = bilinear_interpolate(corners, morph_target, 0.0)
    interp_q1 = bilinear_interpolate(corners, morph_target, 1.0)
    q_runaway = []
    for s_idx in range(NUM_ACTIVE):
        c_q0 = interp_q0[s_idx]; c_q1 = interp_q1[s_idx]
        ps_q0 = pole_roots(c_q0[3], c_q0[4])
        ps_q1 = pole_roots(c_q1[3], c_q1[4])
        pm_q0 = max(abs(p) for p in ps_q0) if ps_q0 else 0.0
        pm_q1 = max(abs(p) for p in ps_q1) if ps_q1 else 0.0
        q_runaway.append((s_idx, pm_q0, pm_q1, pm_q1 - pm_q0))

    # ── Band stacking check ──────────────────────────────────────────────────
    # Find stages whose pole_freq falls within 1 octave of another
    stacking_groups: list[list[int]] = []
    for i, ri in enumerate(rows):
        group = [ri["stage"]]
        for j, rj in enumerate(rows):
            if i == j: continue
            flo = min(ri["p_freq"], rj["p_freq"])
            fhi = max(ri["p_freq"], rj["p_freq"])
            if flo > 0 and fhi / flo <= 2.0:
                group.append(rj["stage"])
        if len(group) > 1 and sorted(group) not in [sorted(g) for g in stacking_groups]:
            stacking_groups.append(sorted(list(set(group))))

    # ─── Format report ────────────────────────────────────────────────────────
    lines: list[str] = []
    sep = "=" * 100

    lines.append(sep)
    lines.append(f"  RCA: {body_name}  |  analysis point: morph={morph_target:.2f}, q={q_target:.2f}")
    lines.append(sep)
    lines.append(f"  Cascade peak:  {full_peak:.6e}  at  {full_peak_freq:.1f} Hz")
    lines.append(f"  Cascade dB:    {20.0*math.log10(max(full_peak, 1e-30)):.1f} dB")
    lines.append("")

    # Per-stage compiled coefficients
    lines.append("── Interpolated coefficients at analysis point ──")
    lines.append(f"  {'stage':>5}  {'c0':>12}  {'c1':>12}  {'c2':>12}  {'c3':>12}  {'c4':>12}")
    for r in rows:
        lines.append(f"  {r['stage']:>5}  {r['c0']:>12.6f}  {r['c1']:>12.6f}  {r['c2']:>12.6f}  {r['c3']:>12.6f}  {r['c4']:>12.6f}")

    lines.append("")
    lines.append("── Pole / zero analysis ──")
    for r in rows:
        s    = r["stage"]
        rp   = "[REAL-POLES]" if r["real_poles"] else ""
        ia   = "[INTERP-AMPLIFIES]" if r["interp_amplifies"] else ""
        cpmax = r["corner_p_max"]
        cp_detail = "  ".join(f"{lbl}={v:.6f}" for lbl, v in r["corner_pm"].items())
        lines.append(
            f"  stage {s}: |p|={r['p_mag']:.6f} pole_f={r['p_freq']:.1f} Hz  "
            f"|z|={r['z_mag']:.6f} zero_f={r['z_freq']:.1f} Hz  "
            f"nearest_z_dist={r['nzd']:.4f}  {rp}{ia}"
        )
        lines.append(f"           corner |p|: {cp_detail}  max_corner={cpmax:.6f}")

    lines.append("")
    lines.append("── Q-axis radius runaway: |p| at q=0 vs q=1 (morph fixed) ──")
    lines.append(f"  {'stage':>5}  {'|p|@q=0':>10}  {'|p|@q=1':>10}  {'delta':>10}  flag")
    for s_idx, pm_q0, pm_q1, delta in q_runaway:
        flag = "<-- RUNAWAY" if pm_q1 > 0.998 and delta > 0.05 else (
               "<-- q pushes past 0.995" if pm_q1 > 0.995 else ""
        )
        lines.append(f"  {s_idx:>5}  {pm_q0:>10.6f}  {pm_q1:>10.6f}  {delta:>+10.6f}  {flag}")

    lines.append("")
    lines.append("── Band stacking groups (pole freq within 1 octave) ──")
    if stacking_groups:
        for g in stacking_groups:
            freqs_in_g = [rows[i]["p_freq"] for i in g]
            lines.append(f"  stages {g}  pole_freqs={[f'{f:.1f}' for f in freqs_in_g]} Hz")
    else:
        lines.append("  (none)")

    lines.append("")
    lines.append("── Per-stage peak |H| and leave-one-out cascade peak ──")
    lines.append(f"  {'stage':>5}  {'peak|H|':>12}  {'peak|H|dB':>10}  {'peak_f Hz':>10}  {'loo_peak':>14}  {'loo_dB':>8}  culprit?")
    for r in rows_sorted:
        loo_dB = 20.0*math.log10(max(r["loo_peak"], 1e-30))
        is_culprit = r["loo_peak"] <= 5.0
        flag = " <-- CULPRIT" if is_culprit else ""
        lines.append(
            f"  {r['stage']:>5}  {r['peak_H']:>12.4e}  {20*math.log10(max(r['peak_H'],1e-30)):>10.1f}  "
            f"{r['peak_H_hz']:>10.1f}  {r['loo_peak']:>14.4e}  {loo_dB:>8.1f}  {flag}"
        )

    # c0 attenuation vs heritage
    lines.append("")
    lines.append("── c0 values (all stages, uniform check vs heritage c0=0.5619 Talking_Hedz) ──")
    c0_vals = [r["c0"] for r in rows]
    lines.append(f"  c0 values: {[f'{v:.4f}' for v in c0_vals]}")
    lines.append(f"  mean={sum(c0_vals)/len(c0_vals):.4f}  min={min(c0_vals):.4f}  max={max(c0_vals):.4f}")
    lines.append(f"  heritage Talking_Hedz c0=0.5619  (val1=-0.43811 => 1+val1=0.56189)")
    lines.append(f"  these cartridges: 1+val1 = {1.0 + (-0.43811):.5f} for all stages (val1 uniform)")

    # Proposed repairs — emit analytical estimates
    lines.append("")
    lines.append("── Repair proposals (per culprit stage, smallest-viable) ──")
    culprits = [r for r in rows if r["loo_peak"] <= 5.0]
    if not culprits:
        # report stages that individually bring cascade below 1000 if removed
        culprits = [r for r in rows if r["loo_peak"] < full_peak * 0.01]
    if not culprits:
        culprits = rows_sorted[:3]  # top 3 by per-stage peak as fallback

    for r in culprits:
        s_idx = r["stage"]
        cur_r_q1_sq = interp_q1[s_idx][4]  # c4 at q=1 => r^2
        cur_r_q1    = math.sqrt(max(cur_r_q1_sq, 0.0))
        target_r    = 0.985  # brings |p| below 0.990 with conservative margin
        if cur_r_q1 > target_r:
            delta_r  = target_r - cur_r_q1
            # New c4 = target_r^2
            new_c4   = target_r ** 2
            # Estimate effect: per-stage peak scales roughly as 1/(1-|p|^2) for high-Q
            # old peak approx: 1/(1-cur_r_q1^2), new: 1/(1-target_r^2)
            old_approx = 1.0/(1.0 - cur_r_q1**2) if cur_r_q1 < 1.0 else float("inf")
            new_approx = 1.0/(1.0 - target_r**2) if target_r < 1.0 else float("inf")
            ratio      = new_approx / old_approx if old_approx > 0 else 1.0
            est_new_cascade = full_peak * ratio
            lines.append(
                f"  stage {s_idx}: lower r at M0_Q100 (and M100_Q100) from {cur_r_q1:.4f} -> {target_r:.4f}  "
                f"(delta_r={delta_r:+.4f})"
            )
            lines.append(
                f"    field change: c4 corner M0_Q100 = {cur_r_q1_sq:.6f} -> {new_c4:.6f}"
            )
            lines.append(
                f"    Equivalent RawStage: set r={target_r:.3f} in M0_Q100 (and M100_Q100) corner for this stage"
            )
            lines.append(
                f"    Estimated cascade peak after: ~{est_new_cascade:.3e}  "
                f"({20.0*math.log10(max(est_new_cascade,1e-30)):.1f} dB) — analytical approximation only"
            )
        else:
            # pole already below target_r; propose c0 reduction instead
            cur_c0  = r["c0"]
            target_c0 = 0.2
            ratio   = target_c0 / cur_c0 if cur_c0 > 0 else 1.0
            est_new_cascade = full_peak * (ratio ** 1)  # linear, single stage
            lines.append(
                f"  stage {s_idx}: reduce c0 from {cur_c0:.4f} -> {target_c0:.4f}  "
                f"(change val1 from {cur_c0-1.0:.5f} -> {target_c0-1.0:.5f})"
            )
            lines.append(
                f"    Estimated cascade peak after: ~{est_new_cascade:.3e} — analytical approximation only"
            )

    return "\n".join(lines)


# ─── Entry point ──────────────────────────────────────────────────────────────

def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    targets = [
        ("speaker_knockerz.json", 0.0, 1.0),
        ("cul_de_sac.json",       0.0, 1.0),
    ]

    for fname, morph, q in targets:
        path = CART_DIR / fname
        if not path.exists():
            print(f"[ERROR] not found: {path}", file=sys.stderr)
            continue
        report = analyze_body(path, morph, q)
        print(report)
        print()

        out_file = OUT_DIR / f"{path.stem}_rca.txt"
        out_file.write_text(report, encoding="utf-8")
        print(f"  -> saved: {out_file}\n")


if __name__ == "__main__":
    main()
