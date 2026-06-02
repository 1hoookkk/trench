"""Author the 7 Speaker Knockerz pills via the v2 pill template math.

Reads the Speaker Knockerz station designs below, computes direct
pole-zero coefficients (Klatt DC-normalized resonators + corrected
notches + 1-pole HPF/LPF), runs a cascade frequency-response audit
against the BODIES.md rubric, and emits one compiled-v1 JSON file per
station to:

    cartridges/engine/_source/shapes/speaker_knockerz/<key>.json

Also prints inventory entries (for hand-pasting into
cartridges/engine/_source/token_inventory_unified_v2.json) and pill
labels (for hand-pasting into tools/bake_phoneme_pills.py).

Audit criteria (BODIES.md lines 36-42):
  - Sub at 50 Hz must be within 3 dB of the median sub level across
    all 7 pills (invariant: "fundamental below 60 Hz stays anchored").
  - sk_choke must have notch depth >= 12 dB at 200 Hz.
  - sk_cry must have a peak at 1.2 kHz that exceeds flanking bands by
    >= 10 dB.

Stdlib only. No numpy.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

REPO = Path("/home/user/trench")
SHAPES_DIR = REPO / "cartridges" / "engine" / "_source" / "shapes" / "speaker_knockerz"
SR = 39062.5
N_BINS = 4096


# ---------------------------------------------------------------------------
# v2 template math — direct pole-zero placement, Klatt DC normalization
# ---------------------------------------------------------------------------

def passthrough() -> list:
    return [1.0, 0.0, 0.0, 0.0, 0.0]


def resonator(f: float, bw: float) -> list:
    """2-pole / 0-zero Klatt DC-normalized resonator.

    Unity gain at DC, peak > 1 at `f`. Cascade-safe: each stage passes
    DC through unchanged so cascading multiple resonators does not
    compound off-peak rolloffs.
    """
    theta = 2.0 * math.pi * f / SR
    r = math.exp(-math.pi * bw / SR)
    c3 = -2.0 * r * math.cos(theta)
    c4 = r * r
    c0 = 1.0 + c3 + c4
    return [c0, 0.0, 0.0, c3, c4]


def notch(f: float, bw: float, d: float) -> list:
    """2-pole / 2-zero corrected notch.

    Zeros placed BETWEEN pole radius and unit circle (rz = r + d·(1-r))
    to produce a genuine dip rather than a resonant bump. `d` in (0,1]:
    0 = pole-zero cancel (flat), 1 = zero on unit circle (infinite
    notch). Use 0.8-0.95 for musical -18 to -30 dB notches.
    """
    theta = 2.0 * math.pi * f / SR
    r = math.exp(-math.pi * bw / SR)
    rz = r + d * (1.0 - r)
    cos_t = math.cos(theta)
    c0 = 1.0
    c1 = -2.0 * rz * cos_t
    c2 = rz * rz
    c3 = -2.0 * r * cos_t
    c4 = r * r
    return [c0, c1, c2, c3, c4]


def lpf1(fc: float) -> list:
    """1-pole lowpass, unity DC gain, -3 dB at fc."""
    omega = 2.0 * math.pi * fc / SR
    alpha = math.exp(-omega)
    return [1.0 - alpha, 0.0, 0.0, -alpha, 0.0]


def hpf1(fc: float) -> list:
    """1-pole highpass, zero DC gain, -3 dB at fc."""
    omega = 2.0 * math.pi * fc / SR
    alpha = math.exp(-omega)
    half = 0.5 * (1.0 + alpha)
    return [half, -half, 0.0, -alpha, 0.0]


# ---------------------------------------------------------------------------
# Speaker Knockerz station designs
# Each pill is 6 active stages (padded to 12 with passthrough) matching
# one of the seven morph-path stations in BODIES.md line 28-34.
# ---------------------------------------------------------------------------

def build_sk_vault() -> list:
    """The Vault — clean sub, heavy weight, clamped highs."""
    return [
        resonator(80, 60),    # low weight
        resonator(180, 120),  # low-mid body
        resonator(300, 200),  # mid warmth
        passthrough(),
        lpf1(500),            # hard clamp
        lpf1(500),            # second 1-pole for 12 dB/oct total
    ]


def build_sk_chest() -> list:
    """Chest Resonance — sub blooms into a warm 150 Hz boost."""
    return [
        resonator(150, 40),   # narrow chest boost (the defining feature)
        resonator(400, 200),  # mid warmth
        passthrough(),
        passthrough(),
        lpf1(1000),           # mild clamp
        passthrough(),
    ]


def build_sk_choke() -> list:
    """The Choke — steep notch at 200 Hz, fundamental separates from
    overtones. Rubric: notch must be >= 12 dB at its center."""
    return [
        notch(200, 80, 0.92), # THE choke (~-22 dB by design)
        resonator(120, 60),   # fill below the notch
        passthrough(),
        passthrough(),
        passthrough(),
        lpf1(1500),           # cap
    ]


def build_sk_rip() -> list:
    """Cardboard Rip — peaky unstable resonance flares at 400 Hz."""
    return [
        resonator(400, 35),   # THE rip — Q ~11, strong flare
        resonator(180, 150),  # low support
        passthrough(),
        passthrough(),
        passthrough(),
        lpf1(1800),           # cap
    ]


def build_sk_rattle() -> list:
    """The Rattle — comb-like flutter in the 800 Hz range.

    Three resonators at 700/820/940 Hz create an interference
    pattern that reads as rattle rather than a single formant. The
    bandwidths are deliberately wide (Q ~5) so the individual peaks
    are modest — the rattle character comes from their interference
    pattern, not from each one screaming. Violates the v2 template's
    "poles > 100 Hz apart" rule on purpose, as the interference IS
    the design."""
    return [
        resonator(700, 160),  # comb tooth 1 — Q ~4
        resonator(820, 160),  # comb tooth 2
        resonator(940, 160),  # comb tooth 3
        passthrough(),
        passthrough(),
        lpf1(2000),           # cap
    ]


def build_sk_cry() -> list:
    """Cone Cry — narrow screaming peak at 1.2 kHz. Rubric: must
    exceed flanking bands by >= 10 dB."""
    return [
        resonator(1200, 80),  # THE cry — Q ~15, narrow but physical
        passthrough(),
        passthrough(),
        passthrough(),
        passthrough(),
        lpf1(2500),           # cap above cry
    ]


def build_sk_fracture() -> list:
    """Total Fracture — mid-band phase-collapses, pure sub + violent top."""
    return [
        notch(350, 250, 0.90),   # kill low-mids
        notch(700, 250, 0.90),   # kill high-mids
        resonator(3200, 250),    # violent top
        passthrough(),
        passthrough(),
        lpf1(5000),              # soft cap
    ]


STATIONS = [
    ("sk_vault",    "The Vault",       "clean sub, heavy weight, clamped highs",                      "VLT", build_sk_vault),
    ("sk_chest",    "Chest Resonance", "sub blooms into a warm 150 Hz boost",                         "CST", build_sk_chest),
    ("sk_choke",    "The Choke",       "steep notch at 200 Hz; fundamental separates from overtones", "CHK", build_sk_choke),
    ("sk_rip",      "Cardboard Rip",   "peaky unstable resonance flares at 400 Hz",                   "RIP", build_sk_rip),
    ("sk_rattle",   "The Rattle",      "comb-like flutter in the 800 Hz range",                       "RAT", build_sk_rattle),
    ("sk_cry",      "Cone Cry",        "narrow screaming peak at 1.2 kHz",                            "CRY", build_sk_cry),
    ("sk_fracture", "Total Fracture",  "mid-band phase-collapses; pure sub plus violent top",         "FRC", build_sk_fracture),
]


# ---------------------------------------------------------------------------
# Cascade frequency response + audit
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
        return 1e9
    return num_mag / den_mag


def cascade_mag(stages: list, f: float) -> float:
    omega = 2.0 * math.pi * f / SR
    m = 1.0
    for stage in stages:
        m *= biquad_mag(stage, omega)
    return m


def cascade_response(stages: list, n_bins: int) -> list:
    out = [0.0] * n_bins
    for i in range(n_bins):
        omega = math.pi * i / (n_bins - 1)
        m = 1.0
        for stage in stages:
            m *= biquad_mag(stage, omega)
        out[i] = m
    return out


def db(x: float) -> float:
    if x <= 0:
        return -999.0
    return 20.0 * math.log10(x)


def cascade_peak(stages: list) -> float:
    """Find the maximum magnitude of the cascade over the audible band
    so we can calibrate boost to ~unity output."""
    response = cascade_response(stages, N_BINS)
    peak = 0.0
    for i, m in enumerate(response):
        hz = i * (SR / 2.0) / (N_BINS - 1)
        if 40.0 <= hz <= 8000.0:
            if m > peak:
                peak = m
    return peak


# ---------------------------------------------------------------------------
# Emit compiled-v1 JSON for one pill
# ---------------------------------------------------------------------------

def pill_json(name: str, stages: list, boost: float) -> dict:
    # Pad to 12 stages (6 active + 6 passthrough)
    padded = list(stages)
    while len(padded) < 12:
        padded.append(passthrough())
    kf_stages = [
        {"c0": s[0], "c1": s[1], "c2": s[2], "c3": s[3], "c4": s[4]}
        for s in padded
    ]
    # Morph-invariant: all 4 corners identical
    keyframes = [
        {"label": "M0_Q0",    "morph": 0.0, "q": 0.0, "boost": boost, "stages": kf_stages},
        {"label": "M100_Q0",  "morph": 1.0, "q": 0.0, "boost": boost, "stages": kf_stages},
        {"label": "M0_Q100",  "morph": 0.0, "q": 1.0, "boost": boost, "stages": kf_stages},
        {"label": "M100_Q100","morph": 1.0, "q": 1.0, "boost": boost, "stages": kf_stages},
    ]
    return {
        "format": "compiled-v1",
        "name": name,
        "provenance": "hand-authored via v2 template, Speaker Knockerz station",
        "sampleRate": SR,
        "stages": 12,
        "keyframes": keyframes,
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    SHAPES_DIR.mkdir(parents=True, exist_ok=True)

    built: list[tuple[str, str, str, str, list, float]] = []
    for key, title, desc, pill, builder in STATIONS:
        stages = builder()
        # Calibrate boost so the sub at 50 Hz plays at unity.
        # This enforces the Speaker Knockerz invariant: the 50 Hz
        # fundamental passes through each pill at the same level,
        # and character peaks rise above that flat sub floor as the
        # pill's "personality". Matches real-world speaker-cone
        # behavior: flat sub with one dominant resonance on top.
        sub_50 = cascade_mag(stages, 50.0)
        boost = 1.0 / sub_50 if sub_50 > 0 else 1.0
        built.append((key, title, desc, pill, stages, boost))

    # ---- audit phase ----
    print(f"\n=== Speaker Knockerz pill audit (v2 template, SR={SR} Hz) ===\n")
    print(f"{'key':<13} {'pill':<5} {'boost':>8}   sub@50Hz   peak@?Hz  (dB rel to sub)")
    print(f"{'-'*13} {'-'*5} {'-'*8}   {'-'*8}   {'-'*20}")

    sub_levels_norm = []
    for key, title, desc, pill, stages, boost in built:
        sub_50 = cascade_mag(stages, 50.0) * boost
        sub_levels_norm.append((key, sub_50))

        # Find the strongest peak in the full band
        response = cascade_response(stages, N_BINS)
        peak_idx = 0
        peak_mag = 0.0
        for i, m in enumerate(response):
            hz = i * (SR / 2.0) / (N_BINS - 1)
            if 40.0 <= hz <= 8000.0 and m > peak_mag:
                peak_mag = m
                peak_idx = i
        peak_hz = peak_idx * (SR / 2.0) / (N_BINS - 1)
        peak_out = peak_mag * boost

        rel_sub_to_peak = db(sub_50) - db(peak_out)
        print(f"{key:<13} {pill:<5} {boost:>8.4f}   {db(sub_50):>+7.2f}    {peak_hz:>6.0f} Hz  ({rel_sub_to_peak:+5.1f} dB)")

    # ---- invariant check: sub level consistency ----
    print(f"\nInvariant check (sub at 50 Hz, normalized):")
    sub_dbs = [db(x) for _, x in sub_levels_norm]
    median_sub = sorted(sub_dbs)[len(sub_dbs) // 2]
    print(f"  median sub-50 level: {median_sub:+.2f} dB")
    worst_delta = max(abs(db_val - median_sub) for db_val in sub_dbs)
    print(f"  worst pill delta from median: {worst_delta:.2f} dB")
    if worst_delta > 3.0:
        print(f"  ** FAIL: sub drops more than 3 dB from median (BODIES.md rubric) **")
        worst_key, worst_val = max(sub_levels_norm, key=lambda p: abs(db(p[1]) - median_sub))
        print(f"     worst offender: {worst_key} at {db(worst_val):+.2f} dB")
    else:
        print(f"  OK: all sub-50 levels within 3 dB of median")

    # ---- station-specific rubric checks ----
    print(f"\nStation-specific rubric checks:")

    # sk_choke: notch depth at 200 Hz >= 12 dB
    choke = next(b for b in built if b[0] == "sk_choke")
    _, _, _, _, choke_stages, choke_boost = choke
    choke_200 = cascade_mag(choke_stages, 200.0) * choke_boost
    # Depth is relative to the sub level (or to unity)
    choke_sub = cascade_mag(choke_stages, 50.0) * choke_boost
    choke_depth_db = db(choke_sub) - db(choke_200)
    print(f"  sk_choke notch at 200 Hz:")
    print(f"    cascade gain at 50 Hz  = {db(choke_sub):+.2f} dB")
    print(f"    cascade gain at 200 Hz = {db(choke_200):+.2f} dB")
    print(f"    notch depth (sub - 200) = {choke_depth_db:.2f} dB")
    if choke_depth_db >= 12.0:
        print(f"    OK: depth >= 12 dB (rubric satisfied)")
    else:
        print(f"    ** FAIL: depth < 12 dB (BODIES.md rubric) **")

    # sk_cry: peak at 1200 Hz exceeds flanking bands by >= 10 dB
    cry = next(b for b in built if b[0] == "sk_cry")
    _, _, _, _, cry_stages, cry_boost = cry
    cry_peak = cascade_mag(cry_stages, 1200.0) * cry_boost
    cry_flank_low = cascade_mag(cry_stages, 600.0) * cry_boost
    cry_flank_high = cascade_mag(cry_stages, 2400.0) * cry_boost
    excess_low = db(cry_peak) - db(cry_flank_low)
    excess_high = db(cry_peak) - db(cry_flank_high)
    print(f"  sk_cry peak at 1200 Hz:")
    print(f"    cascade gain at  600 Hz = {db(cry_flank_low):+.2f} dB")
    print(f"    cascade gain at 1200 Hz = {db(cry_peak):+.2f} dB")
    print(f"    cascade gain at 2400 Hz = {db(cry_flank_high):+.2f} dB")
    print(f"    excess over low flank:  {excess_low:+.2f} dB")
    print(f"    excess over high flank: {excess_high:+.2f} dB")
    if excess_low >= 10.0 and excess_high >= 10.0:
        print(f"    OK: peak exceeds both flanking bands by >= 10 dB")
    else:
        print(f"    ** FAIL: peak does not exceed both flanking bands by >= 10 dB **")

    # ---- stability sweep ----
    print(f"\nStability check (every active stage, all 7 pills):")
    bad_stages = 0
    for key, _, _, _, stages, _ in built:
        for i, s in enumerate(stages):
            if abs(s[3]) >= 2.0 or s[4] >= 1.0:
                print(f"  {key} stage {i+1}: |c3|={abs(s[3]):.4f} c4={s[4]:.4f} — UNSTABLE")
                bad_stages += 1
            if s == passthrough():
                continue
    if bad_stages == 0:
        print(f"  OK: all {sum(len(b[4]) for b in built)} active+passthrough stages stable")

    # ---- emit JSON files ----
    print(f"\nEmitting JSON files:")
    for key, title, desc, pill, stages, boost in built:
        pill_name = f"speaker_knockerz_{key}"
        out_path = SHAPES_DIR / f"{key}.json"
        out_path.write_text(
            json.dumps(pill_json(pill_name, stages, boost), indent=2) + "\n",
            encoding="utf-8",
        )
        print(f"  wrote {out_path.relative_to(REPO)}")

    # ---- inventory snippet for copy-paste ----
    print(f"\n=== INVENTORY ENTRIES (paste into token_inventory_unified_v2.json) ===\n")
    entries = {}
    for key, title, desc, pill, _, _ in built:
        entries[f"speaker_knockerz.{key}"] = {
            "tables_ref": f"speaker_knockerz.{key}",
            "category": "speaker_knockerz",
            "key": key,
            "label": f"Speaker Knockerz / {title} ({desc})",
            "shape_path": f"cartridges/engine/_source/shapes/speaker_knockerz/{key}.json",
            "p2k_phonemes": [],
        }
    print(json.dumps(entries, indent=2))

    # ---- pill label snippet for copy-paste ----
    print(f"\n=== PILL_LABELS ENTRIES (paste into tools/bake_phoneme_pills.py) ===\n")
    for key, _, _, pill, _, _ in built:
        print(f'    "speaker_knockerz.{key}": "{pill}",')

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
