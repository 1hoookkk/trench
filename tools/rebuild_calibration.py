"""Rebuild calibration/ from raw RE artifacts + existing per-skin JSONs.

Produces:
  calibration/index.json          — master lookup (forge/MCP queryable)
  calibration/{Name}.json         — consolidated per-skin (coefficients + response + RE provenance)

Deletes:
  calibration/*.md
  calibration/*.csv
  calibration/INDEX.md

Preserves:
  calibration/*.png               — plots
  calibration/{Name}/S*.png       — per-stage plots
  calibration/parity_comparison.* — top-level comparison artifacts
  calibration/per_stage_breakdown.*
"""

import json
import os
import glob
import sys

CALIBRATION_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "calibration")
RE_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                      "..", "trenchwork_clean", "docs", "analysis")
RE_DIR = os.path.normpath(RE_DIR)

# RE lookup tables
SEMITONE_TABLE = None
Q_RADIUS_TABLE = None
RE_CONSTANTS = None

def load_re_artifacts():
    global SEMITONE_TABLE, Q_RADIUS_TABLE, RE_CONSTANTS
    with open(os.path.join(RE_DIR, "semitone_table.json")) as f:
        SEMITONE_TABLE = json.load(f)
    with open(os.path.join(RE_DIR, "q_radius_table.json")) as f:
        Q_RADIUS_TABLE = json.load(f)
    with open(os.path.join(RE_DIR, "constants.json")) as f:
        RE_CONSTANTS = json.load(f)


def find_nearest_semitone_index(freq_hz: float) -> int:
    """Find the nearest semitone table index for a frequency."""
    values = SEMITONE_TABLE["values"]
    best_i = 0
    best_dist = abs(freq_hz - values[0])
    for i, v in enumerate(values):
        d = abs(freq_hz - v)
        if d < best_dist:
            best_dist = d
            best_i = i
    return best_i


def find_nearest_q_index(radius: float) -> int:
    """Find the nearest Q-radius table index."""
    values = Q_RADIUS_TABLE["values"]
    best_i = 0
    best_dist = abs(radius - values[0])
    for i, v in enumerate(values):
        d = abs(radius - v)
        if d < best_dist:
            best_dist = d
            best_i = i
    return best_i


def compute_zero_placement(stage: dict) -> dict:
    """Recover zero positions from val2/val3 and classify the zero family.

    Kernel form:
        val2 = b1 - a1
        val3 = r² - b2
    So:
        b1 = val2 + a1
        b2 = r² - val3
    Zeros of numerator: 1 + b1*z^-1 + b2*z^-2
    """
    import math

    SR = 39062.5
    a1 = stage.get("val2", 0) + (-2.0 * stage.get("radius", 0) *
         math.cos(2.0 * math.pi * stage.get("pole_freq_hz", 1000) / SR))
    # Recompute a1 from pole params directly
    pole_hz = stage.get("pole_freq_hz", 0)
    radius = stage.get("radius", 0)
    theta_pole = 2.0 * math.pi * pole_hz / SR
    a1 = -2.0 * radius * math.cos(theta_pole)

    val2 = stage.get("val2", 0)
    val3 = stage.get("val3", 0)
    r_sq = radius * radius

    # Recover numerator coefficients
    b1 = val2 + a1
    b2 = r_sq - val3

    zero_info = {
        "b1": round(b1, 6),
        "b2": round(b2, 6),
    }

    # Classify zero family
    zero_energy = abs(val2) + abs(val3)
    if zero_energy < 0.05:
        zero_info["family"] = "PURE"
        zero_info["zero_r"] = 0.0
        zero_info["zero_freq_hz"] = 0.0
        return zero_info

    unit_circle_val3 = r_sq - 1.0
    if abs(val3 - unit_circle_val3) < 0.01 and abs(val2) < 0.1:
        zero_info["family"] = "UNIT_CIRCLE"
        # Zeros on unit circle: b2 ≈ 1.0
        if abs(b1) <= 2.0:
            zero_angle = math.acos(max(-1.0, min(1.0, -b1 / 2.0)))
            zero_info["zero_r"] = 1.0
            zero_info["zero_freq_hz"] = round(zero_angle * SR / (2.0 * math.pi), 1)
        else:
            zero_info["zero_r"] = 1.0
            zero_info["zero_freq_hz"] = 0.0
        return zero_info

    if val3 > 0.0:
        zero_r_sq = r_sq - val3
        if zero_r_sq > 0.0:
            zero_r = math.sqrt(zero_r_sq)
            cos_phi = (-b1 / (2.0 * zero_r)) if zero_r > 1e-6 else 1.0
            cos_phi = max(-1.0, min(1.0, cos_phi))
            zero_angle = math.acos(cos_phi)
            zero_freq = zero_angle * SR / (2.0 * math.pi)

            # Compute hue: octave offset from pole to zero
            hue = math.log2(zero_freq / pole_hz) if pole_hz > 1e-6 and zero_freq > 1e-6 else 0.0
            # Compute color: distance of zero_r from pole_r
            color = (1.0 - zero_r / radius) if radius > 1e-6 else 0.0

            zero_info["family"] = "INTERIOR_ZERO"
            zero_info["zero_r"] = round(zero_r, 6)
            zero_info["zero_freq_hz"] = round(zero_freq, 1)
            zero_info["zero_hue_octaves"] = round(hue, 3)
            zero_info["zero_color"] = round(max(0, min(1, color)), 4)
            zero_info["pole_zero_spread_hz"] = round(abs(zero_freq - pole_hz), 1)
            return zero_info

    zero_info["family"] = "NEAR_ALLPASS"
    zero_info["zero_r"] = 0.0
    zero_info["zero_freq_hz"] = 0.0
    return zero_info


def enrich_stage(stage: dict) -> dict:
    """Add RE provenance fields and zero placement to a stage dict."""
    freq = stage.get("pole_freq_hz", 0)
    radius = stage.get("radius", 0)

    st_idx = find_nearest_semitone_index(freq)
    st_freq = SEMITONE_TABLE["values"][st_idx]

    qr_idx = find_nearest_q_index(radius)
    qr_val = Q_RADIUS_TABLE["values"][qr_idx]

    stage["semitone_index"] = st_idx
    stage["semitone_freq_hz"] = round(st_freq, 2)
    stage["semitone_error_hz"] = round(freq - st_freq, 2)
    stage["q_table_index"] = qr_idx
    stage["q_table_radius"] = round(qr_val, 7)
    stage["q_table_error"] = round(radius - qr_val, 7)
    stage["zeros"] = compute_zero_placement(stage)
    return stage


def classify_strategy(name: str, corners: dict) -> str:
    """Classify topological strategy from corner data."""
    strategies = {
        "BassBox_303": "FREQUENCY_SWAP",
        "Talking_Hedz": "RADIUS_ONLY_MORPH",
        "Ear_Bender": "HF_WALL_FOLD",
        "Meaty_Gizmo": "ONE_VS_FIVE",
        "Fuzzi_Face": "RESONANCE_PHALANX",
        "Razor_Blades": "NYQUIST_SHARPENING",
        "Radio_Craze": "CONSTRUCTIVE_RADIO",
        "Early_Rizer": "FIXED_ANCHOR",
        "Freak_Shifta": "BROADBAND_DIFFUSER",
        "Millennium": "HF_WALL_FOLD",
        "Ooh_to_Eee_(approx)": "VOWEL_FORMANT",
    }
    return strategies.get(name, "UNKNOWN")


def compute_freq_spread(corner: dict) -> dict:
    """Compute frequency statistics for a corner."""
    stages = corner.get("stages", [])
    freqs = [s["pole_freq_hz"] for s in stages]
    if not freqs:
        return {}
    return {
        "freq_min_hz": round(min(freqs), 1),
        "freq_max_hz": round(max(freqs), 1),
        "freq_span_octaves": round(
            (max(freqs) / max(min(freqs), 1.0)).bit_length()
            if min(freqs) > 0 else 0, 2
        ) if min(freqs) > 0 else 0,
        "freq_centroid_hz": round(sum(freqs) / len(freqs), 1),
    }


def rebuild_skin(name: str, old_json_path: str) -> dict:
    """Rebuild a single skin JSON with RE enrichment."""
    with open(old_json_path) as f:
        old = json.load(f)

    skin = {
        "name": old["name"],
        "filter_type": old.get("filterType", 0),
        "stage_count": old.get("stageCount", 6),
        "boost": old.get("boost", 4.0),
        "sample_rate": old.get("sampleRate", 39062.5),
        "strategy": classify_strategy(name, old.get("corners", {})),
        "re_provenance": {
            "semitone_table_size": SEMITONE_TABLE["count"],
            "q_radius_table_size": Q_RADIUS_TABLE["count"],
            "decode_constants": {
                "damping": RE_CONSTANTS.get("damping_const", {}).get("value"),
                "coeff_scaler": RE_CONSTANTS.get("coeff_scaler", {}).get("value"),
                "mantissa_denorm": RE_CONSTANTS.get("mantissa_denorm", {}).get("value"),
                "mantissa_norm": RE_CONSTANTS.get("mantissa_norm", {}).get("value"),
            },
        },
        "corners": {},
    }

    for corner_name, corner_data in old.get("corners", {}).items():
        stages = []
        for s in corner_data.get("stages", []):
            stages.append(enrich_stage(dict(s)))

        import math
        freqs = [s["pole_freq_hz"] for s in stages if s["pole_freq_hz"] > 0]
        freq_min = min(freqs) if freqs else 0
        freq_max = max(freqs) if freqs else 0

        skin["corners"][corner_name] = {
            "stages": stages,
            "cascade_peak_db": corner_data.get("cascade_peak_db"),
            "max_single_stage_db": corner_data.get("max_single_stage_db"),
            "cancellation_db": corner_data.get("cancellation_db"),
            "summary": {
                "freq_min_hz": round(freq_min, 1),
                "freq_max_hz": round(freq_max, 1),
                "freq_span_octaves": round(
                    math.log2(freq_max / freq_min) if freq_min > 0 else 0, 2
                ),
                "mean_radius": round(
                    sum(s["radius"] for s in stages) / len(stages) if stages else 0, 6
                ),
                "shared_c4": stages[0].get("c4_b0") if stages else None,
                "zero_families": _zero_family_summary(stages),
            },
        }

    return skin


def _zero_family_summary(stages: list[dict]) -> dict:
    """Summarize zero placement across all stages in a corner."""
    from collections import Counter
    families = Counter()
    zero_freqs = []
    pole_zero_spreads = []
    for s in stages:
        z = s.get("zeros", {})
        fam = z.get("family", "UNKNOWN")
        families[fam] += 1
        if z.get("zero_freq_hz", 0) > 0:
            zero_freqs.append(z["zero_freq_hz"])
        if "pole_zero_spread_hz" in z:
            pole_zero_spreads.append(z["pole_zero_spread_hz"])

    return {
        "family_counts": dict(families),
        "dominant_family": families.most_common(1)[0][0] if families else "NONE",
        "zero_freq_min_hz": round(min(zero_freqs), 1) if zero_freqs else None,
        "zero_freq_max_hz": round(max(zero_freqs), 1) if zero_freqs else None,
        "mean_pole_zero_spread_hz": round(
            sum(pole_zero_spreads) / len(pole_zero_spreads), 1
        ) if pole_zero_spreads else None,
    }


def build_index(skins: dict) -> dict:
    """Build the master index from all rebuilt skins."""
    entries = []
    for name, skin in sorted(skins.items()):
        m0q0 = skin["corners"].get("M0_Q0", {})
        m100q0 = skin["corners"].get("M100_Q0", {})

        m0q0_stages = m0q0.get("stages", [])
        m100q0_stages = m100q0.get("stages", [])

        entry = {
            "name": skin["name"],
            "file": f"{name}.json",
            "filter_type": skin["filter_type"],
            "strategy": skin["strategy"],
            "stage_count": skin["stage_count"],
            "boost": skin["boost"],
            "m0q0_peak_db": m0q0.get("cascade_peak_db"),
            "m0q0_cancellation_db": m0q0.get("cancellation_db"),
            "m0q0_shared_c4": m0q0.get("summary", {}).get("shared_c4"),
            "m0q0_freq_range": [
                m0q0.get("summary", {}).get("freq_min_hz"),
                m0q0.get("summary", {}).get("freq_max_hz"),
            ],
            "morph_freq_delta": None,
        }

        # Compute morph frequency delta (how much poles move from M0 to M100)
        if m0q0_stages and m100q0_stages and len(m0q0_stages) == len(m100q0_stages):
            deltas = []
            for s0, s1 in zip(m0q0_stages, m100q0_stages):
                deltas.append(abs(s1["pole_freq_hz"] - s0["pole_freq_hz"]))
            entry["morph_freq_delta"] = round(sum(deltas) / len(deltas), 1)

        entries.append(entry)

    return {
        "version": "calibration-v2",
        "count": len(entries),
        "missing": [
            {"filter_type": 34, "name": "Boland Bass", "priority": "HIGH",
             "reason": "No bass-synthesis calibration anchor"},
            {"filter_type": 43, "name": "Eeh to Aah", "priority": "LOW",
             "reason": "Reverse vowel morph, adjacent to Ooh_to_Eee"},
            {"filter_type": 50, "name": "Acid Ravage", "priority": "LOW",
             "reason": "Aggressive acid, adjacent to Razor Blades"},
            {"filter_type": 51, "name": "Bassomatic", "priority": "HIGH",
             "reason": "Distinct bass strategy without calibration"},
            {"filter_type": 52, "name": "Lucifer's Q", "priority": "CRITICAL",
             "reason": "Only Q-dominant body, no oracle for Q-axis lead behavior"},
        ],
        "re_tables": {
            "semitone_table": SEMITONE_TABLE,
            "q_radius_table_size": Q_RADIUS_TABLE["count"],
        },
        "skins": entries,
    }


def main():
    load_re_artifacts()

    # Find all existing skin JSONs
    skin_jsons = glob.glob(os.path.join(CALIBRATION_DIR, "*.json"))
    skin_jsons = [p for p in skin_jsons if os.path.basename(p) not in ("index.json",)]

    skins = {}
    for path in sorted(skin_jsons):
        name = os.path.splitext(os.path.basename(path))[0]
        print(f"  Rebuilding {name}...")
        skins[name] = rebuild_skin(name, path)

    # Write rebuilt per-skin JSONs
    for name, skin in skins.items():
        out_path = os.path.join(CALIBRATION_DIR, f"{name}.json")
        with open(out_path, "w") as f:
            json.dump(skin, f, indent=2)
        print(f"  Wrote {out_path}")

    # Write master index
    index = build_index(skins)
    index_path = os.path.join(CALIBRATION_DIR, "index.json")
    with open(index_path, "w") as f:
        json.dump(index, f, indent=2)
    print(f"  Wrote {index_path}")

    # Delete redundant files
    deleted = 0
    for pattern in ("*.md", "*.csv"):
        for path in glob.glob(os.path.join(CALIBRATION_DIR, pattern)):
            os.remove(path)
            deleted += 1
            print(f"  Deleted {os.path.basename(path)}")

    print(f"\nDone. {len(skins)} skins rebuilt, {deleted} redundant files removed.")


if __name__ == "__main__":
    main()
