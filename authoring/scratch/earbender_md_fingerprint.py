"""Does Ear Bender look MD-shaped, or arbitrary?

Loads the RE'd 4-corner calibration and checks structural fingerprints
that an MD Type 1/2/3 authored filter MUST satisfy:

  1. All pole radii fall on the 512-entry Q table (within small epsilon)
  2. c4 (damping) is shared across all 6 stages at a given corner
  3. Low/mid pole frequencies fit the 64-entry semitone table
  4. HF pole frequencies either fit the table or pin at index 63

If all four hold, Ear Bender is MD-expressible. If any fail, it lives
outside the vocab and needs direct coefficient ingest.
"""
import json
from pathlib import Path

CAL = Path("docs/calibration/Ear_Bender.json")

def analyze_corner(label, corner):
    stages = corner["stages"]
    c4s = [s["c4_b0"] for s in stages]
    shared_c4_spread = max(c4s) - min(c4s)

    q_errors = [abs(s["q_table_error"]) for s in stages]
    max_q_err = max(q_errors)
    mean_q_err = sum(q_errors) / len(q_errors)

    pinned = [i for i, s in enumerate(stages) if s["semitone_index"] == 63]
    unpinned = [i for i in range(6) if i not in pinned]
    sem_errs_unpinned = [abs(stages[i]["semitone_error_hz"]) for i in unpinned]
    max_sem_err_unpinned = max(sem_errs_unpinned) if sem_errs_unpinned else 0.0

    return {
        "label": label,
        "shared_c4_spread": shared_c4_spread,
        "max_q_err": max_q_err,
        "mean_q_err": mean_q_err,
        "pinned_stages_at_idx63": pinned,
        "max_semitone_err_unpinned_hz": max_sem_err_unpinned,
        "pole_freqs_hz": [round(s["pole_freq_hz"], 1) for s in stages],
        "radii": [round(s["radius"], 4) for s in stages],
        "c4": stages[0]["c4_b0"],
    }


def verdict(corners):
    all_shared_c4 = all(c["shared_c4_spread"] < 1e-4 for c in corners.values())
    all_q_fit = all(c["max_q_err"] < 0.015 for c in corners.values())
    low_mid_fit = all(c["max_semitone_err_unpinned_hz"] < 60.0 for c in corners.values())
    pinning_rate = sum(len(c["pinned_stages_at_idx63"]) for c in corners.values()) / (4 * 6)

    return {
        "shared_c4_per_corner": all_shared_c4,
        "radii_on_q_table": all_q_fit,
        "low_mid_poles_on_semitone_grid": low_mid_fit,
        "fraction_stages_pinned_at_top_of_semitone_table": pinning_rate,
    }


def main():
    data = json.loads(CAL.read_text())
    corners = {k: analyze_corner(k, v) for k, v in data["corners"].items()}

    print("=== Ear Bender MD-shape fingerprint ===\n")
    for label, c in corners.items():
        print(f"[{label}]")
        print(f"  shared c4 = {c['c4']:.4f} (spread across stages: {c['shared_c4_spread']:.2e})")
        print(f"  Q-table error: max {c['max_q_err']:.4f}, mean {c['mean_q_err']:.4f}")
        print(f"  pole freqs (Hz): {c['pole_freqs_hz']}")
        print(f"  radii: {c['radii']}")
        print(f"  pinned at semitone_index=63: stages {c['pinned_stages_at_idx63']}")
        print(f"  max semitone err (unpinned): {c['max_semitone_err_unpinned_hz']:.1f} Hz")
        print()

    v = verdict(corners)
    print("=== Verdict ===")
    for k, val in v.items():
        print(f"  {k}: {val}")

    print()
    passed = sum(1 for k in ("shared_c4_per_corner", "radii_on_q_table", "low_mid_poles_on_semitone_grid") if v[k])
    if passed == 3:
        if v["fraction_stages_pinned_at_top_of_semitone_table"] > 0.2:
            print(">>> MD-shaped on radii + low/mid freqs, but HF stages pin at index 63.")
            print(">>> MD vocab PROBABLY reaches this, with a caveat: the 64-entry semitone")
            print(">>> table does not encode the HF pole positions exactly. Two possibilities:")
            print(">>>   (a) MD freq-byte (0-127) has a higher-freq decode path beyond the 64")
            print(">>>       semitone table (likely, since byte is 128-valued).")
            print(">>>   (b) Ear Bender's HF stages live outside MD's reach.")
            print(">>> Resolve: hand-author an MD recipe and null-test.")
        else:
            print(">>> Fully MD-shaped. Vocab reaches Ear Bender.")
    else:
        print(">>> NOT MD-shaped on at least one structural dimension.")
        print(">>> Ear Bender is ROM-native, needs direct coefficient ingest path.")


if __name__ == "__main__":
    main()
