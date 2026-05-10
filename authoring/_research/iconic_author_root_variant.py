"""forge/iconic_author.py — produce RawStage cartridge JSONs for the launch bodies.

Usage:
  python forge/iconic_author.py <body_id> [--out PATH]

body_id ∈ {speaker_knockerz, aluminum_siding, small_talk_ah_ee, cul_de_sac}
"""

import sys
import json
import math
import argparse
import numpy as np

AUTHORING_SR = 39062.5

def biquad_mag(c0, c1, c2, c3, c4, freqs, sr=AUTHORING_SR):
    """Evaluate magnitude response of a single DF2T biquad."""
    w = 2.0 * np.pi * freqs / sr
    e1 = np.exp(-1j * w)
    e2 = np.exp(-2j * w)
    num = c0 + c1 * e1 + c2 * e2
    den = 1.0 + c3 * e1 + c4 * e2
    return np.maximum(np.abs(num / den), 1e-12)

def cascade_db(stages, freqs, sr=AUTHORING_SR):
    """Evaluate multiplied cascade magnitude in dB."""
    mag = np.ones_like(freqs, dtype=np.float64)
    for s in stages:
        mag *= biquad_mag(s["c0"], s["c1"], s["c2"], s["c3"], s["c4"], freqs, sr)
    return 20.0 * np.log10(mag)

def make_coefs(pole_hz, pole_r, zero_hz, zero_r, c0, sr=AUTHORING_SR):
    """Generate intermediate c0..c4 for bisection + RawStage conversion."""
    theta_p = 2.0 * math.pi * pole_hz / sr
    theta_z = 2.0 * math.pi * zero_hz / sr
    
    a1 = -2.0 * pole_r * math.cos(theta_p)
    a2 = pole_r * pole_r
    
    b1 = -2.0 * zero_r * math.cos(theta_z)
    b2 = zero_r * zero_r
    
    return {
        "c0": c0,
        "c1": c0 * b1,
        "c2": c0 * b2,
        "c3": a1,
        "c4": a2,
        "a1": a1,
        "a2": a2,
        "r": pole_r,
        "p_hz": pole_hz
    }

def calibrate_c0(stage_blueprints, target_db, freqs, sr=AUTHORING_SR):
    """Bisect c0 in (0, 1.0] to hit the cascade peak target independently per corner."""
    lo, hi = 1e-6, 1.0
    best_c0 = 1.0
    best_peak = 0.0
    
    for _ in range(60):
        mid = math.sqrt(lo * hi)
        stages = [make_coefs(bp["p_hz"], bp["p_r"], bp["z_hz"], bp["z_r"], mid, sr) for bp in stage_blueprints]
        db = cascade_db(stages, freqs, sr)
        peak = float(np.max(db))
        
        best_c0 = mid
        best_peak = peak
        
        if peak > target_db:
            hi = mid
        else:
            lo = mid
            
    return best_c0, best_peak, stages

def to_raw_stage(c_dict):
    """Convert intermediate coefs to EMU RawStage format for host-SR recompilation."""
    c0 = c_dict["c0"]
    a1 = c_dict["a1"]
    a2 = c_dict["a2"]
    r = c_dict["r"]
    c1 = c_dict["c1"]
    c2 = c_dict["c2"]
    
    return {
        "a1": a1,
        "r": min(r, 0.9999), 
        "val1": c0 - 1.0,      # Attenuation
        "val2": c1 - a1,       # Num offset
        "val3": a2 - c2,       # Num offset (r^2 - c2)
        "flag": 1
    }

def get_body_archetype(body_id):
    """Defines the stage layouts (M0, M100) based on heritage archetypes."""
    if body_id == "small_talk_ah_ee":
        # Archetype: Talking_Hedz (Vocal multi-formant, zeros in HF basket ~4.5k-8.5k)
        return {
            "target_peak": 18.0,
            "M0": [ # Yawn (Oo)
                {"p_hz": 300, "p_r": 0.985, "z_hz": 4500, "z_r": 0.70},
                {"p_hz": 800, "p_r": 0.985, "z_hz": 5200, "z_r": 0.70},
                {"p_hz": 1500, "p_r": 0.970, "z_hz": 6000, "z_r": 0.70},
                {"p_hz": 2500, "p_r": 0.960, "z_hz": 7000, "z_r": 0.70},
                {"p_hz": 3500, "p_r": 0.950, "z_hz": 8000, "z_r": 0.70},
                {"p_hz": 5000, "p_r": 0.940, "z_hz": 9000, "z_r": 0.70},
            ],
            "M100": [ # Shriek (Aaa/Ee bite)
                {"p_hz": 600, "p_r": 0.985, "z_hz": 4500, "z_r": 0.70},
                {"p_hz": 1200, "p_r": 0.985, "z_hz": 5200, "z_r": 0.70},
                {"p_hz": 2600, "p_r": 0.990, "z_hz": 6000, "z_r": 0.70},
                {"p_hz": 3200, "p_r": 0.970, "z_hz": 7000, "z_r": 0.70},
                {"p_hz": 4500, "p_r": 0.960, "z_hz": 8000, "z_r": 0.70},
                {"p_hz": 6000, "p_r": 0.950, "z_hz": 9000, "z_r": 0.70},
            ]
        }
    elif body_id == "aluminum_siding":
        # Archetype: Ear_Bender (HF comb, 1kHz scoop)
        return {
            "target_peak": 22.0,
            "M0": [ # Dull Silver
                {"p_hz": 200, "p_r": 0.950, "z_hz": 1000, "z_r": 0.90}, # Deep scoop at 1k
                {"p_hz": 4000, "p_r": 0.980, "z_hz": 1000, "z_r": 0.90},
                {"p_hz": 6000, "p_r": 0.980, "z_hz": 1100, "z_r": 0.90},
                {"p_hz": 8000, "p_r": 0.980, "z_hz": 12000, "z_r": 0.95},
                {"p_hz": 10000, "p_r": 0.980, "z_hz": 14000, "z_r": 0.95},
                {"p_hz": 12000, "p_r": 0.970, "z_hz": 16000, "z_r": 0.95},
            ],
            "M100": [ # Shatter Point
                {"p_hz": 200, "p_r": 0.900, "z_hz": 1000, "z_r": 0.92},
                {"p_hz": 5000, "p_r": 0.990, "z_hz": 1000, "z_r": 0.92},
                {"p_hz": 7000, "p_r": 0.995, "z_hz": 1100, "z_r": 0.92},
                {"p_hz": 8000, "p_r": 0.990, "z_hz": 13000, "z_r": 0.90},
                {"p_hz": 12000, "p_r": 0.995, "z_hz": 14000, "z_r": 0.90}, # Alum Tear
                {"p_hz": 18000, "p_r": 0.990, "z_hz": 16000, "z_r": 0.90},
            ]
        }
    elif body_id == "speaker_knockerz":
        # Archetype: Razor_Blades (Asymmetric HF/LF, sub anchors)
        return {
            "target_peak": 20.0,
            "M0": [ # Vault
                {"p_hz": 55, "p_r": 0.995, "z_hz": 200, "z_r": 0.92}, # Sub anchored, choke notch
                {"p_hz": 150, "p_r": 0.980, "z_hz": 300, "z_r": 0.85},
                {"p_hz": 400, "p_r": 0.970, "z_hz": 600, "z_r": 0.80},
                {"p_hz": 800, "p_r": 0.960, "z_hz": 1000, "z_r": 0.80},
                {"p_hz": 1200, "p_r": 0.950, "z_hz": 2000, "z_r": 0.70},
                {"p_hz": 3000, "p_r": 0.940, "z_hz": 4000, "z_r": 0.70},
            ],
            "M100": [ # Total Fracture
                {"p_hz": 55, "p_r": 0.995, "z_hz": 200, "z_r": 0.95}, 
                {"p_hz": 400, "p_r": 0.985, "z_hz": 600, "z_r": 0.80},
                {"p_hz": 800, "p_r": 0.990, "z_hz": 1000, "z_r": 0.80},
                {"p_hz": 1200, "p_r": 0.995, "z_hz": 2000, "z_r": 0.85}, # Cone cry
                {"p_hz": 4000, "p_r": 0.980, "z_hz": 5000, "z_r": 0.70},
                {"p_hz": 8000, "p_r": 0.970, "z_hz": 10000, "z_r": 0.70},
            ]
        }
    elif body_id == "cul_de_sac":
        # Archetype: Millennium (Harmonic ladder, zero octave below)
        return {
            "target_peak": 25.0,
            "M0": [ # Iron Pipe
                {"p_hz": 100, "p_r": 0.990, "z_hz": 50, "z_r": 0.85},
                {"p_hz": 200, "p_r": 0.985, "z_hz": 100, "z_r": 0.85},
                {"p_hz": 400, "p_r": 0.980, "z_hz": 200, "z_r": 0.85},
                {"p_hz": 800, "p_r": 0.975, "z_hz": 400, "z_r": 0.85},
                {"p_hz": 1600, "p_r": 0.970, "z_hz": 800, "z_r": 0.85},
                {"p_hz": 3200, "p_r": 0.965, "z_hz": 1600, "z_r": 0.85},
            ],
            "M100": [ # Crystal Dust
                {"p_hz": 100, "p_r": 0.990, "z_hz": 50, "z_r": 0.85}, # Anchor hum
                {"p_hz": 500, "p_r": 0.995, "z_hz": 250, "z_r": 0.92}, # Fracture combs
                {"p_hz": 1100, "p_r": 0.995, "z_hz": 550, "z_r": 0.92},
                {"p_hz": 2300, "p_r": 0.995, "z_hz": 1150, "z_r": 0.92},
                {"p_hz": 4700, "p_r": 0.995, "z_hz": 2350, "z_r": 0.92},
                {"p_hz": 9500, "p_r": 0.995, "z_hz": 4750, "z_r": 0.92},
            ]
        }

def count_peaks(db_arr, threshold_db=8.0):
    """Simple peak finding without scipy."""
    peaks = 0
    for i in range(1, len(db_arr) - 1):
        if db_arr[i] > db_arr[i-1] and db_arr[i] > db_arr[i+1]:
            # Prominence approximation
            if db_arr[i] - min(db_arr[max(0, i-50)], db_arr[min(len(db_arr)-1, i+50)]) > threshold_db:
                peaks += 1
    return peaks

def verify_gates(body_id, db_m0, db_m100, freqs):
    """Run falsifiable failure modes based on BODIES.md contracts."""
    print(f"\n--- VERIFICATION GATES: {body_id} ---")
    passed = True
    
    def get_db(db_arr, target_hz):
        idx = (np.abs(freqs - target_hz)).argmin()
        return db_arr[idx]

    if body_id == "speaker_knockerz":
        sub_m0 = get_db(db_m0, 55)
        sub_m100 = get_db(db_m100, 55)
        flank_m100 = get_db(db_m100, 1600)
        cry_m100 = get_db(db_m100, 1200)
        
        pass_sub = (sub_m0 > 0) and (sub_m100 > 0) and (abs(sub_m0 - sub_m100) < 6.0)
        pass_cry = (cry_m100 - flank_m100) >= 10.0
        print(f"Invariant (Sub < 6dB drop): {'PASS' if pass_sub else 'FAIL'} (M0: {sub_m0:.1f}dB, M100: {sub_m100:.1f}dB)")
        print(f"Cone Cry (>10dB flank): {'PASS' if pass_cry else 'FAIL'} (Cry: {cry_m100:.1f}dB, Flank: {flank_m100:.1f}dB)")
        passed = pass_sub and pass_cry

    elif body_id == "aluminum_siding":
        scoop_m0 = get_db(db_m0, 1000)
        scoop_m100 = get_db(db_m100, 1000)
        
        pass_scoop_m0 = (scoop_m0 <= -12.0)
        pass_scoop_m100 = (scoop_m100 <= -12.0)
        print(f"Invariant (1kHz Scoop M0 <= -12dB): {'PASS' if pass_scoop_m0 else 'FAIL'} ({scoop_m0:.1f}dB)")
        print(f"Invariant (1kHz Scoop M100 <= -12dB): {'PASS' if pass_scoop_m100 else 'FAIL'} ({scoop_m100:.1f}dB)")
        passed = pass_scoop_m0 and pass_scoop_m100

    elif body_id == "small_talk_ah_ee":
        mask = (freqs >= 200) & (freqs <= 4000)
        peaks_m0 = count_peaks(db_m0[mask], threshold_db=8.0)
        peaks_m100 = count_peaks(db_m100[mask], threshold_db=8.0)
        
        pass_m0 = peaks_m0 >= 2
        pass_m100 = peaks_m100 >= 2
        print(f"Invariant (2+ Peaks M0): {'PASS' if pass_m0 else 'FAIL'} (Found {peaks_m0})")
        print(f"Invariant (2+ Peaks M100): {'PASS' if pass_m100 else 'FAIL'} (Found {peaks_m100})")
        passed = pass_m0 and pass_m100

    elif body_id == "cul_de_sac":
        root_m0 = get_db(db_m0, 100)
        root_m100 = get_db(db_m100, 100)
        
        peaks_m100 = count_peaks(db_m100, threshold_db=5.0)
        
        pass_root = abs(root_m0 - root_m100) < 5.0
        pass_comb = peaks_m100 >= 6
        print(f"Invariant (Root Hum anchors): {'PASS' if pass_root else 'FAIL'} (Diff: {abs(root_m0 - root_m100):.1f}dB)")
        print(f"Fracture (6+ Comb Peaks M100): {'PASS' if pass_comb else 'FAIL'} (Found {peaks_m100})")
        passed = pass_root and pass_comb

    return passed

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("body_id", choices=["speaker_knockerz", "aluminum_siding", "small_talk_ah_ee", "cul_de_sac"])
    parser.add_argument("--out", default="cartridge.json")
    args = parser.parse_args()

    arch = get_body_archetype(args.body_id)
    freqs = np.geomspace(20.0, 18000.0, 4096)
    
    print(f"Authoring {args.body_id}...")
    
    # 1. Independent Bisection Calibration
    c0_m0, peak_m0, compiled_m0 = calibrate_c0(arch["M0"], arch["target_peak"], freqs)
    c0_m100, peak_m100, compiled_m100 = calibrate_c0(arch["M100"], arch["target_peak"], freqs)
    
    print(f"M0 Calibration   : c0 = {c0_m0:.5f} | Peak = {peak_m0:+.1f} dB")
    print(f"M100 Calibration : c0 = {c0_m100:.5f} | Peak = {peak_m100:+.1f} dB")

    # 2. Gate Verification
    db_m0 = cascade_db(compiled_m0, freqs)
    db_m100 = cascade_db(compiled_m100, freqs)
    passed = verify_gates(args.body_id, db_m0, db_m100, freqs)
    
    if not passed:
        print("\nWARNING: Failed invariant gates! Cartridge will write, but needs tuning.")
    
    # 3. Export to RawStage format
    raw_m0 = [to_raw_stage(c) for c in compiled_m0]
    while len(raw_m0) < 12: raw_m0.append({"a1":0,"r":0,"val1":0,"val2":0,"val3":0,"flag":0})
        
    raw_m100 = [to_raw_stage(c) for c in compiled_m100]
    while len(raw_m100) < 12: raw_m100.append({"a1":0,"r":0,"val1":0,"val2":0,"val3":0,"flag":0})

    cartridge = {
        "name": args.body_id,
        "sampleRate": AUTHORING_SR,
        "stages": 6,
        "keyframes": [
            {"morph": 0.0, "q": 0.0, "label": "M0_Q0", "boost": 4.0, "stages": raw_m0},
            {"morph": 0.0, "q": 1.0, "label": "M0_Q100", "boost": 4.0, "stages": raw_m0},
            {"morph": 1.0, "q": 0.0, "label": "M100_Q0", "boost": 4.0, "stages": raw_m100},
            {"morph": 1.0, "q": 1.0, "label": "M100_Q100", "boost": 4.0, "stages": raw_m100}
        ]
    }
    
    with open(args.out, "w", encoding="utf-8") as f:
        json.dump(cartridge, f, indent=2)
        
    print(f"\nCartridge successfully written to {args.out} using RawStage schema.")

if __name__ == "__main__":
    main()