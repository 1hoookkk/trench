#!/usr/bin/env python3
import json
import argparse
import numpy as np
import sys
from pathlib import Path

SR = 39062.5

def get_biquad_mag(coeffs, freqs):
    # DF2T: H(z) = (c0 + c1*z^-1 + c2*z^-2) / (1 + c3*z^-1 + c4*z^-2)
    w = 2 * np.pi * freqs / SR
    z = np.exp(1j * w)
    z_inv = 1.0 / z
    z_inv2 = z_inv * z_inv
    num = coeffs["c0"] + coeffs["c1"] * z_inv + coeffs["c2"] * z_inv2
    den = 1.0 + coeffs["c3"] * z_inv + coeffs["c4"] * z_inv2
    return num / den

def check_stability(coeffs):
    # Denominator: 1 + c3*z^-1 + c4*z^-2 = 0 -> z^2 + c3*z + c4 = 0
    # Roots: [-c3 +/- sqrt(c3^2 - 4*c4)] / 2
    c3 = coeffs["c3"]
    c4 = coeffs["c4"]
    discr = c3**2 - 4*c4
    if discr >= 0:
        r1 = (-c3 + np.sqrt(discr)) / 2.0
        r2 = (-c3 - np.sqrt(discr)) / 2.0
    else:
        r1 = complex(-c3 / 2.0, np.sqrt(-discr) / 2.0)
        r2 = r1.conjugate()
    
    m1, m2 = abs(r1), abs(r2)
    return m1 < 1.0 and m2 < 1.0, max(m1, m2)

def validate_cartridge(path):
    with open(path) as f:
        doc = json.load(f)
    
    name = doc.get("name", "unknown")
    kfs = doc.get("keyframes", [])
    if not kfs:
        return "KILL", "No keyframes found."

    freqs = np.logspace(np.log10(20), np.log10(19000), 500)
    
    out_lines = [f"=== {name} ==="]
    global_decision = "PASS"
    global_notes = []

    for kf in kfs:
        label = kf.get("label", "unknown")
        stages = kf["stages"]
        boost = kf.get("boost", 1.0)
        
        h_total = np.ones_like(freqs, dtype=complex)
        max_stage_peak = -999.0
        unstable = False
        max_radius = 0.0
        
        for i, s in enumerate(stages):
            # Stage names 1-6
            sid = i + 1
            
            # 1. Stability
            is_stable, radius = check_stability(s)
            max_radius = max(max_radius, radius)
            if not is_stable:
                unstable = True

            # 2. Magnitude
            h_s = get_biquad_mag(s, freqs)
            if np.any(np.isnan(h_s)) or np.any(np.isinf(h_s)):
                unstable = True
                continue
                
            mag_s = 20 * np.log10(np.abs(h_s) + 1e-9)
            max_stage_peak = max(max_stage_peak, np.max(mag_s))
            h_total *= h_s

        mag_total = 20 * np.log10(np.abs(h_total) * boost + 1e-9)
        cascade_peak = np.max(mag_total)
        # Cancellation: max stage peak vs final cascade peak
        cancellation = max_stage_peak - cascade_peak

        out_lines.append(f"  {label}: cascade={cascade_peak:.1f} dB, max_stage={max_stage_peak:.1f} dB, cancel={cancellation:.1f} dB, r_max={max_radius:.3f}")

        # Validation Logic:
        # 1. Poles MUST be stable
        if unstable or max_radius >= 1.0:
            global_decision = "KILL"
            global_notes.append(f"{label} UNSTABLE (r={max_radius:.3f})")
        
        # 2. Composite cascade MUST NOT exceed safety ceiling (arbitrary +30dB for now, but 0dB ideal)
        # Heritage pack typically expects low cascade peak due to boost=0.004
        if cascade_peak > 24.0:
            if global_decision != "KILL": global_decision = "REVIEW"
            global_notes.append(f"{label} HIGH CASCADE ({cascade_peak:.1f} dB)")

    if global_decision == "PASS" and not global_notes:
        global_notes.append("Stable controlled cancellation.")

    print("\n".join(out_lines))
    print(f"DECISION: {global_decision}")
    print(f"NOTES: {' | '.join(global_notes)}")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--compiled", type=Path, required=True)
    args = ap.parse_args()
    validate_cartridge(args.compiled)

if __name__ == "__main__":
    main()
