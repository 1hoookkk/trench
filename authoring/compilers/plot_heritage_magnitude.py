#!/usr/bin/env python3
import json
import argparse
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# SR matches E-mu reference
SR = 39062.5

def get_biquad_mag(coeffs, freqs):
    # DF2T biquad: y = c0*x + w1; w1 = c1*x - c3*y + w2; w2 = c2*x - c4*y
    # H(z) = (c0 + c1*z^-1 + c2*z^-2) / (1 + c3*z^-1 + c4*z^-2)
    w = 2 * np.pi * freqs / SR
    z = np.exp(1j * w)
    
    # Pre-calculated powers of z for performance
    z_inv = 1.0 / z
    z_inv2 = z_inv * z_inv
    
    num = coeffs["c0"] + coeffs["c1"] * z_inv + coeffs["c2"] * z_inv2
    den = 1.0 + coeffs["c3"] * z_inv + coeffs["c4"] * z_inv2
    
    return num / den

def main():
    ap = argparse.ArgumentParser(description="Plot full magnitude of Morph 0 and 100.")
    ap.add_argument("input", type=Path, help="Path to .compiled.json")
    ap.add_argument("--out", "-o", type=Path, required=True, help="Output PNG path")
    args = ap.parse_args()

    with open(args.input) as f:
        doc = json.load(f)

    # We assume format is compiled-v1
    # Keyframes: [0]=M0_Q0, [2]=M100_Q0
    kfs = doc.get("keyframes", [])
    if not kfs:
        print("Error: no keyframes found in JSON.")
        return

    m0_kf = kfs[0]
    m100_kf = kfs[2]
    
    m0_coeffs = m0_kf["stages"]
    m100_coeffs = m100_kf["stages"]
    
    # Use boost from each frame specifically
    boost_m0 = m0_kf.get("boost", 1.0)
    boost_m100 = m100_kf.get("boost", 1.0)

    freqs = np.logspace(np.log10(20), np.log10(19000), 500)
    
    def total_mag(stages, bst):
        h_total = np.ones_like(freqs, dtype=complex)
        for s in stages:
            h_total *= get_biquad_mag(s, freqs)
        return 20 * np.log10(np.abs(h_total) * bst + 1e-9)

    mag_m0 = total_mag(m0_coeffs, boost_m0)
    mag_m100 = total_mag(m100_coeffs, boost_m100)

    plt.figure(figsize=(10, 6))
    plt.semilogx(freqs, mag_m0, label="Morph 0 (M0_Q0)", color="#007acc", lw=2)
    plt.semilogx(freqs, mag_m100, label="Morph 100 (M100_Q0)", color="#cc0044", lw=2)
    
    plt.title(f"TRENCH Heritage Full Magnitude: {doc['name']}", fontsize=14, pad=20)
    plt.xlabel("Frequency (Hz)", fontsize=12)
    plt.ylabel("Magnitude (dB)", fontsize=12)
    plt.grid(True, which="both", ls="-", alpha=0.5)
    plt.legend(loc="upper right")
    
    # Auto-limit Y but keep at least 40dB range
    all_mag = np.concatenate([mag_m0, mag_m100])
    y_min = max(-100, np.min(all_mag) - 5)
    y_max = min(40, np.max(all_mag) + 5)
    plt.ylim(y_min, y_max)
    plt.xlim(20, 20000)
    
    plt.tight_layout()
    plt.savefig(args.out, dpi=120)
    plt.close()
    print(f"Plot saved to {args.out}")

if __name__ == "__main__":
    main()
