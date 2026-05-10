import json
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as signal
from pathlib import Path
import wave
import struct
import sys

SR = 39062.5

def get_biquad_mag(coeffs, freqs):
    w = 2 * np.pi * freqs / SR
    z = np.exp(1j * w)
    z_inv = 1.0 / z
    z_inv2 = z_inv * z_inv
    num = coeffs["c0"] + coeffs["c1"] * z_inv + coeffs["c2"] * z_inv2
    den = 1.0 + coeffs["c3"] * z_inv + coeffs["c4"] * z_inv2
    return num / den

GUARD = 1e-12

def raw_to_encoded(s):
    # Detect if already encoded
    if "c0" in s:
        return s
    
    # Raw format (a1, r, val1, val2, val3, flag)
    flag = s.get("flag", 1.0)
    a1 = s["a1"]
    r = s["r"]
    val1 = s["val1"]
    val2 = s["val2"]
    val3 = s["val3"]

    if flag >= 0.5:
        # Resonator
        c2 = 2.0 + a1
        c3 = 1.0 - r * r
        a2 = r * r
        b0 = 1.0 + val1
        b1 = a1 + val2
        b2 = a2 - val3
        c4 = b0
        c0 = 2.0 + b1 / b0 if abs(b0) > GUARD else 2.0
        c1 = 1.0 - b2 / b0 if abs(b0) > GUARD else 1.0
    else:
        # Lowpass
        a1_bq = a1 - 2.0
        a2_bq = 1.0 - r
        b0_bq = (1.0 + a1_bq + a2_bq) / 4.0
        c0 = 4.0 if abs(b0_bq) > GUARD else 2.0
        c1 = 0.0 if abs(b0_bq) > GUARD else 1.0
        c2 = a1
        c3 = r
        c4 = b0_bq
    
    return {"c0": c0, "c1": c1, "c2": c2, "c3": c3, "c4": c4}

def interpolate_keyframes(k1, k2, t):
    res = []
    # Ensure both are in encoded format for interpolation
    s1_enc = [raw_to_encoded(s) for s in k1["stages"]]
    s2_enc = [raw_to_encoded(s) for s in k2["stages"]]
    
    for s1, s2 in zip(s1_enc, s2_enc):
        interp = {}
        for k in ["c0", "c1", "c2", "c3", "c4"]:
            interp[k] = s1[k] * (1-t) + s2[k] * t
        res.append(interp)
    return res

def find_peaks(mag_db, freqs):
    # threshold -10dB from max
    #peak_indices, _ = signal.find_peaks(mag_db, height=np.max(mag_db)-12, distance=10)
    peak_indices, _ = signal.find_peaks(mag_db, distance=10)
    # Filter for significant peaks
    valid_peaks = []
    for idx in peak_indices:
        if mag_db[idx] > -20: # arbitrary noise floor
            valid_peaks.append(freqs[idx])
    return valid_peaks

def save_wav(filename, samples, sr):
    with wave.open(str(filename), 'w') as f:
        f.setnchannels(1)
        f.setsampwidth(4) # 32-bit float
        f.setframerate(int(sr))
        # Struct pack float
        data = struct.pack('<' + 'f' * len(samples), *samples)
        f.writeframes(data)

def render_audio(cartridge, duration=2.0):
    kfs = cartridge["keyframes"]
    num_samples = int(duration * SR)
    
    # Saw drone
    t_vals = np.linspace(0, 1, num_samples)
    saw = signal.sawtooth(2 * np.pi * 110 * np.arange(num_samples) / SR)
    
    output = np.zeros(num_samples)
    
    # DF2T states (one per stage)
    states = [[0.0, 0.0] for _ in range(len(kfs[0]["stages"]))]
    
    for i in range(num_samples):
        morph = t_vals[i]
        # Simplistic interpolation between first two keyframes (assuming M0 and M100)
        # In TRENCH, we have many keyframes but here we likely have idx 0 and 1
        stages = interpolate_keyframes(kfs[0], kfs[1], morph)
        boost = kfs[0]["boost"] * (1-morph) + kfs[1]["boost"] * morph
        
        x = saw[i]
        for s_idx, s in enumerate(stages):
            # DF2T
            # y[n] = c0*x[n] + s1[n-1]
            # s1[n] = c1*x[n] - c3*y[n] + s2[n-1]
            # s2[n] = c2*x[n] - c4*y[n]
            y = s["c0"] * x + states[s_idx][0]
            states[s_idx][0] = s["c1"] * x - s["c3"] * y + states[s_idx][1]
            states[s_idx][1] = s["c2"] * x - s["c4"] * y
            x = y
        output[i] = x * boost
        
    return output

def audit(path):
    with open(path) as f:
        cartridge = json.load(f)
    
    name = cartridge.get("name", "talking_bender")
    morph_points = [0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0]
    
    # Correctly identify Q=0 keyframes
    # We expect M0_Q0 at index 0 and M100_Q0 at index 2 (usually)
    # Safer: find by label or morph/q values
    kf_m0 = None
    kf_m100 = None
    for kf in cartridge["keyframes"]:
        if kf["q"] == 0.0:
            if kf["morph"] == 0.0: kf_m0 = kf
            if kf["morph"] == 1.0: kf_m100 = kf
            
    if not kf_m0 or not kf_m100:
        print("Error: Could not find M0_Q0 and M100_Q0 keyframes.")
        sys.exit(1)

    freqs = np.logspace(np.log10(20), np.log10(SR/2), 1000)
    
    fig, axes = plt.subplots(len(morph_points), 1, figsize=(10, 20), sharex=True)
    plt.subplots_adjust(hspace=0.4)
    
    report = []
    
    for i, m in enumerate(morph_points):
        stages = interpolate_keyframes(kf_m0, kf_m100, m)
        boost = kf_m0["boost"] * (1-m) + kf_m100["boost"] * m
        
        h_total = np.ones_like(freqs, dtype=complex)
        for s in stages:
            h_total *= get_biquad_mag(s, freqs)
        
        mag_db = 20 * np.log10(np.abs(h_total) * boost + 1e-9)
        
        pk_db = np.max(mag_db)
        nl_db = np.min(mag_db)
        pks = find_peaks(mag_db, freqs)
        
        report.append({
            "morph": m,
            "cascade_peak_db": round(pk_db, 2),
            "deepest_null_db": round(nl_db, 2),
            "peaks": [round(p, 1) for p in pks]
        })
        
        ax = axes[i]
        ax.semilogx(freqs, mag_db, color='#db2777' if m > 0.5 else '#22d3ee', linewidth=1.5)
        ax.set_title(f"Morph {int(m*100)}% - Peak: {pk_db:.1f}dB, Null: {nl_db:.1f}dB", color='white', fontsize=10)
        ax.set_ylim([-60, 80])
        ax.grid(True, which='both', color='#1e293b', alpha=0.5)
        ax.tick_params(colors='white')
        ax.set_facecolor('#0f172a')
        
    fig.patch.set_facecolor('#070a0f')
    plt.xlabel("Frequency (Hz)", color='white')
    plt.savefig(f"authoring/emu_designer/generated/{name}_audit_traj.png")
    
    print(json.dumps(report, indent=2))
    
    # Audio Rendering
    def render_generic(input_signal, kf_a, kf_b, boost_a, boost_b):
        num_samples = len(input_signal)
        output = np.zeros(num_samples)
        states = [[0.0, 0.0] for _ in range(len(kf_a["stages"]))]
        
        for i in range(num_samples):
            t = i / (num_samples - 1)
            stages = interpolate_keyframes(kf_a, kf_b, t)
            b = boost_a * (1-t) + boost_b * t
            
            x = input_signal[i]
            for s_idx, s in enumerate(stages):
                y = s["c0"] * x + states[s_idx][0]
                states[s_idx][0] = s["c1"] * x - s["c3"] * y + states[s_idx][1]
                states[s_idx][1] = s["c2"] * x - s["c4"] * y
                x = y
            output[i] = x * b
        return output

    print("Rendering audio sweeps...")
    dur = 2.0
    ns = int(dur * SR)
    t = np.arange(ns) / SR
    
    # 1. Saw Drone (110 Hz)
    saw = signal.sawtooth(2 * np.pi * 110 * t)
    out_saw = render_generic(saw, kf_m0, kf_m100, kf_m0["boost"], kf_m100["boost"])
    save_wav(f"authoring/emu_designer/generated/{name}_sweep_saw.wav", out_saw, SR)
    
    # 2. 808 (Pitched Sine Sweep)
    # 150Hz -> 40Hz exp decay
    freq_808 = 40 + (150-40) * np.exp(-5 * t)
    phase_808 = np.cumsum(freq_808) / SR
    sine_808 = np.sin(2 * np.pi * phase_808)
    out_808 = render_generic(sine_808, kf_m0, kf_m100, kf_m0["boost"], kf_m100["boost"])
    save_wav(f"authoring/emu_designer/generated/{name}_sweep_808.wav", out_808, SR)
    
    # 3. Reese (3 Detuned Saws)
    reese = (signal.sawtooth(2 * np.pi * 55 * t) + 
             signal.sawtooth(2 * np.pi * 55.3 * t) * 0.5 + 
             signal.sawtooth(2 * np.pi * 54.7 * t) * 0.5)
    out_reese = render_generic(reese, kf_m0, kf_m100, kf_m0["boost"], kf_m100["boost"])
    save_wav(f"authoring/emu_designer/generated/{name}_sweep_reese.wav", out_reese, SR)

    # 4. Vocal (Formant-friendly Noise pulse)
    # Using noise modulated by a 2Hz pulse for "vocal-like" cadence
    vocal_noise = np.random.uniform(-1, 1, ns) * (0.5 + 0.5 * np.cos(2 * np.pi * 2 * t))
    out_vocal = render_generic(vocal_noise, kf_m0, kf_m100, kf_m0["boost"], kf_m100["boost"])
    save_wav(f"authoring/emu_designer/generated/{name}_sweep_vocal.wav", out_vocal, SR)

    print("Done.")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: audit_morph_trajectory.py <path_to_compiled_json>")
        sys.exit(1)
    audit(sys.argv[1])
