import json
import math
import struct
from pathlib import Path

# ---------------------------------------------------------------------------
# Minifloat Utils (Mirroring sealed engine logic)
# ---------------------------------------------------------------------------
def _mf_decode(packed: int) -> float:
    packed &= 0xFFFF
    if packed == 0xFFFF: return 1.0
    if packed == 0x0000: return 0.0
    u = packed + 1
    exp = ((u >> 12) & 0xF) - 15
    mant = u & 0xFFF
    if exp < -14: return mant / 134217728.0
    return (2.0 ** exp) * ((mant | 0x1000) / 8192.0)

def _mf_encode(value: float) -> int:
    best_raw = 0
    best_err = float("inf")
    for raw in range(0x10000):
        err = abs(_mf_decode(raw) - value)
        if err < best_err:
            best_err = err
            best_raw = raw
            if err == 0: break
    return best_raw

# ---------------------------------------------------------------------------
# Stage Builders
# ---------------------------------------------------------------------------
SR = 39062.5

def build_w_from_coeffs(c0, c1, c2, c3, c4):
    """Inverse solve for d0..d4 and encode to w."""
    d1 = c1
    d0 = (c0 - c1) / 4.0
    d3 = c3
    d2 = (c2 - c3) / 4.0
    d4 = c4
    return [_mf_encode(d) for d in (d0, d1, d2, d3, d4)]

def build_w_from_params(a1, r, v1, v2, v3):
    """Mirror pyruntime.encode.raw_to_encoded logic."""
    c0 = 1.0 + v1
    c1 = a1 + v2
    c2 = r*r - v3
    c3 = a1
    c4 = r*r
    return build_w_from_coeffs(c0, c1, c2, c3, c4)

def resonator_with_zero(p_hz, r, v1, z_hz, zr_factor):
    theta = 2.0 * math.pi * p_hz / SR
    phi = 2.0 * math.pi * z_hz / SR
    a1 = -2.0 * r * math.cos(theta)
    zero_r = r * zr_factor
    b1_t = -2.0 * zero_r * math.cos(phi)
    b2_t = zero_r * zero_r
    return build_w_from_params(a1, r, v1, b1_t - a1, r*r - b2_t)

# ---------------------------------------------------------------------------
# Speaker Knockerz Recipe (v6)
# ---------------------------------------------------------------------------
LANDMARKS = [(190, 250), (354, 450), (800, 700), (1200, 1000), (2700, 4900), (4900, 6000)]

CORNER_CONFIGS = {
    "M0_Q0":   None,
    "M0_Q100":  None,
    "M100_Q0":   (0.40, 0.0, 0.72),
    "M100_Q100": (0.15, 0.0, 0.72),
}

def main():
    root = Path(r"c:\Users\hooki\Trench")
    source_path = root / "trench-juce" / "cartridges" / "p2k" / "P2k_006.json"
    dest_path = root / "trench-juce" / "cartridges" / "speaker_knockerz_baked_v6.json"
    
    with open(source_path) as f:
        p2k = json.load(f)
    
    baked_keyframes = []
    for kf in p2k["keyframes"]:
        label = kf["label"]
        cfg = CORNER_CONFIGS[label]
        
        stages = []
        # 1. Heritage Foundation (0-5)
        for s in kf["stages"][:6]:
            w = build_w_from_coeffs(s["c0"], s["c1"], s["c2"], s["c3"], s["c4"])
            stages.append({"w": w})
            
        # 2. Character Layer (6-11)
        if cfg is None:
            # Identity: (0.25, 0, 0.25, 0, 0)
            id_w = [0x1FFF, 0x0FFF, 0x1FFF, 0x0FFF, 0x0FFF]
            for _ in range(6):
                stages.append({"w": id_w})
        else:
            r, v1, zrf = cfg
            for ph, zh in LANDMARKS:
                w = resonator_with_zero(ph, r, v1, zh, zrf)
                stages.append({"w": w})
        
        baked_kf = {
            "label": label,
            "morph": kf["morph"],
            "q": kf["q"],
            "boost": kf["boost"],
            "stages": stages
        }
        baked_keyframes.append(baked_kf)
        
    result = {
        "format": "compiled-v1",
        "name": "Speaker Knockerz v6 (Baked)",
        "provenance": "sealed-pipeline-repro-v6",
        "sampleRate": 39062.5,
        "keyframes": baked_keyframes
    }
    
    with open(dest_path, "w") as f:
        json.dump(result, f, indent=2)
    
    print(f"Baked {dest_path}")

if __name__ == "__main__":
    main()
