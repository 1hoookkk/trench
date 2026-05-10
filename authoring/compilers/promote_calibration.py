import json
import os
from pathlib import Path
import numpy as np
import math

ROOT = Path(r"C:\Users\hooki\Trench")
CAL_DIR = ROOT / "docs" / "calibration"
OUT_DIR = ROOT / "cartridges" / "p2k"
SR_NATIVE = 39062.5
NUM_STAGES = 12

CAL_TO_FT = {
    "Ooh to Eee (approx)": 33,
    "Talking Hedz": 36,
    "Razor Blades": 41,
    "Radio Craze": 42,
    "Freak Shifta": 46,
}

def compute_a1(pole_freq_hz, radius):
    return -2.0 * radius * math.cos(2.0 * math.pi * pole_freq_hz / SR_NATIVE)

def raw_to_compiled_stage(stage_dict):
    a1 = float(np.float32(stage_dict["a1"]))
    r = float(np.float32(stage_dict["r"]))
    val1 = float(np.float32(stage_dict["val1"]))
    val2 = float(np.float32(stage_dict["val2"]))
    val3 = float(np.float32(stage_dict["val3"]))
    
    b0 = 1.0 + val1
    b1 = a1 + val2
    b2 = r * r - val3
    
    return {"c0": b0, "c1": b1, "c2": b2, "c3": a1, "c4": r * r}

def passthrough():
    return {"c0": 1.0, "c1": 0.0, "c2": 0.0, "c3": 0.0, "c4": 0.0}

def process_calibration(cal_file):
    with open(cal_file, 'r') as f:
        cal = json.load(f)
    
    name = cal["name"]
    boost = cal.get("boost", 4.0)
    
    if name not in CAL_TO_FT:
        print(f"  SKIP: {name} not in calibration mapping")
        return None
    
    ft = CAL_TO_FT[name]
    friendly_name = f"P2k_{ft:03d}"
    
    cartridge = {
        "format": "compiled-v1",
        "name": friendly_name,
        "provenance": "calibration-promote",
        "sampleRate": SR_NATIVE,
        "stages": NUM_STAGES,
        "keyframes": []
    }
    
    for corner_label, corner_data in cal["corners"].items():
        stages = []
        
        for s in corner_data["stages"]:
            pole_freq_hz = s["pole_freq_hz"]
            radius = s["radius"]
            
            a1 = compute_a1(pole_freq_hz, radius)
            
            raw_stage = {
                "a1": a1,
                "r": radius,
                "val1": s["val1"],
                "val2": s["val2"],
                "val3": s["val3"],
            }
            
            encoded = raw_to_compiled_stage(raw_stage)
            stages.append(encoded)
        
        while len(stages) < NUM_STAGES:
            stages.append(passthrough())
        
        cartridge["keyframes"].append({
            "label": corner_label,
            "boost": boost,
            "stages": stages
        })
    
    cartridge["keyframes"].sort(key=lambda k: k["label"])
    
    return ft, cartridge

def main():
    if not OUT_DIR.exists():
        OUT_DIR.mkdir(parents=True)
    
    processed = 0
    
    cal_files = [
        "Ooh_to_Eee_(approx).json",
        "Talking_Hedz.json",
        "Razor_Blades.json",
        "Radio_Craze.json",
        "Freak_Shifta.json",
    ]
    
    for cal_filename in cal_files:
        cal_path = CAL_DIR / cal_filename
        if not cal_path.exists():
            print(f"NOT FOUND: {cal_path}")
            continue
        
        print(f"Processing: {cal_filename}")
        result = process_calibration(cal_path)
        
        if result is None:
            continue
        
        ft, cartridge = result
        out_path = OUT_DIR / f"P2k_{ft:03d}.json"
        
        with open(out_path, 'w') as f:
            json.dump(cartridge, f, indent=2)
        
        print(f"  -> {out_path.name}")
        processed += 1
    
    print(f"\nDONE: Promoted {processed} calibration files to compiled-v1 cartridges")

if __name__ == "__main__":
    main()
