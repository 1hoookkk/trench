import json
import os
from pathlib import Path
import numpy as np

# Ground Truth Paths
ROOT = Path(r"C:\Users\hooki\Trench")
SKINS_DIR = Path(r"C:\Users\hooki\trenchwork_clean\datasets\p2k_skins")
NAMES_FILE = Path(r"C:\Users\hooki\trenchwork_clean\datasets\p2k_filter_names.json")
OUT_DIR = ROOT / "cartridges" / "p2k"

# Domain 1 Constants
SR_NATIVE = 39062.5
NUM_STAGES = 12

def raw_to_compiled_stage(stage_dict):
    """Project raw P2K stage values into the current DF2T coefficient domain."""
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

def main():
    if not OUT_DIR.exists():
        OUT_DIR.mkdir(parents=True)
    
    # Load friendly names
    with open(NAMES_FILE, 'r') as f:
        names_data = json.load(f)
        filter_names = names_data.get("names", {})

    processed = 0
    for skin_file in SKINS_DIR.glob("P2k_*.json"):
        with open(skin_file, 'r') as f:
            skin = json.load(f)
        
        id_str = str(skin["filterType"])
        friendly_name = filter_names.get(id_str, skin["name"])
        
        cartridge = {
            "format": "compiled-v1",
            "name": friendly_name,
            "provenance": "heritage-p2k-extraction",
            "sampleRate": SR_NATIVE,
            "stages": NUM_STAGES,
            "keyframes": []
        }
        
        for label, corner_data in skin["corners"].items():
            stages = []
            # Encode active 6 stages
            for s in corner_data["stages"]:
                stages.append(raw_to_compiled_stage(s))
            
            # Pad to 12
            while len(stages) < NUM_STAGES:
                stages.append(passthrough())
            
            cartridge["keyframes"].append({
                "label": label,
                "boost": skin.get("boost", 4.0),
                "stages": stages
            })
            
        out_path = OUT_DIR / f"{skin['name']}.json"
        with open(out_path, 'w') as f:
            json.dump(cartridge, f, indent=2)
        processed += 1

    print(f"DONE: Promoted {processed} P2K skins to playable cartridges in {OUT_DIR}")

if __name__ == "__main__":
    main()
