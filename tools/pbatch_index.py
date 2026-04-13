"""pbatch_index.py — Index ALL P2k and XML/Cubes into static Shapes.

This script sweeps across the TRENCH RE repositories, extracts every corner 
of every filter body, and saves them as individual static "Shapes" in 
vault/_shapes/. It also generates a 6-stage breakdown plot for every shape.

Workflow:
1. Scan cartridges/p2k/ and trenchwork_clean/vault/heritage/
2. Scan trench_re_vault/ (XML targets)
3. Extract corners -> vault/_shapes/S_XXX.json
4. Plot -> vault/_shapes/plots/S_XXX.png
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import numpy as np

# Add project root to path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pyruntime.body import Body
from pyruntime.constants import SR, NUM_BODY_STAGES
from pyruntime.corner import CornerName, CornerState
from pyruntime.freq_response import stage_response, cascade_response_db, freq_points
from pyruntime.designer_compile import parse_xml, compile_designer_to_body

# Directories
SHAPES_DIR = ROOT / "vault" / "_shapes"
PLOTS_DIR = SHAPES_DIR / "plots"
P2K_DIR = ROOT / "cartridges" / "p2k"
HERITAGE_DIR = Path("C:/Users/hooki/trenchwork_clean/vault/heritage")
RE_VAULT_DIR = Path("C:/Users/hooki/trench_re_vault")

FREQS = freq_points(1024)

def plot_shape_breakdown(name: str, stages: list, out_path: Path):
    """Generate a 6-stage breakdown plot + total response."""
    plt.figure(figsize=(10, 6))
    
    # 1. Plot individual stage responses (transparent)
    colors = ['#ff7f0e', '#1f77b4', '#2ca02c', '#d62728', '#9467bd', '#8c564b']
    for i, enc in enumerate(stages[:6]):
        if abs(enc.c0 - 1.0) < 0.001 and abs(enc.c3) < 0.001:
            continue # Skip passthrough
        
        # Complex response
        resp = stage_response(enc, FREQS, SR)
        mag_db = 20.0 * np.log10(np.maximum(np.abs(resp), 1e-10))
        
        plt.plot(FREQS, mag_db, color=colors[i % len(colors)], alpha=0.3, 
                 label=f"Stage {i+1}", linestyle='--')

    # 2. Plot total cascade response (bold)
    total_db = cascade_response_db(stages, FREQS, SR)
    plt.plot(FREQS, total_db, color='white', linewidth=2.5, label="Total Curve")
    
    # Styling
    plt.xscale('log')
    plt.grid(True, which='both', linestyle=':', alpha=0.5)
    plt.xlim(20, SR/2)
    plt.ylim(-60, 60)
    plt.title(f"TRENCH Shape: {name}", color='cyan', fontsize=14)
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Magnitude (dB)")
    
    # Dark Mode
    ax = plt.gca()
    ax.set_facecolor('#1a1a1a')
    plt.gcf().set_facecolor('#121212')
    ax.tick_params(colors='gray')
    ax.xaxis.label.set_color('gray')
    ax.yaxis.label.set_color('gray')
    for spine in ax.spines.values():
        spine.set_color('#333333')
        
    plt.legend(facecolor='#1a1a1a', edgecolor='#333333', labelcolor='gray')
    plt.tight_layout()
    plt.savefig(out_path, dpi=120)
    plt.close()

def extract_from_body(body: Body, source_tag: str) -> list[dict]:
    """Extract 4 corners from a Body as Shapes."""
    shapes = []
    for cn in [CornerName.A, CornerName.B, CornerName.C, CornerName.D]:
        corner = body.corners.corner(cn)
        # Convert stages to serializable dicts
        encoded_stages = []
        for sp in corner.stages:
            from pyruntime.encode import raw_to_encoded
            enc = raw_to_encoded(sp)
            encoded_stages.append({
                "c0": float(enc.c0), "c1": float(enc.c1), "c2": float(enc.c2),
                "c3": float(enc.c3), "c4": float(enc.c4)
            })
            
        # Metric calculation using encoded stages
        from pyruntime.encode import EncodedCoeffs
        enc_list = [EncodedCoeffs(**st) for st in encoded_stages]
        total_db = cascade_response_db(enc_list, FREQS, SR)
        peak_db = float(np.max(total_db))
        
        # Spectral Centroid
        linear = 10 ** (total_db / 20)
        total_lin = linear.sum()
        centroid = float((FREQS * linear).sum() / total_lin) if total_lin > 0 else 0.0

        shape_id = f"{source_tag}_{body.name}_{cn.json_key()}"
        shapes.append({
            "id": shape_id,
            "name": f"{body.name} ({cn.json_key()})",
            "source": source_tag,
            "peak_db": round(peak_db, 2),
            "centroid_hz": round(centroid, 1),
            "stages": encoded_stages
        })
    return shapes

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--limit", type=int, default=None, help="Limit number of shapes to process")
    args = parser.parse_args()

    SHAPES_DIR.mkdir(parents=True, exist_ok=True)
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)

    all_shapes = []

    # 1. Process cartridges/p2k/
    print("Indexing cartridges/p2k/...")
    p2k_files = sorted(P2K_DIR.glob("*.json"))
    for f in p2k_files:
        try:
            body = Body.from_json(str(f))
            shapes = extract_from_body(body, "P2K")
            all_shapes.extend(shapes)
        except Exception as e:
            print(f"  Error loading {f.name}: {e}")

    # 2. Process heritage/
    if HERITAGE_DIR.exists():
        print(f"Indexing {HERITAGE_DIR}...")
        h_files = sorted(HERITAGE_DIR.glob("*.json"))
        for f in h_files:
            try:
                # Heritage files might be Body or DesignerTemplate
                data = json.loads(f.read_text())
                if "corners" in data:
                    body = Body.from_json(str(f))
                else:
                    # Assume designer template
                    from pyruntime.designer_compile import compile_designer_to_body, parse_xml
                    # If XML is in the directory, handled later. If JSON, we might need a converter.
                    # For now just skip non-body JSONs in heritage
                    continue
                shapes = extract_from_body(body, "HERITAGE")
                all_shapes.extend(shapes)
            except Exception as e:
                print(f"  Error loading {f.name}: {e}")

    # 3. Process filter_design_corpus.json (thousands of targets)
    corpus_file = RE_VAULT_DIR / "analysis" / "filter_design_corpus.json"
    if corpus_file.exists():
        print(f"Indexing corpus: {corpus_file}...")
        try:
            from pyruntime.designer_compile import make_template, compile_designer_to_body
            corpus_data = json.loads(corpus_file.read_text())
            for entry in corpus_data:
                name = entry.get("name", "unknown")
                sections = []
                for s in entry.get("stages", []):
                    sections.append((
                        s.get("type", 0),
                        s.get("low_freq", 0),
                        s.get("low_gain", 0),
                        s.get("high_freq", 0),
                        s.get("high_gain", 0)
                    ))
                template = make_template(name, sections)
                template.frequency = entry.get("frequency", 0.0)
                template.gain = entry.get("gain", 0.0)
                
                body = compile_designer_to_body(template)
                shapes = extract_from_body(body, "CORPUS")
                # Only take A (M0) and C (M100) since these are morph-only
                all_shapes.extend([shapes[0], shapes[2]])
        except Exception as e:
            print(f"  Error processing corpus: {e}")

    print(f"\nFinal processing of {len(all_shapes)} shapes...")
    
    # Write and Plot
    manifest = []
    
    to_process = all_shapes
    if args.limit:
        to_process = all_shapes[:args.limit]

    print(f"Processing and plotting {len(to_process)} shapes...")
    for i, s in enumerate(to_process):
        # Save JSON
        safe_id = s['id'].replace(" ", "_").replace("/", "_").replace("\\", "_").replace(":", "_")
        shape_file = SHAPES_DIR / f"{safe_id}.json"
        
        # Don't rewrite if identical
        if not shape_file.exists():
            with open(shape_file, "w") as f:
                json.dump(s, f, indent=2)
            
        # Plot
        plot_file = PLOTS_DIR / f"{safe_id}.png"
        if not plot_file.exists():
            try:
                from pyruntime.encode import EncodedCoeffs
                enc_stages = [EncodedCoeffs(**st) for st in s['stages']]
                plot_shape_breakdown(s['name'], enc_stages, plot_file)
            except Exception as e:
                print(f"  Failed plotting {s['name']}: {e}")
            
        manifest.append({
            "id": s['id'],
            "name": s['name'],
            "peak_db": s['peak_db'],
            "centroid_hz": s['centroid_hz'],
            "plot": str(plot_file.relative_to(ROOT))
        })
        
        if (i+1) % 50 == 0:
            print(f"  Processed {i+1}/{len(all_shapes)}")

    # Save master index
    index_file = SHAPES_DIR / "shapes_index.json"
    with open(index_file, "w") as f:
        json.dump(manifest, f, indent=2)
        
    print(f"\nSUCCESS: Indexed {len(all_shapes)} shapes.")
    print(f"Index: {index_file}")
    print(f"Plots: {PLOTS_DIR}")

if __name__ == "__main__":
    main()
