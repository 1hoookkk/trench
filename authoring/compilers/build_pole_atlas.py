#!/usr/bin/env python3
"""Build Pole Atlas for Stage Workbench.

Scans raw and compiled filter JSONs to build a unified reference atlas.
"""
import os
import json
import csv
import math
import glob
from pathlib import Path

# Constants for analysis
FS = 39062.5
STEPS = 120

def get_magnitude_from_coeffs(hz, c0, c1, c2, c3, c4):
    w = 2 * math.pi * hz / FS
    num_re = c0 + c1 * math.cos(w) + c2 * math.cos(2*w)
    num_im = -c1 * math.sin(w) - c2 * math.sin(2*w)
    den_re = 1 - c3 * math.cos(w) - c4 * math.cos(2*w)
    den_im = c3 * math.sin(w) + c4 * math.sin(2*w)
    num_mag_sq = num_re**2 + num_im**2
    den_mag_sq = den_re**2 + den_im**2
    return math.sqrt(num_mag_sq / max(1e-12, den_mag_sq))

def get_magnitude_from_stage(hz, s):
    w = 2 * math.pi * hz / FS
    cosW = math.cos(w)
    r = s.get('radius', 0)
    theta = 2 * math.pi * s.get('pole_freq_hz', 0) / FS
    cosT = math.cos(theta)
    mag = 1.0
    if s.get('kind') in ['allpole', 'lowpass2', 'explicit_zero']:
        den = 1 + r**4 + 4*r**2 * cosT**2 - 4*r*(1+r**2)*cosT*cosW + 2*r**2 * (2*cosW**2 - 1)
        mag *= 1.0 / math.sqrt(max(1e-12, den))
        if s.get('kind') == 'explicit_zero' and s.get('zero_radius', 0) > 0:
            rZ = s['zero_radius']
            thetaZ = 2 * math.pi * s.get('zero_freq_hz', 0) / FS
            cosTZ = math.cos(thetaZ)
            num = 1 + rZ**4 + 4*rZ**2 * cosTZ**2 - 4*rZ*(1+rZ**2)*cosTZ*cosW + 2*rZ**2 * (2*cosW**2 - 1)
            mag *= math.sqrt(num)
    return mag

def analyze_peak_db(obj, is_compiled=False):
    max_mag = 0
    min_mag = 100.0
    for i in range(STEPS):
        f = 20 * (1000**(i/STEPS))
        if f > FS/2: break
        if is_compiled:
            m = get_magnitude_from_coeffs(f, obj['c0'], obj['c1'], obj['c2'], obj['c3'], obj['c4'])
        else:
            m = get_magnitude_from_stage(f, obj)
        if m > max_mag: max_mag = m
        if m < min_mag: min_mag = m
    peak_db = 20 * math.log10(max_mag + 1e-12)
    cancel_db = 20 * math.log10(min_mag + 1e-12)
    return peak_db, cancel_db

def process_file(path, atlas):
    with open(path, 'r') as f:
        try:
            data = json.load(f)
        except Exception as e:
            return

    if not isinstance(data, dict):
        return

    fmt = data.get('format', '')
    name = data.get('name', path.stem)
    
    if fmt == 'raw-stage-v1':
        frames = data.get('frames', {})
        if isinstance(frames, dict):
            for frame_id, frame_data in frames.items():
                if not isinstance(frame_data, dict): continue
                stages = frame_data.get('stages', [])
                for idx, s in enumerate(stages):
                    if not isinstance(s, dict): continue
                    peak, cancel = analyze_peak_db(s)
                    atlas.append({
                        "source": "raw-v1",
                        "body_name": name,
                        "corner_label": frame_id,
                        "morph": 0 if 'M0' in frame_id else (1 if 'M100' in frame_id else 0.5),
                        "stage_index": idx,
                        "user_stage": f"S{idx+1}",
                        "pole_freq_hz": s.get('pole_freq_hz'),
                        "pole_radius": s.get('radius'),
                        "zero_freq_hz": s.get('zero_freq_hz'),
                        "zero_radius": s.get('zero_radius'),
                        "stage_peak_db": round(peak, 2),
                        "cancellation_db": round(cancel, 2),
                        "role": s.get('role'),
                        "source_file": str(path)
                    })
                
    elif fmt == 'compiled-v1' or 'compiled' in str(path):
        keyframes = data.get('keyframes', [])
        if not keyframes and 'frames' in data:
            frames_node = data.get('frames')
            if isinstance(frames_node, dict):
                for fid, stages in frames_node.items():
                    if not isinstance(stages, list): continue
                    for idx, s in enumerate(stages):
                        if not isinstance(s, dict) or 'c0' not in s: continue
                        peak, cancel = analyze_peak_db(s, is_compiled=True)
                        atlas.append({
                            "source": "compiled",
                            "body_name": name,
                            "corner_label": fid,
                            "morph": 0 if 'M0' in fid else (1 if 'M100' in fid else 0.5),
                            "stage_index": idx, "user_stage": f"S{idx+1}",
                            "c0": s['c0'], "c1": s['c1'], "c2": s['c2'], "c3": s['c3'], "c4": s['c4'],
                            "stage_peak_db": round(peak, 2), "cancellation_db": round(cancel, 2),
                            "source_file": str(path)
                        })
            return

        for kf in keyframes:
            if not isinstance(kf, dict): continue
            label = kf.get('label', 'UNK')
            morph = kf.get('morph', 0.5)
            q = kf.get('q', 0.5)
            stages = kf.get('stages', [])
            for idx, s in enumerate(stages):
                if not isinstance(s, dict) or 'c0' not in s: continue
                peak, cancel = analyze_peak_db(s, is_compiled=True)
                atlas.append({
                    "source": "compiled-v1",
                    "body_name": name,
                    "corner_label": label,
                    "morph": morph,
                    "q": q,
                    "stage_index": idx,
                    "user_stage": f"S{idx+1}",
                    "c0": s['c0'], "c1": s['c1'], "c2": s['c2'], "c3": s['c3'], "c4": s['c4'],
                    "stage_peak_db": round(peak, 2),
                    "cancellation_db": round(cancel, 2),
                    "source_file": str(path)
                })

    elif 'x3_filter' in str(path):
        filt = data.get('filter', {})
        if not isinstance(filt, dict): return
        sections = filt.get('sections', [])
        for idx, s in enumerate(sections):
            if not isinstance(s, dict) or s.get('type', 0) == 0: continue
            atlas.append({
                "source": "x3_filter",
                "body_name": data.get('template', {}).get('name', name),
                "corner_label": "STATIC",
                "stage_index": idx,
                "user_stage": f"S{idx+1}",
                "pole_freq_hz": s.get('low_freq'),
                "stage_peak_db": s.get('low_gain'),
                "source_file": str(path)
            })

def main():
    atlas = []
    search_patterns = [
        "authoring/emu_designer/generated/*.raw.json",
        "authoring/emu_designer/generated/*.compiled.json",
        "authoring/peak_shelf/generated/*.raw.json",
        "authoring/peak_shelf/generated/*.compiled.json",
        "authoring/imported_x3/*.x3_filter.json",
        "tmp/*.raw.json",
        "tmp/*.compiled.json",
        "trench-juce/plugin/assets/cartridges/*.compiled.json",
        "trench-juce/plugin/assets/cartridges/*.raw.compiled.json"
    ]
    for pattern in search_patterns:
        for path in glob.glob(pattern):
            process_file(Path(path), atlas)
    
    print(f"Extraction complete. Found {len(atlas)} atlas entries.")
    output_json = Path("authoring/atlas/pole_atlas.json")
    output_json.parent.mkdir(parents=True, exist_ok=True)
    with open(output_json, 'w') as f:
        json.dump(atlas, f, indent=2)
    
    output_csv = Path("authoring/atlas/pole_atlas_summary.csv")
    if atlas:
        fieldnames = ["source", "body_name", "corner_label", "morph", "q", "stage_index", "user_stage", 
                      "pole_freq_hz", "pole_radius", "zero_freq_hz", "zero_radius", 
                      "c0", "c1", "c2", "c3", "c4", 
                      "stage_peak_db", "cancellation_db", "role", "source_file"]
        with open(output_csv, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore')
            writer.writeheader()
            writer.writerows(atlas)
    print(f"Atlas saved to {output_json} and {output_csv}")

if __name__ == "__main__":
    main()
