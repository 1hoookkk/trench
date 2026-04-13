"""Render canonical reference WAVs for every P2K skin in datasets/p2k_skins.

Output: Trench/ref/canonical/<name>_<corner>.wav (float32) plus MANIFEST.json.
Pipeline matches tools/parity_null.py: raw stage → SOS cascade → AGC → boost.
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np
import soundfile as sf
from scipy.signal import sosfilt

SKINS = Path(r"C:/Users/hooki/trenchwork_clean/datasets/p2k_skins")
DRY = Path(r"C:/Users/hooki/trenchwork_clean/ref/bypassed-pinknoise.wav")
OUT = Path(r"C:/Users/hooki/Trench/ref/canonical")

AGC_TABLE = np.array(
    [1.0001, 1.0001, 0.996, 0.990, 0.920, 0.500, 0.200, 0.160,
     0.120, 0.120, 0.120, 0.120, 0.120, 0.120, 0.120, 0.120],
    dtype=np.float64,
)


def apply_agc(samples: np.ndarray) -> np.ndarray:
    out = np.empty_like(samples)
    g = 1.0
    for i in range(len(samples)):
        s = samples[i]
        idx = int(g * abs(s)) & 0xF
        ng = g * AGC_TABLE[idx]
        g = ng if ng < 1.0 else 1.0
        out[i] = s * g
    return out


def stage_co(s: dict) -> list[float]:
    a1 = float(s["a1"])
    r = min(float(s["r"]), 0.999999)
    a2 = r * r
    if float(s.get("flag", 1.0)) < 0.5:
        b0, b1, b2 = 1.0, 0.0, 0.0
    else:
        b0 = 1.0 + float(s["val1"])
        b1 = a1 + float(s["val2"])
        b2 = a2 - float(s["val3"])
    return [b0, b1, b2, 1.0, a1, a2]


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    dry_data, sr = sf.read(str(DRY))
    if dry_data.ndim > 1:
        dry_data = dry_data[:, 0]
    dry_data = np.ascontiguousarray(dry_data, dtype=np.float64)
    print(f"sr={sr}  dry={len(dry_data)} samples")

    files = sorted(f for f in os.listdir(SKINS) if f.endswith(".json"))
    print(f"{len(files)} skin files")

    manifest: dict = {}
    total_wavs = 0
    skipped: list[tuple[str, str]] = []

    for fn in files:
        name = fn[:-5]
        try:
            body = json.loads((SKINS / fn).read_text())
        except Exception as exc:
            skipped.append((name, f"parse: {exc}"))
            continue

        corners = body.get("corners")
        if not isinstance(corners, dict) or set(corners) != {
            "M0_Q0", "M0_Q100", "M100_Q0", "M100_Q100"
        }:
            skipped.append((name, f"corners shape: {type(corners).__name__}"))
            continue

        boost = float(body.get("boost", 1.0))
        corner_renders: list[str] = []
        ok = True
        for lbl in ("M0_Q0", "M0_Q100", "M100_Q0", "M100_Q100"):
            stages = corners[lbl].get("stages")
            if not isinstance(stages, list) or not stages:
                skipped.append((name, f"{lbl}: missing stages"))
                ok = False
                break
            try:
                sos = np.array([stage_co(s) for s in stages], dtype=np.float64)
                cascaded = sosfilt(sos, dry_data)
                pred = apply_agc(cascaded) * boost
                wav = OUT / f"{name}_{lbl}.wav"
                sf.write(str(wav), pred.astype(np.float32), int(sr), subtype="FLOAT")
                corner_renders.append(lbl)
                total_wavs += 1
            except Exception as exc:
                skipped.append((name, f"{lbl}: render: {exc}"))
                ok = False
                break

        if ok:
            manifest[name] = {
                "source": f"datasets/p2k_skins/{fn}",
                "boost": boost,
                "stage_count": int(body.get("stageCount", len(stages))),
                "corners": corner_renders,
            }

    (OUT / "MANIFEST.json").write_text(json.dumps(manifest, indent=2))
    print(f"\nrendered {total_wavs} wav files for {len(manifest)} skins")
    if skipped:
        print(f"skipped {len(skipped)}:")
        for name, why in skipped:
            print(f"  {name}: {why}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
