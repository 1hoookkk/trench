"""Render canonical reference WAVs for every P2K source TRENCH has data for.

Two source types:
  1. Raw ROM stage data — datasets/p2k_skins/*.json
     Format: {corners: {LABEL: {stages: [{a1, r, val1, val2, val3, flag}]}}}
     33 numbered (P2k_000..P2k_032) + 2 alternate vocal extractions = 35 files.

  2. Reverse-engineered calibration data — docs/calibration/*.json
     Format: {sample_rate, corners: {LABEL: {stages: [{pole_freq_hz, radius,
              val1, val2, val3, c4_b0, ...}]}}}
     11 files. 6 of these cover filter types (33-55) that have no raw extraction
     in p2k_skins; 5 are duplicates of P2k_025/026/027/029/030 and are skipped
     to keep the gate from double-counting the same filter type.

Output: Trench/ref/canonical/<name>_<corner>.wav (float32)
        plus MANIFEST.json (per-skin source path, boost, stage count, corners)
        plus PROVENANCE.json (chain of custody for both source types)

Pipeline (matches tools/parity_null.py):
    raw stage → SOS cascade (scipy.signal.sosfilt, f64)
              → per-sample AGC (16-entry table from agc.rs)
              → × boost
              → write float32 WAV
"""
from __future__ import annotations

import json
import math
import os
import sys
from pathlib import Path

import numpy as np
import soundfile as sf
from scipy.signal import sosfilt

ROOT = Path(__file__).resolve().parent.parent
SKINS = ROOT / "datasets" / "p2k_skins"
CALIBRATION = ROOT / "docs" / "calibration"
DRY = Path(r"C:/Users/hooki/trenchwork_clean/ref/bypassed-pinknoise.wav")
OUT = ROOT / "ref" / "canonical"

# Calibration files whose ft index already has a raw p2k_skin extraction.
# Skip to avoid double-counting the same filter type. The numbers refer to the
# P2k_NNN.json file already covered.
CALIBRATION_DUPES = {
    "BassBox_303",     # P2k_029
    "Early_Rizer",     # P2k_025
    "Fuzzi_Face",      # P2k_030
    "Meaty_Gizmo",     # P2k_027
    "Millennium",      # P2k_026
}

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


def raw_stage_co(s: dict) -> list[float]:
    """raw p2k_skins stage → SOS row [b0, b1, b2, 1, a1, a2]."""
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


def calibration_stage_co(s: dict, sample_rate: float) -> list[float]:
    """Calibration stage → SOS row [b0, b1, b2, 1, a1, a2].

    Calibration stores pole_freq_hz + radius instead of raw a1. Reconstruct
    a1 = -2 * r * cos(2π * f / sr) where sr is the calibration's authoring
    sample rate (typically 39062.5 Hz for P2K).

    flag is implied: any nonzero val1/val2/val3 → resonator (flag=1);
    all zero → all-pole lowpass (flag=0). Matches the raw_stage_co branches.
    """
    r = min(float(s["radius"]), 0.999999)
    f = float(s["pole_freq_hz"])
    a1 = -2.0 * r * math.cos(2.0 * math.pi * f / float(sample_rate))
    a2 = r * r
    val1 = float(s.get("val1", 0.0))
    val2 = float(s.get("val2", 0.0))
    val3 = float(s.get("val3", 0.0))
    if val1 == 0.0 and val2 == 0.0 and val3 == 0.0:
        b0, b1, b2 = 1.0, 0.0, 0.0
    else:
        b0 = 1.0 + val1
        b1 = a1 + val2
        b2 = a2 - val3
    return [b0, b1, b2, 1.0, a1, a2]


def render_corner(sos_rows: list[list[float]], boost: float, dry: np.ndarray) -> np.ndarray:
    sos = np.asarray(sos_rows, dtype=np.float64)
    cascaded = sosfilt(sos, dry)
    return apply_agc(cascaded) * boost


def render_raw_skin(path: Path, dry: np.ndarray, sr: int) -> tuple[dict | None, str | None]:
    """Render a raw p2k_skin file. Returns (manifest_entry, error)."""
    try:
        body = json.loads(path.read_text())
    except Exception as exc:
        return None, f"parse: {exc}"

    corners = body.get("corners")
    if not isinstance(corners, dict) or set(corners) != {
        "M0_Q0", "M0_Q100", "M100_Q0", "M100_Q100"
    }:
        return None, "corners shape"

    boost = float(body.get("boost", 1.0))
    name = path.stem
    rendered: list[str] = []
    last_stages: list = []
    for label in ("M0_Q0", "M0_Q100", "M100_Q0", "M100_Q100"):
        stages = corners[label].get("stages")
        if not isinstance(stages, list) or not stages:
            return None, f"{label}: missing stages"
        sos_rows = [raw_stage_co(s) for s in stages]
        pred = render_corner(sos_rows, boost, dry)
        sf.write(str(OUT / f"{name}_{label}.wav"), pred.astype(np.float32), sr, subtype="FLOAT")
        rendered.append(label)
        last_stages = stages

    rel = path.relative_to(ROOT).as_posix()
    return {
        "name": name,
        "source": rel,
        "source_type": "raw_p2k_skin",
        "boost": boost,
        "stage_count": int(body.get("stageCount", len(last_stages))),
        "sample_rate_authored": 39062.5,
        "corners": rendered,
    }, None


def render_calibration(path: Path, dry: np.ndarray, sr: int) -> tuple[dict | None, str | None]:
    """Render a calibration file. Returns (manifest_entry, error)."""
    try:
        body = json.loads(path.read_text())
    except Exception as exc:
        return None, f"parse: {exc}"

    corners = body.get("corners")
    if not isinstance(corners, dict) or set(corners) != {
        "M0_Q0", "M0_Q100", "M100_Q0", "M100_Q100"
    }:
        return None, "corners shape"

    boost = float(body.get("boost", 1.0))
    cal_sr = float(body.get("sample_rate", 39062.5))

    raw_name = body.get("name") or path.stem
    safe_name = "cal_" + raw_name.replace(" ", "_").replace("(", "").replace(")", "")

    rendered: list[str] = []
    last_stages: list = []
    for label in ("M0_Q0", "M0_Q100", "M100_Q0", "M100_Q100"):
        stages = corners[label].get("stages")
        if not isinstance(stages, list) or not stages:
            return None, f"{label}: missing stages"
        sos_rows = [calibration_stage_co(s, cal_sr) for s in stages]
        pred = render_corner(sos_rows, boost, dry)
        sf.write(str(OUT / f"{safe_name}_{label}.wav"), pred.astype(np.float32), sr, subtype="FLOAT")
        rendered.append(label)
        last_stages = stages

    rel = path.relative_to(ROOT).as_posix()
    return {
        "name": safe_name,
        "source": rel,
        "source_type": "calibration_re",
        "boost": boost,
        "stage_count": int(body.get("stage_count", len(last_stages))),
        "sample_rate_authored": cal_sr,
        "corners": rendered,
    }, None


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    dry_data, sr = sf.read(str(DRY))
    if dry_data.ndim > 1:
        dry_data = dry_data[:, 0]
    dry_data = np.ascontiguousarray(dry_data, dtype=np.float64)
    print(f"sr={sr}  dry={len(dry_data)} samples")

    manifest: dict[str, dict] = {}
    skipped: list[tuple[str, str]] = []
    total_wavs = 0

    raw_files = sorted(f for f in os.listdir(SKINS) if f.endswith(".json"))
    print(f"raw p2k_skins: {len(raw_files)}")
    for fn in raw_files:
        entry, err = render_raw_skin(SKINS / fn, dry_data, int(sr))
        if err:
            skipped.append((fn, err))
            continue
        manifest[entry["name"]] = {k: v for k, v in entry.items() if k != "name"}
        total_wavs += len(entry["corners"])

    cal_files = sorted(
        f for f in os.listdir(CALIBRATION)
        if f.endswith(".json") and f != "index.json"
    )
    print(f"calibration files: {len(cal_files)}  (skipping {len(CALIBRATION_DUPES)} dupes of raw skins)")
    for fn in cal_files:
        stem = fn[:-5]
        if stem in CALIBRATION_DUPES:
            continue
        entry, err = render_calibration(CALIBRATION / fn, dry_data, int(sr))
        if err:
            skipped.append((fn, err))
            continue
        manifest[entry["name"]] = {k: v for k, v in entry.items() if k != "name"}
        total_wavs += len(entry["corners"])

    (OUT / "MANIFEST.json").write_text(json.dumps(manifest, indent=2))

    provenance = {
        "version": "canonical-v2",
        "description": "TRENCH parity gate references. Two source classes.",
        "raw_p2k_skins": {
            "path": "datasets/p2k_skins/",
            "count": sum(1 for v in manifest.values() if v["source_type"] == "raw_p2k_skin"),
            "format": "{corners: {LABEL: {stages: [{a1, r, val1, val2, val3, flag}]}}}",
            "canonical_sha256_p2k_013_prefix": "9fb1bef0d0212980",
            "upstream": "trenchwork_clean/datasets/p2k_skins (byte-identical)",
        },
        "calibration_re": {
            "path": "docs/calibration/",
            "count": sum(1 for v in manifest.values() if v["source_type"] == "calibration_re"),
            "format": "{sample_rate, corners: {LABEL: {stages: [{pole_freq_hz, radius, val1, val2, val3, c4_b0, ...}]}}}",
            "a1_derivation": "a1 = -2 * radius * cos(2*pi * pole_freq_hz / sample_rate_authored)",
            "skipped_dupes_of_raw": sorted(CALIBRATION_DUPES),
        },
        "dry_input": str(DRY).replace("\\", "/"),
        "render_pipeline": "raw stage -> SOS cascade -> AGC -> * boost -> float32 WAV",
        "playback_sample_rate": int(sr),
        "container_subtype": "FLOAT (float32, 24-bit mantissa, ~-150 dB null floor)",
        "supersedes": "trenchwork_clean/ref/_BROKEN_CE_CAPTURE_QUARANTINE/hedz*.wav (broken Cheat Engine capture, quarantined)",
    }
    (OUT / "PROVENANCE.json").write_text(json.dumps(provenance, indent=2))

    print(f"\nrendered {total_wavs} wav files for {len(manifest)} entries")
    raw_ct = sum(1 for v in manifest.values() if v["source_type"] == "raw_p2k_skin")
    cal_ct = sum(1 for v in manifest.values() if v["source_type"] == "calibration_re")
    print(f"  {raw_ct} raw p2k_skins  +  {cal_ct} calibration_re  =  {len(manifest)}")
    if skipped:
        print(f"skipped {len(skipped)}:")
        for name, why in skipped:
            print(f"  {name}: {why}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
