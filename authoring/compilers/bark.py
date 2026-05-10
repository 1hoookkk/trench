#!/usr/bin/env python3
"""Bark-scale helpers for TRENCH authoring surface.

Bark is the HUMAN MAP. Hz remains storage (heritage frames). Runtime and
compiled-v1 are unaffected.

Conversion: Traunmueller (1990).
    bark = 26.81 * f / (1960 + f) - 0.53
    f    = 1960 * (bark + 0.53) / (26.28 - bark)
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

_REPO_ROOT = Path(__file__).resolve().parents[2]
TABLES_DEFAULT = _REPO_ROOT / "docs" / "sonic_tables" / "tables.json"
LANDMARKS_DEFAULT = _REPO_ROOT / "authoring" / "ui" / "bark_landmarks.json"


def hz_to_bark(hz):
    hz = np.asarray(hz, dtype=float)
    return 26.81 * hz / (1960.0 + hz) - 0.53


def bark_to_hz(bark):
    bark = np.asarray(bark, dtype=float)
    return 1960.0 * (bark + 0.53) / (26.28 - bark)


def load_landmarks(tables_path: Path = TABLES_DEFAULT) -> list[dict]:
    """Full landmark pack from tables.json (landmark entries + vowel formants)."""
    data = json.loads(Path(tables_path).read_text(encoding="utf-8"))
    out: list[dict] = []
    for entry in data.get("landmarks", {}).get("entries", []):
        freq = float(entry["freq_hz"])
        out.append({
            "name": entry["name"],
            "freq_hz": freq,
            "bark": float(hz_to_bark(freq)),
            "source": "landmarks.entries",
        })
    for vowel_id, v in data.get("vowels", {}).items():
        for f_key in ("f1", "f2", "f3", "f4"):
            if f_key in v:
                freq = float(v[f_key])
                out.append({
                    "name": f"{vowel_id}.{f_key}",
                    "freq_hz": freq,
                    "bark": float(hz_to_bark(freq)),
                    "source": f"vowels.{vowel_id}",
                })
    return out


def plot_grid(landmarks_path: Path = LANDMARKS_DEFAULT) -> list[dict]:
    """Curated grid for Bark-axis plotting, loaded from the canonical
    landmark map (authoring/ui/bark_landmarks.json).

    Hz values in the file are treated as canonical; Bark is always
    recomputed from Hz via Traunmueller so the plot axis stays
    internally consistent (the file's `bark` field may be rounded).
    """
    data = json.loads(Path(landmarks_path).read_text(encoding="utf-8"))
    return [
        {
            "id": row["id"],
            "label": row["label"],
            "freq_hz": float(row["hz"]),
            "bark": float(hz_to_bark(row["hz"])),
            "role": row.get("role", ""),
        }
        for row in data.get("rows", [])
    ]


def main() -> int:
    """CLI sanity check: print the landmark grid and conversion samples."""
    print(f"landmark map ({LANDMARKS_DEFAULT.name}):")
    for g in plot_grid():
        print(f"  {g['id']:12s}  {g['label']:14s}  {g['freq_hz']:7.1f} Hz   "
              f"bark {g['bark']:6.2f}   -- {g['role']}")
    print()
    print("samples:")
    for hz in (60, 125, 500, 1000, 2000, 5000, 10000, 20000):
        b = float(hz_to_bark(hz))
        print(f"  {hz:5d} Hz -> bark {b:6.2f}   round-trip -> {float(bark_to_hz(b)):7.1f} Hz")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
