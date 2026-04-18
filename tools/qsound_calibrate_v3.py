"""Diagnostic: per-channel RMS per pan, from the raw file."""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
from scipy.io import wavfile

DATA_DIR = Path(
    r"C:/Users/hooki/trench_re_vault/datasets/qsound_spatial_v1"
    r"/2026-03-05_pan_batch_fullgrid_a"
)
PAN_VALUES = list(range(-60, 70, 10))

print(f"{'pan':>5}  {'L_rms':>10}  {'R_rms':>10}  {'R-L dB':>8}  {'total_rms':>10}  {'N':>6}")
for pan in PAN_VALUES:
    sign = "+" if pan >= 0 else "-"
    path = DATA_DIR / f"pan_{sign}{abs(pan):03d}.wav"
    sr, data = wavfile.read(str(path))
    data = data.astype(np.float64)
    # skip transient
    s = int(0.25 * sr)
    e = min(s + int(1.0 * sr), data.shape[0])
    w = data[s:e]
    L = w[:, 0]; R = w[:, 1]
    Lr = float(np.sqrt(np.mean(L * L)))
    Rr = float(np.sqrt(np.mean(R * R)))
    Tr = float(np.sqrt(np.mean((L + R) * (L + R) / 2)))
    drb = 20 * math.log10(Rr / (Lr + 1e-20))
    print(f"{pan:>5}  {Lr:>10.6f}  {Rr:>10.6f}  {drb:>+8.3f}  {Tr:>10.6f}  {len(L):>6}")
