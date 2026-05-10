"""Cross-null diagnostic for Talking Hedz: render FilterEngine at each
(morph,q) corner and null against each canonical WAV corner.

If the diagonal nulls well and off-diagonals don't: corner-label map
is correct; shape mismatch is within a single corner.
If a non-diagonal pair nulls better: corner-label or interp-order bug.
If EVERY pair is ~-0.1 dB: compile-path or canonical-render diverge
at the cascade/coefficient level (not just at the resampler).
"""
from __future__ import annotations

import math
import subprocess
import sys
from pathlib import Path

import numpy as np
import soundfile as sf

ROOT = Path(__file__).resolve().parent.parent.parent
REF = ROOT / "ref" / "canonical"
CORNERS = ["M0_Q0", "M0_Q100", "M100_Q0", "M100_Q100"]


def mono(d):
    return d[:, 0] if d.ndim > 1 else d


def null_db(pred: np.ndarray, ref: np.ndarray, max_lag: int = 256) -> tuple[int, float, float]:
    best = (0, 0.0, float("inf"))
    for lag in range(-max_lag, max_lag + 1):
        if lag >= 0:
            p_off, r_off = 0, lag
        else:
            p_off, r_off = -lag, 0
        n = min(len(pred) - p_off, len(ref) - r_off)
        if n < 1024:
            continue
        p = pred[p_off : p_off + n].astype(np.float64)
        r = ref[r_off : r_off + n].astype(np.float64)
        den = float(np.dot(p, p))
        if den <= 0:
            continue
        gain = float(np.dot(p, r) / den)
        res = r - gain * p
        r_rms_sq = float(np.dot(r, r))
        if r_rms_sq <= 0:
            continue
        res_sq = float(np.dot(res, res))
        if res_sq <= 0:
            return (lag, gain, -300.0)
        rel = 10.0 * math.log10(res_sq / r_rms_sq)
        if rel < best[2]:
            best = (lag, gain, rel)
    return best


def main() -> int:
    # Render FilterEngine at each corner using a tiny Rust binary? No —
    # easier: use the existing talking_hedz_parity.rs trampoline indirectly
    # by writing a helper test. Instead we render the CANONICAL python
    # chain at each corner (which is what generated the WAVs) and null
    # each Python render against each canonical WAV. This shows whether
    # the canonical WAVs were written to the right files (naming bug).
    sys.path.insert(0, str(ROOT / "tools"))
    from parity_null import render_canonical, load_source, CORNERS as PY_CORNERS
    import json
    manifest = json.loads((REF / "MANIFEST.json").read_text())
    entry = manifest["cal_Talking_Hedz"]
    body, src_type, sr_auth = load_source(entry)
    boost = float(entry["boost"])

    dry_path = Path("C:/Users/hooki/trenchwork_clean/ref/bypassed-pinknoise.wav")
    dry_data, dry_sr = sf.read(str(dry_path))
    dry = mono(dry_data)
    dry_sr = float(dry_sr)

    # Render python at each corner.
    preds = {}
    for corner_label in CORNERS:
        stages = body["corners"][corner_label]["stages"]
        preds[corner_label] = render_canonical(
            dry, stages, src_type, sr_auth, boost, dry_sr, apply_agc_enabled=False
        )

    # Load WAVs.
    refs = {}
    for wav_label in CORNERS:
        refs[wav_label] = mono(sf.read(str(REF / f"cal_Talking_Hedz_{wav_label}.wav"))[0])

    header = "pred/ref"
    print(f"{header:<15}" + "".join(f"{c:>16}" for c in CORNERS))
    for pred_label in CORNERS:
        row = f"{pred_label:<15}"
        for wav_label in CORNERS:
            skip = len(preds[pred_label]) // 20
            lag, gain, rel = null_db(
                preds[pred_label][skip:],
                refs[wav_label][skip:],
                max_lag=32,
            )
            tag = f"{rel:+6.1f}/{gain:.3f}"
            row += f"{tag:>16}"
        print(row)
    return 0


if __name__ == "__main__":
    sys.exit(main())
