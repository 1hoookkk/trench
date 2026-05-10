"""Null: FilterEngine dump vs canonical WAVs vs Python renders.

Three rails:
  - ENG:  target/hedz_engine_{label}.f32       (Rust FilterEngine @ 44100)
  - CAN:  ref/canonical/cal_Talking_Hedz_{label}.wav  (Python canonical)
  - PY:   render_canonical() at that corner    (fresh Python render)

Diagonal nulls CAN<->PY at -300 dB (confirmed in crossnull_diag).
Where does ENG diverge? — this script measures it.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
import soundfile as sf

ROOT = Path(__file__).resolve().parent.parent.parent
TARGET = ROOT / "target"
REF = ROOT / "ref" / "canonical"
sys.path.insert(0, str(ROOT / "tools"))
from parity_null import render_canonical, load_source  # noqa


def load_f32(p: Path) -> np.ndarray:
    return np.frombuffer(p.read_bytes(), dtype="<f4").astype(np.float32, copy=True)


def mono(d):
    return d[:, 0] if d.ndim > 1 else d


def best_null(pred, ref, max_lag=256):
    best = (0, 0.0, float("inf"))
    for lag in range(-max_lag, max_lag + 1):
        if lag >= 0:
            po, ro = 0, lag
        else:
            po, ro = -lag, 0
        n = min(len(pred) - po, len(ref) - ro)
        if n < 1024:
            continue
        p = pred[po:po + n].astype(np.float64)
        r = ref[ro:ro + n].astype(np.float64)
        den = float(np.dot(p, p))
        if den <= 0:
            continue
        gain = float(np.dot(p, r) / den)
        res = r - gain * p
        rms2 = float(np.dot(res, res))
        rr = float(np.dot(r, r))
        if rr <= 0:
            continue
        rel = -300.0 if rms2 <= 0 else 10 * math.log10(rms2 / rr)
        if rel < best[2]:
            best = (lag, gain, rel)
    return best


def main() -> int:
    manifest = json.loads((REF / "MANIFEST.json").read_text())
    entry = manifest["cal_Talking_Hedz"]
    body, src_type, sr_auth = load_source(entry)
    boost = float(entry["boost"])

    dry_path = Path(r"C:/Users/hooki/trenchwork_clean/ref/bypassed-pinknoise.wav")
    dry_data, dry_sr = sf.read(str(dry_path))
    dry = mono(dry_data)
    dry_sr = float(dry_sr)

    skip = len(dry) // 20

    for label in ("M0_Q0", "M100_Q100"):
        eng = load_f32(TARGET / f"hedz_engine_{label}.f32")
        can = mono(sf.read(str(REF / f"cal_Talking_Hedz_{label}.wav"))[0])
        stages = body["corners"][label]["stages"]
        py = render_canonical(dry, stages, src_type, sr_auth, boost, dry_sr, apply_agc_enabled=False)

        print(f"--- {label} ---")
        n_min = min(len(eng), len(can), len(py))
        print(f"   samples: eng={len(eng)} can={len(can)} py={len(py)} min={n_min}")
        # ENG vs CAN
        lag, gain, rel = best_null(eng[skip:n_min], can[skip:n_min])
        print(f"   ENG vs CAN: lag={lag:+3} gain={gain:+.4f} null={rel:+7.2f} dB")
        # ENG vs PY
        lag, gain, rel = best_null(eng[skip:n_min], py[skip:n_min])
        print(f"   ENG vs PY : lag={lag:+3} gain={gain:+.4f} null={rel:+7.2f} dB")
        # CAN vs PY
        lag, gain, rel = best_null(can[skip:n_min], py[skip:n_min])
        print(f"   CAN vs PY : lag={lag:+3} gain={gain:+.4f} null={rel:+7.2f} dB  (should be -300)")

        # Peak sanity
        print(f"   peaks: eng={np.max(np.abs(eng)):.4f} can={np.max(np.abs(can)):.4f} py={np.max(np.abs(py)):.4f}")
        print(f"   rms : eng={np.sqrt(np.mean(eng.astype(np.float64)**2)):.4f} "
              f"can={np.sqrt(np.mean(can.astype(np.float64)**2)):.4f} "
              f"py={np.sqrt(np.mean(py.astype(np.float64)**2)):.4f}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
