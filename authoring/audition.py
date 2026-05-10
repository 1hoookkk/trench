"""forge/audition.py — listen to a cartridge.

The canonical verification pipeline. Cartridges are authored at 39062.5 Hz
(E-MU silicon clock = 10 MHz / 256). Their coefficients only describe the
intended transfer function at THAT sample rate, so any source must be
resampled to 39062.5 Hz before the cascade runs, then the output resamples
back to the playback SR.

Pipeline (matches reference/canonical_audio/PROVENANCE.json):

    source WAV (any SR)
      → resample to 39062.5 Hz
      → SOS cascade  (scipy.signal.sosfilt with rows [c0, c1, c2, 1, c3, c4])
      → AGC          (per-sample 16-entry soft-limit table, from agc.rs)
      → * boost      (cartridge-level scalar; heritage default 4.0)
      → resample to 44100 Hz
      → float32 WAV

Outputs four WAVs per cartridge — M0_Q0, M0_Q100, M100_Q0, M100_Q100 — to
cartridges/factory/generated/qlaw/audition/. Drop those into your DAW or
play them in Foobar/VLC; that's the ground truth your ears compare against
reference/canonical_audio/cal_*.wav.

Usage:
    python forge/audition.py <cartridge.json>
    python forge/audition.py <cartridge.json> --source <path.wav>
    python forge/audition.py <cartridge.json> --source pinknoise --duration 6
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import soundfile as sf
from scipy.signal import resample as fft_resample, sosfilt


REPO = Path(__file__).resolve().parent.parent
DEFAULT_DRY = Path("C:/Users/hooki/trenchwork_clean/ref/bypassed-pinknoise.wav")
OUT_DIR = REPO / "cartridges" / "factory" / "generated" / "qlaw" / "audition"

AUTHORING_SR = 39062.5
PLAYBACK_SR = 44100
CORNERS = ("M0_Q0", "M0_Q100", "M100_Q0", "M100_Q100")

# AGC_TABLE ported from runtime/trench-core/src/agc.rs (16-entry soft-limit).
AGC_TABLE = np.array(
    [1.0001, 1.0001, 0.996, 0.990, 0.920, 0.500, 0.200, 0.160,
     0.120, 0.120, 0.120, 0.120, 0.120, 0.120, 0.120, 0.120],
    dtype=np.float64,
)


def to_mono(x: np.ndarray) -> np.ndarray:
    return x[:, 0] if x.ndim > 1 else x


def resample_to(sig: np.ndarray, sr_in: float, sr_out: float) -> np.ndarray:
    if sr_in == sr_out:
        return sig.astype(np.float64, copy=False)
    n_out = int(round(len(sig) * sr_out / sr_in))
    return fft_resample(sig.astype(np.float64, copy=False), n_out)


def synth_pink_noise(duration_s: float, sr: float, seed: int = 42) -> np.ndarray:
    """Voss-McCartney pink noise, normalised to 0.3 RMS so the cartridge
    has headroom before AGC slams it."""
    rng = np.random.default_rng(seed)
    n = int(round(duration_s * sr))
    # 16 octaves of summed white noise, each updated half as often as the previous.
    n_rows = 16
    cols = np.zeros((n_rows, n), dtype=np.float64)
    cols[0] = rng.standard_normal(n)
    for k in range(1, n_rows):
        period = 1 << k
        idx = np.arange(n) // period
        unique = rng.standard_normal(idx.max() + 1)
        cols[k] = unique[idx]
    pink = cols.sum(axis=0)
    rms = float(np.sqrt(np.mean(pink ** 2)))
    if rms > 0:
        pink *= 0.3 / rms
    return pink


def load_source(path_or_keyword: str | Path, duration_s: float) -> tuple[np.ndarray, float]:
    """Load source WAV or synthesize one. Returns (mono samples, sr)."""
    if str(path_or_keyword).lower() == "pinknoise":
        return synth_pink_noise(duration_s, AUTHORING_SR), AUTHORING_SR
    p = Path(path_or_keyword)
    if not p.exists():
        # Fallback: synthesize if canonical pinknoise is missing.
        print(f"[warn] source {p} not found — synthesizing pink noise instead.", file=sys.stderr)
        return synth_pink_noise(duration_s, AUTHORING_SR), AUTHORING_SR
    data, sr = sf.read(str(p), always_2d=False)
    return to_mono(np.asarray(data, dtype=np.float64)), float(sr)


def apply_agc(samples: np.ndarray) -> np.ndarray:
    """Per-sample soft-limit AGC (matches runtime/trench-core/src/agc.rs).
    Slow O(n) Python loop; fine for audition lengths."""
    out = np.empty_like(samples)
    gain = 1.0
    for i in range(samples.size):
        s = samples[i]
        idx = int(gain * abs(s)) & 0xF
        new_gain = gain * AGC_TABLE[idx]
        gain = new_gain if new_gain < 1.0 else 1.0
        out[i] = s * gain
    return out


def assemble_sos(stages: list[dict]) -> np.ndarray | None:
    """compiled-v1 stages → scipy SOS rows. Drop passthroughs (c3=c4=0)."""
    rows = []
    for s in stages:
        c0 = float(s["c0"]); c1 = float(s["c1"]); c2 = float(s["c2"])
        c3 = float(s["c3"]); c4 = float(s["c4"])
        if c3 == 0.0 and c4 == 0.0:
            continue
        rows.append([c0, c1, c2, 1.0, c3, c4])
    if not rows:
        return None
    return np.asarray(rows, dtype=np.float64)


def heritage_to_compiled(cart: dict) -> dict:
    """Convert reference/calibration/*.json heritage schema → compiled-v1.

    Heritage stages carry pole_freq_hz/radius/val1/val2/val3; we reconstruct
    a1 = -2r·cos(2π·f/sr) and produce (c0, c1, c2, c3, c4) per the canonical
    SOS pipeline (parity_null.py): b0 = 1 + val1, b1 = a1 + val2, b2 = a2 - val3.
    """
    import math
    sr = float(cart.get("sample_rate", AUTHORING_SR))
    boost = float(cart.get("boost", 1.0))
    keyframes = []
    for label, corner in cart.get("corners", {}).items():
        stages_in = corner.get("stages", [])
        stages_out = []
        for s in stages_in:
            r = min(float(s.get("radius", 0.0)), 0.999999)
            f = float(s.get("pole_freq_hz", 0.0))
            if r == 0.0 or f == 0.0:
                stages_out.append({"c0": 1.0, "c1": 0.0, "c2": 0.0, "c3": 0.0, "c4": 0.0})
                continue
            a1 = -2.0 * r * math.cos(2.0 * math.pi * f / sr)
            a2 = r * r
            val1 = float(s.get("val1", 0.0))
            val2 = float(s.get("val2", 0.0))
            val3 = float(s.get("val3", 0.0))
            stages_out.append({
                "c0": 1.0 + val1,
                "c1": a1 + val2,
                "c2": a2 - val3,
                "c3": a1,
                "c4": a2,
            })
        while len(stages_out) < 12:
            stages_out.append({"c0": 1.0, "c1": 0.0, "c2": 0.0, "c3": 0.0, "c4": 0.0})
        morph = 0.0 if "M0" in label else 1.0
        q = 0.0 if "Q0" in label else 1.0
        keyframes.append({
            "label": label, "morph": morph, "q": q,
            "boost": boost, "stages": stages_out,
        })
    return {
        "format": "compiled-v1",
        "name": cart.get("name", "heritage"),
        "sampleRate": sr,
        "boost": boost,
        "keyframes": keyframes,
    }


def normalize_cartridge(cart: dict) -> dict:
    """Auto-detect schema. Heritage has 'corners'; compiled-v1 has 'keyframes'."""
    if "keyframes" in cart:
        return cart
    if "corners" in cart:
        return heritage_to_compiled(cart)
    raise ValueError("cartridge has neither 'keyframes' (compiled-v1) nor 'corners' (heritage)")


def find_corner(cart: dict, label: str) -> dict | None:
    for kf in cart.get("keyframes", []):
        if kf.get("label") == label:
            return kf
    return None


def render_corner(cart: dict, label: str, dry_at_authoring: np.ndarray) -> np.ndarray | None:
    kf = find_corner(cart, label)
    if kf is None:
        return None
    sos = assemble_sos(kf.get("stages", []))
    if sos is None:
        return None
    boost = float(kf.get("boost", cart.get("boost", 1.0)))
    cascaded = sosfilt(sos, dry_at_authoring)
    agc = apply_agc(cascaded)
    return agc * boost


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("cartridge", type=Path, help="path to compiled-v1 cartridge JSON")
    ap.add_argument("--source", type=str, default=str(DEFAULT_DRY),
                    help="source WAV path, or 'pinknoise' to synthesize. "
                         f"Default: {DEFAULT_DRY}")
    ap.add_argument("--duration", type=float, default=6.0,
                    help="duration of synthesized source (seconds). default 6.0")
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR,
                    help=f"output directory. default {OUT_DIR}")
    args = ap.parse_args(argv)

    if not args.cartridge.exists():
        print(f"cartridge not found: {args.cartridge}", file=sys.stderr)
        return 1
    cart_raw = json.loads(args.cartridge.read_text(encoding="utf-8"))
    cart = normalize_cartridge(cart_raw)
    cart_name = cart.get("name", args.cartridge.stem)

    # 1. Load source at its native SR.
    src, src_sr = load_source(args.source, args.duration)
    print(f"source: {args.source}  ({len(src)} samples @ {src_sr:.1f} Hz)")

    # 2. Resample to authoring SR (39062.5 Hz).
    dry_authored = resample_to(src, src_sr, AUTHORING_SR)
    print(f"resampled to authoring SR: {len(dry_authored)} samples @ {AUTHORING_SR:.1f} Hz")

    # 3. Render each corner.
    args.out_dir.mkdir(parents=True, exist_ok=True)
    written = []
    for label in CORNERS:
        out_authored = render_corner(cart, label, dry_authored)
        if out_authored is None:
            print(f"  {label}: corner missing in cartridge — skipped")
            continue
        # 4. Resample to playback SR (44100 Hz).
        out_play = resample_to(out_authored, AUTHORING_SR, PLAYBACK_SR)
        # 5. Write float32 WAV.
        out_path = args.out_dir / f"{cart_name}_{label}.wav"
        sf.write(str(out_path), out_play.astype(np.float32),
                 PLAYBACK_SR, subtype="FLOAT")
        peak_db = 20.0 * np.log10(max(float(np.max(np.abs(out_play))), 1e-12))
        rms_db = 20.0 * np.log10(max(float(np.sqrt(np.mean(out_play ** 2))), 1e-12))
        print(f"  {label}: peak {peak_db:+5.1f} dB   rms {rms_db:+5.1f} dB   → {out_path.name}")
        written.append(out_path)

    if not written:
        print("nothing written.")
        return 1
    print(f"\n{len(written)} WAV(s) written to {args.out_dir}")
    print("listen and decide.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
