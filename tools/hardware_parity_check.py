"""Hardware parity check — compare a hardware capture of a P2K filter
against the Rust/Python engine's predicted magnitude response.

This is the first check in the gate pyramid that compares against
ground truth (real hardware) rather than Rust-vs-Rust or
Python-vs-Python internal consistency. Per the project conversation:
every prior gate is tautological — they verify our code agrees with
our code, not that our code agrees with E-mu hardware.

Method:

  1. Load the hardware capture WAV (mono, assumed stationary input of
     known spectral shape — default pink noise per user guidance).
  2. Compute its magnitude spectrum via Welch PSD.
  3. De-slope the assumed input spectrum out, yielding an estimate
     of |H(f)| as heard from hardware.
  4. Load the matching compiled-v1 cartridge. Interpolate the 4
     corners at the captured (morph, Q) using the same Q-first then
     morph-axis bilinear as the runtime.
  5. Evaluate the predicted cascade magnitude |H(e^jω)| at the same
     frequencies.
  6. Normalize both curves to their value at a fixed mid-band anchor
     (default 1 kHz) so we compare SHAPE, not absolute gain — the
     hardware capture's absolute level depends on the test rig's
     input trim and amp gain we don't know.
  7. Report: mean dB error, max dB error, Pearson correlation in dB
     space. Optionally write a PNG overlay plot if matplotlib is
     available.

Exit code is 0 on success (metrics within tolerance), non-zero on
failure. Tolerances are conservative defaults that pass a
reasonable-quality match and fail obvious divergence.

Usage:

    python tools/hardware_parity_check.py \\
        --capture <hardware.wav> \\
        --cartridge <compiled-v1.json> \\
        --morph 0.5 --q 0.5 \\
        [--input pink|white] [--plot out.png] [--tol-mean-db 6]
"""
from __future__ import annotations

import argparse
import json
import math
import struct
import sys
import wave
from pathlib import Path
from typing import Sequence

import numpy as np


# ── Anchor frequency and default tolerances ──────────────────────────

ANCHOR_HZ = 1000.0
DEFAULT_TOL_MEAN_DB = 8.0
DEFAULT_TOL_MAX_DB = 18.0
DEFAULT_TOL_CORR = 0.70

# Band we actually compare over. Below 40 Hz the test-rig noise floor
# dominates; above 16 kHz the hardware's input AA filter and capture
# interface colour the signal. Stick to the meat.
COMPARE_LO_HZ = 60.0
COMPARE_HI_HZ = 16000.0


def load_wav_mono(path: Path) -> tuple[np.ndarray, int]:
    # Stdlib wave handles PCM16 but chokes on IEEE float (fmt tag 3), which
    # Trench's hardware captures ship as. Fall back to soundfile/scipy for
    # anything it rejects — keeps the fast path allocation-free on normal
    # PCM captures.
    try:
        with wave.open(str(path), "rb") as w:
            n = w.getnframes()
            sr = w.getframerate()
            ch = w.getnchannels()
            sw = w.getsampwidth()
            raw = w.readframes(n)
        if sw != 2:
            raise SystemExit(f"unsupported sample width {sw} bytes in {path}")
        fmt = "<" + "h" * (n * ch)
        flat = np.array(struct.unpack(fmt, raw), dtype=np.float32) / 32768.0
        if ch > 1:
            flat = flat.reshape(n, ch).mean(axis=1)
        return flat, sr
    except wave.Error:
        try:
            import soundfile as sf
        except ImportError:
            try:
                from scipy.io import wavfile as _wavfile
            except ImportError as exc:
                raise SystemExit(
                    f"non-PCM WAV at {path}; install soundfile or scipy"
                ) from exc
            sr, data = _wavfile.read(str(path))
            flat = np.asarray(data, dtype=np.float32)
            if flat.dtype.kind in "iu":
                flat = flat / float(np.iinfo(flat.dtype).max)
            if flat.ndim > 1:
                flat = flat.mean(axis=1)
            return flat.astype(np.float32), int(sr)
        data, sr = sf.read(str(path), always_2d=False)
        flat = np.asarray(data, dtype=np.float32)
        if flat.ndim > 1:
            flat = flat.mean(axis=1)
        return flat, int(sr)


def welch_psd(xs: np.ndarray, sr: int, nfft: int = 4096, hop: int = 2048) -> tuple[np.ndarray, np.ndarray]:
    """Segmented Hann-windowed PSD, averaged across segments. Returns
    (freqs_hz, psd) where psd is power per bin (not per Hz — we're
    only comparing shapes)."""
    window = np.hanning(nfft).astype(np.float32)
    norm = float(np.sum(window * window))
    freqs = np.fft.rfftfreq(nfft, d=1.0 / sr)
    acc = np.zeros(freqs.size, dtype=np.float64)
    segs = 0
    for start in range(0, xs.size - nfft + 1, hop):
        frame = xs[start : start + nfft] * window
        spec = np.fft.rfft(frame)
        acc += (spec.real * spec.real + spec.imag * spec.imag) / norm
        segs += 1
    if segs == 0:
        raise SystemExit("capture too short for Welch PSD")
    return freqs, (acc / segs).astype(np.float32)


def input_slope_db(freqs_hz: np.ndarray, model: str) -> np.ndarray:
    """dB correction to subtract from the output spectrum to recover
    |H(f)|^2 from |X(f)|^2 * |N(f)|^2. Pink = -3 dB/octave slope;
    white = flat; impulse = flat (an impulse has flat spectrum)."""
    safe = np.maximum(freqs_hz, 1.0)
    if model == "pink":
        return -10.0 * np.log10(safe)
    if model in ("white", "impulse"):
        return np.zeros_like(freqs_hz)
    raise SystemExit(f"unknown input model: {model}")


def impulse_response_spectrum(xs: np.ndarray, sr: int,
                              trim_silence_db: float = -80.0
                              ) -> tuple[np.ndarray, np.ndarray]:
    """For impulse-response captures: trim trailing silence, take the
    FFT magnitude of the useful prefix. The filter excitation is short
    (< ~200 ms) and the capture's tail is ~0, so Welch-averaging over
    fixed segments drowns the IR in zeros. Direct FFT of the whole
    non-silent portion recovers |H(f)| up to a flat input spectrum."""
    abs_xs = np.abs(xs)
    threshold = np.max(abs_xs) * (10.0 ** (trim_silence_db / 20.0))
    # Last index above threshold (treat everything past that as silence).
    nonzero = np.where(abs_xs > threshold)[0]
    if nonzero.size == 0:
        raise SystemExit("capture is silent")
    last = int(nonzero[-1]) + 1
    useful = xs[:last]
    # Next power of 2 for FFT friendliness
    n = 1 << (useful.size - 1).bit_length()
    padded = np.zeros(n, dtype=np.float32)
    padded[: useful.size] = useful
    spec = np.fft.rfft(padded)
    freqs = np.fft.rfftfreq(n, d=1.0 / sr)
    return freqs, np.abs(spec).astype(np.float32)


# ── Cartridge load + corner interpolation ────────────────────────────

def load_compiled_v1(path: Path) -> dict:
    d = json.loads(path.read_text(encoding="utf-8"))
    if d.get("format") != "compiled-v1":
        raise SystemExit(f"not compiled-v1: {d.get('format')!r}")
    kfs = d.get("keyframes") or []
    if len(kfs) != 4:
        raise SystemExit(f"compiled-v1 needs 4 keyframes, got {len(kfs)}")
    by_label = {kf["label"]: kf for kf in kfs}
    for lab in ("M0_Q0", "M100_Q0", "M0_Q100", "M100_Q100"):
        if lab not in by_label:
            raise SystemExit(f"keyframe {lab} missing")
    return {
        "boost": float(by_label["M0_Q0"].get("boost", 1.0)),
        "corners": {
            lab: [
                (s["c0"], s["c1"], s["c2"], s["c3"], s["c4"])
                for s in by_label[lab]["stages"]
            ]
            for lab in ("M0_Q0", "M100_Q0", "M0_Q100", "M100_Q100")
        },
    }


def interpolate_corners(cart: dict, morph: float, q: float) -> list[tuple[float, ...]]:
    """Q-first, then morph-axis bilinear on the 12-stage kernel
    coefficients. Mirrors cartridge.rs:370-440."""
    corners = cart["corners"]
    n_stages = len(corners["M0_Q0"])
    out: list[tuple[float, ...]] = []
    for s in range(n_stages):
        c00 = np.array(corners["M0_Q0"][s])     # M0, Q0
        c10 = np.array(corners["M100_Q0"][s])   # M100, Q0
        c01 = np.array(corners["M0_Q100"][s])   # M0, Q100
        c11 = np.array(corners["M100_Q100"][s]) # M100, Q100
        # Q first
        q0_blend = c00 + morph * (c10 - c00)
        q1_blend = c01 + morph * (c11 - c01)
        final = q0_blend + q * (q1_blend - q0_blend)
        out.append(tuple(float(v) for v in final))
    return out


def cascade_magnitude(stages: Sequence[tuple[float, ...]],
                      freqs_hz: np.ndarray, sr: float,
                      boost: float = 1.0) -> np.ndarray:
    """|H(e^jω)| for the 12-stage kernel cascade. Each stage's kernel
    coefficients (c0..c4) unpack to a DF2T second-order section
    whose transfer function is:

        H(z) = (c0 + c1 z^-1 + c2 z^-2) / (1 + c3 z^-1 + c4 z^-2)

    See cascade.rs DF2T recurrence."""
    w = 2.0 * np.pi * freqs_hz / sr
    z1 = np.exp(-1j * w)
    z2 = np.exp(-2j * w)
    mag = np.full(freqs_hz.size, boost, dtype=np.float64)
    for c0, c1, c2, c3, c4 in stages:
        num = c0 + c1 * z1 + c2 * z2
        den = 1.0 + c3 * z1 + c4 * z2
        stage_mag = np.abs(num / den)
        mag = mag * stage_mag
    return mag.astype(np.float32)


# ── Compare ──────────────────────────────────────────────────────────

def compare_shapes(
    freqs: np.ndarray,
    hw_db: np.ndarray,
    engine_db: np.ndarray,
    lo_hz: float,
    hi_hz: float,
    anchor_hz: float,
) -> dict:
    band = (freqs >= lo_hz) & (freqs <= hi_hz)
    fb = freqs[band]
    hb = hw_db[band]
    eb = engine_db[band]

    anchor_idx = int(np.argmin(np.abs(fb - anchor_hz)))
    hb_n = hb - hb[anchor_idx]
    eb_n = eb - eb[anchor_idx]

    err = hb_n - eb_n
    abs_err = np.abs(err)
    corr = float(np.corrcoef(hb_n, eb_n)[0, 1])

    return {
        "mean_err_db": float(np.mean(abs_err)),
        "max_err_db": float(np.max(abs_err)),
        "corr": corr,
        "n_points": int(fb.size),
        "freqs_compared": (float(fb[0]), float(fb[-1])),
        "_arrays": {"freqs": fb, "hw_db": hb_n, "eng_db": eb_n},
    }


def maybe_plot(result: dict, out_path: Path, title: str) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return
    fb = result["_arrays"]["freqs"]
    hb = result["_arrays"]["hw_db"]
    eb = result["_arrays"]["eng_db"]
    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.semilogx(fb, hb, label="hardware", alpha=0.85, linewidth=1.3)
    ax.semilogx(fb, eb, label="engine",   alpha=0.85, linewidth=1.3)
    ax.set_xlabel("Hz")
    ax.set_ylabel("dB (relative to 1 kHz)")
    ax.set_title(title)
    ax.legend()
    ax.grid(True, which="both", alpha=0.3)
    ax.set_xlim (fb[0], fb[-1])
    fig.tight_layout()
    fig.savefig(out_path, dpi=110)
    plt.close(fig)


# ── CLI ──────────────────────────────────────────────────────────────

def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--capture",   type=Path, required=True)
    ap.add_argument("--cartridge", type=Path, required=True)
    ap.add_argument("--morph", type=float, required=True)
    ap.add_argument("--q",     type=float, required=True)
    ap.add_argument("--input", choices=("pink", "white", "impulse"),
                    default="impulse",
                    help="assumed input signal spectrum (default: impulse — "
                         "for filter impulse-response captures that ring out "
                         "and go silent)")
    ap.add_argument("--plot", type=Path, default=None,
                    help="optional PNG output path for spectrum overlay")
    ap.add_argument("--tol-mean-db", type=float, default=DEFAULT_TOL_MEAN_DB)
    ap.add_argument("--tol-max-db",  type=float, default=DEFAULT_TOL_MAX_DB)
    ap.add_argument("--tol-corr",    type=float, default=DEFAULT_TOL_CORR)
    ap.add_argument("--coeff-sr",    type=float, default=None,
                    help="sample rate to evaluate the filter at. Defaults to "
                         "the capture SR (which matches what hardware produces "
                         "and what the plugin hears at host SR). Do NOT trust "
                         "the cartridge's declared sampleRate field — it's "
                         "stale authoring metadata, not the execution rate.")
    args = ap.parse_args(argv)

    xs, sr_cap = load_wav_mono(args.capture)
    if args.input == "impulse":
        freqs, mag = impulse_response_spectrum(xs, sr_cap)
        hw_db = 20.0 * np.log10(np.maximum(mag, 1e-20))
    else:
        freqs, psd = welch_psd(xs, sr_cap)
        hw_db = 10.0 * np.log10(np.maximum(psd, 1e-20))
    # Divide out assumed input spectrum (dB).
    hw_db = hw_db - input_slope_db(freqs, args.input)

    cart = load_compiled_v1(args.cartridge)
    stages = interpolate_corners(cart, args.morph, args.q)

    # Use capture SR by default. The JSON's `sampleRate` field is stale
    # authoring metadata; hardware parity tests against Talking Hedz at
    # all 4 captured positions show near-unity correlation when evaluating
    # at capture SR, and ~0.7 when using the declared 39062.5 — see
    # the memory entry on this.
    coeff_sr = args.coeff_sr or float(sr_cap)
    eng_mag = cascade_magnitude(stages, freqs, coeff_sr, cart.get("boost", 1.0))
    eng_db = 20.0 * np.log10(np.maximum(eng_mag, 1e-20))

    result = compare_shapes(freqs, hw_db, eng_db,
                            COMPARE_LO_HZ, COMPARE_HI_HZ, ANCHOR_HZ)

    if args.plot is not None:
        title = (f"{args.capture.name}  |  M={args.morph:.2f} Q={args.q:.2f}"
                 f"  |  mean={result['mean_err_db']:.1f} dB"
                 f"  corr={result['corr']:.2f}")
        maybe_plot(result, args.plot, title)

    print(json.dumps({
        "capture":   str(args.capture),
        "cartridge": str(args.cartridge),
        "morph": args.morph, "q": args.q,
        "input_model": args.input,
        "sr_capture": sr_cap,
        "compare_band_hz": [COMPARE_LO_HZ, COMPARE_HI_HZ],
        "mean_err_db": round(result["mean_err_db"], 3),
        "max_err_db":  round(result["max_err_db"], 3),
        "correlation": round(result["corr"], 4),
        "n_points": result["n_points"],
    }, indent=2))

    ok = (result["mean_err_db"] <= args.tol_mean_db
          and result["max_err_db"]  <= args.tol_max_db
          and result["corr"]        >= args.tol_corr)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
