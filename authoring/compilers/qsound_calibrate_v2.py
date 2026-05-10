"""
qsound_calibrate_v2.py

Refines v1 calibration findings:
  * sub-sample ITD via zero-padded (upsampled) cross-correlation
  * ILD using multiple definitions (RMS, peak, per-harmonic)
  * per-harmonic L/R magnitude ratios (bypasses crossover choice)
  * sanity probes: L-channel stability across pans, peak sample inspection

Goal: determine what QSound is *actually* doing at 220 Hz stimulus before
forcing a band model onto measurements that may be dominated by a few
harmonics.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
from scipy.io import wavfile
from scipy.signal import correlate


DATA_DIR = Path(
    r"C:/Users/hooki/trench_re_vault/datasets/qsound_spatial_v1"
    r"/2026-03-05_pan_batch_fullgrid_a"
)
SR = 48000
PAN_VALUES = list(range(-60, 70, 10))
STIM = 220.0
WIN_START = 0.25
WIN_DUR = 1.0


def load_window(pan: int):
    sign = "+" if pan >= 0 else "-"
    path = DATA_DIR / f"pan_{sign}{abs(pan):03d}.wav"
    sr, data = wavfile.read(str(path))
    assert sr == SR
    data = data.astype(np.float64)
    if np.issubdtype(wavfile.read(str(path))[1].dtype, np.integer):
        data /= np.iinfo(wavfile.read(str(path))[1].dtype).max
    start = int(WIN_START * sr)
    n = int(WIN_DUR * sr)
    data = data[start:start + n]
    return data[:, 0], data[:, 1]


def upsampled_xcorr_lag(L: np.ndarray, R: np.ndarray, upsample: int = 16,
                        max_lag_samples: int = 128) -> float:
    """Measure L/R lag to sub-sample precision via FFT upsampling of xcorr."""
    L = L - L.mean()
    R = R - R.mean()
    L /= np.linalg.norm(L) + 1e-20
    R /= np.linalg.norm(R) + 1e-20
    xc = correlate(R, L, mode="full")
    # upsample the xcorr via zero-pad FFT interpolation
    n = len(xc)
    X = np.fft.fft(xc)
    X_up = np.zeros(n * upsample, dtype=complex)
    half = n // 2
    X_up[:half] = X[:half]
    X_up[-half:] = X[-half:]
    xc_up = np.real(np.fft.ifft(X_up)) * upsample
    center = (len(L) - 1) * upsample
    lo = center - max_lag_samples * upsample
    hi = center + max_lag_samples * upsample + 1
    window = xc_up[lo:hi]
    peak_idx = int(np.argmax(np.abs(window))) - max_lag_samples * upsample
    return peak_idx / upsample  # in original samples, fractional


def phase_lag_at_harmonic(L: np.ndarray, R: np.ndarray, f: float,
                          sr: int = SR) -> tuple[float, float]:
    """Return (phase_deg, lag_samples) at frequency f, with bin-centered DFT."""
    n = len(L)
    # windowed DFT at specific frequency
    k = np.arange(n)
    win = np.hanning(n)
    twid = np.exp(-2j * np.pi * f * k / sr)
    XL = np.sum(L * win * twid)
    XR = np.sum(R * win * twid)
    # phase difference: R relative to L
    dphi = np.angle(XR / XL)  # radians, in (-pi, pi]
    phase_deg = math.degrees(dphi)
    # convert to group-delay-style lag at this freq
    lag = dphi / (2 * math.pi * f) * sr  # samples
    return phase_deg, lag


def harmonic_magnitudes(x: np.ndarray, sr: int = SR, base: float = STIM,
                        n_harmonics: int = 20) -> list[float]:
    n = len(x)
    X = np.fft.rfft(x * np.hanning(n))
    mag = np.abs(X)
    freqs = np.fft.rfftfreq(n, 1.0 / sr)
    out = []
    for k in range(1, n_harmonics + 1):
        f = base * k
        if f > sr / 2 - 50:
            break
        idx = int(round(f * n / sr))
        lo = max(0, idx - 3)
        hi = min(len(mag), idx + 4)
        # power = sum of squared mags in main lobe
        power = float((mag[lo:hi] ** 2).sum())
        out.append((f, power))
    return out


def rms(x): return float(np.sqrt(np.mean(x * x)))
def db(x): return 10.0 * math.log10(x + 1e-30)


def main():
    print("=== QSound calibration v2 ===\n")

    # Center reference
    Lc, Rc = load_window(0)
    refs = {
        "L_rms": rms(Lc),
        "R_rms": rms(Rc),
        "L_peak": float(np.max(np.abs(Lc))),
        "R_peak": float(np.max(np.abs(Rc))),
    }
    print(f"Center: L_rms={refs['L_rms']:.6f}  R_rms={refs['R_rms']:.6f}  "
          f"L/R ratio dB = {20*math.log10(refs['R_rms']/refs['L_rms']):+.3f}")
    print()

    # per-pan measurement
    print("Sub-sample ITD measurements (R leads when positive):")
    print(f"{'pan':>5}  {'xc_lag':>8}  {'phase_lag@220':>14}  "
          f"{'phase_lag@440':>14}  {'phase_lag@660':>14}")
    itd_table = []
    for pan in PAN_VALUES:
        L, R = load_window(pan)
        lag_xc = upsampled_xcorr_lag(L, R)
        p220, l220 = phase_lag_at_harmonic(L, R, 220)
        p440, l440 = phase_lag_at_harmonic(L, R, 440)
        p660, l660 = phase_lag_at_harmonic(L, R, 660)
        itd_table.append((pan, lag_xc, l220, l440, l660, p220, p440, p660))
        print(f"{pan:>5}  {lag_xc:>8.3f}  "
              f"{l220:>+8.3f}({p220:>+7.1f}°)  "
              f"{l440:>+8.3f}({p440:>+7.1f}°)  "
              f"{l660:>+8.3f}({p660:>+7.1f}°)")
    print()

    # per-pan per-harmonic L and R magnitudes (to see spectral shaping)
    print("Per-harmonic |R| - |L| (dB) per pan:")
    header = "  pan  "
    harms_to_show = [220, 440, 660, 880, 1100, 1760, 2640, 4400, 8800]
    header += "  ".join(f"{int(f)}Hz".rjust(8) for f in harms_to_show)
    print(header)
    per_harm_table = {}
    for pan in PAN_VALUES:
        L, R = load_window(pan)
        hL = dict((f, e) for f, e in harmonic_magnitudes(L))
        hR = dict((f, e) for f, e in harmonic_magnitudes(R))
        row = f"  {pan:>4}  "
        deltas = {}
        for f in harms_to_show:
            if f in hL and f in hR and hL[f] > 1e-20 and hR[f] > 1e-20:
                d = db(hR[f]) - db(hL[f])
                row += f"{d:>+7.2f}  "
                deltas[f] = d
            else:
                row += f"{'n/a':>8}  "
        per_harm_table[pan] = deltas
        print(row)
    print()

    # Per-channel per-harmonic dB, referenced to center case
    print("L-channel harmonic levels (dB) rel to center:")
    header = "  pan  " + "  ".join(f"{int(f)}Hz".rjust(8) for f in harms_to_show)
    print(header)
    Lc_harm = dict((f, e) for f, e in harmonic_magnitudes(Lc))
    Rc_harm = dict((f, e) for f, e in harmonic_magnitudes(Rc))
    L_dbtable = {}
    R_dbtable = {}
    for pan in PAN_VALUES:
        L, R = load_window(pan)
        hL = dict((f, e) for f, e in harmonic_magnitudes(L))
        row = f"  {pan:>4}  "
        row_l = {}
        for f in harms_to_show:
            if f in hL and hL[f] > 1e-20 and Lc_harm.get(f, 0) > 1e-20:
                d = db(hL[f]) - db(Lc_harm[f])
                row += f"{d:>+7.2f}  "
                row_l[f] = d
            else:
                row += f"{'n/a':>8}  "
        L_dbtable[pan] = row_l
        print(row)
    print()

    print("R-channel harmonic levels (dB) rel to center:")
    print(header)
    for pan in PAN_VALUES:
        L, R = load_window(pan)
        hR = dict((f, e) for f, e in harmonic_magnitudes(R))
        row = f"  {pan:>4}  "
        row_r = {}
        for f in harms_to_show:
            if f in hR and hR[f] > 1e-20 and Rc_harm.get(f, 0) > 1e-20:
                d = db(hR[f]) - db(Rc_harm[f])
                row += f"{d:>+7.2f}  "
                row_r[f] = d
            else:
                row += f"{'n/a':>8}  "
        R_dbtable[pan] = row_r
        print(row)
    print()

    # dump
    out = {
        "itd": [
            {
                "pan": r[0],
                "xc_lag_samples": r[1],
                "phase_lag_220_samples": r[2],
                "phase_lag_440_samples": r[3],
                "phase_lag_660_samples": r[4],
                "phase_deg_220": r[5],
                "phase_deg_440": r[6],
                "phase_deg_660": r[7],
            }
            for r in itd_table
        ],
        "per_harm_R_minus_L_db": {str(k): v for k, v in per_harm_table.items()},
        "L_db_rel_center": {str(k): v for k, v in L_dbtable.items()},
        "R_db_rel_center": {str(k): v for k, v in R_dbtable.items()},
    }
    Path(__file__).with_suffix(".json").write_text(json.dumps(out, indent=2, default=str))


if __name__ == "__main__":
    main()
