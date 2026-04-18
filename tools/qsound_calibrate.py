"""
qsound_calibrate.py
===================

Calibrates the parametric QSound spatial model against measured WAV output
from QMixer.dll pan-batch captures. Resolves UNKNOWNs in
docs/archive/qsound_spatial.md:

  * ITD output unit (coef output -> samples @ 48 kHz conversion factor)
  * az input unit (radians vs degrees when feeding sin(az))
  * ITD ear-mapping convention (which channel is delayed at +az)
  * ILD dB output unit
  * 3-Band topology: empirical band-gain fit + best crossover pair

Produces:

  * Per-pan ITD / ILD / per-band measurements
  * Comparison against the parametric model in qsound_spatial.md
  * Markdown addendum written by a separate step (writer, not this script)

Usage:
    python tools/qsound_calibrate.py

Data source:
    C:/Users/hooki/trench_re_vault/datasets/qsound_spatial_v1/
        2026-03-05_pan_batch_fullgrid_a/
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
from scipy.io import wavfile
from scipy.signal import correlate


# --------------------------------------------------------------------------- #
# Config
# --------------------------------------------------------------------------- #

DATA_DIR = Path(
    r"C:/Users/hooki/trench_re_vault/datasets/qsound_spatial_v1"
    r"/2026-03-05_pan_batch_fullgrid_a"
)
EXPECTED_SR = 48000
STIMULUS_FREQ = 220.0  # saw fundamental, Hz
PAN_VALUES = list(range(-60, 70, 10))  # -60..+60 step 10

# steady-state window (skip 0.25s transient, take 1.0s)
WIN_START_SEC = 0.25
WIN_DUR_SEC = 1.0

# crossover candidates for band-law fit
CROSSOVER_CANDIDATES = [
    (200.0, 2000.0),
    (500.0, 5000.0),
    (300.0, 3000.0),
    (500.0, 2000.0),
    (1000.0, 4000.0),
]


# --------------------------------------------------------------------------- #
# Parametric model (from docs/archive/qsound_spatial.md)
# --------------------------------------------------------------------------- #

ITD_COEFS = [
    3578.7646232504208,
    -99.11300089026597,
    -960.7096166779854,
    631.0360877559084,
    -229.78233226297414,
    -173.11867492761178,
]

ILD_COEFS = [
    6.819731516349544,
    -2.5008130981129657,
    0.8210199212077876,
    -0.198121089137829,
    0.021401263894609834,
    8.819077411548193e-05,
]

# band_law features: 1, log2(dist/0.25), log2(dist/0.25)^2,
# sin(az), cos(az), sin(2*az), cos(2*az),
# sin(3*az), cos(3*az), el/30, sin(az)*el/30, cos(az)*el/30
BAND_LAW = {
    "l": {
        "low":  [-72.211726, 0.386310, -1.609965, -2.942365, 4.905630,
                 0.894425, -1.686776, -0.082666, 0.833466,
                 -5.542e-17, 2.501e-16, 1.788e-16],
        "mid":  [-74.892487, 0.460364, -1.645200, -2.794224, 9.547364,
                 1.100615, -2.989885, -0.168168, 1.485241,
                 1.372e-17, 2.241e-17, 3.424e-16],
        "high": [-73.911507, 0.438043, -1.634712, -2.835356, 8.046734,
                 1.041892, -2.546882, -0.132837, 1.302472,
                 -3.615e-17, 3.924e-17, -1.875e-16],
    },
    "r": {
        "low":  [-72.211726, 0.386310, -1.609965, 2.942365, 4.905630,
                 -0.894425, -1.686776, 0.082666, 0.833466,
                 6.880e-17, 9.784e-17, -3.152e-16],
        "mid":  [-74.892487, 0.460364, -1.645200, 2.794224, 9.547364,
                 -1.100615, -2.989885, 0.168168, 1.485241,
                 -3.922e-16, 1.316e-16, 9.925e-17],
        "high": [-73.911507, 0.438043, -1.634712, 2.835356, 8.046734,
                 -1.041892, -2.546882, 0.132837, 1.302472,
                 -1.004e-16, 1.863e-16, -2.686e-16],
    },
}


def itd_parametric(az_deg: float, el_deg: float = 0.0, az_as: str = "radians") -> float:
    """Raw parametric ITD output (unit unknown). `az_as` picks az unit."""
    if az_as == "radians":
        az = math.radians(az_deg)
    elif az_as == "degrees":
        az = az_deg
    else:
        raise ValueError(az_as)
    feats = [
        math.sin(az),
        math.sin(2 * az),
        math.sin(3 * az),
        math.sin(4 * az),
        math.sin(5 * az),
        math.sin(az) * abs(el_deg) / 30.0,
    ]
    return sum(c * f for c, f in zip(ITD_COEFS, feats))


def ild_parametric(az_deg: float, el_deg: float = 0.0, az_as: str = "radians") -> float:
    if az_as == "radians":
        az = math.radians(az_deg)
    elif az_as == "degrees":
        az = az_deg
    else:
        raise ValueError(az_as)
    feats = [
        math.sin(az),
        math.sin(2 * az),
        math.sin(3 * az),
        math.sin(4 * az),
        math.sin(5 * az),
        math.sin(az) * abs(el_deg) / 30.0,
    ]
    return sum(c * f for c, f in zip(ILD_COEFS, feats))


def band_parametric(
    channel: str,
    band: str,
    az_deg: float,
    dist: float = 1.0,
    el_deg: float = 0.0,
) -> float:
    """Parametric per-band gain in dB."""
    az = math.radians(az_deg)
    ldist = math.log2(dist / 0.25)
    feats = [
        1.0,
        ldist,
        ldist ** 2,
        math.sin(az),
        math.cos(az),
        math.sin(2 * az),
        math.cos(2 * az),
        math.sin(3 * az),
        math.cos(3 * az),
        el_deg / 30.0,
        math.sin(az) * el_deg / 30.0,
        math.cos(az) * el_deg / 30.0,
    ]
    coefs = BAND_LAW[channel][band]
    return sum(c * f for c, f in zip(coefs, feats))


# --------------------------------------------------------------------------- #
# Measurement
# --------------------------------------------------------------------------- #

def load_wav(path: Path):
    sr, data = wavfile.read(str(path))
    if data.dtype != np.float32:
        # promote to float32 and normalize if integer encoding
        data = data.astype(np.float32)
        if np.issubdtype(wavfile.read(str(path))[1].dtype, np.integer):
            data /= np.iinfo(wavfile.read(str(path))[1].dtype).max
    if data.ndim != 2 or data.shape[1] != 2:
        raise RuntimeError(f"expected stereo, got shape {data.shape} from {path}")
    return sr, data


def steady_window(x: np.ndarray, sr: int) -> np.ndarray:
    start = int(WIN_START_SEC * sr)
    n = int(WIN_DUR_SEC * sr)
    end = min(start + n, x.shape[0])
    return x[start:end]


def measure_itd_samples(L: np.ndarray, R: np.ndarray, max_lag: int = 256) -> int:
    """Cross-correlate L vs R. Returns lag where R = L shifted by lag samples.

    Sign convention: positive lag means R *leads* L (R appears earlier).
    """
    # normalize to reduce bias from amplitude differences
    L = L - L.mean()
    R = R - R.mean()
    Ln = L / (np.linalg.norm(L) + 1e-20)
    Rn = R / (np.linalg.norm(R) + 1e-20)
    xc = correlate(Rn, Ln, mode="full")
    center = len(L) - 1
    lo = center - max_lag
    hi = center + max_lag + 1
    xc_window = xc[lo:hi]
    peak = int(np.argmax(np.abs(xc_window))) - max_lag
    return peak


def measure_ild_db(L: np.ndarray, R: np.ndarray) -> float:
    lrms = float(np.sqrt(np.mean(L * L)))
    rrms = float(np.sqrt(np.mean(R * R)))
    if lrms < 1e-12 or rrms < 1e-12:
        return float("nan")
    return 20.0 * math.log10(rrms / lrms)


def band_energies(x: np.ndarray, sr: int, crossovers: tuple[float, float]) -> tuple[float, float, float]:
    """Sum harmonic bin magnitudes^2 within each band. Returns (low, mid, high) linear energy."""
    n = len(x)
    X = np.fft.rfft(x * np.hanning(n))
    freqs = np.fft.rfftfreq(n, 1.0 / sr)
    mag2 = (np.abs(X) ** 2)

    # find harmonic bins of 220 Hz up to Nyquist
    harmonics = []
    k = 1
    while STIMULUS_FREQ * k < sr / 2 - 100:
        f = STIMULUS_FREQ * k
        # find bin closest to harmonic
        bin_idx = int(round(f * n / sr))
        # sum a small neighborhood (main lobe of hann ~ +-2 bins)
        lo = max(0, bin_idx - 3)
        hi = min(len(mag2), bin_idx + 4)
        harmonics.append((f, mag2[lo:hi].sum()))
        k += 1

    xo_low, xo_hi = crossovers
    low = mid = high = 0.0
    for f, e in harmonics:
        if f < xo_low:
            low += e
        elif f < xo_hi:
            mid += e
        else:
            high += e
    return low, mid, high


def db(x: float) -> float:
    return 10.0 * math.log10(x + 1e-30)


# --------------------------------------------------------------------------- #
# Driver
# --------------------------------------------------------------------------- #

def fmt_float(v, w=8, d=3):
    return f"{v:{w}.{d}f}"


def main() -> None:
    print(f"\n=== QSound calibration ===\nData: {DATA_DIR}\n")

    # ---- Load and window every pan -------------------------------------- #
    panwavs = {}
    for pan in PAN_VALUES:
        sign = "+" if pan >= 0 else "-"
        fname = f"pan_{sign}{abs(pan):03d}.wav"
        path = DATA_DIR / fname
        sr, data = load_wav(path)
        if sr != EXPECTED_SR:
            raise RuntimeError(f"unexpected sr {sr} in {path}")
        w = steady_window(data, sr)
        panwavs[pan] = (sr, w)

    # center reference ----------------------------------------------------- #
    sr_c, wc = panwavs[0]
    Lc = wc[:, 0].astype(np.float64)
    Rc = wc[:, 1].astype(np.float64)

    # reference band energies per channel at center (dB)
    refs_by_xo = {}
    for xo in CROSSOVER_CANDIDATES:
        l_low, l_mid, l_high = band_energies(Lc, sr_c, xo)
        r_low, r_mid, r_high = band_energies(Rc, sr_c, xo)
        refs_by_xo[xo] = {
            "l": (db(l_low), db(l_mid), db(l_high)),
            "r": (db(r_low), db(r_mid), db(r_high)),
        }

    # ---- Per-pan measurements ------------------------------------------- #
    itd_measured = {}   # pan -> lag samples (R leads L when positive)
    ild_measured = {}   # pan -> R_dB - L_dB
    band_measured = {xo: {} for xo in CROSSOVER_CANDIDATES}

    for pan, (sr, w) in panwavs.items():
        L = w[:, 0].astype(np.float64)
        R = w[:, 1].astype(np.float64)
        lag = measure_itd_samples(L, R)
        itd_measured[pan] = lag
        ild_measured[pan] = measure_ild_db(L, R)
        for xo in CROSSOVER_CANDIDATES:
            l_low, l_mid, l_high = band_energies(L, sr, xo)
            r_low, r_mid, r_high = band_energies(R, sr, xo)
            band_measured[xo][pan] = {
                "l": (db(l_low), db(l_mid), db(l_high)),
                "r": (db(r_low), db(r_mid), db(r_high)),
            }

    # ---- Report: ITD ---------------------------------------------------- #
    print("ITD measurement (R leads L when lag > 0):")
    print(
        f"{'pan':>5}  {'lag_samp':>8}  {'lag_us':>8}  "
        f"{'param_rad':>10}  {'param_deg':>10}  {'ratio_rad':>9}"
    )

    rad_ratios = []
    deg_ratios = []
    itd_rows = []
    for pan in PAN_VALUES:
        lag = itd_measured[pan]
        us = lag * 1e6 / EXPECTED_SR
        p_rad = itd_parametric(pan, 0.0, "radians")
        p_deg = itd_parametric(pan, 0.0, "degrees")
        ratio_rad = (p_rad / lag) if lag != 0 else float("inf")
        ratio_deg = (p_deg / lag) if lag != 0 else float("inf")
        if pan != 0 and abs(lag) > 0:
            rad_ratios.append(p_rad / lag)
            deg_ratios.append(p_deg / lag)
        itd_rows.append((pan, lag, us, p_rad, p_deg, ratio_rad))
        print(
            f"{pan:>5}  {lag:>8d}  {us:>8.1f}  "
            f"{p_rad:>10.2f}  {p_deg:>10.2f}  {ratio_rad:>9.2f}"
        )

    rad_ratio_mean = float(np.mean(rad_ratios)) if rad_ratios else float("nan")
    rad_ratio_std = float(np.std(rad_ratios)) if rad_ratios else float("nan")
    deg_ratio_mean = float(np.mean(deg_ratios)) if deg_ratios else float("nan")
    deg_ratio_std = float(np.std(deg_ratios)) if deg_ratios else float("nan")

    print(f"\nratio (radians interp) mean = {rad_ratio_mean:.3f}  std = {rad_ratio_std:.3f}")
    print(f"ratio (degrees interp) mean = {deg_ratio_mean:.3f}  std = {deg_ratio_std:.3f}")

    # ear mapping
    # positive lag (R leads L) at pan +60 means R channel is EARLY.
    # If pan +60 = source on right, L should be DELAYED (i.e. L comes later -> R leads).
    # So lag>0 at pan +60 -> L is delayed ear.
    lag_pos60 = itd_measured[60]
    lag_neg60 = itd_measured[-60]
    print(f"\near mapping probe: lag(+60) = {lag_pos60}  lag(-60) = {lag_neg60}")

    # ---- Report: ILD ---------------------------------------------------- #
    print("\nILD measurement (20log10(R_rms/L_rms)):")
    print(f"{'pan':>5}  {'meas_dB':>9}  {'param_rad':>10}  {'param_deg':>10}  {'res_rad':>9}")
    ild_rows = []
    rad_res = []
    for pan in PAN_VALUES:
        m = ild_measured[pan]
        p_rad = ild_parametric(pan, 0.0, "radians")
        p_deg = ild_parametric(pan, 0.0, "degrees")
        ild_rows.append((pan, m, p_rad, p_deg, m - p_rad))
        if pan != 0:
            rad_res.append(abs(m - p_rad))
        print(f"{pan:>5}  {m:>9.3f}  {p_rad:>10.3f}  {p_deg:>10.3f}  {m - p_rad:>9.3f}")
    ild_rad_mae = float(np.mean(rad_res)) if rad_res else float("nan")
    print(f"\nILD MAE (radians interp, pan!=0) = {ild_rad_mae:.3f} dB")

    # ---- Report: band-law ---------------------------------------------- #
    print("\n3-Band gain measurements (per pan, vs center reference):")

    best_xo = None
    best_err = float("inf")
    xo_results = {}
    for xo in CROSSOVER_CANDIDATES:
        residuals = []
        print(f"\n  crossovers = {xo} Hz")
        print(
            f"  {'pan':>5}  "
            f"{'L_low_m':>8} {'L_low_p':>8}  "
            f"{'L_mid_m':>8} {'L_mid_p':>8}  "
            f"{'L_hi_m':>7} {'L_hi_p':>7}   |   "
            f"{'R_low_m':>8} {'R_low_p':>8}  "
            f"{'R_mid_m':>8} {'R_mid_p':>8}  "
            f"{'R_hi_m':>7} {'R_hi_p':>7}"
        )
        per_pan_xo = {}
        for pan in PAN_VALUES:
            ref_l = refs_by_xo[xo]["l"]
            ref_r = refs_by_xo[xo]["r"]
            meas = band_measured[xo][pan]
            l_deltas = tuple(meas["l"][i] - ref_l[i] for i in range(3))
            r_deltas = tuple(meas["r"][i] - ref_r[i] for i in range(3))

            # parametric - also take delta from pan=0
            pl = tuple(band_parametric("l", b, pan) for b in ("low", "mid", "high"))
            pr = tuple(band_parametric("r", b, pan) for b in ("low", "mid", "high"))
            pl0 = tuple(band_parametric("l", b, 0.0) for b in ("low", "mid", "high"))
            pr0 = tuple(band_parametric("r", b, 0.0) for b in ("low", "mid", "high"))
            pl_delta = tuple(pl[i] - pl0[i] for i in range(3))
            pr_delta = tuple(pr[i] - pr0[i] for i in range(3))

            per_pan_xo[pan] = {
                "l_meas": l_deltas, "l_param": pl_delta,
                "r_meas": r_deltas, "r_param": pr_delta,
            }

            if pan != 0:
                for a, b in zip(l_deltas, pl_delta):
                    residuals.append(abs(a - b))
                for a, b in zip(r_deltas, pr_delta):
                    residuals.append(abs(a - b))

            print(
                f"  {pan:>5}  "
                f"{l_deltas[0]:>8.2f} {pl_delta[0]:>8.2f}  "
                f"{l_deltas[1]:>8.2f} {pl_delta[1]:>8.2f}  "
                f"{l_deltas[2]:>7.2f} {pl_delta[2]:>7.2f}   |   "
                f"{r_deltas[0]:>8.2f} {pr_delta[0]:>8.2f}  "
                f"{r_deltas[1]:>8.2f} {pr_delta[1]:>8.2f}  "
                f"{r_deltas[2]:>7.2f} {pr_delta[2]:>7.2f}"
            )
        mae = float(np.mean(residuals)) if residuals else float("inf")
        xo_results[xo] = {"mae": mae, "per_pan": per_pan_xo}
        print(f"  MAE vs parametric = {mae:.3f} dB")
        if mae < best_err:
            best_err = mae
            best_xo = xo

    print(f"\nBest crossover: {best_xo} with MAE {best_err:.3f} dB")

    # ---- Dump JSON summary for addendum consumption ---------------------- #
    out = {
        "data_dir": str(DATA_DIR),
        "sr": EXPECTED_SR,
        "stimulus_hz": STIMULUS_FREQ,
        "window_sec": [WIN_START_SEC, WIN_START_SEC + WIN_DUR_SEC],
        "itd": {
            "per_pan": [
                {
                    "pan": pan,
                    "lag_samples": itd_measured[pan],
                    "lag_us": itd_measured[pan] * 1e6 / EXPECTED_SR,
                    "parametric_rad": itd_parametric(pan, 0.0, "radians"),
                    "parametric_deg": itd_parametric(pan, 0.0, "degrees"),
                }
                for pan in PAN_VALUES
            ],
            "ratio_radians_mean": rad_ratio_mean,
            "ratio_radians_std": rad_ratio_std,
            "ratio_degrees_mean": deg_ratio_mean,
            "ratio_degrees_std": deg_ratio_std,
            "lag_pos60": lag_pos60,
            "lag_neg60": lag_neg60,
        },
        "ild": {
            "per_pan": [
                {
                    "pan": pan,
                    "measured_dB": ild_measured[pan],
                    "parametric_rad_dB": ild_parametric(pan, 0.0, "radians"),
                    "parametric_deg_dB": ild_parametric(pan, 0.0, "degrees"),
                }
                for pan in PAN_VALUES
            ],
            "mae_radians_dB": ild_rad_mae,
        },
        "bands": {
            "best_crossovers_hz": list(best_xo),
            "best_mae_dB": best_err,
            "all": {
                f"{xo[0]:.0f}_{xo[1]:.0f}": {
                    "mae_dB": xo_results[xo]["mae"],
                }
                for xo in CROSSOVER_CANDIDATES
            },
        },
    }
    out_path = Path(__file__).with_suffix(".json")
    out_path.write_text(json.dumps(out, indent=2))
    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
