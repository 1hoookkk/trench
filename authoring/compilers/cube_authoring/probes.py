"""Probe suite — the 7 sources CUBE_GATE.md §6 runs every cube through.

Five probes are synthesized in code (sub, sine_sweep, saw_bass, reese, noise).
Two require real audio samples (vocal, drum_loop); those load from
`cartridges/factory/_source/probes/` when present and fail loudly otherwise —
no silent substitution.

All probes return mono float32 arrays at the runtime sample rate
(CUBE_GATE.md §0 = 44100 Hz). Each is one second long, level-normalized
so RMS equals `PROBE_RMS`.
"""
from __future__ import annotations

import wave
from pathlib import Path

import numpy as np

RUNTIME_SR = 44100
PROBE_SECONDS = 1.0
PROBE_RMS = 0.1          # -20 dB nominal — matches CUBE_GATE.md blow-up math
FADE_MS = 5.0            # short anti-click fade on deterministic probes

PROBES_DIR_DEFAULT = Path(__file__).resolve().parents[2] / "cartridges" / "engine" \
    / "_source" / "probes"

ProbeID = str
PROBE_IDS: tuple[ProbeID, ...] = (
    "sub",
    "sine_sweep",
    "saw_bass",
    "reese",
    "vocal",
    "drum_loop",
    "noise",
)


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

def _apply_fade(x: np.ndarray, sr: int, ms: float = FADE_MS) -> np.ndarray:
    n = int(round(sr * ms / 1000.0))
    if n <= 0 or 2 * n >= x.size:
        return x
    ramp = np.linspace(0.0, 1.0, n, dtype=np.float32)
    x = x.copy()
    x[:n] *= ramp
    x[-n:] *= ramp[::-1]
    return x


def _normalize_rms(x: np.ndarray, target_rms: float = PROBE_RMS) -> np.ndarray:
    rms = float(np.sqrt(np.mean(x.astype(np.float64) ** 2)))
    if rms <= 1e-12:
        return x.astype(np.float32)
    return (x * (target_rms / rms)).astype(np.float32)


def _t(sr: int, seconds: float = PROBE_SECONDS) -> np.ndarray:
    n = int(round(sr * seconds))
    return np.arange(n, dtype=np.float64) / sr


# -----------------------------------------------------------------------------
# Synthesized probes
# -----------------------------------------------------------------------------

def probe_sub(sr: int = RUNTIME_SR) -> np.ndarray:
    """60 Hz sine. DC/LF probe."""
    t = _t(sr)
    x = np.sin(2.0 * np.pi * 60.0 * t).astype(np.float32)
    return _normalize_rms(_apply_fade(x, sr))


def probe_sine_sweep(sr: int = RUNTIME_SR) -> np.ndarray:
    """Exponential (log) sweep 20 Hz -> 20 kHz. Broadband pole tracing."""
    t = _t(sr)
    f0, f1 = 20.0, 20000.0
    dur = float(t[-1] - t[0]) if t.size >= 2 else PROBE_SECONDS
    k = (f1 / f0) ** (1.0 / dur)
    # Phase is integral of instantaneous frequency:
    # phi(t) = 2 pi f0 * (k^t - 1) / ln(k)
    phase = 2.0 * np.pi * f0 * (np.power(k, t) - 1.0) / np.log(k)
    x = np.sin(phase).astype(np.float32)
    return _normalize_rms(_apply_fade(x, sr))


def probe_saw_bass(sr: int = RUNTIME_SR) -> np.ndarray:
    """55 Hz bandlimited-ish sawtooth. Harmonic-rich LF."""
    t = _t(sr)
    f0 = 55.0
    nyq = sr / 2.0
    n_harm = int(nyq // f0)
    x = np.zeros_like(t)
    for k in range(1, n_harm + 1):
        x += (np.sin(2.0 * np.pi * f0 * k * t)) / k
    x *= (2.0 / np.pi)
    x = x.astype(np.float32)
    return _normalize_rms(_apply_fade(x, sr))


def probe_reese(sr: int = RUNTIME_SR) -> np.ndarray:
    """Two detuned saws near C1 (~32.7 Hz), 7 Hz apart. Midrange stress."""
    t = _t(sr)
    f_lo, f_hi = 32.7, 39.7
    nyq = sr / 2.0
    x = np.zeros_like(t)
    for f0 in (f_lo, f_hi):
        n_harm = int(nyq // f0)
        for k in range(1, n_harm + 1):
            x += (np.sin(2.0 * np.pi * f0 * k * t)) / k
    x *= (1.0 / np.pi)
    x = x.astype(np.float32)
    return _normalize_rms(_apply_fade(x, sr))


def probe_noise(sr: int = RUNTIME_SR, seed: int = 0xC0B0) -> np.ndarray:
    """Gaussian white noise. Flat excitation."""
    rng = np.random.default_rng(seed)
    n = int(round(sr * PROBE_SECONDS))
    x = rng.standard_normal(n).astype(np.float32)
    return _normalize_rms(_apply_fade(x, sr))


# -----------------------------------------------------------------------------
# Sample-backed probes
# -----------------------------------------------------------------------------

class MissingProbeSample(FileNotFoundError):
    """Raised when a sample-backed probe is requested and the file is absent.

    The gate is expected to either curate the sample or declare the probe
    out-of-scope for the cube's declared scope (CUBE_GATE.md §7) — not to
    fall back to a substitute source.
    """


def _load_wav_mono(path: Path, sr: int) -> np.ndarray:
    with wave.open(str(path), "rb") as w:
        file_sr = w.getframerate()
        nch = w.getnchannels()
        width = w.getsampwidth()
        raw = w.readframes(w.getnframes())
    if width == 2:
        data = np.frombuffer(raw, dtype=np.int16).astype(np.float32) / 32768.0
    elif width == 3:
        # 24-bit PCM
        b = np.frombuffer(raw, dtype=np.uint8).reshape(-1, 3)
        vals = b[:, 0].astype(np.int32) | (b[:, 1].astype(np.int32) << 8) \
             | (b[:, 2].astype(np.int32) << 16)
        vals = np.where(vals & 0x800000, vals - 0x1000000, vals)
        data = vals.astype(np.float32) / (2 ** 23)
    elif width == 4:
        data = np.frombuffer(raw, dtype=np.int32).astype(np.float32) / (2 ** 31)
    else:
        raise ValueError(f"unsupported wav sampwidth {width}")
    if nch > 1:
        data = data.reshape(-1, nch).mean(axis=1)
    if file_sr != sr:
        raise ValueError(
            f"probe sample {path} is {file_sr} Hz; gate requires {sr} Hz. "
            f"Resample upstream (clean-room) before placing under probes/."
        )
    target_n = int(round(sr * PROBE_SECONDS))
    if data.size >= target_n:
        data = data[:target_n]
    else:
        data = np.pad(data, (0, target_n - data.size))
    return _normalize_rms(data.astype(np.float32))


def probe_vocal(sr: int = RUNTIME_SR,
                 probes_dir: Path | None = None) -> np.ndarray:
    """Dry vocal 'ah', ~1 s, mono, at the runtime sample rate.

    Loaded from probes/vocal.wav. If absent, raises MissingProbeSample — the
    caller must either curate the file or declare 'vocal' out-of-scope for
    the cube per CUBE_GATE.md §7. No synthesis substitute is provided; a
    formant-synthesis vocal would itself be a Z-plane output, which is what
    the gate is *testing* — using it as an input is circular.
    """
    d = probes_dir or PROBES_DIR_DEFAULT
    path = d / "vocal.wav"
    if not path.exists():
        raise MissingProbeSample(
            f"vocal probe missing: {path}. Either place a dry 1 s 44.1 kHz "
            f"mono vocal WAV there, or declare 'vocal' out-of-scope for the "
            f"cube (CUBE_GATE.md §7)."
        )
    return _load_wav_mono(path, sr)


def probe_drum_loop(sr: int = RUNTIME_SR,
                     probes_dir: Path | None = None) -> np.ndarray:
    """Normalized breakbeat, ~1 s, mono, at the runtime sample rate.

    Loaded from probes/drum_loop.wav. Raises MissingProbeSample if absent.
    Not synthesizable cleanly — transient density / stereo mix / swing are
    what this probe tests, and a stub would hide exactly what we need to
    measure.
    """
    d = probes_dir or PROBES_DIR_DEFAULT
    path = d / "drum_loop.wav"
    if not path.exists():
        raise MissingProbeSample(
            f"drum_loop probe missing: {path}. Place a 1 s 44.1 kHz mono "
            f"breakbeat there, or declare 'drum_loop' out-of-scope (§7)."
        )
    return _load_wav_mono(path, sr)


# -----------------------------------------------------------------------------
# Registry + scope filtering
# -----------------------------------------------------------------------------

_PROBE_FACTORIES = {
    "sub": probe_sub,
    "sine_sweep": probe_sine_sweep,
    "saw_bass": probe_saw_bass,
    "reese": probe_reese,
    "vocal": probe_vocal,
    "drum_loop": probe_drum_loop,
    "noise": probe_noise,
}


def get_probe(probe_id: ProbeID, sr: int = RUNTIME_SR) -> np.ndarray:
    if probe_id not in _PROBE_FACTORIES:
        raise ValueError(f"unknown probe_id {probe_id!r}")
    return _PROBE_FACTORIES[probe_id](sr)


def available_probes(sr: int = RUNTIME_SR,
                      probes_dir: Path | None = None
                      ) -> dict[ProbeID, bool]:
    """Return which probes are currently materializable without scope-gating.

    Synthesized probes are always True. Sample-backed probes are True iff
    their WAV exists at the expected runtime SR.
    """
    d = probes_dir or PROBES_DIR_DEFAULT
    out = {p: True for p in ("sub", "sine_sweep", "saw_bass", "reese", "noise")}
    for p, fname in (("vocal", "vocal.wav"), ("drum_loop", "drum_loop.wav")):
        out[p] = (d / fname).exists()
    return out


def in_scope_probes(scope: str) -> tuple[ProbeID, ...]:
    """Which probes are in-scope for a declared cube scope.

    hf_resonance           : broadband HF behavior; skips LF-only probes.
    broadband              : all 7.
    formant                : vocal/spectral focus; skips pure-LF.
    """
    if scope == "hf_resonance":
        return ("sine_sweep", "reese", "drum_loop", "noise")
    if scope == "formant":
        return ("sine_sweep", "saw_bass", "reese", "vocal", "noise")
    if scope == "broadband":
        return PROBE_IDS
    # Unknown scope — be conservative, return all. Gate will report per-probe.
    return PROBE_IDS
