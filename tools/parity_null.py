"""Parity null: TRENCH DSP pipeline vs canonical references.

Walks ref/canonical/MANIFEST.json. Each entry carries a per-skin source
path (relative to the Trench root) and a source_type:

  - raw_p2k_skin  : datasets/p2k_skins/*.json (33 numbered + 2 vocal)
                    stage form {a1, r, val1, val2, val3, flag}
  - calibration_re: docs/calibration/*.json (6 ft=33-55 extractions)
                    stage form {pole_freq_hz, radius, val1, val2, val3, ...}
                    a1 is reconstructed: a1 = -2r*cos(2pi*f/sample_rate_authored)

Pipeline:
    source stage -> (b0, b1, b2, a1, a2) -> scipy.signal.sosfilt
                 -> [per-sample AGC (16-entry table from agc.rs) — bypassed with --no-agc]
                 -> * boost
                 -> best-fit gain null vs reference WAV (lag=0)

Hard fail at rel_null > FAIL_THRESHOLD_DB.

Flags:
    --no-agc                Skip the per-sample AGC step on both rebake
                            and null paths. Canonical Talking Hedz WAVs
                            are baked AGC-off because the bit-masked AGC
                            state machine (`(gain * abs_sample) as u32 & 0xF`)
                            diverges between f32 (Rust) and f64 (Python);
                            AGC correctness is covered by
                            `trench-core/src/agc.rs` unit tests.
    --rebake NAME           Regenerate manifest entry NAME's 4 corner WAVs.
                            Combine with --no-agc for AGC-off rebakes.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
import soundfile as sf
from scipy.signal import sosfilt

# Port of trench-core/src/resampler.rs — a 32-tap, 512-phase Blackman-
# windowed sinc polyphase FIR with per-phase DC normalisation. Canonical
# renders use this instead of scipy.signal.resample_poly so both sides
# of the parity null use the SAME FIR at the host<->native boundary.
# See tools/rust_resampler.py for the port; xlang_null against a Rust
# dump measures -100 dB (bit-close).
from rust_resampler import resample as rust_resample

ROOT = Path(__file__).resolve().parent.parent
DRY = Path(r"C:/Users/hooki/trenchwork_clean/ref/bypassed-pinknoise.wav")
REF_DIR = ROOT / "ref" / "canonical"

CORNERS = ("M0_Q0", "M0_Q100", "M100_Q0", "M100_Q100")
# Null worse than this triggers a failure exit code.
FAIL_THRESHOLD_DB = -140.0

# Cascade coefficients are authored at the E-mu native sample rate. The
# runtime FilterEngine resamples host -> native -> cascade -> host; parity
# has to mirror that or every pole shifts in analog frequency.
NATIVE_SR = 39062.5

# AGC_TABLE ported from trenchwork_clean/trench-core/src/agc.rs
AGC_TABLE = np.array(
    [1.0001, 1.0001, 0.996, 0.990, 0.920, 0.500, 0.200, 0.160,
     0.120, 0.120, 0.120, 0.120, 0.120, 0.120, 0.120, 0.120],
    dtype=np.float64,
)


def db(x: float) -> float:
    x = float(x)
    return -300.0 if x <= 0.0 else 20.0 * np.log10(x)


def mono(d: np.ndarray) -> np.ndarray:
    return d[:, 0] if d.ndim > 1 else d


def stage_coefficients_raw(stage: dict) -> tuple[float, float, float, float, float]:
    """Raw p2k_skin stage → (b0, b1, b2, a1, a2)."""
    a1 = float(stage["a1"])
    r = min(float(stage["r"]), 0.999999)
    a2 = r * r
    if float(stage.get("flag", 1.0)) < 0.5:
        b0, b1, b2 = 1.0, 0.0, 0.0
    else:
        b0 = 1.0 + float(stage["val1"])
        b1 = a1 + float(stage["val2"])
        b2 = a2 - float(stage["val3"])
    return b0, b1, b2, a1, a2


def stage_coefficients_compiled_v1(
    stage: dict,
) -> tuple[float, float, float, float, float]:
    """Compiled-v1 stage → (b0, b1, b2, a1, a2). Direct mapping, no
    recomputation: the compiled-v1 kernel already holds DF2T coefficients
    with `c0=b0, c1=b1, c2=b2, c3=a1, c4=a2` (see
    `trench-juce/forge/pyruntime/render.py:47`). Use this when the parity
    test wants bit-close coefficient identity with the Rust runtime,
    instead of the calibration_re path that re-derives `a1` from
    `pole_freq_hz` and accumulates compile-path drift."""
    return (
        float(stage["c0"]),
        float(stage["c1"]),
        float(stage["c2"]),
        float(stage["c3"]),
        float(stage["c4"]),
    )


def stage_coefficients_calibration(
    stage: dict, sample_rate_authored: float
) -> tuple[float, float, float, float, float]:
    """Calibration stage → (b0, b1, b2, a1, a2). Reconstructs a1 from
    pole_freq_hz + radius at the calibration's authoring sample rate."""
    r = min(float(stage["radius"]), 0.999999)
    f = float(stage["pole_freq_hz"])
    a1 = -2.0 * r * math.cos(2.0 * math.pi * f / float(sample_rate_authored))
    a2 = r * r
    val1 = float(stage.get("val1", 0.0))
    val2 = float(stage.get("val2", 0.0))
    val3 = float(stage.get("val3", 0.0))
    if val1 == 0.0 and val2 == 0.0 and val3 == 0.0:
        b0, b1, b2 = 1.0, 0.0, 0.0
    else:
        b0 = 1.0 + val1
        b1 = a1 + val2
        b2 = a2 - val3
    return b0, b1, b2, a1, a2


def build_sos(stages: list[dict], source_type: str, sample_rate_authored: float) -> np.ndarray:
    rows = []
    for st in stages:
        if source_type == "calibration_re":
            b0, b1, b2, a1, a2 = stage_coefficients_calibration(st, sample_rate_authored)
        elif source_type == "compiled_v1":
            b0, b1, b2, a1, a2 = stage_coefficients_compiled_v1(st)
        else:
            b0, b1, b2, a1, a2 = stage_coefficients_raw(st)
        rows.append([b0, b1, b2, 1.0, a1, a2])
    return np.asarray(rows, dtype=np.float64)


def get_corner_stages(body: dict, source_type: str, label: str) -> list[dict] | None:
    """Resolve a corner's stage list across source formats.

    - raw_p2k_skin / calibration_re: `body['corners'][label]['stages']`
    - compiled_v1: `body['keyframes'][i]['stages']` where keyframes[i].label == label
    """
    if source_type == "compiled_v1":
        for kf in body.get("keyframes", []):
            if kf.get("label") == label:
                return kf.get("stages")
        return None
    corner = body.get("corners", {}).get(label)
    return corner.get("stages") if corner else None


def get_body_boost(body: dict, source_type: str) -> float:
    """Top-level boost across source formats. Compiled-v1 carries per-keyframe
    boost but heritage cartridges use one scalar for all corners — take the
    first keyframe's boost as the canonical scalar. Parity tests only run
    cartridges where every corner shares one boost (verified at bake time
    in tools/bake_hedz_from_p2k.py)."""
    if source_type == "compiled_v1":
        kfs = body.get("keyframes", [])
        if kfs:
            return float(kfs[0].get("boost", 1.0))
        return 1.0
    return float(body.get("boost", 1.0))


def apply_agc(samples: np.ndarray) -> np.ndarray:
    """Per-sample post-cascade soft limiter (EmulatorX.dll AGC)."""
    out = np.empty_like(samples)
    gain = 1.0
    for i, s in enumerate(samples):
        idx = int(gain * abs(s)) & 0xF
        new_gain = gain * AGC_TABLE[idx]
        gain = new_gain if new_gain < 1.0 else 1.0
        out[i] = s * gain
    return out


def apply_dc_block(samples: np.ndarray, sample_rate: float = 44100.0) -> np.ndarray:
    """One-pole DC blocker matching trench-core engine.rs DcBlocker."""
    r = 1.0 - (2.0 * np.pi * 20.0 / sample_rate)
    out = np.empty_like(samples)
    x_prev = 0.0
    y_prev = 0.0
    for i, x in enumerate(samples):
        y = x - x_prev + r * y_prev
        x_prev = x
        y_prev = y if np.isfinite(y) else 0.0
        out[i] = y_prev
    return out


def best_fit_null(pred: np.ndarray, ref: np.ndarray, lag_search: int = 0):
    """Lag-aligned best-fit-gain null. Canonical references are rendered
    through the exact same pipeline as `pred`, so lag=0 is correct and we
    skip the search by default."""
    n = min(len(pred), len(ref))
    p = pred[:n].astype(np.float64, copy=False)
    r = ref[:n].astype(np.float64, copy=False)
    den = float(np.dot(p, p))
    if den <= 0.0:
        return (0, 0.0, 0.0, 0.0, 0.0)
    gain = float(np.dot(p, r) / den)
    res = r - gain * p
    rms = float(np.sqrt(float(np.dot(res, res)) / n))
    peak = float(np.max(np.abs(res)))
    ref_rms = float(np.sqrt(float(np.dot(r, r)) / n))
    return (0, gain, rms, peak, ref_rms)


def load_source(entry: dict) -> tuple[dict, str, float]:
    """Resolve a manifest entry's source file. Returns (body, source_type,
    sample_rate_authored)."""
    src = ROOT / entry["source"]
    body = json.loads(src.read_text())
    src_type = entry.get("source_type", "raw_p2k_skin")
    sr_auth = float(entry.get("sample_rate_authored", 39062.5))
    return body, src_type, sr_auth


def render_canonical(
    dry: np.ndarray,
    stages: list[dict],
    src_type: str,
    sr_auth: float,
    boost: float,
    dry_sr: float,
    cascade_sr: float = NATIVE_SR,
    apply_agc_enabled: bool = True,
) -> np.ndarray:
    """Canonical render chain that mirrors trench-core's FilterEngine:

        host SR -> resample down to native SR -> DF2T cascade (SOS)
                -> [per-sample AGC] -> * boost -> resample up to host SR

    Cascade coefficients are authored at `sr_auth` (normally `cascade_sr`),
    so running sosfilt at `dry_sr` instead shifts every pole in analog
    frequency. Pulling the cascade back down to `cascade_sr` realigns it
    with what the runtime engine sees.

    The host<->native resampling uses `rust_resample` — the Python port of
    trench-core/src/resampler.rs (32-tap 512-phase Blackman-windowed sinc
    polyphase FIR with per-phase DC normalisation). scipy's
    `resample_poly` uses a Kaiser window with ~23 taps per subfilter, so
    its magnitude response near Nyquist differs from the runtime
    FilterEngine; using the ported FIR means both sides of the parity
    null share the exact same resampler.

    If `apply_agc_enabled` is False the AGC step is skipped. The AGC is a
    bit-masked nonlinear state machine (`(gain * abs_sample) as u32 & 0xF`
    indexing a 16-entry table with cliff-like values); run at f32 vs f64
    the index-flip trajectories diverge after thousands of samples, so
    cascade+resampler parity tests disable AGC on both sides and rely on
    `trench-core/src/agc.rs` unit tests for AGC correctness.

    Returns a float64 ndarray the same length as `dry` (truncated or
    zero-padded as needed, since polyphase resampler length is approximate).
    """
    if dry_sr == cascade_sr:
        dry_native = dry.astype(np.float64, copy=False)
    else:
        dry_native = rust_resample(
            dry.astype(np.float32, copy=False),
            input_sr=float(dry_sr),
            output_sr=float(cascade_sr),
        ).astype(np.float64, copy=False)
    sos = build_sos(stages, src_type, sr_auth)
    cascaded = sosfilt(sos, dry_native)
    if apply_agc_enabled:
        limited = apply_agc(cascaded)
    else:
        limited = cascaded
    scaled = limited * boost
    if dry_sr == cascade_sr:
        out = scaled
    else:
        out = rust_resample(
            scaled.astype(np.float32, copy=False),
            input_sr=float(cascade_sr),
            output_sr=float(dry_sr),
        ).astype(np.float64, copy=False)
    n = len(dry)
    if len(out) == n:
        return out
    if len(out) > n:
        return out[:n]
    padded = np.zeros(n, dtype=out.dtype)
    padded[: len(out)] = out
    return padded


def rebake(name: str, apply_agc_enabled: bool = True) -> int:
    """Regenerate a manifest entry's 4 corner WAVs through render_canonical.

    Writes to ref/canonical/{name}_{corner}.wav as 32-bit float mono @ 44100 Hz.
    """
    manifest_path = REF_DIR / "MANIFEST.json"
    if not manifest_path.exists():
        print(f"rebake: manifest missing at {manifest_path}")
        return 2
    manifest = json.loads(manifest_path.read_text())
    if name not in manifest:
        print(f"rebake: manifest has no entry named {name!r}")
        return 2
    if not DRY.exists():
        print(f"rebake: dry input missing at {DRY}")
        return 2

    entry = manifest[name]
    body, src_type, sr_auth = load_source(entry)
    boost = float(entry.get("boost", get_body_boost(body, src_type)))
    dry_data, dry_sr = sf.read(str(DRY))
    dry = mono(dry_data)
    dry_sr = float(dry_sr)

    agc_status = "on" if apply_agc_enabled else "OFF (bypassed)"
    print(f"rebake {name}: src_type={src_type} sr_auth={sr_auth} boost={boost} agc={agc_status}")
    print(f"  dry: {DRY.name} {len(dry)} samples @ {dry_sr:.1f} Hz")

    for label in CORNERS:
        stages = get_corner_stages(body, src_type, label)
        if stages is None:
            print(f"  {label}: missing in source body; skipped")
            continue
        out = render_canonical(
            dry,
            stages,
            src_type,
            sr_auth,
            boost,
            dry_sr,
            apply_agc_enabled=apply_agc_enabled,
        )
        out_path = REF_DIR / f"{name}_{label}.wav"
        sf.write(str(out_path), out.astype(np.float32), 44100, subtype="FLOAT")
        peak = float(np.max(np.abs(out))) if len(out) else 0.0
        print(f"  {label}: wrote {out_path.name}  len={len(out)}  peak={peak:.4f}")
    return 0


def main() -> int:
    # Flag parsing: --no-agc may appear anywhere after the subcommand.
    argv = sys.argv[1:]
    apply_agc_enabled = True
    if "--no-agc" in argv:
        apply_agc_enabled = False
        argv = [a for a in argv if a != "--no-agc"]

    # --rebake NAME: regenerate a manifest entry's 4 corner WAVs.
    if argv and argv[0] == "--rebake":
        if len(argv) != 2:
            print("usage: parity_null.py --rebake NAME [--no-agc]")
            return 2
        return rebake(argv[1], apply_agc_enabled=apply_agc_enabled)

    manifest_path = REF_DIR / "MANIFEST.json"
    if not REF_DIR.exists() or not manifest_path.exists() or not DRY.exists():
        print(
            f"parity inputs missing; skip "
            f"({REF_DIR.exists()=} {manifest_path.exists()=} {DRY.exists()=})"
        )
        return 0

    manifest = json.loads(manifest_path.read_text())
    dry_data, dry_sr = sf.read(str(DRY))
    dry = mono(dry_data)
    dry_sr = float(dry_sr)

    agc_status = "AGC" if apply_agc_enabled else "(AGC bypassed)"
    print(f"dry input:  {DRY.name}  {len(dry)} samples  @ {dry_sr:.1f} Hz")
    print(f"refs:       {REF_DIR}  ({len(manifest)} entries)")
    print(
        f"pipeline:   host->native ({dry_sr:.1f}->{NATIVE_SR:.1f}) -> SOS cascade"
        f" -> {agc_status} -> *boost -> native->host -> gain null"
    )
    print(f"fail at:    rel_null > {FAIL_THRESHOLD_DB:.0f} dB")
    print()
    print(f"{'name':28} {'src':4} {'boost':>6} {'worst_corner':>12} {'gain':>8} {'rel_null':>11}")

    failures: list[tuple[str, str, float]] = []
    import gc
    for name in sorted(manifest):
        entry = manifest[name]
        body, src_type, sr_auth = load_source(entry)
        worst = (None, 0.0, 1.0)  # corner, rel_db, gain
        for label in CORNERS:
            wav = REF_DIR / f"{name}_{label}.wav"
            if not wav.exists():
                continue
            ref = mono(sf.read(str(wav))[0])
            stages = get_corner_stages(body, src_type, label)
            if stages is None:
                continue
            boost = float(entry.get("boost", get_body_boost(body, src_type)))
            pred = render_canonical(
                dry,
                stages,
                src_type,
                sr_auth,
                boost,
                dry_sr,
                apply_agc_enabled=apply_agc_enabled,
            )
            _lag, gain, rms, _peak, ref_rms = best_fit_null(pred, ref)
            rel = db(rms) - db(ref_rms)
            if worst[0] is None or rel > worst[1]:
                worst = (label, rel, gain)
            del ref, pred
        gc.collect()
        if worst[0] is None:
            print(f"{name:28} no refs")
            continue
        boost = float(entry.get("boost", 1.0))
        marker = " " if worst[1] <= FAIL_THRESHOLD_DB else "!"
        tag = {"raw_p2k_skin": "raw", "calibration_re": "cal", "compiled_v1": "c-v1"}.get(src_type, src_type[:4])
        print(
            f"{marker}{name:27} {tag:>4} {boost:6.2f} {worst[0]:>12} "
            f"{worst[2]:8.4f} {worst[1]:+8.1f} dB"
        )
        if worst[1] > FAIL_THRESHOLD_DB:
            failures.append((name, worst[0], worst[1]))

    print()
    if failures:
        print(f"FAIL: {len(failures)}/{len(manifest)} entries exceeded {FAIL_THRESHOLD_DB:.0f} dB")
        for name, corner, rel in failures:
            print(f"  {name}  worst corner {corner}  rel_null={rel:+.1f}dB")
        return 1
    print(f"OK: {len(manifest)}/{len(manifest)} entries null at <= {FAIL_THRESHOLD_DB:.0f} dB")
    return 0


if __name__ == "__main__":
    sys.exit(main())
