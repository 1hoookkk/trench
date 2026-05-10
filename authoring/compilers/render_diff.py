#!/usr/bin/env python3
"""render_diff.py — audibility-tolerance WAV diff for TRENCH render-diff.

Compares two mono IEEE-float32 WAV files. The gate is scoped to "do the
two renders sound the same," not "are the bits identical." All
thresholds are expressed **relative to the reference signal's peak**,
which makes the gate invariant to overall output level — a 10x louder
render with the same error pattern reports the same numbers.

Why relative thresholds:
    The Rust reference runs an f64 cascade internally; the shipping C++
    plugin runs an f32 cascade under `/fp:fast`. Precision and reorder
    differences produce ~-40 to -60 dB error relative to ref-peak —
    inaudible, and the product of the precision choice, not a bug. An
    absolute-dBFS gate would either reject benign drift or miss real
    bugs depending on how loud the test signal happens to be.

Metrics reported every run:
    samples=<N>  sr=<SR>
    ref  peak=<linear> (<dBFS>) rms=<dBFS>
    test peak=<linear> (<dBFS>) rms=<dBFS>
    diff peak=<linear> (<dBFS abs>, <dB rel-ref-peak>)
    diff rms =<linear> (<dBFS abs>, <dB rel-ref-peak>)
    first_drift@<sample>            (only if any sample fails)

Pass requires BOTH:
    diff_peak / ref_peak <= `--max-peak-db`  (default -40 dB rel ref peak)
    diff_rms  / ref_peak <= `--max-rms-db`   (default -60 dB rel ref peak)

Exit codes:
    0  both metrics within tolerance
    1  peak or RMS tolerance exceeded
    2  shape mismatch (channels / rate / length)
    3  NaN in either input

Usage:
    python tools/render_diff.py <ref.wav> <test.wav>
    python tools/render_diff.py <ref.wav> <test.wav> --max-peak-db -30
    python tools/render_diff.py <ref.wav> <test.wav> --max-rms-db -50
"""

import argparse
import math
import struct
import sys
from pathlib import Path


def read_f32_wav(path: Path):
    """Return (samples: list[float], sample_rate: int, channels: int)."""
    data = path.read_bytes()
    if data[0:4] != b"RIFF" or data[8:12] != b"WAVE":
        raise ValueError(f"{path}: not a RIFF/WAVE file")

    fmt = None
    payload = None
    cursor = 12
    while cursor < len(data):
        chunk_id = data[cursor : cursor + 4]
        (chunk_size,) = struct.unpack("<I", data[cursor + 4 : cursor + 8])
        body = data[cursor + 8 : cursor + 8 + chunk_size]
        if chunk_id == b"fmt ":
            fmt = body
        elif chunk_id == b"data":
            payload = body
        cursor += 8 + chunk_size + (chunk_size & 1)  # pad byte if odd

    if fmt is None or payload is None:
        raise ValueError(f"{path}: missing fmt or data chunk")

    (audio_format,) = struct.unpack("<H", fmt[0:2])
    (channels,) = struct.unpack("<H", fmt[2:4])
    (sample_rate,) = struct.unpack("<I", fmt[4:8])
    (bits_per_sample,) = struct.unpack("<H", fmt[14:16])

    EXT = 0xFFFE
    IEEE_FLOAT = 3
    if audio_format == EXT:
        if len(fmt) < 40:
            raise ValueError(f"{path}: EXTENSIBLE fmt chunk truncated")
        (sub_format,) = struct.unpack("<I", fmt[24:28])
        effective_format = sub_format
    else:
        effective_format = audio_format

    if effective_format != IEEE_FLOAT:
        raise ValueError(
            f"{path}: expected IEEE float (fmt=3 or EXTENSIBLE w/ float GUID), "
            f"got fmt={audio_format} effective={effective_format}"
        )
    if bits_per_sample != 32:
        raise ValueError(
            f"{path}: expected 32-bit samples, got {bits_per_sample}-bit"
        )

    n_frames = len(payload) // 4 // channels
    samples = list(struct.unpack(f"<{n_frames * channels}f", payload))
    return samples, sample_rate, channels


def linear_to_db(value: float) -> float:
    """Convert a non-negative linear amplitude to dB. -inf for 0."""
    if value <= 0.0:
        return float("-inf")
    return 20.0 * math.log10(value)


def format_db(db: float) -> str:
    if db == float("-inf"):
        return "-inf"
    return f"{db:+.2f}"


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__.strip().splitlines()[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("ref", type=Path, help="reference WAV (Rust)")
    parser.add_argument("test", type=Path, help="WAV to compare (C++)")
    parser.add_argument(
        "--max-peak-db",
        type=float,
        default=-40.0,
        help="peak |diff| limit in dB relative to ref peak (default -40)",
    )
    parser.add_argument(
        "--max-rms-db",
        type=float,
        default=-60.0,
        help="RMS(diff) limit in dB relative to ref peak (default -60)",
    )
    args = parser.parse_args()

    ref, ref_sr, ref_ch = read_f32_wav(args.ref)
    test, test_sr, test_ch = read_f32_wav(args.test)

    if ref_sr != test_sr:
        print(f"FAIL: sample-rate mismatch {ref_sr} vs {test_sr}", file=sys.stderr)
        return 2
    if ref_ch != test_ch:
        print(f"FAIL: channel-count mismatch {ref_ch} vs {test_ch}", file=sys.stderr)
        return 2
    if len(ref) != len(test):
        print(
            f"FAIL: length mismatch ref={len(ref)} samples, test={len(test)} samples",
            file=sys.stderr,
        )
        return 2

    n = len(ref)
    ref_peak = 0.0
    ref_sum_sq = 0.0
    test_peak = 0.0
    test_sum_sq = 0.0
    max_err = 0.0
    sum_sq_err = 0.0
    first_drift_sample = -1

    for i, (a, b) in enumerate(zip(ref, test)):
        if math.isnan(a) or math.isnan(b):
            print(f"FAIL: NaN at sample {i} (ref={a}, test={b})", file=sys.stderr)
            return 3
        abs_a = abs(a)
        abs_b = abs(b)
        if abs_a > ref_peak:
            ref_peak = abs_a
        if abs_b > test_peak:
            test_peak = abs_b
        ref_sum_sq += a * a
        test_sum_sq += b * b
        d = a - b
        sum_sq_err += d * d
        err = abs(d)
        if err > max_err:
            max_err = err

    ref_rms = math.sqrt(ref_sum_sq / n) if n else 0.0
    test_rms = math.sqrt(test_sum_sq / n) if n else 0.0
    diff_rms = math.sqrt(sum_sq_err / n) if n else 0.0

    # Relative-to-ref-peak thresholds.
    peak_rel_db = linear_to_db(max_err / ref_peak) if ref_peak > 0 else float("-inf")
    rms_rel_db = linear_to_db(diff_rms / ref_peak) if ref_peak > 0 else float("-inf")

    # Locate first drift using the tightest of the two relative bounds
    # (so the sample index matches a real criterion, not a stale one).
    max_peak_linear = ref_peak * (10.0 ** (args.max_peak_db / 20.0)) if ref_peak > 0 else 0.0
    if max_peak_linear > 0:
        for i, (a, b) in enumerate(zip(ref, test)):
            if abs(a - b) > max_peak_linear:
                first_drift_sample = i
                break

    print(f"samples={n}  sr={ref_sr}")
    print(
        f"ref  peak={ref_peak:.9e} ({format_db(linear_to_db(ref_peak))} dBFS) "
        f"rms={format_db(linear_to_db(ref_rms))} dBFS"
    )
    print(
        f"test peak={test_peak:.9e} ({format_db(linear_to_db(test_peak))} dBFS) "
        f"rms={format_db(linear_to_db(test_rms))} dBFS"
    )
    print(
        f"diff peak={max_err:.9e} "
        f"({format_db(linear_to_db(max_err))} dBFS abs, "
        f"{format_db(peak_rel_db)} dB rel-ref-peak)"
    )
    print(
        f"diff rms ={diff_rms:.9e} "
        f"({format_db(linear_to_db(diff_rms))} dBFS abs, "
        f"{format_db(rms_rel_db)} dB rel-ref-peak)"
    )
    print(
        f"gate: peak<={args.max_peak_db:+.2f} dB rel  |  rms<={args.max_rms_db:+.2f} dB rel"
    )

    peak_fail = peak_rel_db > args.max_peak_db
    rms_fail = rms_rel_db > args.max_rms_db

    if peak_fail or rms_fail:
        reasons = []
        if peak_fail:
            reasons.append(f"peak {format_db(peak_rel_db)} > {args.max_peak_db:+.2f}")
        if rms_fail:
            reasons.append(f"rms {format_db(rms_rel_db)} > {args.max_rms_db:+.2f}")
        print(f"FAIL: {' ; '.join(reasons)}", file=sys.stderr)
        if first_drift_sample >= 0:
            print(
                f"      first sample exceeding peak threshold: "
                f"{first_drift_sample} (ref={ref[first_drift_sample]} "
                f"test={test[first_drift_sample]})",
                file=sys.stderr,
            )
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
