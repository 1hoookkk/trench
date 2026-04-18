#!/usr/bin/env python3
"""render_diff.py — bit-accurate WAV diff for the TRENCH render-diff harness.

Compares two mono IEEE-float32 WAV files sample-by-sample and reports the
maximum absolute error. Exit code 0 if max error <= --tolerance (default
0.0 = bit-identical); non-zero otherwise.

Intended pairing:
    Rust reference WAV from `cargo test -p trench-core render_diff_harness`
    C++ offline-render WAV from the JUCE plugin's TRENCH_OFFLINE_RENDER
    mode (upcoming).

Usage:
    python tools/render_diff.py <ref.wav> <test.wav> [--tolerance N]

No external dependencies — parses the RIFF/WAVE container manually so
IEEE-float32 samples (which the stdlib `wave` module does not handle) are
read correctly.
"""

import argparse
import math
import struct
import sys
from pathlib import Path


def read_f32_wav(path: Path):
    """Return (samples: list[float], sample_rate: int, channels: int).

    Accepts WAVE_FORMAT_IEEE_FLOAT (fmt=3) with 32 bits/sample, mono or
    stereo. Rejects anything else — this tool is for the render-diff
    contract, not a general WAV reader.
    """
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

    # WAVE_FORMAT_EXTENSIBLE (0xFFFE) stashes the real format code in the
    # first 4 bytes of the SubFormat GUID at fmt[24:28]. `hound` writes
    # float WAVs as EXTENSIBLE by default; plain IEEE-float (fmt=3) is
    # also accepted.
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


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    parser.add_argument("ref", type=Path, help="reference WAV (Rust)")
    parser.add_argument("test", type=Path, help="WAV to compare (C++)")
    parser.add_argument(
        "--tolerance",
        type=float,
        default=0.0,
        help="max absolute error allowed; default 0.0 (bit-identical)",
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

    max_err = 0.0
    first_drift = None
    drift_count = 0
    for i, (a, b) in enumerate(zip(ref, test)):
        # NaN in either side is an immediate fail; NaN != NaN so abs(a-b) is NaN.
        if math.isnan(a) or math.isnan(b):
            print(
                f"FAIL: NaN at sample {i} (ref={a}, test={b})",
                file=sys.stderr,
            )
            return 3
        err = abs(a - b)
        if err > max_err:
            max_err = err
        if err > args.tolerance:
            drift_count += 1
            if first_drift is None:
                first_drift = (i, a, b, err)

    print(f"samples={len(ref)} max_abs_err={max_err:.9e} tolerance={args.tolerance:.9e}")
    if drift_count:
        i, a, b, err = first_drift
        print(
            f"FAIL: {drift_count} sample(s) exceed tolerance. First at {i}: "
            f"ref={a} test={b} err={err:.9e}",
            file=sys.stderr,
        )
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
